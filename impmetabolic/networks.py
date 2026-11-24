"""Networks of a mapped cohort (imp-metabolic map): community network per sample, network per MAG, pan-reactome
network per taxon; topology of each (global properties, hubs, scope from seed metabolites)."""
import concurrent.futures as cf
import logging
import multiprocessing as mp
import os
import re

import pandas as pd

from . import network, pangenome, topology
from .store import read_table, write_manifest, write_table

log = logging.getLogger(__name__)
SEEDS = os.path.join(os.path.dirname(__file__), "data", "seeds_minimal.tsv")


def resolve_seeds(path, metabolites, currency):
    names = [ln.strip() for ln in open(path) if ln.strip() and not ln.startswith("#")]
    by_name = metabolites.drop_duplicates("name").set_index("name")["mnx_id"]
    ids, missing = set(), []
    for x in names:
        if x.startswith("MNXM") or x in ("WATER",):
            ids.add(x)
        elif x in by_name.index:
            ids.add(by_name[x])
        else:
            missing.append(x)
    if missing:
        log.warning("seed names not found in MetaNetX: %s", ", ".join(missing))
    return ids | set(currency)


_W = {}


def _init(reactions, currency, seed_ids, out, graphml):
    _W.update(eq=reactions, currency=currency, seeds=seed_ids, out=out, graphml=graphml)


def analyse(job):
    """one network: (level, id, sample, reaction ids, hubs?, set, relative expression per reaction or None)
    -> (topology row, hubs table, scope table, keystone table)"""
    level, nid, sample, rxns, hubs_too, rset, mta = job
    net = network.build(rxns, _W["eq"])
    light = level == "bin"                             # many MAG networks: fewer path-length roots
    s = topology.summary(net, _W["currency"], path_sources=50 if light else 200)
    s.update(level=level, id=nid, sample_id=sample, set=rset)
    sc_tab = None
    if _W["seeds"]:
        sc, fired = topology.scope(net, _W["seeds"], currency=_W["currency"])
        s["scope_metabolites"] = len(sc - _W["seeds"])
        s["scope_reactions"] = len(fired)
        s["scope_reaction_share"] = len(fired) / max(1, len(net.reactions()))
        sc_tab = pd.DataFrame({"set": rset, "level": level, "id": nid, "metabolite": sorted(sc - _W["seeds"])})
    h = None
    if hubs_too:
        h = topology.hubs(net, _W["currency"])
        h.insert(0, "id", nid)
        h.insert(0, "level", level)
        h.insert(0, "set", rset)
    if _W["graphml"]:
        d = os.path.join(_W["out"], "graphml", rset, level)
        os.makedirs(d, exist_ok=True)
        net.to_graphml(os.path.join(d, "%s.graphml" % re.sub(r"[^\w.-]+", "_", str(nid))))
    ks = None
    if level == "community":
        ks = topology.keystones(net, _W["currency"], pd.Series(mta) if mta else None)
        ks.insert(0, "id", nid)
        ks.insert(0, "set", rset)
        s["keystones"] = int(ks["keystone"].sum())
    return s, h, sc_tab, ks


def run(mapping, refdb, out, rank=("species", "genus", "family"), tiers=("high", "medium"), seeds=SEEDS, graphml=False, fmt="auto",
        bin_hubs=False, threads=4, sets=("specific", "all", "active"), active_tpm=1.0, extra_currency=(),
        keep_currency=(), currency_mode="pairs", pair_min_reactions=100, extra_pairs=()):
    reactions = read_table(refdb, "reactions", columns=["mnx_id", "equation"])
    mets = read_table(refdb, "metabolites", columns=["mnx_id", "name"], required=False)
    # currency metabolites by their role per reaction (currency.py); mode 'global' = the old fixed list everywhere
    from . import currency as cx
    if currency_mode == "global":
        from .refdb.build import CURRENCY_NAMES
        cur, cur_tab = cx.build(mets, None, "global", extra_singles=tuple(CURRENCY_NAMES) + tuple(extra_currency or ()),
                                keep=keep_currency or ())
    else:
        cur, cur_tab = cx.build(mets, reactions, "pairs", pair_min_reactions, extra_currency or (),
                                extra_pairs or (), keep_currency or ())
    log.info("currency (%s): %d single metabolites, %d cofactor pairs (%d curated, %d frequent in the reference "
             "reactions)", currency_mode, len(cur.singles), len(cur.pairs), int((cur_tab["source"] == "curated").sum()),
             int((cur_tab["source"] != "curated").sum()))
    currency = cur
    seed_ids = resolve_seeds(seeds, mets, cur.singles) if seeds else set()
    com = read_table(mapping, "community_reactions")
    br = read_table(mapping, "bin_reactions")
    bins = read_table(mapping, "bins")
    for c in ("completeness", "contamination"):
        bins[c] = pd.to_numeric(bins[c], errors="coerce")
    q = bins[bins["quality"].isin(tiers)]
    if "confidence" not in com.columns:                          # mappings of 0.1.0.dev0 before the confidence
        com["confidence"] = br["confidence"] = "specific"
        sets = tuple(x for x in sets if x != "specific")
    for df in (com, br):
        for c in ("mt_tpm", "mta"):
            df[c] = pd.to_numeric(df[c], errors="coerce") if c in df.columns else float("nan")
    if "active" in sets and not com["mt_tpm"].notna().any():
        log.info("no metaT: reaction set 'active' left out")
        sets = tuple(x for x in sets if x != "active")
    used = set(com["mnx_id"]) | set(br["mnx_id"])
    reactions = reactions[reactions["mnx_id"].isin(used)]          # workers get only the reactions needed
    ranks = (rank,) if isinstance(rank, str) else tuple(rank)
    jobs, pans, pan_sums, cov = [], [], [], []
    for rset in sets:
        if rset == "active":                                     # expressed reaction groups (metaT TPM)
            c, b = com[com["mt_tpm"] >= active_tpm], br[br["mt_tpm"] >= active_tpm]
        else:
            c = com if rset == "all" else com[com["confidence"] == "specific"]
            b = br if rset == "all" else br[br["confidence"] == "specific"]
        bq = b[b["bin_id"].isin(q["bin_id"])]
        jobs += [("community", sid, sid, list(g["mnx_id"]), True, rset,
                  {k: v for k, v in zip(g["mnx_id"], g["mta"]) if v == v} or None) for sid, g in c.groupby("sample_id")]
        jobs += [("bin", bid, bid.split(":")[0], list(g["mnx_id"]), bin_hubs, rset, None)
                 for bid, g in bq.groupby("bin_id")]
        npan = []
        for rk in ranks:                       # pan-reactomes per species and per higher rank (as AGORA2)
            pan, pan_sum = pangenome.pan_reactome(b, bins, rk, tiers)
            pans.append(pan.assign(set=rset))
            pan_sums.append(pan_sum.assign(set=rset))
            jobs += [("pan_" + rk, tax, "", list(g["mnx_id"]), True, rset, None) for tax, g in pan.groupby("taxon")]
            npan.append("%d %s" % (pan["taxon"].nunique(), rk))
        # share of the community's reactions found in its MAGs
        cv = c.groupby("sample_id")["mnx_id"].nunique().rename("community_reactions").to_frame()
        cv["in_mags"] = bq.groupby("sample_id")["mnx_id"].nunique()
        cv["in_mags"] = cv["in_mags"].fillna(0).astype(int)
        cv["share_in_mags"] = cv["in_mags"] / cv["community_reactions"]
        cov.append(cv.reset_index().assign(set=rset))
        log.info("%s reactions: %d communities, %d MAGs (%s), pan-reactomes of taxa with >= 2 MAGs: %s", rset,
                 c["sample_id"].nunique(), bq["bin_id"].nunique(), ", ".join(tiers), ", ".join(npan))
    log.info("%d networks, %d threads", len(jobs), threads)
    topo, hub_rows, scopes, kst = [], [], [], []
    with cf.ProcessPoolExecutor(max_workers=max(1, threads), mp_context=mp.get_context("spawn"),
                                initializer=_init, initargs=(reactions, currency, seed_ids, out, graphml)) as ex:
        step = max(1, len(jobs) // 20)
        for i, (s, h, sc, ks) in enumerate(ex.map(analyse, jobs, chunksize=4), 1):
            topo.append(s)
            if ks is not None:
                kst.append(ks)
            if h is not None:
                hub_rows.append(h)
            if sc is not None:
                scopes.append(sc)
            if i % step == 0 or i == len(jobs):
                log.info("%d/%d networks (%s, %s %s: %d reactions)", i, len(jobs), s["set"], s["level"], s["id"],
                         s["reactions"])
    hubs = pd.concat(hub_rows, ignore_index=True) if hub_rows else pd.DataFrame()
    if len(hubs) and len(mets):                                   # metabolite names next to the MNXM@compartment ids
        names = mets.drop_duplicates("mnx_id").set_index("mnx_id")["name"]
        hubs.insert(hubs.columns.get_loc("metabolite") + 1, "name",
                    hubs["metabolite"].str.split("@").str[0].map(names).fillna(""))
    tabs = {"topology": pd.DataFrame(topo), "hubs": hubs, "pan_reactome": pd.concat(pans, ignore_index=True),
            "pan_summary": pd.concat(pan_sums, ignore_index=True), "mag_coverage": pd.concat(cov, ignore_index=True),
            "scope": pd.concat(scopes, ignore_index=True) if scopes else pd.DataFrame(),
            "keystones": pd.concat(kst, ignore_index=True) if kst else pd.DataFrame(),
            "currency": cur_tab}
    first = ["set", "level", "id", "sample_id"]
    tabs["topology"] = tabs["topology"][first + [c for c in tabs["topology"].columns if c not in first]]
    rows = {}
    for k, v in tabs.items():
        write_table(v, out, k, fmt)
        rows[k] = int(len(v))
    write_manifest(out, "networks", rows, {"mapping": os.path.abspath(mapping), "refdb": os.path.abspath(refdb),
                                           "ranks": list(ranks), "tiers": list(tiers), "seeds": seeds or "",
                                           "sets": list(sets), "active_tpm": active_tpm,
                                           "currency_mode": currency_mode,
                                           "pair_min_reactions": pair_min_reactions})
    return tabs
