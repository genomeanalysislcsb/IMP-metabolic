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
    """one network: (level, id, sample, reaction ids, hubs?) -> (topology row, hubs table, scope table)"""
    level, nid, sample, rxns, hubs_too = job
    net = network.build(rxns, _W["eq"])
    light = level == "bin"                             # many MAG networks: fewer path-length roots
    s = topology.summary(net, _W["currency"], path_sources=50 if light else 200)
    s.update(level=level, id=nid, sample_id=sample)
    sc_tab = None
    if _W["seeds"]:
        sc, fired = topology.scope(net, _W["seeds"])
        s["scope_metabolites"] = len(sc - _W["seeds"])
        s["scope_reactions"] = len(fired)
        s["scope_reaction_share"] = len(fired) / max(1, len(net.reactions()))
        sc_tab = pd.DataFrame({"level": level, "id": nid, "metabolite": sorted(sc - _W["seeds"])})
    h = None
    if hubs_too:
        h = topology.hubs(net, _W["currency"])
        h.insert(0, "id", nid)
        h.insert(0, "level", level)
    if _W["graphml"]:
        d = os.path.join(_W["out"], "graphml", level)
        os.makedirs(d, exist_ok=True)
        net.to_graphml(os.path.join(d, "%s.graphml" % re.sub(r"[^\w.-]+", "_", str(nid))))
    return s, h, sc_tab


def run(mapping, refdb, out, rank="species", tiers=("high", "medium"), seeds=SEEDS, graphml=False, fmt="auto",
        bin_hubs=False, threads=4):
    reactions = read_table(refdb, "reactions", columns=["mnx_id", "equation"])
    currency = set(read_table(refdb, "currency", required=False).get("mnx_id", pd.Series(dtype=str)))
    mets = read_table(refdb, "metabolites", columns=["mnx_id", "name"], required=False)
    seed_ids = resolve_seeds(seeds, mets, currency) if seeds else set()
    com = read_table(mapping, "community_reactions")
    br = read_table(mapping, "bin_reactions")
    bins = read_table(mapping, "bins")
    for c in ("completeness", "contamination"):
        bins[c] = pd.to_numeric(bins[c], errors="coerce")
    q = bins[bins["quality"].isin(tiers)]
    brq = br[br["bin_id"].isin(q["bin_id"])]
    pan, pan_sum = pangenome.pan_reactome(br, bins, rank, tiers)
    used = set(com["mnx_id"]) | set(br["mnx_id"])
    reactions = reactions[reactions["mnx_id"].isin(used)]          # workers get only the reactions needed
    jobs = [("community", sid, sid, list(g["mnx_id"]), True) for sid, g in com.groupby("sample_id")]
    jobs += [("bin", bid, bid.split(":")[0], list(g["mnx_id"]), bin_hubs) for bid, g in brq.groupby("bin_id")]
    jobs += [("pan_" + rank, tax, "", list(g["mnx_id"]), True) for tax, g in pan.groupby("taxon")]
    log.info("%d networks: %d communities, %d MAGs (%s), %d pan-reactomes (%s with >= 2 MAGs); %d threads",
             len(jobs), com["sample_id"].nunique(), brq["bin_id"].nunique(), ", ".join(tiers),
             pan["taxon"].nunique(), rank, threads)
    topo, hub_rows, scopes = [], [], []
    with cf.ProcessPoolExecutor(max_workers=max(1, threads), mp_context=mp.get_context("spawn"),
                                initializer=_init, initargs=(reactions, currency, seed_ids, out, graphml)) as ex:
        step = max(1, len(jobs) // 20)
        for i, (s, h, sc) in enumerate(ex.map(analyse, jobs, chunksize=4), 1):
            topo.append(s)
            if h is not None:
                hub_rows.append(h)
            if sc is not None:
                scopes.append(sc)
            if i % step == 0 or i == len(jobs):
                log.info("%d/%d networks (%s %s: %d reactions)", i, len(jobs), s["level"], s["id"], s["reactions"])
    # share of the community's reactions found in its MAGs
    cov = com.groupby("sample_id")["mnx_id"].nunique().rename("community_reactions").to_frame()
    cov["in_mags"] = brq.groupby("sample_id")["mnx_id"].nunique()
    cov["in_mags"] = cov["in_mags"].fillna(0).astype(int)
    cov["share_in_mags"] = cov["in_mags"] / cov["community_reactions"]
    tabs = {"topology": pd.DataFrame(topo), "hubs": pd.concat(hub_rows, ignore_index=True) if hub_rows else
            pd.DataFrame(), "pan_reactome": pan, "pan_summary": pan_sum, "mag_coverage": cov.reset_index(),
            "scope": pd.concat(scopes, ignore_index=True) if scopes else pd.DataFrame()}
    first = ["level", "id", "sample_id"]
    tabs["topology"] = tabs["topology"][first + [c for c in tabs["topology"].columns if c not in first]]
    rows = {}
    for k, v in tabs.items():
        write_table(v, out, k, fmt)
        rows[k] = int(len(v))
    write_manifest(out, "networks", rows, {"mapping": os.path.abspath(mapping), "refdb": os.path.abspath(refdb),
                                           "rank": rank, "tiers": list(tiers), "seeds": seeds or "",
                                           "currency": sorted(currency)})
    return tabs
