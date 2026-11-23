"""Networks of a mapped cohort (imp-metabolic map): community network per sample, network per MAG, pan-reactome
network per taxon; topology of each (global properties, hubs, scope from seed metabolites)."""
import logging
import os

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


def run(mapping, refdb, out, rank="species", tiers=("high", "medium"), seeds=SEEDS, graphml=False, fmt="auto",
        bin_hubs=False):
    reactions = read_table(refdb, "reactions", columns=["mnx_id", "equation"])
    currency = set(read_table(refdb, "currency", required=False).get("mnx_id", pd.Series(dtype=str)))
    mets = read_table(refdb, "metabolites", columns=["mnx_id", "name"], required=False)
    seed_ids = resolve_seeds(seeds, mets, currency) if seeds else set()
    com = read_table(mapping, "community_reactions")
    br = read_table(mapping, "bin_reactions")
    bins = read_table(mapping, "bins")
    for c in ("completeness", "contamination"):
        bins[c] = pd.to_numeric(bins[c], errors="coerce")
    topo, hub_rows, scopes = [], [], []

    def analyse(level, nid, sample, rxns, hubs_too):
        net = network.build(rxns, reactions)
        s = topology.summary(net, currency)
        s.update(level=level, id=nid, sample_id=sample)
        if seed_ids:
            sc, fired = topology.scope(net, seed_ids)
            s["scope_metabolites"] = len(sc - seed_ids)
            s["scope_reactions"] = len(fired)
            s["scope_reaction_share"] = len(fired) / max(1, len(net.reactions()))
            scopes.append(pd.DataFrame({"level": level, "id": nid, "metabolite": sorted(sc - seed_ids)}))
        topo.append(s)
        if hubs_too:
            h = topology.hubs(net, currency)
            h.insert(0, "id", nid)
            h.insert(0, "level", level)
            hub_rows.append(h)
        if graphml:
            os.makedirs(os.path.join(out, "graphml", level), exist_ok=True)
            net.to_graphml(os.path.join(out, "graphml", level, "%s.graphml" % str(nid).replace(":", "_")
                                        .replace("/", "_").replace(" ", "_")))

    for sid, g in com.groupby("sample_id"):
        analyse("community", sid, sid, g["mnx_id"], True)
        log.info("community %s: %d reactions", sid, len(g))
    q = bins[bins["quality"].isin(tiers)]
    brq = br[br["bin_id"].isin(q["bin_id"])]
    for bid, g in brq.groupby("bin_id"):
        analyse("bin", bid, bid.split(":")[0], g["mnx_id"], bin_hubs)
    log.info("%d MAG networks (%s)", brq["bin_id"].nunique(), ", ".join(tiers))
    pan, pan_sum = pangenome.pan_reactome(br, bins, rank, tiers)
    for tax, g in pan.groupby("taxon"):
        analyse("pan_" + rank, tax, "", g["mnx_id"], True)
    log.info("%d pan-reactomes (%s with >= 2 MAGs)", pan["taxon"].nunique(), rank)
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
