"""Annotation of reaction sets (co-expression modules): KEGG pathway enrichment, ECs and the genera carrying the
reactions. KEGG pathways of a MetaNetX reaction: through its KEGG reaction ids (MetaNetX cross-references) and
its ECs (KEGG pathway - EC links); global and overview maps (map01xxx) left out."""
import os

import numpy as np
import pandas as pd
from scipy import stats

from .store import read_table


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if not n:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(q, 0, 1)
    return out


def reaction_pathways(refdb):
    """DataFrame mnx_id, pathway (mapNNNNN), name"""
    raw = os.path.join(refdb, "raw")
    parts = []
    x = read_table(refdb, "kegg_to_mnx_reaction", required=False)
    lr = read_table(raw, "kegg_link_pathway_reaction", required=False)
    if len(x) and len(lr) and {"path", "rn"} <= set(lr.columns):
        parts.append(x.rename(columns={"source_id": "rn"}).merge(lr, on="rn")[["mnx_id", "path"]])
    re_ = read_table(refdb, "reaction_ec", required=False)
    le = read_table(raw, "kegg_link_pathway_enzyme", required=False)
    if len(re_) and len(le) and {"path", "ec"} <= set(le.columns):
        parts.append(re_.merge(le, on="ec")[["mnx_id", "path"]])
    if not parts:
        return pd.DataFrame(columns=["mnx_id", "pathway", "name"])
    d = pd.concat(parts, ignore_index=True).drop_duplicates().rename(columns={"path": "pathway"})
    d = d[d["pathway"].str.match(r"map0[0-9]") & ~d["pathway"].str.match(r"map01[0-9]{3}")]
    names = read_table(raw, "kegg_list_pathway", required=False)
    if len(names):
        nm = dict(zip(names["kegg_id"].str.replace("path:", "", regex=False), names["description"]))
        d["name"] = d["pathway"].map(nm).fillna("")
    else:
        d["name"] = ""
    return d


def enrichment(module_sets, universe, rp, min_k=2):
    """module_sets: {key: set of mnx ids}; universe: set -> one row per module and pathway (hypergeometric)"""
    rp = rp[rp["mnx_id"].isin(universe)]
    U = set(rp["mnx_id"])
    by_path = rp.groupby("pathway")["mnx_id"].apply(set)
    names = rp.drop_duplicates("pathway").set_index("pathway")["name"]
    rows = []
    for key, s in module_sets.items():
        s = s & U
        if not s:
            continue
        for pw, members in by_path.items():
            k = len(s & members)
            if k >= min_k:
                rows.append((key, pw, names.get(pw, ""), k, len(s), len(members), len(U),
                             float(stats.hypergeom.sf(k - 1, len(U), len(members), len(s)))))
    return pd.DataFrame(rows, columns=["key", "pathway", "name", "k", "module_size", "pathway_size", "universe", "p"])


def annotate_modules(sets, active, refdb, bin_reactions, bins, cond_samples, top=5):
    """sets: {(condition, module): mnx ids}; active: {condition: mnx ids}; -> (module_pathways, module_annotation)"""
    rp = reaction_pathways(refdb)
    rows = []
    for cn in sorted(active):
        ms = {m: s for (c, m), s in sets.items() if c == cn}
        e = enrichment(ms, active[cn], rp)
        if len(e):
            e["q"] = bh(e["p"].values)
            e.insert(0, "condition", cn)
            rows.append(e.rename(columns={"key": "module"}))
    paths = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["condition", "module", "pathway", "name", "k", "module_size", "pathway_size", "universe", "p", "q"])
    rec = read_table(refdb, "reaction_ec", required=False)
    ecs = rec.groupby("mnx_id")["ec"].apply(set) if len(rec) else pd.Series(dtype=object)
    genus = {}
    if len(bin_reactions) and len(bins) and "genus" in bins.columns:
        g = bins.set_index("bin_id")["genus"].replace("", "unclassified").fillna("unclassified")
        br = bin_reactions[["sample_id", "bin_id", "mnx_id"]].copy()
        br["genus"] = br["bin_id"].map(g).fillna("unclassified")
        genus = br
    ann = []
    for (cn, m), s in sorted(sets.items()):
        p = paths[(paths["condition"] == cn) & (paths["module"] == m) & (paths["q"] < 0.05)].sort_values("p")
        e = pd.Series([x for r in s for x in ecs.get(r, ())]).value_counts()
        row = {"condition": cn, "module": m, "reaction_groups": len(s),
               "pathways": "; ".join("%s (%d)" % (n or pw, k) for pw, n, k in
                                     list(zip(p["pathway"], p["name"], p["k"]))[:top]) if len(p) else "",
               "ecs": " ".join("%s(%d)" % (k, v) for k, v in e.head(top).items()), "distinct_ecs": len(e)}
        if len(genus):
            x = genus[genus["mnx_id"].isin(s) & genus["sample_id"].isin(cond_samples.get(cn, ()))]
            c = x.groupby("genus").size().sort_values(ascending=False)
            row["genera"] = " ".join("%s:%.0f%%" % (k, 100 * v / c.sum()) for k, v in c.head(top).items()) \
                if c.sum() else ""
        ann.append(row)
    return paths, pd.DataFrame(ann)
