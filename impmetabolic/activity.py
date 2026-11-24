"""Taxon-resolved abundance and activity: per sample, the metaG and metaT reads / TPM of the genes of each taxon
(GTDB genus, family or species of their MAG) summed per KEGG pathway (or per reaction group), and the relative
expression (metaT TPM / metaG TPM), as the species- and gene-level activity of Villette et al. 2025.

Each gene counts once per pathway (genes linked to several reactions of a pathway are not counted twice); only
genes on MAGs of the chosen quality tiers, and only gene -> reaction links of the chosen confidence.
"""
import logging
import os

import numpy as np
import pandas as pd

from .annotate import reaction_pathways
from .store import read_table, write_manifest, write_table

log = logging.getLogger(__name__)
COLUMNS = ["sample_id", "rank", "taxon", "unit", "feature", "name", "n_genes", "mg_reads", "mt_reads", "mg_tpm",
           "mt_tpm", "mta"]


def _sum(s):
    return s.sum(min_count=1)


def sample_activity(gr, bins, rp, rank="genus", unit="pathway", confidence="specific", tiers=("high", "medium")):
    """gene_reactions of one sample -> rows in COLUMNS"""
    b = bins[bins["quality"].isin(tiers)].set_index("bin_id")[rank].fillna("").replace("", np.nan)
    g = gr[gr["bin_id"].isin(b.index)].copy()
    if confidence != "all":
        g = g[g["confidence"] == confidence]
    g["taxon"] = g["bin_id"].map(b)
    g = g.dropna(subset=["taxon"])
    for c in ("mg_reads", "mt_reads", "mg_tpm", "mt_tpm"):
        g[c] = pd.to_numeric(g[c], errors="coerce") if c in g.columns else np.nan
    if unit == "pathway":
        g = g.merge(rp[["mnx_id", "pathway", "name"]], left_on="mnx_id", right_on="mnx_id")
        g = g.rename(columns={"pathway": "feature"})
    else:
        g["feature"] = g["group_id"] if "group_id" in g.columns else g["mnx_id"]
        g["name"] = ""
    g = g.drop_duplicates(["gene", "taxon", "feature"])             # a gene counts once per feature
    a = g.groupby(["sample_id", "taxon", "feature"], as_index=False).agg(
        name=("name", "first"), n_genes=("gene", "nunique"), mg_reads=("mg_reads", _sum),
        mt_reads=("mt_reads", _sum), mg_tpm=("mg_tpm", _sum), mt_tpm=("mt_tpm", _sum))
    a["mta"] = a["mt_tpm"] / a["mg_tpm"].where(a["mg_tpm"] > 0)
    a["rank"], a["unit"] = rank, unit
    return a[COLUMNS]


def run(mapping, refdb, out, rank="genus", unit="pathway", confidence="specific", tiers=("high", "medium"),
        fmt="auto"):
    bins = read_table(mapping, "bins")
    rp = reaction_pathways(refdb) if unit == "pathway" else None
    if unit == "pathway" and rp.empty:
        raise SystemExit("no KEGG pathway links in the refdb (raw/kegg_link_pathway_*)")
    sdir = os.path.join(mapping, "samples")
    parts = []
    for sid in sorted(os.listdir(sdir)):
        gr = read_table(os.path.join(sdir, sid), "gene_reactions", required=False)
        if gr.empty:
            continue
        a = sample_activity(gr, bins, rp, rank, unit, confidence, tiers)
        parts.append(a)
        log.info("%s: %d %s x %s features (%d taxa)", sid, len(a), rank, unit, a["taxon"].nunique())
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLUMNS)
    write_table(df, out, "taxon_activity", fmt)
    write_manifest(out, "activity", {"taxon_activity": int(len(df))},
                   {"mapping": os.path.abspath(mapping), "refdb": os.path.abspath(refdb), "rank": rank, "unit": unit,
                    "confidence": confidence, "tiers": list(tiers)})
    return df
