"""Genes -> MetaNetX reactions, per IMP3 sample, then reactions of the community and of each MAG.

Two evidence paths (both by default):
  ec  the cascade of the metabolic-network prototype: the gene's EC numbers; without them ECs via its KOs; without
      those ECs via its KEGG reactions; EC -> MetaNetX reactions (MetaNetX reac_prop classifs). Partial ECs
      (1.3.99.-) are left out unless --partial-ec.
  rn  KEGG reactions of the gene (Mantis kegg_reaction) and of its KOs (KEGG ko -> rn) -> MetaNetX reactions
      through MetaNetX's KEGG cross-references; more specific than ECs (one EC often covers many reactions).
"""
import concurrent.futures as cf
import logging
import multiprocessing as mp
import os

import pandas as pd

from . import imp3
from .store import read_table, write_manifest, write_table

log = logging.getLogger(__name__)
REF_TABLES = ("reaction_ec", "ko_ec", "rn_ec", "ko_rn", "kegg_to_mnx_reaction")


def load_ref(refdb):
    r = {k: read_table(refdb, k, required=k != "kegg_to_mnx_reaction") for k in REF_TABLES}
    r["kegg_to_mnx_reaction"] = r["kegg_to_mnx_reaction"].rename(columns={"source_id": "rn"})
    return r


def _ids(genes, col, new):
    x = genes[["gene", col]]
    x = x[x[col].fillna("") != ""].assign(**{new: lambda d: d[col].str.split(",")}).explode(new)
    x[new] = x[new].str.strip().str.replace(r"^(ko|ec|rn|EC):", "", regex=True)
    return x[x[new] != ""][["gene", new]].drop_duplicates()


def map_genes(genes, ref, methods=("ec", "rn"), partial_ec=False):
    """genes (imp3.genes) -> DataFrame gene, mnx_id, evidence (path:source, e.g. ec:KO), via (EC or RN id)"""
    parts = []
    ec = _ids(genes, "enzyme_ec", "ec")
    ko = _ids(genes, "kegg_ko", "ko")
    rn = _ids(genes, "kegg_reaction", "rn")
    if "ec" in methods:
        direct = ec.assign(source="EC")
        via_ko = ko[~ko["gene"].isin(direct["gene"])].merge(ref["ko_ec"], on="ko")[["gene", "ec"]] \
            .drop_duplicates().assign(source="KO")
        done = set(direct["gene"]) | set(via_ko["gene"])
        via_rn = rn[~rn["gene"].isin(done)].merge(ref["rn_ec"], on="rn")[["gene", "ec"]] \
            .drop_duplicates().assign(source="RN")
        g_ec = pd.concat([direct, via_ko, via_rn], ignore_index=True)
        if not partial_ec:
            g_ec = g_ec[~g_ec["ec"].str.contains("-", regex=False)]
        m = g_ec.merge(ref["reaction_ec"], on="ec")
        parts.append(pd.DataFrame({"gene": m["gene"], "mnx_id": m["mnx_id"], "evidence": "ec:" + m["source"],
                                   "via": m["ec"]}))
    if "rn" in methods and len(ref["kegg_to_mnx_reaction"]):
        g_rn = pd.concat([rn.assign(source="RN"),
                          ko.merge(ref["ko_rn"], on="ko")[["gene", "rn"]].drop_duplicates().assign(source="KO")],
                         ignore_index=True).drop_duplicates(["gene", "rn"])
        m = g_rn.merge(ref["kegg_to_mnx_reaction"], on="rn")
        parts.append(pd.DataFrame({"gene": m["gene"], "mnx_id": m["mnx_id"], "evidence": "rn:" + m["source"],
                                   "via": m["rn"]}))
    if not parts:
        return pd.DataFrame(columns=["gene", "mnx_id", "evidence", "via"])
    out = pd.concat(parts, ignore_index=True)
    out = out.groupby(["gene", "mnx_id"], as_index=False).agg(
        evidence=("evidence", lambda s: ",".join(sorted(set(s)))), via=("via", lambda s: ",".join(sorted(set(s)))))
    return out


def map_sample(base, sample, ref, methods=("ec", "rn"), partial_ec=False):
    """-> dict of tables: gene_reactions, community_reactions, bin_reactions, bins, stats"""
    g = imp3.genes(base)
    gr = map_genes(g, ref, methods, partial_ec)
    c2b = imp3.contig2bin(base)
    loc = g[["gene", "contig"]].assign(bin=lambda d: d["contig"].map(c2b))
    gr = gr.merge(loc, on="gene", how="left")
    gr.insert(0, "sample_id", sample)
    gr["bin_id"] = [sample + ":" + b if isinstance(b, str) else "" for b in gr["bin"]]
    gr = gr.drop(columns=["bin"])
    reads = imp3.gene_reads(base)
    if len(reads):
        w = reads.pivot_table(index="gene", columns="omic", values="reads", aggfunc="sum")
        for o in ("mg", "mt"):
            gr[o + "_reads"] = gr["gene"].map(w[o]) if o in w.columns else float("nan")
    else:
        gr["mg_reads"] = gr["mt_reads"] = float("nan")

    def agg(df, keys):
        # reads of a gene count once per reaction it maps to
        a = df.groupby(keys + ["mnx_id"], as_index=False).agg(n_genes=("gene", "nunique"),
                                                               mg_reads=("mg_reads", "sum"),
                                                               mt_reads=("mt_reads", "sum"))
        return a
    com = agg(gr, ["sample_id"])
    binr = agg(gr[gr["bin_id"] != ""], ["sample_id", "bin_id"])
    st = {"sample_id": sample, "genes": len(g), "genes_with_reaction": int(gr["gene"].nunique()),
          "genes_binned_with_reaction": int(gr.loc[gr["bin_id"] != "", "gene"].nunique()),
          "reactions": int(gr["mnx_id"].nunique()), "bins": int(binr["bin_id"].nunique())}
    for e in ("ec:EC", "ec:KO", "ec:RN", "rn:RN", "rn:KO"):
        st["genes_" + e.replace(":", "_")] = int(gr.loc[gr["evidence"].str.contains(e, regex=False), "gene"].nunique())
    return {"gene_reactions": gr, "community_reactions": com, "bin_reactions": binr,
            "bins": imp3.bins(base, sample), "stats": pd.DataFrame([st])}


def _job(args):
    sid, base, refdb, out, methods, partial, fmt = args
    try:
        ref = load_ref(refdb)
        res = map_sample(base, sid, ref, methods, partial)
        write_table(res.pop("gene_reactions"), os.path.join(out, "samples", sid), "gene_reactions", fmt)
        return sid, res, None
    except Exception as e:                                  # noqa: BLE001 (one broken sample must not stop the run)
        return sid, None, "%s: %s" % (type(e).__name__, e)


def map_samples(samples, refdb, out, methods=("ec", "rn"), partial_ec=False, threads=4, fmt="auto"):
    """samples: list of (sample_id, imp3_dir) -> cohort tables in out/ and gene tables in out/samples/<id>/"""
    parts = {k: [] for k in ("community_reactions", "bin_reactions", "bins", "stats")}
    errors = {}
    jobs = [(s, d, refdb, out, tuple(methods), partial_ec, fmt) for s, d in samples]
    with cf.ProcessPoolExecutor(max_workers=max(1, threads), mp_context=mp.get_context("spawn")) as ex:
        for i, (sid, res, err) in enumerate(ex.map(_job, jobs), 1):
            if err:
                errors[sid] = err
                log.error("%d/%d %s FAILED %s", i, len(jobs), sid, err)
                continue
            for k in parts:
                parts[k].append(res[k])
            s = res["stats"].iloc[0]
            log.info("%d/%d %s: %d of %d genes -> %d reactions (%d binned genes, %d bins)", i, len(jobs), sid,
                     s["genes_with_reaction"], s["genes"], s["reactions"], s["genes_binned_with_reaction"], s["bins"])
    rows = {}
    for k, v in parts.items():
        df = pd.concat(v, ignore_index=True) if v else pd.DataFrame()
        write_table(df, out, k, fmt)
        rows[k] = int(len(df))
    write_manifest(out, "mapping", rows, {"refdb": os.path.abspath(refdb), "methods": list(methods),
                                          "partial_ec": partial_ec, "failed_samples": errors})
    return errors
