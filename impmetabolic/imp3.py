"""Readers for one IMP3 sample directory: gene annotations (Mantis attributes in the annotation GFF), gene -> contig,
contig -> bin (MetaWrap), bin quality and GTDB taxonomy, reads per gene (metaG / metaT)."""
import glob
import gzip
import io
import os
import re

import pandas as pd

GFF = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
ANNOT_KEYS = ("enzyme_ec", "kegg_ko", "kegg_reaction", "kegg_module", "pfam", "tcdb")
GTDB_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")


def find(base, *patterns):
    for p in patterns:
        for cand in ((p, p + ".gz") if not p.endswith(".gz") else (p, p[:-3])):
            hits = [h for h in sorted(glob.glob(os.path.join(base, cand))) if os.path.isfile(h) and os.path.getsize(h)]
            if hits:
                return hits[0]
    return None


def _open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return io.TextIOWrapper(gzip.open(path)) if magic == b"\x1f\x8b" else open(path)


def table(path, header="infer", names=None, comment=None):
    if path is None:
        return pd.DataFrame()
    with _open(path) as fh:
        t = fh.read()
    if not t.strip():
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO(t), sep="\t", header=0 if header == "infer" else header, names=names, dtype=str,
                     keep_default_na=False, comment=comment)
    df.columns = [str(c).strip() for c in df.columns]
    return df.apply(lambda c: c.str.strip())


def genes(base, keys=ANNOT_KEYS):
    """CDS of the annotation GFF -> DataFrame gene, contig, start, end, strand, <keys> (comma-separated ids)"""
    p = find(base, GFF)
    if p is None:
        raise FileNotFoundError("no %s in %s" % (GFF, base))
    rx = {k: re.compile(r"(?:^|;)%s=([^;]*)" % re.escape(k)) for k in ("ID",) + tuple(keys)}
    rows = []
    with _open(p) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = f[8]
            m = rx["ID"].search(a)
            if not m:
                continue
            r = [m.group(1), f[0], int(f[3]), int(f[4]), f[6]]
            for k in keys:
                mk = rx[k].search(a)
                r.append(mk.group(1) if mk else "")
            rows.append(r)
    return pd.DataFrame(rows, columns=["gene", "contig", "start", "end", "strand"] + list(keys))


def contig2bin(base):
    t = table(find(base, "Binning/MetaWrap/bins.contigs"), header=None, names=["contig", "bin"])
    return dict(zip(t["contig"], t["bin"])) if len(t) else {}


def bins(base, sample):
    """bin_id, bin, completeness, contamination, quality (MIMAG), GTDB ranks"""
    t = table(find(base, "Stats/all_bin_stats.tsv", "Stats/mgmt/mgmt.bins.tsv"))
    cols = ["bin_id", "sample_id", "bin", "completeness", "contamination", "quality", "classification"] + \
        list(GTDB_RANKS)
    if t.empty or "Bin" not in t.columns:
        return pd.DataFrame(columns=cols)
    comp = pd.to_numeric(t.get("Completeness"), errors="coerce")
    cont = pd.to_numeric(t.get("Contamination"), errors="coerce")
    q = ["high" if c >= 90 and x <= 5 else "medium" if c >= 50 and x <= 10 else "low" for c, x in zip(comp, cont)]
    df = pd.DataFrame({"bin_id": sample + ":" + t["Bin"], "sample_id": sample, "bin": t["Bin"], "completeness": comp,
                       "contamination": cont, "quality": q, "classification": t.get("classification", "")})
    ranks = df["classification"].fillna("").str.split(";")
    for i, r in enumerate(GTDB_RANKS):
        df[r] = [(x[i][3:] if len(x) > i and len(x[i]) > 3 else "") for x in ranks]
    return df[cols]


def gene_reads(base, omics=("mg", "mt")):
    """reads per gene and omic: DataFrame gene, omic, reads (IMP3 Stats/<omic>/gene_read_count.txt)"""
    out = []
    for o in omics:
        p = find(base, "Stats/%s/gene_read_count.txt" % o)
        t = table(p, header=None)
        if t.empty or t.shape[1] < 5:
            continue
        num = t.apply(lambda c: pd.to_numeric(c, errors="coerce"))
        out.append(pd.DataFrame({"gene": t.iloc[:, 3], "omic": o, "reads": num.iloc[:, -1].fillna(0)}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["gene", "omic", "reads"])
