"""imp-metabolic command line

  imp-metabolic refdb download --data DIR [--sources metanetx,kegg,...]   reference files (needs internet)
  imp-metabolic refdb build --data DIR --out REFDB                         reference tables (Parquet)
  imp-metabolic map --refdb REFDB --imp-root ROOT | --samples SHEET --out MAPPING
                                                                           genes -> reactions, community, MAGs
  imp-metabolic networks --mapping MAPPING --refdb REFDB --out NETS        networks + topology
  imp-metabolic activity --mapping MAPPING --refdb REFDB --out ACT        genus x pathway activity (MT/MG)
  imp-metabolic coexpression --mapping MAPPING --out COEX [--metadata M --traits group]
                                                                           co-expression modules (WGCNA)
"""
import argparse
import logging
import os
import sys

from . import __version__


def _samples(a):
    import pandas as pd
    if a.samples:
        s = pd.read_csv(a.samples, sep="," if a.samples.endswith(".csv") else "\t", dtype=str, comment="#")
        if "sample_id" not in s.columns or "imp_dir" not in s.columns:
            sys.exit("--samples needs the columns sample_id and imp_dir")
        base = os.path.dirname(os.path.abspath(a.samples))
        return [(r.sample_id, r.imp_dir if os.path.isabs(r.imp_dir) else os.path.join(base, r.imp_dir))
                for r in s.itertuples()]
    out = []
    for d in sorted(os.listdir(a.imp_root)):
        p = os.path.join(a.imp_root, d)
        if os.path.isdir(p) and os.path.exists(os.path.join(p, "Analysis")):
            out.append((d, p))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-metabolic", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    rd = sub.add_parser("refdb", help="reference database")
    rsub = rd.add_subparsers(dest="action", required=True)
    d = rsub.add_parser("download", help="download the reference files")
    d.add_argument("--data", required=True, help="directory for the downloaded files")
    d.add_argument("--sources", default="", help="comma-separated databases [all of the catalogue]")
    d.add_argument("--workers", type=int, default=5)
    d.add_argument("--force", action="store_true", help="download files that exist already again")
    b = rsub.add_parser("build", help="build the reference tables")
    b.add_argument("--data", required=True, help="directory of the downloaded files")
    b.add_argument("--out", required=True, help="reference database directory")
    b.add_argument("--sources", default="")
    b.add_argument("--no-raw", action="store_true", help="only the derived tables, not every source table")
    b.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")

    m = sub.add_parser("map", help="genes -> MetaNetX reactions per IMP3 sample; community and MAG reactions")
    src = m.add_mutually_exclusive_group(required=True)
    src.add_argument("--imp-root", help="directory whose sub-directories are IMP3 sample outputs")
    src.add_argument("--samples", help="TSV/CSV with sample_id, imp_dir")
    m.add_argument("--refdb", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--methods", default="ec,rn", help="evidence paths: ec, rn or both [ec,rn]")
    m.add_argument("--partial-ec", action="store_true", help="also use partial ECs (1.3.99.-)")
    m.add_argument("--max-ec-reactions", type=int, default=10, help="an EC link is 'specific' when the EC stands "
                   "for at most this many MetaNetX reactions, else 'broad' [10]")
    m.add_argument("--threads", type=int, default=4)
    m.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")

    n = sub.add_parser("networks", help="community, MAG and pan-reactome networks and their topology")
    n.add_argument("--mapping", required=True, help="output of imp-metabolic map")
    n.add_argument("--refdb", required=True)
    n.add_argument("--out", required=True)
    n.add_argument("--ranks", default="species,genus,family", help="GTDB ranks of the pan-reactomes (taxa with >= 2 "
                   "MAGs) [species,genus,family]")
    n.add_argument("--quality", default="high,medium", help="MAG tiers used [high,medium]")
    n.add_argument("--seeds", default="default", help="seed metabolites for the scope (file of MetaNetX names / "
                   "MNXM ids; 'default' = minimal glucose medium; 'none' = no scope)")
    n.add_argument("--graphml", action="store_true", help="also write every network as GraphML (Cytoscape)")
    n.add_argument("--bin-hubs", action="store_true", help="hub metabolites also for every MAG network")
    n.add_argument("--sets", default="specific,all,active", help="reaction sets analysed: specific (KEGG-reaction "
                   "evidence or specific ECs), all, active (reaction groups expressed: metaT TPM >= --active-tpm; "
                   "left out without metaT) [specific,all,active]")
    n.add_argument("--currency-mode", choices=("pairs", "global"), default="pairs", help="pairs: cofactors are "
                   "currency only where a reaction uses them (cofactor pairs, carriers), not in their own synthesis "
                   "or degradation; global: a fixed list removed everywhere [pairs]")
    n.add_argument("--pair-min-reactions", type=int, default=100, help="metabolite pairs on opposite sides of at least "
                   "this many reference reactions are cofactor pairs, next to the curated ones [100]")
    n.add_argument("--currency-pairs", help="file of further cofactor pairs: two MetaNetX names or MNXM ids per line, "
                   "tab-separated")
    n.add_argument("--currency", help="file of further single currency metabolites (MetaNetX names or MNXM ids)")
    n.add_argument("--keep-currency", default="", help="names to take out of the currency list (comma-separated), "
                   "e.g. 'S-adenosyl-L-methionine,S-adenosyl-L-homocysteine'")
    n.add_argument("--active-tpm", type=float, default=1.0, help="metaT TPM of a reaction group to count as active [1]")
    n.add_argument("--threads", type=int, default=4)
    n.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    c = sub.add_parser("coexpression", help="co-expression network (WGCNA) of reaction groups across the cohort")
    c.add_argument("--mapping", required=True, help="output of imp-metabolic map")
    c.add_argument("--out", required=True)
    c.add_argument("--value", choices=("auto", "mta", "mt", "mg"), default="auto", help="auto: metaT / metaG ratio "
                   "with metaG and metaT, metaT TPM with metaT only, metaG TPM with metaG only [auto]")
    c.add_argument("--reactions", choices=("all", "specific"), default="all", help="reaction groups used [all]")
    c.add_argument("--min-prevalence", type=float, default=0.5, help="share of samples a reaction group must be "
                   "present in [0.5]")
    c.add_argument("--pseudocount", type=float, help="TPM added to metaT and metaG for the ratio [half the smallest "
                   "positive TPM]")
    c.add_argument("--transform", choices=("yeo-johnson", "log2", "none"), default="yeo-johnson")
    c.add_argument("--correlation", choices=("pearson", "spearman"), default="pearson")
    c.add_argument("--power", type=int, help="soft-threshold power [lowest power with scale-free R2 >= --r2]")
    c.add_argument("--r2", type=float, default=0.8)
    c.add_argument("--min-module-size", type=int, default=20)
    c.add_argument("--tree-cut", choices=("dynamic", "static"), default="dynamic")
    c.add_argument("--reassign-kme", type=float, default=0.5, help="after the tree cut, unassigned reaction groups "
                   "join the module with the highest kME if >= this (0 = off) [0.5]")
    c.add_argument("--deep-split", type=int, choices=range(5), default=2, help="0 (few large modules) .. 4 [2]")
    c.add_argument("--merge-diss", type=float, default=0.18, help="merge modules whose eigengenes correlate above "
                   "1 - this [0.18]")
    c.add_argument("--metadata", help="TSV/CSV with sample_id and trait columns (ids as in the mapping)")
    c.add_argument("--block", help="metadata column to test the traits within as well (e.g. site, when groups and "
                   "sites are confounded)")
    c.add_argument("--adjust", help="metadata column whose effect is removed before the network (each reaction group "
                   "centred within each level), e.g. site")
    c.add_argument("--by", help="condition-specific modules: metadata column of the conditions (e.g. group); one "
                   "network per condition (within --block levels), module overlap between conditions")
    c.add_argument("--active-tpm", type=float, default=1.0, help="with --by: metaT TPM for a reaction group to be "
                   "active in a sample [1]")
    c.add_argument("--permutations", type=int, default=0, help="with --by: shuffle the conditions within each block "
                   "this many times -> null distribution of the module classes [0]")
    c.add_argument("--threads", type=int, default=1, help="with --permutations: worker processes [1]")
    c.add_argument("--refdb", help="with --by: annotate the modules (KEGG pathway enrichment, ECs, genera)")
    c.add_argument("--min-jaccard", type=float, default=0.2, help="with --by: Jaccard of mutual best partners to call "
                   "a module conserved [0.2]")
    c.add_argument("--traits", default="", help="metadata columns to associate with the modules (comma-separated)")
    c.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    t = sub.add_parser("activity", help="taxon-resolved abundance and activity (genus x KEGG pathway, metaT / metaG)")
    t.add_argument("--mapping", required=True)
    t.add_argument("--refdb", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--rank", default="genus", help="GTDB rank of the MAGs [genus]")
    t.add_argument("--unit", choices=("pathway", "reaction"), default="pathway", help="feature per taxon [pathway]")
    t.add_argument("--reactions", choices=("specific", "all"), default="specific")
    t.add_argument("--quality", default="high,medium")
    t.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(name)s] %(message)s", datefmt="%H:%M:%S")

    if a.cmd == "refdb" and a.action == "download":
        from .refdb.download import download
        res = download(a.data, [x for x in a.sources.split(",") if x] or None, a.workers, not a.force)
        bad = [r for r in res if not r["ok"]]
        for r in bad:
            logging.error("failed %s: %s", r["url"], r["error"])
        logging.info("%d files downloaded, %d failed", len(res) - len(bad), len(bad))
        return 1 if bad else 0
    if a.cmd == "refdb" and a.action == "build":
        from .refdb.build import build
        mf = build(a.data, a.out, [x for x in a.sources.split(",") if x] or None, a.format, not a.no_raw)
        ver = mf.get("metanetx_version") or "version not found; header: " + " | ".join(mf.get("metanetx_header", [])[:3])
        logging.info("refdb %s: %d tables, MetaNetX %s", a.out, len(mf["tables"]), ver)
        return 0
    if a.cmd == "map":
        from .mapping import map_samples
        s = _samples(a)
        if not s:
            sys.exit("no samples found")
        err = map_samples(s, a.refdb, a.out, [x for x in a.methods.split(",") if x], a.partial_ec, a.threads,
                          a.format, a.max_ec_reactions)
        return 2 if err else 0
    if a.cmd == "networks":
        from .networks import SEEDS, run
        seeds = SEEDS if a.seeds == "default" else (None if a.seeds == "none" else a.seeds)
        run(a.mapping, a.refdb, a.out, tuple(x for x in a.ranks.split(",") if x), tuple(a.quality.split(",")), seeds, a.graphml, a.format, a.bin_hubs,
            a.threads, tuple(x for x in a.sets.split(",") if x), a.active_tpm,
            [ln.strip() for ln in open(a.currency) if ln.strip() and not ln.startswith("#")] if a.currency else (),
            [x.strip() for x in a.keep_currency.split(",") if x.strip()], a.currency_mode, a.pair_min_reactions,
            [tuple(ln.rstrip("\n").split("\t")[:2]) for ln in open(a.currency_pairs)
             if ln.strip() and not ln.startswith("#") and "\t" in ln] if a.currency_pairs else ())
        return 0
    if a.cmd == "activity":
        from .activity import run as act
        act(a.mapping, a.refdb, a.out, a.rank, a.unit, a.reactions, tuple(a.quality.split(",")), a.format)
        return 0
    if a.cmd == "coexpression":
        import pandas as pd
        from .coexpression import run as coex
        meta = None
        if a.metadata:
            meta = pd.read_csv(a.metadata, sep="," if a.metadata.endswith(".csv") else "\t", dtype=str)
            if "sample_id" not in meta.columns:
                meta = meta.rename(columns={meta.columns[0]: "sample_id"})
        if a.by:
            if meta is None:
                sys.exit("--by needs --metadata")
            from .coexpression import run_by
            run_by(a.mapping, a.out, a.by, a.block, meta, a.value, a.reactions, a.min_prevalence, a.active_tpm,
                   a.pseudocount, a.transform, a.correlation, a.power, a.r2, a.min_module_size, a.merge_diss,
                   a.tree_cut, a.deep_split, a.reassign_kme, 0.3, a.min_jaccard, a.format, a.permutations,
                   a.threads, a.refdb)
            return 0
        coex(a.mapping, a.out, a.value, a.reactions, a.min_prevalence, a.pseudocount, a.transform, a.correlation,
             a.power, a.r2, a.min_module_size, a.merge_diss, meta, [x for x in a.traits.split(",") if x], a.format,
             a.tree_cut, a.deep_split, a.block, a.reassign_kme, adjust=a.adjust)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
