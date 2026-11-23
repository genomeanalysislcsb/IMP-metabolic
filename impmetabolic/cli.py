"""imp-metabolic command line

  imp-metabolic refdb download --data DIR [--sources metanetx,kegg,...]   reference files (needs internet)
  imp-metabolic refdb build --data DIR --out REFDB                         reference tables (Parquet)
  imp-metabolic map --refdb REFDB --imp-root ROOT | --samples SHEET --out MAPPING
                                                                           genes -> reactions, community, MAGs
  imp-metabolic networks --mapping MAPPING --refdb REFDB --out NETS        networks + topology
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
    n.add_argument("--rank", default="species", help="GTDB rank of the pan-reactomes [species]")
    n.add_argument("--quality", default="high,medium", help="MAG tiers used [high,medium]")
    n.add_argument("--seeds", default="default", help="seed metabolites for the scope (file of MetaNetX names / "
                   "MNXM ids; 'default' = minimal glucose medium; 'none' = no scope)")
    n.add_argument("--graphml", action="store_true", help="also write every network as GraphML (Cytoscape)")
    n.add_argument("--bin-hubs", action="store_true", help="hub metabolites also for every MAG network")
    n.add_argument("--sets", default="specific,all", help="reaction sets analysed: specific (links with KEGG-reaction "
                   "evidence or specific ECs) and / or all [specific,all]")
    n.add_argument("--threads", type=int, default=4)
    n.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
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
        run(a.mapping, a.refdb, a.out, a.rank, tuple(a.quality.split(",")), seeds, a.graphml, a.format, a.bin_hubs,
            a.threads, tuple(x for x in a.sets.split(",") if x))
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
