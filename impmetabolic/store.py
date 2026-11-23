"""Tables on disk: Parquet when pyarrow is installed (else TSV.gz), a manifest.json per directory, and DuckDB views
over a directory of tables for fast joins (optional; pandas is used without DuckDB)."""
import datetime
import json
import os

import pandas as pd

from . import __version__

EXTS = (".parquet", ".tsv.gz", ".tsv")


def parquet_ok():
    try:
        import pyarrow  # noqa: F401
        return True
    except ImportError:
        return False


def write_table(df, outdir, name, fmt="auto"):
    os.makedirs(outdir, exist_ok=True)
    fmt = ("parquet" if parquet_ok() else "tsv") if fmt == "auto" else fmt
    path = os.path.join(outdir, name + (".parquet" if fmt == "parquet" else ".tsv.gz"))
    for e in EXTS:                                   # one file per table: remove a copy in the other format
        other = os.path.join(outdir, name + e)
        if other != path and os.path.exists(other):
            os.remove(other)
    if fmt == "parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def table_path(d, name):
    for e in EXTS:
        p = os.path.join(d, name + e)
        if os.path.exists(p):
            return p
    return None


def read_table(d, name, columns=None, required=True):
    p = table_path(d, name)
    if p is None:
        if required:
            raise FileNotFoundError("table %s not found in %s" % (name, d))
        return pd.DataFrame(columns=columns or [])
    if p.endswith(".parquet"):
        return pd.read_parquet(p, columns=columns)
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, usecols=columns)


def tables(d):
    out = {}
    for f in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        for e in EXTS:
            if f.endswith(e):
                out[f[:-len(e)]] = os.path.join(d, f)
    return out


def write_manifest(d, kind, tables_rows, extra=None):
    m = {"kind": kind, "created": datetime.datetime.now().isoformat(timespec="seconds"),
         "impmetabolic": __version__, "tables": tables_rows}
    m.update(extra or {})
    with open(os.path.join(d, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2)
    return m


def read_manifest(d):
    p = os.path.join(d, "manifest.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def duck(d):
    """DuckDB connection with one view per table of directory d (None without duckdb)"""
    try:
        import duckdb
    except ImportError:
        return None
    con = duckdb.connect()
    for name, p in tables(d).items():
        q = p.replace("'", "''")
        src = ("read_parquet('%s')" % q) if p.endswith(".parquet") else \
            ("read_csv_auto('%s', delim='\t', header=true, all_varchar=true)" % q)
        con.execute('CREATE VIEW "%s" AS SELECT * FROM %s' % (name, src))
    return con
