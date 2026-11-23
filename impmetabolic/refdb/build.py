"""Build the reference database: every catalogue file with a loader -> one table (Parquet), plus the derived
tables used for the mapping and the networks:

  reactions        mnx_id, equation, reference, ec_list, is_balanced, is_transport   (MetaNetX reac_prop)
  reaction_ec      mnx_id, ec                                                        (one row per EC)
  metabolites      mnx_id, name, formula, charge, mass, inchikey, ...                (MetaNetX chem_prop)
  ko_ec, rn_ec, ko_rn                                                                (KEGG links)
  <db>_to_mnx_<reaction|metabolite|compartment>  mnx_id, source_id                   (MetaNetX xrefs)
  currency         mnx_id, name: currency metabolites (cofactors, inorganic ions) by name
"""
import importlib
import logging
import os
import re

import pandas as pd

from ..store import write_manifest, write_table
from .catalogue import DATA_CATALOGUE, OPTIONAL_LOADERS

log = logging.getLogger(__name__)

# MetaNetX names of currency metabolites (removed from the metabolite graph before topology, kept in the network)
CURRENCY_NAMES = ("H2O", "H(+)", "ATP", "ADP", "AMP", "GTP", "GDP", "UTP", "UDP", "CTP", "CDP", "NAD(+)", "NADH",
                  "NADP(+)", "NADPH", "FAD", "FADH2", "phosphate", "diphosphate", "CO2", "O2", "NH4(+)",
                  "hydrogencarbonate", "hydrogen peroxide", "coenzyme A", "H2", "Na(+)", "K(+)", "Fe(2+)", "Fe(3+)")
CURRENCY_IDS = ("WATER",)

XREF_SPECS = {
    "metanetx_chem_xref": (("rhea", "metacyc", "bigg", "seed", "kegg", "chebi", "mnx"), "metabolite",
                           ("metabolite", "biomass", "water")),
    "metanetx_reac_xref": (("rhea", "metacyc", "bigg", "seed", "kegg", "mnx"), "reaction", ("reaction", "empty")),
    "metanetx_comp_xref": (("rhea", "bigg", "seed", "mnx"), "compartment", ("compartment", "in", "out", "boundary")),
}


def _loader_ok(loader):
    name = getattr(getattr(loader, "func", loader), "__name__", "")
    lib = OPTIONAL_LOADERS.get(name)
    if not lib:
        return True
    try:
        importlib.import_module(lib)
        return True
    except ImportError:
        log.warning("skipping %s: needs %s", name, lib)
        return False


def load_all(data_dir, sources=None):
    out, files = {}, {}
    for group in DATA_CATALOGUE:
        if sources and group.database not in sources:
            continue
        for fg in group.file_groups:
            d = os.path.join(data_dir, group.database, fg.local_subdir) if fg.local_subdir else \
                os.path.join(data_dir, group.database)
            for f in fg.files:
                p = os.path.join(d, f.name)
                if f.loader is None or not _loader_ok(f.loader):
                    continue
                if not os.path.exists(p):
                    log.warning("missing %s (not downloaded?)", p)
                    continue
                log.info("loading %s", p)
                res = f.loader(__import__("pathlib").Path(p))
                if f.key not in files:
                    files[f.key] = {"file": os.path.relpath(p, data_dir), "bytes": 0, "mtime": 0}
                files[f.key]["bytes"] += os.path.getsize(p)
                files[f.key]["mtime"] = max(files[f.key]["mtime"], int(os.path.getmtime(p)))
                if not isinstance(res, pd.DataFrame):
                    out.update(res)
                elif f.key in out:              # parts of one table (ModelSEED compound_00.tsv, ...)
                    out[f.key] = pd.concat([out[f.key], res], ignore_index=True)
                    files[f.key]["file"] += "," + os.path.basename(p)
                    continue
                else:
                    out[f.key] = res
    return out, files


def metanetx_version(data_dir):
    p = os.path.join(data_dir, "metanetx", "latest", "reac_prop.tsv")
    if not os.path.exists(p):
        return ""
    with open(p) as fh:
        for ln in fh:
            if not ln.startswith("#"):
                break
            m = re.search(r"(MNXref|MetaNetX)\s+[Vv]ersion\s*:?\s*([\w.]+)", ln)
            if m:
                return m.group(2)
    return ""


def explode(df, col, new, sep):
    x = df[[c for c in df.columns if c != col]].assign(**{new: df[col].fillna("").str.split(sep)}).explode(new)
    x[new] = x[new].str.strip()
    return x[x[new].fillna("") != ""]


def derived(t):
    out = {}
    for key, (srcs, main, valid) in XREF_SPECS.items():
        if key not in t:
            continue
        df = t[key]
        for s in srcs:
            b = df[(df["source_db"] == s) & df["type"].isin(valid)][["mnx_id", "source_id"]].drop_duplicates()
            if len(b):
                out["%s_to_mnx_%s" % (s, main)] = b.reset_index(drop=True)
    if "metanetx_reac_prop" in t:
        r = t["metanetx_reac_prop"]
        r = r[r["type"] == "reaction"].rename(columns={"mnx_equation": "equation"})
        out["reactions"] = r[["mnx_id", "equation", "reference", "ec_list", "is_balanced", "is_transport"]] \
            .reset_index(drop=True)
        out["reaction_ec"] = explode(out["reactions"][["mnx_id", "ec_list"]], "ec_list", "ec", ";")[["mnx_id", "ec"]] \
            .drop_duplicates().reset_index(drop=True)
    if "metanetx_chem_prop" in t:
        m = t["metanetx_chem_prop"]
        out["metabolites"] = m[[c for c in ("mnx_id", "name", "reference", "formula", "charge", "mass", "inchikey",
                                           "smiles") if c in m.columns]].reset_index(drop=True)
        cur = out["metabolites"][out["metabolites"]["name"].isin(CURRENCY_NAMES) |
                                 out["metabolites"]["mnx_id"].isin(CURRENCY_IDS)][["mnx_id", "name"]]
        out["currency"] = cur.drop_duplicates().reset_index(drop=True)
    for key, a, b in (("kegg_link_orthology_enzyme", "ko", "ec"), ("kegg_link_reaction_enzyme", "rn", "ec"),
                      ("kegg_link_orthology_reaction", "ko", "rn")):
        if key in t and {a, b} <= set(t[key].columns):
            out["%s_%s" % (a, b)] = t[key][[a, b]].drop_duplicates().reset_index(drop=True)
    return out


def build(data_dir, out_dir, sources=None, fmt="auto", raw=True):
    t, files = load_all(data_dir, sources)
    d = derived(t)
    rows = {}
    for name, df in (list(t.items()) if raw else []) + list(d.items()):
        write_table(df, os.path.join(out_dir, "raw") if name in t and name not in d else out_dir, name, fmt)
        rows[name] = int(len(df))
        log.info("%-40s %10d rows", name, len(df))
    return write_manifest(out_dir, "refdb", rows, {"metanetx_version": metanetx_version(data_dir), "files": files,
                                                   "data_dir": os.path.abspath(data_dir)})
