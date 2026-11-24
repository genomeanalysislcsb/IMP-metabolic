"""co-expression (WGCNA) on a synthetic cohort with three planted modules and a group effect on one of them"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from impmetabolic import coexpression as cx  # noqa: E402
from impmetabolic.store import read_table, write_table  # noqa: E402


def cohort(tmp, n=40, with_mt=True, with_mg=True, seed=3):
    rng = np.random.default_rng(seed)
    samples = ["S%02d" % i for i in range(n)]
    group = np.array(["case"] * (n // 2) + ["control"] * (n - n // 2))
    f = rng.normal(size=(3, n))
    f[0] += np.where(group == "case", 1.5, -1.5)                 # module A follows the group
    rows = []
    for j in range(150):
        mod = j // 30 if j < 90 else None                            # 3 x 30 planted, 60 noise
        act = (1.2 * f[mod] if mod is not None else 0) + rng.normal(scale=0.5, size=n)
        mg = rng.lognormal(3, 0.3, size=n)
        mt = mg * np.exp(act)                                        # relative expression carries the signal
        for i, s in enumerate(samples):
            if j < 5:                                                # substrate variants: same genes as MNXR000
                rows.append({"sample_id": s, "mnx_id": "MNXV%03d" % j, "confidence": "broad",
                             "mg_tpm": mg[i] if with_mg else np.nan, "mt_tpm": mt[i] if with_mt else np.nan})
            if j >= 140 and i % 3:                                   # rare reaction groups: filtered out
                continue
            rows.append({"sample_id": s, "mnx_id": "MNXR%03d" % j, "confidence": "specific",
                         "mg_tpm": mg[i] if with_mg else np.nan, "mt_tpm": mt[i] if with_mt else np.nan})
    d = tmp / "map"
    write_table(pd.DataFrame(rows), str(d), "community_reactions", "tsv")
    return str(d), pd.DataFrame({"sample_id": samples, "group": group})


def test_values_auto(tmp_path):
    for mg, mt, want in ((True, True, "mta"), (False, True, "mt"), (True, False, "mg")):
        d, _ = cohort(tmp_path / want, with_mg=mg, with_mt=mt)
        g, t = cx.matrices(read_table(d, "community_reactions"))
        X, v, pc, prev = cx.values(g, t, "auto")
        assert v == want and X.shape[1] == 145                       # 10 rare groups below 50 % prevalence
        Xc, mem = cx.collapse_duplicates(X)
        assert Xc.shape[1] == 140 and mem["MNXR000"] == "MNXR000 MNXV000"


def test_planted_modules(tmp_path):
    d, meta = cohort(tmp_path)
    tabs = cx.run(d, str(tmp_path / "coex"), metadata=meta, trait_cols=["group"], min_module_size=10, fmt="tsv")
    m = tabs["modules"].set_index("mnx_id")["module"]
    for k in range(3):                                               # each planted module mostly in one module
        top = m[["MNXR%03d" % j for j in range(30 * k, 30 * k + 30)]].value_counts()
        assert top.index[0] != "M0" and top.iloc[0] >= 25
    mods = [m[["MNXR%03d" % j for j in range(30 * k, 30 * k + 30)]].value_counts().index[0] for k in range(3)]
    assert len(set(mods)) == 3
    assert (m[["MNXR%03d" % j for j in range(90, 140)]] == "M0").mean() > 0.8     # noise unassigned
    tr = tabs["module_traits"]
    sig = tr[(tr["test"] == "mann-whitney") & (tr["q"] < 0.01) & (tr["block"] == "all")]["module"]
    assert set(sig) == {mods[0]}                                     # only module A differs between the groups
    assert tabs["modules"]["intramodular_hub"].sum() >= 9
    assert set(read_table(str(tmp_path / "coex"), "module_summary")["module"]) == set(mods)
    assert not tabs["modules"]["mnx_id"].str.startswith("MNXV").any()             # collapsed into MNXR0xx
    assert (tabs["modules"].set_index("mnx_id").loc["MNXR001", "n_same_profile"]) == 2


def test_tree_cuts(tmp_path):
    """dynamic and static cut both separate the three planted modules"""
    d, meta = cohort(tmp_path)
    for cut in ("dynamic", "static"):
        tabs = cx.run(d, str(tmp_path / cut), min_module_size=10, fmt="tsv", tree_cut=cut)
        m = tabs["modules"].set_index("mnx_id")["module"]
        top = [m[["MNXR%03d" % j for j in range(30 * k, 30 * k + 30)]].value_counts().index[0] for k in range(3)]
        assert len(set(top)) == 3 and "M0" not in top, cut


def test_traits_within_block(tmp_path):
    """group confounded with site: within-site tests are reported next to the pooled ones"""
    d, meta = cohort(tmp_path)
    meta["site"] = np.where(np.arange(len(meta)) % 2 == 0, "A", "B")
    tabs = cx.run(d, str(tmp_path / "b"), metadata=meta, trait_cols=["group"], block="site", min_module_size=10,
                  fmt="tsv")
    tr = tabs["module_traits"]
    assert set(tr["block"]) == {"all", "site=A", "site=B"}
    m = tabs["modules"].set_index("mnx_id")["module"]
    mod_a = m[["MNXR%03d" % j for j in range(30)]].value_counts().index[0]
    w = tr[(tr["module"] == mod_a) & (tr["test"] == "mann-whitney")]
    assert (w["p"] < 0.01).all() and len(w) == 3                  # pooled and in each site


def test_kme_reassignment(tmp_path):
    """members of a planted module left out by the cut are reassigned by kME; noise stays unassigned"""
    d, meta = cohort(tmp_path)
    a = cx.run(d, str(tmp_path / "a"), min_module_size=10, fmt="tsv", reassign_kme=0)
    b = cx.run(d, str(tmp_path / "b"), min_module_size=10, fmt="tsv", reassign_kme=0.5)
    planted = ["MNXR%03d" % j for j in range(90)]
    ma, mb = (x["modules"].set_index("mnx_id")["module"] for x in (a, b))
    assert (mb[planted] != "M0").sum() >= (ma[planted] != "M0").sum()
    assert (mb[["MNXR%03d" % j for j in range(90, 140)]] == "M0").mean() > 0.8


def test_adjust_block(tmp_path):
    """a site shift on all reaction groups forms no module after --adjust site"""
    d, meta = cohort(tmp_path)
    meta["site"] = np.where(np.arange(len(meta)) < 20, "A", "B")
    com = read_table(d, "community_reactions")
    shift = com["sample_id"].map(dict(zip(meta["sample_id"], np.where(meta["site"] == "A", 8.0, 1.0))))
    com["mt_tpm"] = com["mt_tpm"].astype(float) * shift
    write_table(com, d, "community_reactions", "tsv")
    tabs = cx.run(d, str(tmp_path / "adj"), metadata=meta, trait_cols=["site"], adjust="site", min_module_size=10,
                  fmt="tsv")
    tr = tabs["module_traits"]
    assert not (tr[(tr["test"] == "mann-whitney")]["p"] < 0.01).any()


def test_condition_specific_modules(tmp_path):
    """module A co-expressed only in 'case', B and C in both conditions, 20 reaction groups active only in case"""
    rng = np.random.default_rng(5)
    n = 40
    samples = ["S%02d" % i for i in range(2 * n)]
    group = np.array(["case"] * n + ["control"] * n)
    f = rng.normal(size=(3, 2 * n))
    rows = []
    for j in range(170):
        mod = j // 30 if j < 90 else None
        if mod == 0:                                   # A: shared factor in case only
            act = np.where(group == "case", 1.2 * f[0], 0) + rng.normal(scale=0.5, size=2 * n)
        else:
            act = (1.2 * f[mod] if mod is not None else 0) + rng.normal(scale=0.5, size=2 * n)
        mg = rng.lognormal(3, 0.3, size=2 * n)
        mt = mg * np.exp(act)
        if j >= 150:                                   # active only in case
            mt = np.where(group == "case", mt, 0.0)
        for i, s in enumerate(samples):
            rows.append({"sample_id": s, "mnx_id": "MNXR%03d" % j, "confidence": "specific", "mg_tpm": mg[i],
                         "mt_tpm": mt[i]})
    d = tmp_path / "map"
    write_table(pd.DataFrame(rows), str(d), "community_reactions", "tsv")
    meta = pd.DataFrame({"sample_id": samples, "group": group})
    tabs = cx.run_by(str(d), str(tmp_path / "by"), "group", metadata=meta, min_module_size=10, fmt="tsv",
                     permutations=2, threads=1)
    ct = tabs["class_tests"].set_index("class")
    assert {"observed", "null_mean", "p_greater"} <= set(ct.columns) and ct.loc["specific", "observed"] >= 1
    act = tabs["condition_active"]
    assert act[(act["condition"] == "control") & act["mnx_id"].isin(["MNXR160"])]["active"].sum() == 0
    mods = tabs["modules"]
    case = mods[mods["condition"] == "case"].set_index("mnx_id")["module"]
    ctrl = mods[mods["condition"] == "control"].set_index("mnx_id")["module"]
    a_case = case[["MNXR%03d" % j for j in range(30)]].value_counts().index[0]
    assert a_case != "M0"
    assert (ctrl.reindex(["MNXR%03d" % j for j in range(30)]) == "M0").mean() > 0.7   # no module A in controls
    m = tabs["module_matches"].set_index(["condition", "module"])["class"]
    assert m[("case", a_case)] == "specific"
    for k in (1, 2):
        mb = case[["MNXR%03d" % j for j in range(30 * k, 30 * k + 30)]].value_counts().index[0]
        assert m[("case", mb)] == "conserved"
