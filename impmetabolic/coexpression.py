"""Co-expression networks (weighted gene co-expression network analysis, WGCNA) of reaction groups across a cohort,
after Villette et al. 2025 (npj Biofilms Microbiomes 11:142; PyWGCNA) with the MT/MG ratio of Roume et al. 2015.

Value per sample and reaction group (--value):
  auto  relative expression (metaT TPM / metaG TPM) when the cohort has metaG and metaT, metaT TPM with metaT only,
        metaG TPM with metaG only (co-abundance instead of co-expression)
  mta, mt, mg  one of them explicitly

Steps: reaction groups present in >= --min-prevalence of the samples (present = metaG TPM > 0 for mta / mg, metaT
TPM > 0 for mt); values of absent samples are missing, not 0; reaction groups with identical profiles (the same
genes, e.g. substrate variants of one EC) are collapsed into one feature (listed in same_profile); ratio with a pseudocount (--pseudocount, default half
the smallest positive TPM); Yeo-Johnson power transform per reaction group (as sklearn's PowerTransformer in the
paper; --transform log2 / none); correlations from the pairs of samples where both are present; signed adjacency
((1 + r) / 2) ** power with the power chosen by the scale-free topology fit (lowest power with R2 >= --r2, or
--power); topological overlap (TOM); average-linkage clustering of 1 - TOM; modules by an adaptive
branch cut after the Dynamic Tree method (dynamic_cut, deepSplit 2 as in cutreeHybrid, without its PAM stage;
--tree-cut static: the static cut with most modules of at least --min-module-size); modules merged when their eigengenes correlate above 1 - --merge-diss (0.18 as in the
paper); eigengenes = first principal component; module membership kME, intramodular connectivity kWithin, hubs
(top 5 % kTotal) and intramodular hubs (top 10 % kWithin per module); module diversity (Shannon index of the summed
values of the module's members, as in the paper); module-trait association with sample metadata (categorical: point-
biserial correlation per level and Kruskal-Wallis; numeric: Spearman), BH-FDR over all modules and traits.
"""
import logging
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from .store import read_table, write_manifest, write_table

log = logging.getLogger(__name__)


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


# ----------------------------------------------------------------------------- input matrix
def matrices(com, confidence=None):
    """community_reactions -> (mg, mt) TPM matrices samples x reaction groups (NaN where the omic is missing)"""
    d = com if confidence in (None, "all") else com[com["confidence"] == confidence]
    out = []
    for c in ("mg_tpm", "mt_tpm"):
        v = pd.to_numeric(d[c], errors="coerce") if c in d.columns else pd.Series(np.nan, index=d.index)
        out.append(d.assign(v=v).pivot_table(index="sample_id", columns="mnx_id", values="v", aggfunc="sum")
                   if v.notna().any() else None)
    return out


def values(mg, mt, value="auto", min_prev=0.5, pseudocount=None):
    """-> (samples x features matrix with NaN for absent, value used, pseudocount, prevalence)"""
    if value == "auto":
        value = "mta" if mg is not None and mt is not None else ("mt" if mt is not None else "mg")
    if value in ("mta", "mg") and mg is None or value in ("mta", "mt") and mt is None:
        raise ValueError("value %s needs %s" % (value, "metaG and metaT" if value == "mta" else value))
    base = mg if value in ("mta", "mg") else mt
    base = base.fillna(0)
    present = base > 0
    prev = present.mean()
    keep = prev[prev >= min_prev].index
    if value == "mta":
        g, t = mg.reindex(columns=keep).fillna(0), mt.reindex(index=mg.index, columns=keep).fillna(0)
        if pseudocount is None:
            pos = np.concatenate([g.values[g.values > 0], t.values[t.values > 0]])
            pseudocount = float(pos.min()) / 2 if len(pos) else 1e-3
        X = (t + pseudocount) / (g + pseudocount)
    else:
        X = base[keep].astype(float)
    X = X.where(present[keep])
    return X, value, pseudocount, prev[keep]


def collapse_duplicates(X):
    """reaction groups with identical profiles (the same genes: substrate variants of one EC, subunits) are one
    feature: they would otherwise form modules of their own and inflate the connectivity. -> (X with one
    representative per profile, Series representative -> members separated by spaces)"""
    key = X.round(6).fillna(-1.0).T.apply(lambda r: hash(tuple(r.values)), axis=1)
    groups = pd.Series(X.columns, index=X.columns).groupby(key.values)
    rep, members = [], {}
    for _, cols in groups:
        cols = sorted(cols)
        rep.append(cols[0])
        members[cols[0]] = " ".join(cols)
    rep = [c for c in X.columns if c in set(rep)]
    return X[rep], pd.Series(members)


def transform(X, how="yeo-johnson"):
    """per feature, on the present values; standardised (as PowerTransformer(standardize=True))"""
    if how == "none":
        return X
    out = X.copy()
    for c in X.columns:
        v = X[c].dropna()
        if len(v) < 3 or v.nunique() < 2:
            out[c] = np.nan
            continue
        if how == "log2":
            y = np.log2(v.values)
        else:
            y, _ = stats.yeojohnson(v.values)
        y = (y - y.mean()) / (y.std() or 1.0)
        out.loc[v.index, c] = y
    return out.dropna(axis=1, how="all")


def correlation(X, method="pearson", min_periods=5):
    C = X.corr(method=method, min_periods=min_periods).fillna(0).values.astype(np.float32)
    np.fill_diagonal(C, 1.0)
    return C


# ----------------------------------------------------------------------------- WGCNA
def scale_free_fit(C, powers=range(1, 31), bins=10):
    """-> DataFrame power, r2 (signed: -sign(slope) * R2), slope, mean_k, median_k"""
    rows = []
    for b in powers:
        A = ((1 + C) / 2) ** b
        np.fill_diagonal(A, 0)
        k = A.sum(axis=1)
        h, edges = np.histogram(k, bins=bins)
        mid = (edges[:-1] + edges[1:]) / 2
        ok = h > 0
        if ok.sum() < 3:
            rows.append((b, np.nan, np.nan, k.mean(), np.median(k)))
            continue
        x, y = np.log10(mid[ok]), np.log10(h[ok] / h.sum())
        sl, ic, r, _, _ = stats.linregress(x, y)
        rows.append((b, -np.sign(sl) * r * r, sl, k.mean(), np.median(k)))
    return pd.DataFrame(rows, columns=["power", "r2", "slope", "mean_k", "median_k"])


def pick_power(fit, r2=0.8):
    ok = fit[fit["r2"] >= r2]
    if len(ok):
        return int(ok["power"].iloc[0]), True
    return int(fit.loc[fit["r2"].idxmax(), "power"]) if fit["r2"].notna().any() else 12, False


def tom(A):
    """signed topological overlap of an adjacency matrix (diagonal set to 1)"""
    A = A.astype(np.float32, copy=True)
    np.fill_diagonal(A, 0)
    k = A.sum(axis=1)
    L = A @ A
    T = (L + A) / (np.minimum.outer(k, k) + 1 - A)
    np.fill_diagonal(T, 1)
    return T


MAX_CORE_SCATTER = (0.64, 0.73, 0.82, 0.91, 0.95)            # per deepSplit 0..4, as in WGCNA's cutreeHybrid


def dynamic_cut(Z, min_size=20, deep_split=2, cut_height=0.99):
    """adaptive branch cut of a hierarchical tree, after the 'Dynamic Tree' idea of Langfelder et al. 2008
    (Bioinformatics 24:719): walking down from the root, a branch is split into its two sub-branches when both are
    large enough (>= min_size) and each is separated from the joining height by a gap of at least min_gap; otherwise
    the branch is one module (or unassigned if smaller than min_size). Heights are scaled to the tree's range;
    branches joining above cut_height (default 99 % of the range) are always split; deep_split 0..4 sets the core
    scatter and gap as WGCNA's cutreeHybrid (maxCoreScatter 0.64..0.95, minGap = (1 - maxCoreScatter) * 3/4).
    Unlike cutreeHybrid there is no second (PAM) stage: members of no branch stay unassigned (M0)."""
    n = Z.shape[0] + 1
    h = Z[:, 2]
    lo, hi = float(h.min()), float(h.max())
    rng = (hi - lo) or 1.0
    scaled = (h - lo) / rng
    max_scatter = MAX_CORE_SCATTER[max(0, min(4, int(deep_split)))]
    min_gap = (1 - max_scatter) * 3 / 4
    size = np.ones(2 * n - 1, dtype=int)
    for i in range(n - 1):
        size[n + i] = size[int(Z[i, 0])] + size[int(Z[i, 1])]
    # core height of a branch: mean of the merge heights inside it (leaves: 0)
    tot = np.zeros(2 * n - 1)
    cnt = np.zeros(2 * n - 1)
    for i in range(n - 1):
        a, b = int(Z[i, 0]), int(Z[i, 1])
        tot[n + i] = tot[a] + tot[b] + scaled[i]
        cnt[n + i] = cnt[a] + cnt[b] + 1
    core = np.divide(tot, cnt, out=np.zeros_like(tot), where=cnt > 0)

    def leaves(node):
        out, stack = [], [node]
        while stack:
            x = stack.pop()
            if x < n:
                out.append(x)
            else:
                stack += [int(Z[x - n, 0]), int(Z[x - n, 1])]
        return out
    labels = np.zeros(n, dtype=int)
    k = 0
    stack = [2 * n - 2]
    while stack:
        node = stack.pop()
        if node < n:
            continue
        i = node - n
        a, b = int(Z[i, 0]), int(Z[i, 1])
        if scaled[i] > cut_height:
            stack += [a, b]
            continue
        split = all(size[c] >= min_size and scaled[i] - core[c] >= min_gap for c in (a, b))
        if split:
            stack += [a, b]
        elif size[node] >= min_size and core[node] <= max_scatter:
            k += 1
            labels[leaves(node)] = k
        else:                               # a loose branch: look for modules inside it
            stack += [c for c in (a, b) if size[c] >= min_size]
    return labels


def cut_tree(Z, D, min_size=20, deep_split=2, method="dynamic"):
    """module labels (0 = unassigned): 'dynamic' = dynamic_cut; 'static' = the static cut of the tree with most
    clusters of >= min_size members"""
    if method == "dynamic":
        return dynamic_cut(Z, min_size, deep_split), "dynamic branch cut (deepSplit %d)" % deep_split
    best, best_n = None, -1
    for t in np.unique(np.round(Z[:, 2], 4))[::-1]:
        lab = fcluster(Z, t=t, criterion="distance")
        sizes = np.bincount(lab)
        n = int((sizes >= min_size).sum())
        if n > best_n:
            best, best_n = lab, n
    sizes = np.bincount(best)
    lab = np.where(sizes[best] >= min_size, best, 0)
    return lab, "static cut"


def eigengenes(Xt, labels):
    """first principal component per module (sign: positive correlation with the members' mean)"""
    Z = Xt.fillna(0)
    out = {}
    for m in sorted(set(labels) - {0}):
        cols = Z.columns[labels == m]
        Y = Z[cols].values
        Y = Y - Y.mean(axis=0)
        u, s, vt = np.linalg.svd(Y, full_matrices=False)
        pc = u[:, 0] * s[0]
        if np.corrcoef(pc, Y.mean(axis=1))[0, 1] < 0:
            pc = -pc
        out[m] = (pc - pc.mean()) / (pc.std() or 1.0)
    return pd.DataFrame(out, index=Xt.index)


def merge_modules(Xt, labels, diss=0.18):
    ME = eigengenes(Xt, labels)
    if ME.shape[1] < 2:
        return labels
    D = 1 - np.corrcoef(ME.values.T)
    cl = fcluster(linkage(squareform(D, checks=False), "average"), t=diss, criterion="distance")
    m2c = dict(zip(ME.columns, cl))
    return np.array([m2c.get(m, 0) if m else 0 for m in labels])


def relabel(labels):
    """M1 = largest module, ...; 0 -> M0 (unassigned)"""
    s = pd.Series(labels)
    order = s[s > 0].value_counts().index
    new = {m: i + 1 for i, m in enumerate(order)}
    return np.array([new.get(m, 0) for m in labels])


def shannon(v):
    v = np.asarray(v, float)
    v = v[v > 0]
    if not len(v):
        return np.nan
    p = v / v.sum()
    return float(-(p * np.log(p)).sum())


def traits(ME, meta, columns):
    rows = []
    for col in columns:
        if col not in meta.columns:
            continue
        x = meta[col].reindex(ME.index)
        num = pd.to_numeric(x, errors="coerce")
        if num.notna().sum() >= max(5, 0.8 * x.notna().sum()) and num.nunique() > 2:
            for m in ME.columns:
                ok = num.notna()
                r, p = stats.spearmanr(ME.loc[ok, m], num[ok])
                rows.append((m, col, "", "spearman", r, p))
            continue
        lv = [v for v in sorted(x.dropna().astype(str).unique()) if (x.astype(str) == v).sum() >= 2]
        if len(lv) < 2:
            continue
        for m in ME.columns:
            groups = [ME.loc[x.astype(str) == v, m].values for v in lv]
            p = stats.kruskal(*groups).pvalue if len(lv) > 2 else stats.mannwhitneyu(*groups).pvalue
            rows.append((m, col, "", "kruskal" if len(lv) > 2 else "mann-whitney", np.nan, p))
            for v in lv:
                ind = (x.astype(str) == v).astype(float)
                ok = x.notna()
                r, p = stats.pearsonr(ME.loc[ok, m], ind[ok])
                rows.append((m, col, v, "point-biserial", r, p))
    df = pd.DataFrame(rows, columns=["module", "trait", "level", "test", "r", "p"])
    if len(df):
        df["q"] = bh(df["p"].fillna(1).values)
    return df


def run(mapping, out, value="auto", confidence="all", min_prev=0.5, pseudocount=None, how="yeo-johnson",
        method="pearson", power=None, r2=0.8, min_module_size=20, merge_diss=0.18, metadata=None,
        trait_cols=(), fmt="auto", tree_cut="dynamic", deep_split=2):
    com = read_table(mapping, "community_reactions")
    mg, mt = matrices(com, confidence)
    X, value, pc, prev = values(mg, mt, value, min_prev, pseudocount)
    n0 = X.shape[1]
    X, members = collapse_duplicates(X)
    log.info("value %s%s; %d samples, %d reaction groups present in >= %.0f %% of them, %d distinct profiles", value,
             " (pseudocount %.3g TPM)" % pc if value == "mta" else "", X.shape[0], n0, 100 * min_prev, X.shape[1])
    Xt = transform(X, how)
    C = correlation(Xt, method, min_periods=max(5, int(0.5 * min_prev * len(Xt))))
    fit = scale_free_fit(C)
    if power is None:
        power, ok = pick_power(fit, r2)
        log.info("soft-threshold power %d (%s)", power, "scale-free R2 >= %.2f" % r2 if ok else
                 "no power reached R2 %.2f: highest R2" % r2)
    A = ((1 + C) / 2) ** power
    np.fill_diagonal(A, 0)
    T = tom(A)
    D = 1 - T
    np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), "average")
    labels, how_cut = cut_tree(Z, D, min_module_size, deep_split, tree_cut)
    labels = relabel(merge_modules(Xt, labels, merge_diss))
    ME = eigengenes(Xt, labels)
    log.info("%d modules (%s), %d of %d reaction groups unassigned", ME.shape[1], how_cut, int((labels == 0).sum()),
             len(labels))
    kTotal = A.sum(axis=1)
    kWithin = np.array([A[i, labels == labels[i]].sum() if labels[i] else np.nan for i in range(len(labels))])
    Zf = Xt.fillna(0)
    kME = np.array([np.corrcoef(Zf.iloc[:, i], ME[labels[i]])[0, 1] if labels[i] else np.nan
                    for i in range(len(labels))])
    mods = pd.DataFrame({"mnx_id": Xt.columns, "module": ["M%d" % m for m in labels], "kME": kME,
                         "kWithin": kWithin, "kTotal": kTotal, "prevalence": prev.reindex(Xt.columns).values,
                         "mean_value": X[Xt.columns].mean().values})
    mods["same_profile"] = mods["mnx_id"].map(members).fillna(mods["mnx_id"])
    mods["n_same_profile"] = mods["same_profile"].str.count(" ") + 1
    asg = mods["module"] != "M0"
    mods["hub"] = asg & (mods["kTotal"] >= mods.loc[asg, "kTotal"].quantile(0.95))
    mods["intramodular_hub"] = asg & (mods["kWithin"] >= mods.groupby("module")["kWithin"].transform(
        lambda s: s.quantile(0.9)))
    summ = []
    for m in sorted(set(labels) - {0}):
        idx = np.where(labels == m)[0]
        sub = A[np.ix_(idx, idx)]
        summ.append({"module": "M%d" % m, "size": len(idx), "density": float(sub.sum() / max(1, len(idx) *
                                                                                            (len(idx) - 1))),
                     "mean_kWithin": float(np.nanmean(kWithin[idx])), "sum_kWithin": float(np.nansum(kWithin[idx])),
                     "mean_kME": float(np.nanmean(kME[idx])),
                     "diversity": shannon(X.iloc[:, idx].sum(axis=0, skipna=True).values)})
    ME.columns = ["M%d" % m for m in ME.columns]
    tabs = {"modules": mods, "module_summary": pd.DataFrame(summ), "soft_threshold": fit,
            "eigengenes": ME.reset_index().melt(id_vars="sample_id", var_name="module", value_name="eigengene")}
    if metadata is not None and len(trait_cols):
        meta = metadata.set_index("sample_id")
        meta.index = meta.index.astype(str)
        tabs["module_traits"] = traits(ME, meta, trait_cols)
    rows = {}
    for k, v in tabs.items():
        write_table(v, out, k, fmt)
        rows[k] = int(len(v))
    write_manifest(out, "coexpression", rows, {"mapping": os.path.abspath(mapping), "value": value,
                                               "confidence": confidence, "min_prevalence": min_prev,
                                               "pseudocount": pc, "transform": how, "correlation": method,
                                               "power": power, "min_module_size": min_module_size,
                                               "merge_diss": merge_diss, "tree_cut": how_cut})
    return tabs
