"""Pan-reactomes: the reactions of all MAGs of one taxon (GTDB species by default) across the samples of a cohort;
core reactions are found in most MAGs of the taxon, accessory ones in some."""
import pandas as pd


def pan_reactome(bin_reactions, bins, rank="species", tiers=("high", "medium"), min_mags=2, core=0.9, rare=0.15):
    """-> (taxon x reaction table: n_mags, share, class), per-taxon summary"""
    b = bins[bins["quality"].isin(tiers) & (bins[rank].fillna("") != "")][["bin_id", rank]]
    n = b.groupby(rank)["bin_id"].nunique()
    keep = n[n >= min_mags].index
    br = bin_reactions[["bin_id", "mnx_id"]].drop_duplicates().merge(b[b[rank].isin(keep)], on="bin_id")
    t = br.groupby([rank, "mnx_id"], as_index=False)["bin_id"].nunique().rename(columns={"bin_id": "n_mags"})
    t["share"] = t["n_mags"] / t[rank].map(n)
    t["class"] = pd.cut(t["share"], [-1, rare, core - 1e-9, 1.01], labels=["rare", "accessory", "core"]).astype(str)
    t = t.rename(columns={rank: "taxon"})
    s = t.groupby("taxon").agg(pan=("mnx_id", "nunique"), core=("class", lambda c: (c == "core").sum()),
                               accessory=("class", lambda c: (c == "accessory").sum()),
                               rare=("class", lambda c: (c == "rare").sum())).reset_index()
    s.insert(1, "n_mags", s["taxon"].map(n).astype(int))
    s.insert(1, "rank", rank)
    t.insert(0, "rank", rank)
    return t, s
