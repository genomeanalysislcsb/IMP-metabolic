"""Currency metabolites by their role in each reaction, not by identity.

A cofactor is currency in a reaction that only uses it (ATP -> ADP in a kinase, NAD+ -> NADH in a dehydrogenase,
SAM -> SAH in a methyltransferase, acetyl-CoA -> CoA in an acetyltransferase), but not in its own synthesis,
degradation or interconversion (ADP + Pi -> ATP, AMP -> IMP -> ... -> uric acid, methionine -> SAM). Per reaction:

  1. cofactor pairs (one on each side: ATP / ADP, NAD+ / NADH, L-glutamate / 2-oxoglutarate, ...) and single
     carriers / small inorganic molecules (H2O, H+, phosphate, CO2, NH4+, CoA, ...) are taken out one after the
     other, most frequent first,
  2. but a pair only while both sides keep a metabolite that is not a single (water, phosphate, ...), and a single
     only while its side keeps another metabolite: what is left is the reaction's main conversion. If removing a
     pair would leave only water, phosphate etc., the reaction converts the cofactor itself and the pair stays
     (ADP + phosphate -> ATP + H2O keeps ADP -> ATP).

Pairs = a curated list (by MetaNetX name) plus the pairs that occur together on opposite sides of many reactions of
the reference (MetaNetX; --pair-min-reactions), i.e. the shuttles of group transfer reactions.
"""
from collections import Counter

import pandas as pd

from .network import PROTON, parse_equation

SINGLES = ("H2O", "H(+)", "phosphate", "hydrogenphosphate", "diphosphate", "CO2", "O2", "NH4(+)", "ammonia",
           "hydrogencarbonate", "hydrogen peroxide", "H2", "Na(+)", "K(+)", "Fe(2+)", "Fe(3+)", "Mg(2+)", "Zn(2+)",
           "Mn(2+)", "Cl(-)", "Ca(2+)", "Co(2+)", "Cu(2+)", "CoA", "coenzyme A", "holo-[acyl-carrier protein]",
           "holo-[ACP]", "(6S)-5,6,7,8-tetrahydrofolate", "tetrahydrofolate")
PAIRS = (("ATP", "ADP"), ("ATP", "AMP"), ("ADP", "AMP"), ("GTP", "GDP"), ("GTP", "GMP"), ("UTP", "UDP"),
         ("UTP", "UMP"), ("CTP", "CDP"), ("CTP", "CMP"), ("ITP", "IDP"), ("dATP", "dADP"), ("dTTP", "dTDP"),
         ("NAD(+)", "NADH"), ("NADP(+)", "NADPH"), ("FAD", "FADH2"), ("FMN", "FMNH2"),
         ("S-adenosyl-L-methionine", "S-adenosyl-L-homocysteine"), ("acetyl-CoA", "CoA"),
         ("L-glutamate", "2-oxoglutarate"), ("L-glutamine", "L-glutamate"), ("L-aspartate", "oxaloacetate"),
         ("pyruvate", "L-alanine"), ("ubiquinone-8", "ubiquinol-8"), ("menaquinone-8", "menaquinol-8"),
         ("demethylmenaquinone-8", "demethylmenaquinol-8"), ("oxidized ferredoxin", "reduced ferredoxin"),
         ("thioredoxin disulfide", "thioredoxin"), ("glutathione disulfide", "glutathione"),
         ("oxidized flavodoxin", "reduced flavodoxin"), ("UDP-alpha-D-glucose", "UDP"),
         ("5,10-methylenetetrahydrofolate", "(6S)-5,6,7,8-tetrahydrofolate"),
         ("10-formyltetrahydrofolate", "(6S)-5,6,7,8-tetrahydrofolate"))
IGNORE = (PROTON, "WATER")


def sides(eq):
    """metabolite ids of the two sides of a MetaNetX equation (compartments dropped)"""
    a, b = parse_equation(eq)
    return {m for _, m, _ in a}, {m for _, m, _ in b}


def frequent_pairs(reactions, min_reactions=100, ignore=IGNORE):
    """unordered metabolite pairs on opposite sides of >= min_reactions reactions -> Series (a, b) -> count"""
    c = Counter()
    for eq in reactions["equation"]:
        L, R = sides(eq)
        L, R = L - set(ignore), R - set(ignore)
        if len(L) * len(R) > 64:                      # very large generic equations
            continue
        for a in L:
            for b in R:
                if a != b:
                    c[tuple(sorted((a, b)))] += 1
    s = pd.Series(c, dtype=float)
    return s[s >= min_reactions].sort_values(ascending=False)


class Currency:
    """mode 'pairs': per-reaction removal as described above; mode 'global': the singles are removed from every
    reaction (older behaviour). ids are MetaNetX metabolite ids without compartment."""

    def __init__(self, singles=(), pairs=(), weights=None, mode="pairs"):
        self.singles, self.mode = set(singles), mode
        self.pairs = {frozenset(p) for p in pairs if len(set(p)) == 2}
        w = weights or {}
        self.weight = {p: w.get(p, 0) for p in self.pairs}
        self.by_met = {}
        for p in self.pairs:
            for m in p:
                self.by_met.setdefault(m, []).append(p)

    def __contains__(self, m):
        return m in self.singles

    def __iter__(self):
        return iter(self.singles)

    def __len__(self):
        return len(self.singles) + len(self.pairs)

    def strip(self, L, R):
        """(L, R) metabolite ids -> the main conversion (L', R')"""
        if self.mode == "global":
            return L - self.singles, R - self.singles
        L, R = set(L), set(R)
        found = {p for m in L for p in self.by_met.get(m, ()) if len(p & L) == 1 and len(p & R) == 1}
        for p in sorted(found, key=lambda p: -self.weight.get(p, 0)):
            if (L - p - self.singles) and (R - p - self.singles):      # something other than water etc. remains
                L, R = L - p, R - p
        for m in sorted(self.singles & (L | R)):
            if m in L and len(L) > 1:
                L = L - {m}
            elif m in R and len(R) > 1:
                R = R - {m}
        return L, R


def build(metabolites, reactions=None, mode="pairs", min_reactions=100, extra_singles=(), extra_pairs=(),
          keep=()):
    """-> (Currency, table of the singles and pairs used with names, counts and source)"""
    m = metabolites[["mnx_id", "name"]].fillna("")
    by_name = {}
    for i, n in zip(m["mnx_id"], m["name"]):
        by_name.setdefault(n.lower(), set()).add(i)

    def ids(x):
        return {x} if x.startswith("MNXM") or x == "WATER" else by_name.get(x.lower(), set())
    keep = {k.lower() for k in keep}
    singles = set()
    for n in tuple(SINGLES) + tuple(extra_singles):
        if n.lower() not in keep:
            singles |= ids(n)
    singles |= {"WATER"}
    rows = [("single", n, "", "curated", None) for n in SINGLES if ids(n) and n.lower() not in keep]
    pairs, weights = set(), {}
    for a, b in tuple(PAIRS) + tuple(extra_pairs):
        if a.lower() in keep or b.lower() in keep:
            continue
        for x in ids(a):
            for y in ids(b):
                pairs.add(tuple(sorted((x, y))))
        if ids(a) and ids(b):
            rows.append(("pair", a, b, "curated", None))
    if mode == "pairs" and reactions is not None and len(reactions):
        fp = frequent_pairs(reactions, min_reactions)
        name = dict(zip(m["mnx_id"], m["name"]))
        for (x, y), n in fp.items():
            if name.get(x, x).lower() in keep or name.get(y, y).lower() in keep:
                continue
            pairs.add((x, y))
            weights[frozenset((x, y))] = n
            rows.append(("pair", name.get(x, x), name.get(y, y), "frequent in the reference", int(n)))
    tab = pd.DataFrame(rows, columns=["kind", "a", "b", "source", "reactions"])
    return Currency(singles, pairs, weights, mode), tab
