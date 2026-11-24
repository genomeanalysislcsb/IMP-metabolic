"""Bipartite metabolite-reaction network (ported from the metabolic-network prototype, A. Kokorin and M. Mahieux).

Nodes: reactions (MetaNetX MNXR ids) and metabolites (MET@COMPARTMENT, e.g. MNXM3@MNXD1); edges substrate ->
reaction -> product with the stoichiometric coefficient. MetaNetX equations are undirected: the direction of the
edges follows the written equation; reversible = None until thermodynamics / curation are added."""
import re

import networkx as nx
import pandas as pd

_TERM = re.compile(r"^\s*(\S+)\s+(\S+?)@(\S+)\s*$")


def parse_side(side):
    out = []
    for term in side.split(" + "):
        m = _TERM.match(term)
        if not m:
            continue
        try:
            coeff = float(m.group(1))
        except ValueError:
            coeff = m.group(1)                       # generic coefficients (n, 2n, ...) kept as text
        out.append((coeff, m.group(2), m.group(3)))
    return out


def parse_equation(eq):
    """'1 MNXM1@MNXD1 + 2 MNXM2@MNXD1 = 1 MNXM3@MNXD1' -> (reactants, products) as (coeff, met, comp)"""
    if "=" not in eq:
        return [], []
    left, right = eq.split("=", 1)
    left, right = left.replace("<", ""), right.replace(">", "")
    return parse_side(left), parse_side(right)


PROTON = "MNXM1"


def reaction_groups(reactions, ignore=(PROTON,)):
    """MetaNetX reactions that are the same conversion: the same metabolites on the two sides, ignoring
    compartments, coefficients, direction and protons (transports keep their compartments).
    -> Series mnx_id -> group id (the smallest MNXR id of the group)"""
    keys = {}
    for r, eq in zip(reactions["mnx_id"], reactions["equation"]):
        a, b = parse_equation(eq)
        sa = tuple(sorted({m for _, m, _ in a if m not in ignore}))
        sb = tuple(sorted({m for _, m, _ in b if m not in ignore}))
        if sa == sb:                                   # transport: keep the compartments
            sa = tuple(sorted("%s@%s" % (m, c) for _, m, c in a))
            sb = tuple(sorted("%s@%s" % (m, c) for _, m, c in b))
        keys[r] = (sa, sb) if sa <= sb else (sb, sa)
    k = pd.Series(keys)
    if k.empty:
        return pd.Series(dtype=str)
    order = pd.Series(k.index, index=k.index).str.extract(r"(\d+)$", expand=False).astype(float)
    df = pd.DataFrame({"key": k.astype(str), "order": order}).sort_values("order")
    first = df[~df["key"].duplicated()]
    rep = pd.Series(first.index, index=first["key"].values)
    return df["key"].map(rep)


class MetabolicNetwork(nx.DiGraph):
    def add_reaction(self, rxn_id, equation, **attrs):
        reactants, products = parse_equation(equation)
        self.add_node(rxn_id, bipartite="reaction", equation=equation, reversible=attrs.pop("reversible", None),
                      **attrs)
        for side, d in ((reactants, 1), (products, -1)):
            for coeff, met, comp in side:
                mid = "%s@%s" % (met, comp)
                if mid not in self:
                    self.add_node(mid, bipartite="metabolite", metabolite=met, compartment=comp)
                self.add_edge(*((mid, rxn_id) if d == 1 else (rxn_id, mid)), stoichiometry=coeff)

    def reactions(self):
        return [n for n, d in self.nodes(data=True) if d.get("bipartite") == "reaction"]

    def metabolites(self):
        return [n for n, d in self.nodes(data=True) if d.get("bipartite") == "metabolite"]

    def metabolite_graph(self, currency=(), merge_compartments=False, directed=False):
        """metabolite projection: A - B when A is a substrate and B a product of the same reaction; currency
        metabolites (ids without compartment, e.g. MNXM3) are left out so they do not form fake hubs"""
        cur = set(currency)
        G = nx.DiGraph() if directed else nx.Graph()

        def key(n):
            return self.nodes[n]["metabolite"] if merge_compartments else n
        for r in self.reactions():
            subs = [key(m) for m in self.predecessors(r) if self.nodes[m]["metabolite"] not in cur]
            prods = [key(m) for m in self.successors(r) if self.nodes[m]["metabolite"] not in cur]
            for a in subs:
                for b in prods:
                    if a != b:
                        if G.has_edge(a, b):
                            G[a][b]["reactions"].add(r)
                        else:
                            G.add_edge(a, b, reactions={r})
        return G

    def reaction_graph(self, currency=()):
        """reaction projection: r1 - r2 when they share a metabolite that is a product of one and a substrate of the
        other (or, as MetaNetX equations are undirected, on any side); currency metabolites left out"""
        cur = set(currency)
        G = nx.Graph()
        G.add_nodes_from(self.reactions())
        for m in self.metabolites():
            if self.nodes[m]["metabolite"] in cur:
                continue
            rs = list(set(self.predecessors(m)) | set(self.successors(m)))
            for i in range(len(rs)):
                for j in range(i + 1, len(rs)):
                    G.add_edge(rs[i], rs[j])
        return G

    def to_graphml(self, path):
        H = nx.DiGraph()
        for n, d in self.nodes(data=True):
            H.add_node(n, **{k: (v if isinstance(v, (str, int, float)) else str(v)) for k, v in d.items()
                             if v is not None})
        for u, v, d in self.edges(data=True):
            H.add_edge(u, v, stoichiometry=str(d.get("stoichiometry", "")))
        nx.write_graphml(H, path)


def build(reaction_ids, reactions, attrs=None):
    """network of the given MetaNetX reactions; reactions: refdb table (mnx_id, equation, ...); attrs: optional
    DataFrame indexed by mnx_id with node attributes (n_genes, mg_reads, ...)"""
    eq = reactions.set_index("mnx_id")["equation"]
    net = MetabolicNetwork()
    for r in pd.unique(pd.Series(list(reaction_ids))):
        if r in eq.index:
            a = attrs.loc[r].to_dict() if attrs is not None and r in attrs.index else {}
            net.add_reaction(r, eq[r], **a)
    return net
