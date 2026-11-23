"""Topological analyses of metabolic networks: global properties, hub metabolites, scope (network expansion from
a set of seed metabolites), Jaccard similarity of reaction sets."""
import random

import networkx as nx
import pandas as pd


def summary(net, currency=()):
    """global properties of a MetabolicNetwork and of its metabolite graph (without currency metabolites)"""
    G = net.metabolite_graph(currency)
    comps = sorted(nx.connected_components(G), key=len, reverse=True) if G.number_of_nodes() else []
    giant = G.subgraph(comps[0]) if comps else G
    out = {"reactions": len(net.reactions()), "metabolites": len(net.metabolites()),
           "graph_nodes": G.number_of_nodes(), "graph_edges": G.number_of_edges(),
           "components": len(comps), "giant_share": len(comps[0]) / G.number_of_nodes() if comps else 0.0,
           "mean_degree": 2 * G.number_of_edges() / G.number_of_nodes() if G.number_of_nodes() else 0.0,
           "clustering": nx.average_clustering(G) if G.number_of_nodes() else 0.0,
           "mean_path_length": mean_path_length(giant)}
    return out


def mean_path_length(G, sources=200, seed=1):
    """mean shortest-path length of a connected graph, estimated from up to `sources` BFS roots"""
    n = G.number_of_nodes()
    if n < 2:
        return 0.0
    nodes = list(G.nodes)
    random.Random(seed).shuffle(nodes)
    tot = cnt = 0
    for s in nodes[:sources]:
        d = nx.single_source_shortest_path_length(G, s)
        tot += sum(d.values())
        cnt += len(d) - 1
    return tot / cnt if cnt else 0.0


def hubs(net, currency=(), top=20, k=200, seed=1):
    """metabolites with the highest betweenness (approximated from k sources) and their degree"""
    G = net.metabolite_graph(currency)
    if G.number_of_nodes() < 3:
        return pd.DataFrame(columns=["metabolite", "degree", "betweenness"])
    bc = nx.betweenness_centrality(G, k=min(k, G.number_of_nodes()), seed=seed)
    df = pd.DataFrame({"metabolite": list(bc), "betweenness": list(bc.values())})
    df["degree"] = df["metabolite"].map(dict(G.degree()))
    return df.sort_values("betweenness", ascending=False).head(top)[["metabolite", "degree", "betweenness"]]


def scope(net, seeds, reversible=True, merge_compartments=True):
    """network expansion: a reaction fires when all its substrates are in the scope (with reversible=True also
    from right to left); its products join the scope. -> (scope metabolites, fired reactions)"""
    def key(m):
        return net.nodes[m]["metabolite"] if merge_compartments else m
    rx = []
    for r in net.reactions():
        s = {key(m) for m in net.predecessors(r)}
        p = {key(m) for m in net.successors(r)}
        rx.append((r, s, p))
        if reversible:
            rx.append((r, p, s))
    have = set(seeds)
    fired = set()
    changed = True
    while changed:
        changed = False
        for r, s, p in rx:
            if s <= have and not p <= have:
                have |= p
                fired.add(r)
                changed = True
            elif s <= have:
                fired.add(r)
    return have, fired


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else float("nan")
