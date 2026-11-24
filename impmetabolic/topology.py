"""Topological analyses of metabolic networks: global properties, hub metabolites, scope (network expansion from
a set of seed metabolites), Jaccard similarity of reaction sets."""
import random

import networkx as nx
import pandas as pd


def summary(net, currency=(), path_sources=200, clustering=True):
    """global properties of a MetabolicNetwork and of its metabolite graph (without currency metabolites)"""
    G = net.metabolite_graph(currency)
    comps = sorted(nx.connected_components(G), key=len, reverse=True) if G.number_of_nodes() else []
    giant = G.subgraph(comps[0]) if comps else G
    out = {"reactions": len(net.reactions()), "metabolites": len(net.metabolites()),
           "graph_nodes": G.number_of_nodes(), "graph_edges": G.number_of_edges(),
           "components": len(comps), "giant_share": len(comps[0]) / G.number_of_nodes() if comps else 0.0,
           "mean_degree": 2 * G.number_of_edges() / G.number_of_nodes() if G.number_of_nodes() else 0.0,
           "clustering": nx.average_clustering(G) if clustering and G.number_of_nodes() else float("nan"),
           "mean_path_length": mean_path_length(giant, path_sources) if path_sources else float("nan")}
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
    from right to left); its products join the scope. -> (scope metabolites, fired reactions).
    Linear in the size of the network: each (reaction, direction) counts its substrates not yet reached."""
    from collections import defaultdict, deque

    def key(m):
        return net.nodes[m]["metabolite"] if merge_compartments else m
    sides, missing, waiting = [], [], defaultdict(list)
    for r in net.reactions():
        s = {key(m) for m in net.predecessors(r)}
        p = {key(m) for m in net.successors(r)}
        for a, b in ((s, p), (p, s)) if reversible else ((s, p),):
            i = len(sides)
            sides.append((r, b))
            missing.append(len(a))
            for m in a:
                waiting[m].append(i)
    have, fired = set(), set()
    queue = deque()

    def reach(m):
        if m not in have:
            have.add(m)
            queue.append(m)
    for m in seeds:
        reach(m)
    for i, n in enumerate(missing):                 # reactions without substrates
        if n == 0:
            fired.add(sides[i][0])
            for m in sides[i][1]:
                reach(m)
    while queue:
        m = queue.popleft()
        for i in waiting.get(m, ()):
            missing[i] -= 1
            if missing[i] == 0:
                fired.add(sides[i][0])
                for p in sides[i][1]:
                    reach(p)
    return have, fired


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else float("nan")


def keystones(net, currency=(), activity=None, k=300, seed=1, load_q=0.95, mta_q=0.75):
    """reactions with a high betweenness relative to their degree in the reaction graph ('load points', Rahman &
    Schomburg 2006) and a high relative expression (metaT / metaG), the keystone criterion of Roume et al. 2015.
    activity: Series mnx_id -> mta (optional). -> DataFrame mnx_id, degree, betweenness, load, mta, keystone"""
    G = net.reaction_graph(currency)
    G.remove_nodes_from([n for n, d in dict(G.degree()).items() if d == 0])
    if G.number_of_nodes() < 3:
        return pd.DataFrame(columns=["mnx_id", "degree", "betweenness", "load", "mta", "keystone"])
    bc = nx.betweenness_centrality(G, k=min(k, G.number_of_nodes()), seed=seed)
    df = pd.DataFrame({"mnx_id": list(bc), "betweenness": list(bc.values())})
    df["degree"] = df["mnx_id"].map(dict(G.degree()))
    df["load"] = df["betweenness"] / df["degree"]
    df["mta"] = df["mnx_id"].map(activity) if activity is not None else float("nan")
    high_load = df["load"] >= df["load"].quantile(load_q)
    if df["mta"].notna().any():
        df["keystone"] = high_load & (df["mta"] >= df["mta"].quantile(mta_q))
    else:
        df["keystone"] = high_load
    return df.sort_values("load", ascending=False)[["mnx_id", "degree", "betweenness", "load", "mta", "keystone"]]
