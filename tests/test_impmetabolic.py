import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import fixture  # noqa: E402
from impmetabolic import network, topology  # noqa: E402
from impmetabolic.cli import main as cli  # noqa: E402
from impmetabolic.store import read_manifest, read_table  # noqa: E402


def setup(tmp):
    data, refdb, root = str(tmp / "data"), str(tmp / "refdb"), str(tmp / "imp")
    fixture.reference(data)
    fixture.sample(root, "S1")
    fixture.sample(root, "S2")
    assert cli(["refdb", "build", "--data", data, "--out", refdb, "--format", "tsv"]) == 0
    return data, refdb, root


def test_refdb_build(tmp_path):
    _, refdb, _ = setup(tmp_path)
    r = read_table(refdb, "reactions")
    assert len(r) == 8 and set(r.columns) >= {"mnx_id", "equation", "ec_list", "is_transport"}
    assert set(read_table(refdb, "reaction_ec")["ec"]) == {"2.7.1.1", "5.3.1.9", "2.7.1.11", "4.1.2.13", "1.1.1.1"}
    assert set(read_table(refdb, "kegg_to_mnx_reaction")["source_id"]) == {"R00299", "R02740", "R04779", "R01068",
                                                                           "R00754"}
    assert {"ko_ec", "rn_ec", "ko_rn"} <= set(read_manifest(refdb)["tables"])
    cur = set(read_table(refdb, "currency")["name"])
    assert {"ATP", "H(+)", "H2O", "NAD(+)"} <= cur and "D-glucose" not in cur
    assert read_manifest(refdb)["metanetx_version"] == "4.4"
    from impmetabolic.refdb.build import metanetx_version
    assert metanetx_version(["MetaNetX/MNXref reconciliation", "MNXref Version 2025/03/03"]) == "2025/03/03"
    assert metanetx_version(["Licence: CC BY 4.0", "release: 4.5"]) == "4.5"
    sc = read_table(os.path.join(refdb, "raw"), "seed_compounds")          # parts concatenated
    assert list(sc["id"]) == ["cpd00001", "cpd00002"] and len(read_table(os.path.join(refdb, "raw"), "seed_reactions")) == 1
    from impmetabolic.refdb.download import plan
    urls = [u for u, _ in plan(str(tmp_path / "x"), ["seed"])]
    assert len(urls) == 111 and all("/refs/tags/v2.1.0/Biochemistry/" in u for u in urls)


def test_equation_parser():
    re_, pr = network.parse_equation("1 MNXM1@MNXD1 + 2 MNXM2@MNXD1 = n MNXM3@MNXD2")
    assert re_ == [(1.0, "MNXM1", "MNXD1"), (2.0, "MNXM2", "MNXD1")] and pr == [("n", "MNXM3", "MNXD2")]


def test_map_and_networks(tmp_path):
    _, refdb, root = setup(tmp_path)
    mp_ = str(tmp_path / "map")
    assert cli(["map", "--refdb", refdb, "--imp-root", root, "--out", mp_, "--threads", "1", "--format", "tsv"]) == 0
    gr = read_table(os.path.join(mp_, "samples", "S1"), "gene_reactions")
    ev = dict(zip(gr["gene"] + ">" + gr["mnx_id"], gr["evidence"]))
    assert ev["S1_G1>MNXR1"] == "ec:EC,rn:KO"            # direct EC, and its KO's KEGG reaction
    assert ev["S1_G2>MNXR2"] == "ec:KO,rn:KO"
    assert ev["S1_G4>MNXR4"] == "ec:RN,rn:RN"            # only the KEGG reaction
    assert ev["S1_G5>MNXR5"] == "ec:KO,rn:KO"
    assert not gr["gene"].isin(["S1_G6", "S1_G7"]).any()   # partial EC, no annotation
    assert set(gr.loc[gr["gene"] == "S1_G1", "bin_id"]) == {"S1:bin.1"}
    assert set(gr.loc[gr["mnx_id"] == "MNXR6", "group_id"]) == {"MNXR2"}            # same conversion
    com = read_table(mp_, "community_reactions")
    c1 = com[com["sample_id"] == "S1"].set_index("mnx_id")
    assert set(c1.index) == {"MNXR1", "MNXR2", "MNXR3", "MNXR4", "MNXR5", "MNXR7"}   # MNXR6 in the MNXR2 group
    assert c1.loc["MNXR2", "variants"] == "2" and (c1["confidence"] == "specific").all()
    r1 = com[(com["sample_id"] == "S1") & (com["mnx_id"] == "MNXR1")].iloc[0]
    assert float(r1["mg_reads"]) == 10 and float(r1["mt_reads"]) == 60
    br = read_table(mp_, "bin_reactions")
    assert set(br.loc[br["bin_id"] == "S1:bin.1", "mnx_id"]) == {"MNXR1", "MNXR2", "MNXR3"}
    assert set(br.loc[br["bin_id"] == "S1:bin.2", "mnx_id"]) == {"MNXR4", "MNXR5", "MNXR7"}
    # TPM and relative expression (metaT TPM / metaG TPM): metaG reads 10, 20, ... per gene (length 300), metaT
    # 3x with G1 doubled and G5 silent -> mta of MNXR1 (G1) = 2 x mta of MNXR2 (G2); MNXR5 not expressed
    for c in ("mg_tpm", "mt_tpm", "mta"):
        c1[c] = c1[c].astype(float)
    assert abs(c1["mg_tpm"].sum() - 1e6 * (1 + 2 + 3 + 4 + 5 + 5) / 28) < 1    # G1-G5 of 7 genes; G5 in 2 groups
    assert abs(c1.loc["MNXR1", "mta"] / c1.loc["MNXR2", "mta"] - 2) < 1e-6
    assert c1.loc["MNXR5", "mt_tpm"] == 0 and c1.loc["MNXR5", "mta"] == 0
    st = read_table(mp_, "stats")
    assert int(st.loc[st["sample_id"] == "S1", "genes_with_reaction"].iloc[0]) == 5
    assert int(st.loc[st["sample_id"] == "S1", "reactions_ec_only"].iloc[0]) == 2          # MNXR6, MNXR7
    # partial ECs on request: G6 then maps nowhere still (no 1.1.1.- reaction), but nothing breaks
    assert cli(["map", "--refdb", refdb, "--imp-root", root, "--out", str(tmp_path / "map2"), "--threads", "1",
                "--format", "tsv", "--methods", "ec", "--partial-ec"]) == 0
    # ECs with more than one reaction count as broad: MNXR7 (only via EC 1.1.1.1) broad, MNXR5 (KEGG) specific
    m3 = str(tmp_path / "map3")
    assert cli(["map", "--refdb", refdb, "--imp-root", root, "--out", m3, "--threads", "1", "--format", "tsv",
                "--max-ec-reactions", "1"]) == 0
    c3 = read_table(m3, "community_reactions")
    c3 = c3[c3["sample_id"] == "S1"].set_index("mnx_id")["confidence"]
    assert c3["MNXR7"] == "broad" and c3["MNXR5"] == "specific" and c3["MNXR2"] == "specific"
    n3 = str(tmp_path / "nets3")
    assert cli(["networks", "--mapping", m3, "--refdb", refdb, "--out", n3, "--format", "tsv", "--threads", "1"]) == 0
    t3 = read_table(n3, "topology")
    t3 = dict(zip(t3["set"] + ">" + t3["level"] + ">" + t3["id"], t3["reactions"]))
    assert t3["specific>community>S1"] == "5" and t3["all>community>S1"] == "6"

    nets = str(tmp_path / "nets")
    assert cli(["networks", "--mapping", mp_, "--refdb", refdb, "--out", nets, "--format", "tsv",
                "--graphml", "--threads", "2"]) == 0
    t0 = read_table(nets, "topology")
    t = t0[t0["set"] == "all"]
    lv = dict(zip(t["level"] + ">" + t["id"], t["reactions"]))
    assert lv["community>S1"] == "6" and lv["bin>S1:bin.1"] == "3"
    assert "bin>S1:bin.2" in lv                            # medium-quality MAG
    assert lv["pan_species>Bacteroides uniformis"] == "3"   # 2 MAGs (S1, S2); bin.2 has no species name
    assert lv["pan_genus>Blautia_A"] == "3"                  # ... but a genus: R4, R5, R7
    assert lv["pan_family>Bacteroidaceae"] == "3"
    pan = read_table(nets, "pan_reactome")
    pan = pan[pan["set"] == "all"]
    assert (pan["class"] == "core").all() and set(pan["n_mags"]) == {"2"}
    sc = read_table(nets, "scope")
    s1 = set(sc.loc[(sc["set"] == "all") & (sc["level"] == "community") & (sc["id"] == "S1"), "metabolite"])
    assert {"MNXM100", "MNXM101", "MNXM102", "MNXM103", "MNXM104"} <= s1   # glucose -> triose phosphates
    assert "MNXM106" not in s1                              # no ethanol in the seeds
    ta = t0[t0["set"] == "active"]
    assert dict(zip(ta["level"] + ">" + ta["id"], ta["reactions"]))["community>S1"] == "4"   # without MNXR5, MNXR7
    ks = read_table(nets, "keystones")
    assert set(ks["set"]) == {"specific", "all", "active"} and {"mnx_id", "load", "mta", "keystone"} <= set(ks.columns)
    hb = read_table(nets, "hubs")
    assert "name" in hb.columns and set(hb.loc[hb["metabolite"] == "MNXM100@MNXD1", "name"]) <= {"D-glucose 6-phosphate"}
    cov = read_table(nets, "mag_coverage")
    assert float(cov.loc[(cov["sample_id"] == "S1") & (cov["set"] == "all"), "share_in_mags"].iloc[0]) == 1.0
    assert os.path.exists(os.path.join(nets, "graphml", "all", "community", "S1.graphml"))


def test_topology_and_currency():
    rx = pd.DataFrame({"mnx_id": [r[0] for r in fixture.RXNS], "equation": [r[1] for r in fixture.RXNS]})
    net = network.build(rx["mnx_id"], rx)
    G = net.metabolite_graph(currency={"MNXM3", "MNXM7", "MNXM1", "MNXM8", "MNXM10"})
    assert "MNXM3@MNXD1" not in G and G.has_edge("MNXM99@MNXD1", "MNXM100@MNXD1")
    s = topology.summary(net, {"MNXM3", "MNXM7", "MNXM1", "MNXM8", "MNXM10"})
    assert s["reactions"] == 5 and s["components"] == 2      # glycolysis fragment, alcohol dehydrogenase
    have, fired = topology.scope(net, {"MNXM99", "MNXM3"})
    assert fired == {"MNXR1", "MNXR2", "MNXR3", "MNXR4"}
    assert topology.jaccard({"a", "b"}, {"b", "c"}) == 1 / 3
