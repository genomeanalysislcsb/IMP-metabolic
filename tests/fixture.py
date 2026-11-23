"""Synthetic reference files (in the formats of MetaNetX and KEGG REST) and IMP3 sample folders for the tests.

Reactions: a glycolysis fragment R1-R4 and an alcohol dehydrogenase R5. Two samples; bin.1 (Bacteroides uniformis,
in both samples) carries R1-R3, bin.2 R4 (only through its KEGG reaction) and R5; an unbinned gene has only a
partial EC (1.1.1.-) and maps to nothing."""
import os

MNX_HEAD = "### MetaNetX/MNXref reconciliation ###\n#Licence: CC BY 4.0\n#MNXref Version: 4.4\n"
METS = [("MNXM1", "H(+)"), ("WATER", "H2O"), ("MNXM3", "ATP"), ("MNXM7", "ADP"), ("MNXM8", "NAD(+)"),
        ("MNXM10", "NADH"), ("MNXM99", "D-glucose"), ("MNXM100", "D-glucose 6-phosphate"),
        ("MNXM101", "D-fructose 6-phosphate"), ("MNXM102", "D-fructose 1,6-bisphosphate"),
        ("MNXM103", "dihydroxyacetone phosphate"), ("MNXM104", "D-glyceraldehyde 3-phosphate"),
        ("MNXM105", "ethanol"), ("MNXM106", "acetaldehyde")]
RXNS = [("MNXR1", "1 MNXM99@MNXD1 + 1 MNXM3@MNXD1 = 1 MNXM100@MNXD1 + 1 MNXM7@MNXD1 + 1 MNXM1@MNXD1", "2.7.1.1",
         "R00299", "K00844"),
        ("MNXR2", "1 MNXM100@MNXD1 = 1 MNXM101@MNXD1", "5.3.1.9", "R02740", "K01810"),
        ("MNXR3", "1 MNXM101@MNXD1 + 1 MNXM3@MNXD1 = 1 MNXM102@MNXD1 + 1 MNXM7@MNXD1", "2.7.1.11", "R04779",
         "K00850"),
        ("MNXR4", "1 MNXM102@MNXD1 = 1 MNXM103@MNXD1 + 1 MNXM104@MNXD1", "4.1.2.13", "R01068", "K01623"),
        ("MNXR5", "1 MNXM105@MNXD1 + 1 MNXM8@MNXD1 = 1 MNXM106@MNXD1 + 1 MNXM10@MNXD1 + 1 MNXM1@MNXD1",
         "1.1.1.1", "R00754", "K00001")]


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def reference(data):
    m = os.path.join(data, "metanetx", "latest")
    w(os.path.join(m, "reac_prop.tsv"), MNX_HEAD + "#ID\tmnx_equation\treference\tclassifs\tis_balanced\tis_transport\n" +
      "".join("%s\t%s\tkeggR:%s\t%s\tB\t\n" % (r, eq, rn, ec) for r, eq, ec, rn, ko in RXNS) +
      "MNXR9\t1 MNXM99@MNXD2 = 1 MNXM99@MNXD1\tmnx:MNXR9\t\tB\tT\n")
    w(os.path.join(m, "reac_xref.tsv"), MNX_HEAD + "#source\tID\tdescription\n" +
      "".join("keggR:%s\t%s\tx\nkegg.reaction:%s\t%s\tx\n" % (rn, r, rn, r) for r, eq, ec, rn, ko in RXNS))
    w(os.path.join(m, "chem_prop.tsv"), MNX_HEAD + "#ID\tname\treference\tformula\tcharge\tmass\tInChI\tInChIKey\tSMILES\n"
      + "".join("%s\t%s\tmnx:%s\tX\t0\t1\t\t\t\n" % (i, n, i) for i, n in METS))
    w(os.path.join(m, "chem_xref.tsv"), MNX_HEAD + "#source\tID\tdescription\n" +
      "".join("mnx:%s\t%s\t%s\n" % (i, i, n) for i, n in METS))
    k = os.path.join(data, "kegg")
    w(os.path.join(k, "link_enzyme", "orthology"), "".join("ko:%s\tec:%s\n" % (ko, ec) for _, _, ec, _, ko in RXNS))
    w(os.path.join(k, "link_enzyme", "reaction"), "".join("rn:%s\tec:%s\n" % (rn, ec) for _, _, ec, rn, _ in RXNS))
    w(os.path.join(k, "link_reaction", "orthology"), "".join("ko:%s\trn:%s\n" % (ko, rn) for _, _, _, rn, ko in RXNS))


GENES = [  # gene, contig, attributes
    ("G1", "c1", "enzyme_ec=2.7.1.1;kegg_ko=K00844"),
    ("G2", "c1", "kegg_ko=K01810"),
    ("G3", "c2", "kegg_ko=K00850"),
    ("G4", "c3", "kegg_reaction=R01068"),
    ("G5", "c3", "kegg_ko=K00001,K99999"),
    ("G6", "c4", "enzyme_ec=1.1.1.-"),
    ("G7", "c4", "product=hypothetical protein"),
]


def sample(root, sid, with_bin2=True):
    b = os.path.join(root, sid)
    gff = "##gff-version 3\n" + "".join("%s\tBakta\tCDS\t%d\t%d\t.\t+\t0\tID=%s_%s;Name=x;%s\n" % (
        c, 10 * i + 1, 10 * i + 300, sid, g, a) for i, (g, c, a) in enumerate(GENES))
    w(os.path.join(b, "Analysis", "annotation", "annotation_CDS_RNA_hmms.gff"), gff)
    w(os.path.join(b, "Binning", "MetaWrap", "bins.contigs"),
      "c1\tbin.1\nc2\tbin.1\n" + ("c3\tbin.2\n" if with_bin2 else ""))
    w(os.path.join(b, "Stats", "all_bin_stats.tsv"),
      "Bin \tCompleteness \tContamination \tclassification\n"
      "bin.1 \t95.0 \t1.0 \td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;"
      "g__Bacteroides;s__Bacteroides uniformis\n" +
      ("bin.2 \t60.0 \t2.0 \td__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;"
       "g__Blautia_A;s__\n" if with_bin2 else ""))
    for o, f in (("mg", 1), ("mt", 3)):
        w(os.path.join(b, "Stats", o, "gene_read_count.txt"), "".join(
            "%s\t%d\t%d\t%s_%s\t%d\n" % (c, 10 * i, 10 * i + 300, sid, g, (i + 1) * 10 * f)
            for i, (g, c, a) in enumerate(GENES)))
    w(os.path.join(b, "Stats", "mg", "mg.read_stats.tsv"), "x\n")
    return b
