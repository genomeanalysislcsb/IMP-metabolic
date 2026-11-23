# Changelog

## 0.1.0.dev0 (unreleased)

First version, built on the metabolic-network prototype (Artem Kokorin, LCSB).

- **Reference database** (`imp-metabolic refdb download | build`): the prototype's catalogue (MetaNetX, KEGG REST,
  ModelSEED, BiGG, MetaCyc via gapseq, Rhea, ECDomainMiner, TCDB, ChEBI); one Parquet table per source table
  instead of a SQLite file, plus derived tables (`reactions`, `reaction_ec`, `metabolites`, KEGG links, MetaNetX
  cross-reference bridges, `currency`) and a manifest with the MetaNetX version and the source files. MetaCyc .rda
  and ChEBI SDF loaders are optional (pyreadr, RDKit). ModelSEED is pinned to release v2.1.0, whose compound and
  reaction tables come in parts (`compound_00.tsv` …), concatenated by the build.
- **Mapping** (`imp-metabolic map`): genes of IMP3 samples -> MetaNetX reactions through the prototype's EC cascade
  (EC, else KO -> EC, else KEGG reaction -> EC; partial ECs left out) and through KEGG reactions (the gene's and its
  KOs') with MetaNetX's KEGG cross-references; evidence per gene and reaction; contig, MAG and metaG / metaT reads
  per gene; reactions of the community and of each MAG.
- **Networks** (`imp-metabolic networks`): bipartite metabolite-reaction networks (the prototype's
  `MetabolicNetwork`) of each community, each MAG and each GTDB species' pan-reactome (core / accessory / rare
  reactions); topology without currency metabolites (size, components, degree, clustering, path length, hub
  metabolites by betweenness), scope expansion from seed metabolites, share of the community's reactions in its
  MAGs; GraphML export.
