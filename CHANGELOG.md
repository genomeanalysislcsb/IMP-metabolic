# Changelog

## 0.1.0.dev0 (unreleased)

First version, built on the metabolic-network prototype (Artem Kokorin and Marie Mahieux, LCSB).

- **Reference database** (`imp-metabolic refdb download | build`): the prototype's catalogue (MetaNetX, KEGG REST,
  ModelSEED, BiGG, MetaCyc via gapseq, Rhea, ECDomainMiner, TCDB, ChEBI); one Parquet table per source table
  instead of a SQLite file, plus derived tables (`reactions`, `reaction_ec`, `metabolites`, KEGG links, MetaNetX
  cross-reference bridges, `currency`) and a manifest with the MetaNetX version and the source files. MetaCyc .rda
  and ChEBI SDF loaders are optional (pyreadr, RDKit). ModelSEED is pinned to release v2.1.0, whose compound and
  reaction tables come in parts (`compound_00.tsv` …), concatenated by the build.
- **Mapping** (`imp-metabolic map`): genes of IMP3 samples -> MetaNetX reactions through the prototype's EC cascade
  (EC, else KO -> EC, else KEGG reaction -> EC; partial ECs left out) and through KEGG reactions (the gene's and its
  KOs') with MetaNetX's KEGG cross-references; evidence per gene and reaction; contig, MAG and metaG / metaT reads
  per gene; reactions of the community and of each MAG. Each gene–reaction link is *specific* (KEGG-reaction
  evidence, or an EC standing for at most `--max-ec-reactions` MetaNetX reactions, default 10) or *broad* (only
  through ECs with many substrate variants); MetaNetX reactions that are the same conversion (other compartment,
  protons, coefficients) form one reaction group in the community and MAG tables.
- **Networks** (`imp-metabolic networks`): bipartite metabolite-reaction networks (the prototype's
  `MetabolicNetwork`) of each community, each MAG and the pan-reactome of each GTDB species, genus and family
  (`--ranks`) (core / accessory / rare
  reactions); topology without currency metabolites (size, components, degree, clustering, path length, hub
  metabolites by betweenness), scope expansion from seed metabolites, share of the community's reactions in its
  MAGs; GraphML export. Analysed for the specific and for all reactions (`--sets`), in parallel (`--threads`).
- **Abundance and activity**: metaG and metaT TPM per gene and reaction group, relative expression (metaT TPM /
  metaG TPM, Roume et al. 2015; Villette et al. 2025); reaction set *active* (metaT TPM ≥ `--active-tpm`); keystone
  reactions per community (load points with high relative expression, Roume et al. 2015), table `keystones`.
- coexpression `--by <condition> [--block site]`: active reaction groups, network and modules per condition;
  module overlap between the conditions of a site (Jaccard, hypergeometric test, BH) and per module its best
  partner and class (conserved / partial / specific), with the share of its reactions active only in its condition.
- coexpression: `--adjust site`: each transformed reaction group centred within each site before the network, so
  site differences do not form modules (Expobiome: the site module M8).
- coexpression: kME-based (re)assignment after the tree cut as WGCNA's blockwiseModules (`--reassign-kme` 0.5,
  members below kME 0.3 dropped): the branch cut has no PAM stage and left 2,188 of 3,190 Expobiome profiles
  unassigned.
- coexpression: `--block` (e.g. site): module–trait tests also within each block level, for designs where groups
  and sites are confounded (Expobiome: PD only in Kassel, RA only in Berlin, HC in both); column `block` in
  `module_traits` (`all` = pooled).
- coexpression: own adaptive branch cut after WGCNA's Dynamic Tree method (`--deep-split`, `--tree-cut static`);
  the package dynamicTreeCut does not run with NumPy 2 (np.in1d, float ranks) and is no longer used.
- coexpression: reaction groups with identical profiles (same genes: substrate variants of one EC, subunits) are one
  feature (`same_profile`, `n_same_profile` in `modules`); networks: `keystones` keeps the keystones and the top
  10 % load per community (5.5 M rows for Expobiome before).
- **Currency metabolites by role** (`--currency-mode pairs`, default): cofactor pairs (curated + frequent in the
  reference reactions, `--pair-min-reactions`) and carriers / inorganic molecules are left out of a reaction only
  where it uses them; their own synthesis and degradation stay in the metabolite and reaction graphs and the scope;
  table `currency`; `--currency-pairs FILE` for further pairs; `--currency-mode global` = the earlier fixed list.
- **Currency metabolites** matched case-insensitively and extended (CoA, which MetaNetX 4.4 calls "CoA"; SAM / SAH;
  nucleoside mono-/diphosphates; quinones, ferredoxin, thioredoxin, glutathione; metal ions); `networks` applies the
  current list also to older refdb builds and takes `--currency FILE` / `--keep-currency NAMES`; the names used
  are logged and stored in the manifest.
- **Co-expression** (`imp-metabolic coexpression`): WGCNA of reaction groups across a cohort in Python, after
  Villette et al. 2025: metaT / metaG ratio by default, metaT or metaG alone when only one omic is available;
  prevalence filter with absent values kept missing, one pseudocount, Yeo-Johnson transform, signed adjacency with
  scale-free power, TOM, dynamic tree cut (optional dependency dynamicTreeCut), eigengene merging, kME / kWithin,
  hubs, module diversity, module–trait tests with BH-FDR over all modules and traits.
