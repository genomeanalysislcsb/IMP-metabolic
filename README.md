# IMP-metabolic

Metabolic networks of microbial communities and their MAGs, reconstructed from the gene annotations of
[IMP3](https://github.com/LCSB-BioCore/IMP3) runs, for whole cohorts:

- a **community network** per sample (all genes of the sample's assembly),
- a **network per MAG** (the genes on the MAG's contigs),
- a **pan-reactome network per taxon** (all MAGs of one GTDB species across the cohort: core and accessory reactions),

with **topological analyses** first (size, connectivity, hub metabolites, scope from seed metabolites, comparison of
reaction sets) and **flux-based analyses** (gap-filling, FBA) as the next step. Results are tables that
[IMP-explorer](https://github.com/LCSB-BioCore/imp-explorer) shows next to the other layers of a study.

IMP-metabolic is built on the *metabolic-network* prototype by Artem Kokorin and Marie Mahieux (LCSB,
2025–2026): its reference-data catalogue, the EC cascade from Mantis annotations to MetaNetX reactions and the
bipartite `MetabolicNetwork`.

## Install

```bash
pip install .            # pandas, networkx, pyarrow
pip install .[fast]      # + DuckDB (fast queries over the reference tables)
```

## Usage

```bash
# 1. reference data (once; needs internet, e.g. on the HPC access node; ~2 GB)
imp-metabolic refdb download --data refdata/
imp-metabolic refdb build --data refdata/ --out refdb/

# 2. genes -> reactions for every sample of an IMP3 cohort (sub-directories of --imp-root, or --samples sheet)
imp-metabolic map --refdb refdb/ --imp-root /path/to/imp3/ --out mapping/ --threads 8

# 3. networks and their topology
imp-metabolic networks --mapping mapping/ --refdb refdb/ --out networks/ [--graphml]
```

## How it works

**Reference database.** The files of the catalogue (MetaNetX, KEGG REST, ModelSEED, BiGG, MetaCyc via gapseq,
Rhea, ECDomainMiner, TCDB, ChEBI) are loaded into one table each (`refdb/raw/`) and combined into the tables used
downstream (`refdb/`): `reactions` (MetaNetX equations, ECs, balance and transport flags), `reaction_ec`,
`metabolites`, KEGG links `ko_ec`, `rn_ec`, `ko_rn`, the MetaNetX cross-references `<db>_to_mnx_<type>` and
`currency` (cofactors and inorganic ions by name). The tables are Parquet files with a `manifest.json` (MetaNetX
version, source files and dates) instead of a SQLite database: columnar and compressed, read directly into pandas,
and queried in place with DuckDB when installed (the same engine as IMP-explorer); the build is a plain directory of
files that can be versioned and shared.

**Mapping.** Per gene (CDS of `Analysis/annotation/annotation_CDS_RNA_hmms.gff`, Mantis attributes):

| path | evidence | route |
|---|---|---|
| `ec` | `ec:EC`, `ec:KO`, `ec:RN` | the gene's ECs; without them the ECs of its KOs; without those the ECs of its KEGG reactions → MetaNetX reactions with that EC (partial ECs such as 1.3.99.- left out unless `--partial-ec`) |
| `rn` | `rn:RN`, `rn:KO` | the gene's KEGG reactions and those of its KOs → MetaNetX reactions through MetaNetX's KEGG cross-references (more specific: one EC often stands for many reactions) |

Contig → MAG from `Binning/MetaWrap/bins.contigs`, MAG quality and GTDB taxonomy from `Stats/all_bin_stats.tsv`,
reads per gene (metaG, metaT) from `Stats/<omic>/gene_read_count.txt`.

**Networks.** Bipartite graphs of reactions and metabolites (`MNXM…@compartment`), edges substrate → reaction →
product with the stoichiometry. MetaNetX equations are undirected, so reversibility is unknown until curation /
thermodynamics are added; the scope treats reactions as reversible by default. For topology, the metabolite graph
(substrate – product of the same reaction) leaves out the currency metabolites, which otherwise become hubs.

## Outputs

| directory / table | content |
|---|---|
| `mapping/samples/<id>/gene_reactions` | gene, contig, bin, MetaNetX reaction, evidence, via (EC / KEGG reaction), metaG and metaT reads |
| `mapping/community_reactions` | sample × reaction: genes, metaG and metaT reads |
| `mapping/bin_reactions` | sample × MAG × reaction: the same per MAG |
| `mapping/bins`, `mapping/stats` | MAG quality and taxonomy; genes mapped per evidence path |
| `networks/topology` | one row per network (community, MAG, pan-reactome): reactions, metabolites, components, share in the largest component, mean degree, clustering, mean path length, scope size |
| `networks/hubs` | metabolites with the highest betweenness per community / pan-reactome network |
| `networks/pan_reactome`, `pan_summary` | taxon × reaction: MAGs carrying it, share, core / accessory / rare; per taxon pan, core and accessory sizes |
| `networks/scope` | metabolites reachable from the seeds (default: minimal glucose medium, `impmetabolic/data/seeds_minimal.tsv`) |
| `networks/mag_coverage` | share of the community's reactions found in its MAGs |
| `networks/graphml/` | the networks as GraphML (Cytoscape), with `--graphml` |

## Roadmap

- Rhea reactions and ECDomainMiner (Pfam → EC) / TCDB (transporters) evidence for more coverage
- currency metabolites by graph structure (GraphCM) next to the name list
- reversibility and thermodynamics (ModelSEED ΔG), spontaneous and exchange reactions
- gap-filling and FBA (COBRApy; community models), SBML export
- IMP-explorer page: community and pan-reactome networks of a cohort

## Tests

```bash
python3 -m pytest tests     # synthetic MetaNetX / KEGG files and IMP3 samples (tests/fixture.py)
```

## Licence

MIT (see LICENSE).
