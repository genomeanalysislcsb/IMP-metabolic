# -*- coding: utf-8 -*-
"""Catalogue of the reference data: where each file is downloaded from and how it is loaded into a table.

Ported from metabolic-network (Artem Kokorin, LCSB) data_catalogue.py; optional loaders (MetaCyc .rda via
pyreadr, ChEBI SDF via RDKit) import their library only when used.
"""
import json

import pandas as pd

from pathlib import Path
from urllib.parse import urljoin

from functools import partial
from collections.abc import Iterator, Callable, Sequence, Hashable

__all__ = ['DATA_CATALOGUE', 'OPTIONAL_LOADERS']

# tables whose loader needs an optional library; skipped (with a warning) when it is missing
OPTIONAL_LOADERS = {'load_metacyc_rda_file': 'pyreadr', 'load_chebi_file': 'rdkit'}


class File:
    __slots__ = ('name', 'key', 'loader')

    def __init__(
            self,
            name: str,
            key: str,
            loader: Callable[[Path], pd.DataFrame | dict[str, pd.DataFrame]] | None = None
    ):
        self.name = name
        self.key = key
        self.loader = loader


class FileGroup:  # one same remote path, multiple files
    __slots__ = ('local_subdir', 'remote_path', 'files')

    def __init__(
            self,
            local_subdir: str,  # output_dir/<DB_NAME>/local_subdir(s)
            remote_path: str,
            files: tuple[File, ...] = ()
    ):
        self.local_subdir = local_subdir
        self.remote_path = remote_path
        self.files = files


class SourceGroup:  # one same base URL
    __slots__ = ('database', 'base_url', 'file_groups')

    def __init__(
            self,
            database: str,
            base_url: str,
            file_groups: tuple[FileGroup, ...] = ()
    ):
        self.database = database
        self.base_url = base_url
        self.file_groups = file_groups

    def iter_downloads(self, output_dir: Path) -> Iterator[tuple[str, Path]]:
        """
        yields tuples (full_url, local_file_path)
        """
        output_dir_path = Path(output_dir).resolve()

        for group in self.file_groups:
            base = self.base_url.rstrip('/') + '/'
            remote_path = group.remote_path.strip('/')
            remote_parent = urljoin(base, remote_path + '/') if remote_path else base

            local_dir = output_dir_path / self.database
            if group.local_subdir:
                local_dir = local_dir / group.local_subdir

            for file in group.files:
                yield urljoin(remote_parent, file.name), local_dir / file.name


def make_csv_loader(
    sep: str | None = '\t',
    comment: str | None = None,
    dtype: pd._typing.DtypeArg | None = str,
    names: Sequence[Hashable] | None = None,
    fillna_value: str | None = ''  # filling NA values with empty strings
) -> Callable[[Path], pd.DataFrame]:
    # if no names are passed the behavior is identical to header=0
    # and column names are inferred from the first line of the file,
    # if column names are passed explicitly to names then the behavior
    # is identical to header=None (https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html):
    def _loader(file_path: Path) -> pd.DataFrame:
        df = pd.read_csv(
            file_path,
            sep=sep, comment=comment,
            dtype=dtype, names=names
        )
        if fillna_value is not None:
            df = df.fillna(fillna_value)
        return df

    return _loader


def load_metanetx_file(file_path: Path, names: tuple[str, ...]) -> pd.DataFrame:
    def _get_source_db_and_id(raw_source: str) -> tuple[str, str]:
        source_info, source_id = raw_source, raw_source  # MNXR01 -> (MNXR01, MNXR01)
        if ':' in raw_source:
            source_info, source_id = raw_source.split(':', 1)  # seed.reaction:1000 -> (seed.reaction, 1000)

        # "MetaNetX special identifiers": https://www.metanetx.org/mnxdoc/mnxref.html:
        special = ('biomass', 'water', 'empty', 'boundary', 'in', 'out', 'spontaneous')
        # IDs from common DBs
        # BiGG / SEED / KEGG / MetaCyc /
        # SABIO-RK (https://sabiork.h-its.org/) / Reactome /
        # Virtual Metabolic Human (https://www.vmh.life/) /
        # enviPath (https://envipath.org/) /
        # LIPID MAPS (https://www.lipidmaps.org/)
        common = ('bigg', 'seed', 'kegg', 'metacyc', 'sabiork', 'reactome', 'vmh', 'envipath', 'lipidmaps')
        # self-consistent source names (SLM is for Swiss Lipids; https://www.swisslipids.org/#/) and compartment names:
        # CL: https://obofoundry.org/ontology/cl.html; CCO: https://brg.ai.sri.com/CCO/; GO: https://geneontology.org/?
        other = ('chebi', 'hmdb', 'glycosphingo', 'slm', 'cco', 'go', 'cl')

        tmp = source_info.lower()

        if tmp in special:
            return 'special', source_id
        if tmp == 'mnx' or tmp.startswith(('mnxm', 'mnxr', 'mnxc', 'mnxd')):  # MNXD1/2 are generic compartments
            return 'mnx', source_id  # MetaNetX internal id
        if tmp == 'rh' or tmp.startswith('rhea'):
            return 'rhea', source_id  # Rhea id
        for prefix in common:
            if tmp.startswith(prefix):
                return prefix, source_id
        if tmp in other:
            return tmp, source_id

        return 'unknown', source_id

    def _type_from_mnx_id(mnx_id: str) -> str:
        m = mnx_id.lower()
        if m.startswith('mnxm'):
            return 'metabolite'
        if m.startswith('mnxr'):
            return 'reaction'
        if m.startswith(('mnxc', 'mnxd')):
            return 'compartment'  # MNXD1/2 are generic compartments
        return m

    df = make_csv_loader(
        sep='\t', comment='#',
        dtype=str, names=names,
        fillna_value='')(file_path)
    # taking only valid MNX compounds (M) / reactions (R) / compartments (C)...
    mask_mnx = df['mnx_id'].str.startswith('MNX')
    # ...or special identifiers ("MetaNetX special identifiers" on https://www.metanetx.org/mnxdoc/mnxref.html)
    mask_special = df['mnx_id'].isin(
        ['BIOMASS', 'WATER',  # for compounds
         'EMPTY',  # for reactions
         'BOUNDARY', 'IN', 'OUT',  # for compartments
         'SPONTANEOUS']  # for peptides/genes
    )
    df = df[mask_mnx | mask_special]
    # classifying entities into metabolites/reactions/compartments/...
    df['type'] = df['mnx_id'].map(_type_from_mnx_id)
    # adding records' source DB and corresponding ID
    df['source_db'], df['source_id'] = zip(*df['reference'].map(_get_source_db_and_id))

    return df


def load_kegg_list_file(file_path: Path) -> pd.DataFrame:
    return make_csv_loader(
        sep='\t', dtype=str,
        names=('kegg_id', 'description'),
        fillna_value=None
    )(file_path)


def load_kegg_link_file(file_path: Path) -> pd.DataFrame:
    df = make_csv_loader(
        sep='\t', dtype=str,
        names=('a', 'b'),
        fillna_value=None)(file_path)  # 'a' and 'b' are temporary names
    df.columns = [x.split(':')[0] for x in df.iloc[0]]  # getting actual col names from the 1st line prefixes
    df = df.map(lambda x: x.split(':', 1)[1])  # removing 'id:' for every cell
    return df


def load_kegg_pathway_link_file(file_path: Path) -> pd.DataFrame:
    tmp = load_kegg_link_file(file_path)
    # only retrieve pathway IDs starting with 'map' (the ones starting with 'ec'/'rn' are duplicates):
    df = tmp[tmp['path'].str.startswith('map')]
    return df


def load_json_file(file_path: Path) -> dict:
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def load_bigg_compartments_file(file_path: Path) -> pd.DataFrame:
    """
    main reference: http://bigg.ucsd.edu/data_access
    """
    data_comp = load_json_file(file_path)
    return pd.DataFrame(data_comp['compartments'])


def load_bigg_universal_file(file_path: Path) -> dict[str, pd.DataFrame]:
    def _convert_iterables(df):
        # small fix for processing iterables from universal data
        f = lambda x: json.dumps(x, ensure_ascii=False) if isinstance(x, (dict, list)) else x
        for col in df.columns:
            df[col] = df[col].map(f)
        return df

    data_univ = load_json_file(file_path)
    df_univ_met = _convert_iterables(pd.DataFrame(data_univ['metabolites']))
    df_univ_rxn = _convert_iterables(pd.DataFrame(data_univ['reactions']))

    return {
        'bigg_universal_metabolites': df_univ_met,  # metabolites from the universal model
        'bigg_universal_reactions': df_univ_rxn  # (NOTE: All reactions in the universal model are irreversible)
    }


def load_metacyc_rda_file(file_path: Path) -> pd.DataFrame:
    import pyreadr                      # optional dependency
    return pyreadr.read_r(str(file_path))[file_path.stem]  # DataFrame


def load_rhea_file(file_path: Path) -> pd.DataFrame:
    # TODO: Implement loading of .rdf/.owl files. Potentially, we could load Rhea data using web interface (via SPARQL).
    raise NotImplementedError


def load_ecdm_main_file(file_path: Path) -> pd.DataFrame:
    df_main = make_csv_loader()(file_path)
    # 1. filling NA values with empty strings (by default within the make_csv_loader func)
    # 2. sorting all rows by confidence scores in descending order
    # 3. removing duplicates from EC-Number column, keeping EC-Number <-> Pfam-Domain mapping with the highest score
    df_main = (
        df_main
        .sort_values('Confidence-Score', ascending=False)
        .drop_duplicates(subset=['EC-Number'], keep='first')
    )

    return df_main


def load_tcdb_substrates_file(file_path: Path) -> pd.DataFrame:
    return make_csv_loader(
        sep='\t', dtype=str,
        names=('tc_id', 'description'),
        fillna_value=''
    )(file_path)


def load_chebi_file(file_path: Path) -> pd.DataFrame:
    from rdkit.Chem import PandasTools  # optional dependency
    from rdkit.rdBase import DisableLog
    DisableLog('rdApp.*')  # removing RDKit-related logs
    return PandasTools.LoadSDF(
        str(file_path),  # first argument of LoadSDF must be of str type
        molColName=None, smilesName=None,  # no molecular preprocessing, raw data only
        removeHs=False, strictParsing=False
    )


# All specified files are downloaded but only those with
# available loader are used for the SQL DB construction:
DATA_CATALOGUE = (
    # MetaNetX:
    #
    # main reference: https://www.metanetx.org/mnxdoc/mnxref.html
    #
    #
    # ##############################################################################################################
    # PROP FILES:
    # ##############################################################################################################
    #
    # 1)
    # chem_prop.tsv file has multiple comment lines, then
    # #ID	name	reference	formula	charge	mass	InChI	InChIKey	SMILES
    # line, then actual tab-separated data columns.
    #
    # 1. ID ('mnx_id' will be used): The identifier of the chemical compound in the MNXref namespace [MNX_ID]
    # 2. name: The common name of the compound [STRING]
    # 3. reference: The reference compound, i.e. a compound selected in an external resource that best represent this entry [REFERENCE]
    # 4. formula: The formula of the compound [STRING]
    # 5. charge: The charge of the compound [INTEGER]
    # 6. mass: The mass of the compound [REAL]
    # 7. InChI: The standard InChI of the compound [STRING]
    # 8. InChIKey: The standard InChIKey of the compound [STRING]
    # 9. SMILES: The SMILES of the compound [STRING]
    #
    # 2)
    # reac_prop.tsv file has multiple comment lines, then
    # #ID	mnx_equation	reference	classifs	is_balanced	is_transport
    # line, then actual tab-separated data columns.
    #
    # 1. ID ('mnx_id' will be used): The identifier of the reaction in the MNXref namespace [MNX_ID]
    # 2. mnx_equation: Equation of the reaction in the MNXref namespace (compartmentalized and undirected) [EQUA]
    # 3. reference: The original best resource from where this reaction comes from [REFERENCE]
    # 4. classifs ('ec_list' will be used): The EC(s) associated to this reaction [STRING]
    # 5. is_balanced: Is the equation balanced with respect to elemental composition and charge [BOOLEAN]
    # 6. is_transport: Is this a transport reaction [BOOLEAN]
    #
    # 3)
    # comp_prop.tsv file has multiple comment lines, then
    # #ID	name	reference
    # line, then actual tab-separated data columns.
    #
    # 1. ID ('mnx_id' will be used): The identifier of the subcellular compartment in the namespace [MNX_ID]
    # 2. name: The common name of that compartment [STRING]
    # 3. reference: The original best resource from where this compartment comes from [REFERENCE]
    #
    #
    # ##############################################################################################################
    # XREF FILES:
    # ##############################################################################################################
    #
    # (chem / reac / comp)_xref.tsv files have multiple comment lines, then
    # #source	ID	description
    # line, then actual tab-separated data columns.
    #
    # 1. source ('reference' will be used): The identifier of a (chemical compound / biochemical reaction / subcellular compartment) in an external resource [XREF]
    # 2. ID ('mnx_id' will be used): The corresponding identifier in the MNXref namespace [MNX_ID]
    # 3. description: The (chemical compound / biochemical reaction / subcellular compartment) description given by the external resource [STRING]
    #
    #
    # Following source DBs are found in the files (df['reference'].apply(lambda x: x.split(':', 1)[0]).unique())
    # for metabolites:
    # ['BIOMASS' 'mnx' 'seed.compound' 'seedM' 'MNXM01' 'CHEBI' 'MNXM02'
    #  'bigg.metabolite' 'biggM' 'chebi' 'envipath' 'envipathM' 'hmdb'
    #  'kegg.compound' 'keggC' 'metacyc.compound' 'metacycM' 'reactome'
    #  'reactomeM' 'vmhM' 'vmhmetabolite' 'MNXM03' 'MNXM1' 'sabiork.compound'
    #  'sabiorkM' 'SLM' 'lipidmaps' 'lipidmapsM' 'slm' 'kegg.drug' 'keggD'
    #  'glycosphingo' 'kegg.glycan' 'keggG' 'rheaP' 'rheaG' 'WATER']
    #
    #  for reactions:
    #  ['EMPTY' 'bigg.reaction' 'biggR' 'metacyc.reaction' 'metacycR' 'mnx' 'rh'
    #  'rhea' 'seed.reaction' 'seedR' 'vmhR' 'vmhreaction' 'MNXR01' 'MNXR02'
    #  'MNXR03' 'sabiork.reaction' 'sabiorkR' 'kegg.reaction' 'keggR']
    #
    #  for compartments:
    #  ['mnx' 'cco' 'rh' 'rhea' 'GO' 'go' 'rheaC' 'seedC' 'CL' 'bigg.compartment'
    #  'biggC' 'cl']
    SourceGroup(
        database='metanetx',
        base_url='https://www.metanetx.org',
        file_groups=(
            FileGroup(
                local_subdir='latest',
                remote_path='ftp/latest/',
                files=(
                    File(name='chem_xref.tsv', key='metanetx_chem_xref',
                         loader=partial(load_metanetx_file, names=('reference', 'mnx_id', 'description'))),
                    File(name='chem_prop.tsv', key='metanetx_chem_prop',
                         loader=partial(load_metanetx_file, names=(
                             'mnx_id', 'name', 'reference', 'formula', 'charge', 'mass', 'inchi', 'inchikey', 'smiles'
                         ))),
                    File(name='chem_depr.tsv', key='metanetx_chem_depr', loader=None),
                    File(name='chem_isom.tsv', key='metanetx_chem_isom', loader=None),
                    File(name='comp_xref.tsv', key='metanetx_comp_xref',
                         loader=partial(load_metanetx_file, names=('reference', 'mnx_id', 'description'))),
                    File(name='comp_prop.tsv', key='metanetx_comp_prop',
                         loader=partial(load_metanetx_file, names=('mnx_id', 'name', 'reference'))),
                    File(name='comp_depr.tsv', key='metanetx_comp_depr', loader=None),
                    File(name='reac_xref.tsv', key='metanetx_reac_xref',
                         loader=partial(load_metanetx_file, names=('reference', 'mnx_id', 'description'))),
                    File(name='reac_prop.tsv', key='metanetx_reac_prop',
                         loader=partial(load_metanetx_file, names=(
                             'mnx_id', 'mnx_equation', 'reference', 'ec_list', 'is_balanced', 'is_transport'
                         ))),
                    File(name='reac_depr.tsv', key='metanetx_reac_depr', loader=None),
                    File(name='reac_isom.tsv', key='metanetx_reac_isom', loader=None)
                )  # TODO: add loaders for "isom" files, make corrections in load_metanetx_file func afterwards
            ),
        )
    ),

    # KEGG:
    #
    # main reference: www.kegg.jp/kegg/rest/keggapi.html
    #
    # KEGG list file has a following per-line structure: KEGG_ID\tDescription (no header).
    # For compound list:
    # C00001\tH2O; Water
    # C00002\tATP; Adenosine 5'-triphosphate
    #
    # KEGG link file (A-to-B) has a following per-line structure:
    # KEGG_ID_B\tKEGG_ID_A (no header), i.e., the order is reversed.
    # For enzyme-to-kegg_orthology link:
    # ko:K00001	ec:1.1.1.1
    # ko:K00002	ec:1.1.1.2
    SourceGroup(
        database='kegg',
        base_url='https://rest.kegg.jp',
        file_groups=(
            FileGroup(
                local_subdir='link_enzyme',
                remote_path='link/enzyme/',
                files=(
                    File(name='orthology', key='kegg_link_orthology_enzyme', loader=load_kegg_link_file),
                    File(name='pathway', key='kegg_link_pathway_enzyme', loader=load_kegg_pathway_link_file),
                    File(name='reaction', key='kegg_link_reaction_enzyme', loader=load_kegg_link_file)
                )
            ),  # ko/path/rn <-> ec mapping
            FileGroup(
                local_subdir='link_reaction',
                remote_path='link/reaction/',
                files=(
                    File(name='orthology', key='kegg_link_orthology_reaction', loader=load_kegg_link_file),
                    File(name='pathway', key='kegg_link_pathway_reaction', loader=load_kegg_pathway_link_file),
                    File(name='compound', key='kegg_link_compound_reaction', loader=load_kegg_link_file)
                )
            ),  # ko/path/cpd <-> rn mapping
            FileGroup(
                local_subdir='list',
                remote_path='list/',
                files=(
                    File(name='enzyme', key='kegg_list_enzyme', loader=load_kegg_list_file),
                    File(name='orthology', key='kegg_list_orthology', loader=load_kegg_list_file),
                    File(name='compound', key='kegg_list_compound', loader=load_kegg_list_file),
                    File(name='reaction', key='kegg_list_reaction', loader=load_kegg_list_file),
                    File(name='pathway', key='kegg_list_pathway', loader=load_kegg_list_file)
                )
            )  # information tables
        )
    ),

    # ModelSEED:
    #
    # main reference: https://github.com/ModelSEED/ModelSEEDDatabase/tree/master/Biochemistry
    #
    # SEED files for reaction and compound data are standard .tsv files.
    #
    # This function returns two pandas.DataFrame (for compounds and reactions) in a dict.
    # DataFrame's columns for compounds:
    #
    # [
    # 'id', 'abbreviation', 'name', 'formula', 'mass', 'source', 'inchikey',
    # 'charge', 'is_core', 'is_obsolete', 'linked_compound', 'is_cofactor',
    # 'deltag', 'deltagerr', 'pka', 'pkb', 'abstract_compound',
    # 'comprised_of', 'aliases', 'smiles', 'notes'
    # ]
    #
    # DataFrame's columns for reactions:
    #
    # [
    # 'id', 'abbreviation', 'name', 'code', 'stoichiometry', 'is_transport',
    # 'equation', 'definition', 'reversibility', 'direction',
    # 'abstract_reaction', 'pathways', 'aliases', 'ec_numbers', 'deltag',
    # 'deltagerr', 'compound_ids', 'status', 'is_obsolete', 'linked_reaction',
    # 'notes', 'source'
    # ]
    SourceGroup(
        database='seed',
        base_url='https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/refs/heads/master',
        file_groups=(
            FileGroup(
                local_subdir='',
                remote_path='Biochemistry/',
                files=(
                    File(name='compounds.tsv', key='seed_compounds', loader=make_csv_loader()),
                    File(name='reactions.tsv', key='seed_reactions', loader=make_csv_loader())
                )
            ),
        )
    ),

    # BiGG (http://bigg.ucsd.edu/data_access for reference):
    #
    # main reference: http://bigg.ucsd.edu/data_access
    SourceGroup(
        database='bigg',
        base_url='https://raw.githubusercontent.com/cdanielmachado/carveme/refs/heads/master/carveme/data',
        file_groups=(
            FileGroup(
                local_subdir='from_carveme/generated',
                remote_path='generated/',
                files=(
                    File(name='bigg_gprs.csv.gz', key='bigg_gprs', loader=None),
                    File(name='model_specific_data.csv.gz', key='model_specific_data', loader=None),
                    File(name='gene_annotations.tsv.gz', key='gene_annotations', loader=None)
                )
            ),
            FileGroup(
                local_subdir='from_carveme/input',
                remote_path='input/',
                files=(
                    File(name='bigg_models.tsv', key='bigg_models', loader=None),
                    File(name='biomass_db.tsv', key='biomass_db', loader=None),
                    File(name='manually_curated.tsv', key='manually_curated', loader=None),
                    File(name='media_db.tsv', key='media_db', loader=None),
                    File(name='mnx_compounds.tsv', key='mnx_compounds', loader=None)
                )
            )
        )  # this data from both generated/ and input/ can be potentially used for building a draft metabolic network
    ),

    SourceGroup(
        database='bigg',
        base_url='http://bigg.ucsd.edu',
        file_groups=(
            FileGroup(
                local_subdir='original_bigg',
                remote_path='api/v2/',
                files=(
                    File(name='compartments', key='bigg_compartments', loader=load_bigg_compartments_file),
                )
            ),  # compartments data
            FileGroup(
                local_subdir='original_bigg',
                remote_path='static/namespace/',
                files=(
                    File(name='bigg_models_metabolites.txt', key='bigg_metabolites', loader=make_csv_loader()),
                    File(name='bigg_models_reactions.txt', key='bigg_reactions', loader=make_csv_loader()),
                    File(name='universal_model.json', key='bigg_universal', loader=load_bigg_universal_file)
                )
            ),  # ! (NOTE: All reactions in the universal model are irreversible)
        )
    ),

    # MetaCyc data:
    #
    # DataFrames with the following columns are created for gapseq-related data:
    #
    # 'meta_genes.csv': ['rxn', 'genes', 'pwy', 'ncbi', 'uniprot', 'location', 'go', 'org']
    # 'meta_rea.tbl': ['id', 'name', 'altname', 'ec', 'kegg', 'pwy', 'enzymatic', 'gibbs',
    #    'direction']
    # 'meta_pwy.tbl': ['id', 'name', 'altname', 'hierarchy', 'taxrange', 'reaId', 'reaEc',
    #    'keyRea', 'reaName', 'reaNr', 'ecNr', 'superpathway', 'status',
    #    'spont']
    # 'meta_rea_pwy-gapseq.tbl': ['rxn', 'name', 'ec', 'pathway', 'pathway.name', 'subsystem', 'gapseq'] <- rxn-pwy mapping
    #
    #
    # DataFrames with the following columns are created for ftmsRanalysis-related data:
    #
    # 'mc_compounds.rda':  ['COMPOUND', 'URL', 'COMMON-NAME', 'MONOISOTOPIC-MW', 'MOLECULAR-WEIGHT',
    #    'MF', 'CHEMICAL-FORMULA', 'DBLINKS']
    # 'mc_reactions.rda':  ['REACTION', 'URL', 'SYSTEMATIC-NAME', 'EC-NUMBER', 'IN-PATHWAY',
    #    'REACTION-LIST', 'LEFT', 'RIGHT', 'REACTION-DIRECTION', 'COMPOUNDS',
    #    'DBLINKS']
    # 'mc_modules.rda':  ['MODULE', 'URL', 'COMMON-NAME', 'TYPES', 'SUPER-PATHWAYS',
    #    'SUB-PATHWAYS', 'REACTION-LIST', 'REACTION-LAYOUT', 'SPECIES',
    #    'TAXONOMIC-RANGE', 'PREDECESSORS', 'PATHWAY-LINKS', 'COMMENT',
    #    'DBLINKS', 'IS_SUPER_PATHWAY']
    # 'mc_reaction_module_node_map.rda':  ['REACTION', 'MODULE_NODE', 'MODULE']  <- this is reaction-to-pathway mapping
    SourceGroup(
        database='metacyc',
        base_url='https://raw.githubusercontent.com',
        file_groups=(
            FileGroup(
                local_subdir='from_gapseq',  # data from gapseq repo
                remote_path='jotech/gapseq/master/dat/',
                files=(
                    File(name='meta_genes.csv', key='metacyc_gapseq_genes', loader=make_csv_loader(sep=',')),
                    File(name='meta_pwy.tbl', key='metacyc_gapseq_pathways', loader=make_csv_loader()),
                    File(name='meta_rea.tbl', key='metacyc_gapseq_reactions', loader=make_csv_loader()),
                    File(name='meta_rea_pwy-gapseq.tbl', key='metacyc_gapseq_pwy_rxn_map', loader=make_csv_loader())
                )
            ),
            FileGroup(
                local_subdir='from_ftmsRanalysis',  # data used by ftmsRanalysis package
                remote_path='EMSL-Computing/MetaCycData/master/R/',
                files=(
                    File(name='mc_compounds.R', key='metacyc_ftms_metabolites_R', loader=None),
                    File(name='mc_reactions.R', key='metacyc_ftms_reactions_R', loader=None),
                    File(name='mc_modules.R', key='metacyc_ftms_pathways_R', loader=None)
                    # 'mc_compound_reaction_map.R', 'mc_module_reaction_map.R', 'mc_reaction_compound_map.R'
                )  # These files contain only format descriptions for the ones below
            ),
            FileGroup(
                local_subdir='from_ftmsRanalysis',  # data used by ftmsRanalysis package
                remote_path='EMSL-Computing/MetaCycData/master/data/',
                files=(
                    File(name='mc_compounds.rda', key='metacyc_ftms_metabolites', loader=load_metacyc_rda_file),
                    File(name='mc_reactions.rda', key='metacyc_ftms_reactions', loader=load_metacyc_rda_file),
                    File(name='mc_modules.rda', key='metacyc_ftms_pathways', loader=load_metacyc_rda_file),
                    File(name='mc_reaction_module_node_map.rda', key='metacyc_ftms_pwy_rxn_map', loader=load_metacyc_rda_file)
                    # 'mc_compounds_per_formula.rda', 'mc_compound_reaction_map.rda',
                    # 'mc_reaction_compound_map.rda', 'mc_module_reaction_map.rda',
                    # 'mc_module_superpathway_map.rda', 'mc_reaction_module_map.rda',
                )  # These files contain metacyc data stored in .rda files ('modules' == 'pathways')
            )
        )
    ),

    # Rhea (manually curated chemical and transport reactions of biological interest):
    #
    # main reference: https://www.rhea-db.org/help/download
    SourceGroup(
        database='rhea',
        base_url='https://ftp.expasy.org',
        file_groups=(
            FileGroup(
                local_subdir='',
                remote_path='databases/rhea/rdf/',
                files=(
                    File(name='rhea.rdf.gz', key='rhea_rdf_gz', loader=None),
                    File(name='chebi.owl.gz', key='chebi_owl_gz', loader=None),
                    File(name='rhea_rdf_documentation.pdf', key='rhea_docs', loader=None)
                )  # TODO: implement load_rhea_file function depending on the chosen strategy (rdf vs SPARQL)
            ),
        )
    ),

    # ECDomainMiner data: (1) associations of ECDomainMiner with confidence score above threshold
    # and (2) one-to-one associations between EC numbers and Pfam domains.
    #
    # main reference: https://ecdm.loria.fr
    #
    # ECDomainMiner data file contains EC-number <-> Pfam-domain mapping with confidence score above internal threshold:
    #
    # "EC-Number	Pfam-Domain	Confidence-Score	Type	In-InterPro"
    # "2.6.99.2	PF03740	1	GOLD	Yes"
    # "3.5.4.27	PF02289	1	GOLD	Yes"
    # "3.5.4.29	PF05165	1	GOLD	Yes"
    # "2.1.1.242	PF04445	1	GOLD	Yes"
    # "2.7.1.161	PF01982	1	GOLD	Yes"
    # "2.4.99.17	PF02547	0.9999643333	GOLD	Yes"
    #
    # This file is used to identify some proteins without the exact EC, but with Pfam ID reported in the annotation.
    #
    # ###############################################################################################################
    #
    # Additionally, we use file with one-to-one EC-number <-> Pfam-domain associations for more strict search:
    #
    # "EC number	Pfam domain	Score	EC name	Pfam name	Label	InterPro ?	P-value SIFTS	P-value Swissprot	P-value Trembl	P-value UniRule"
    # "2.7.7.86	PF03281	0.5827015	Cyclic GMP-AMP synthase.	Mab-21 protein	Gold	1.86219866885e-60	2.96446857355e-22	-	-	SD	NotInterpro"
    # "4.1.1.55	PF09084	0.580242	4,5-dihydroxyphthalate decarboxylase.	NMT1/THI5 like	Gold	-	6.81480172922e-13	1.62084356843e-123	-	SD	NotInterpro"
    SourceGroup(
        database='ecdm',
        base_url='https://ecdm.loria.fr',
        file_groups=(
            FileGroup(
                local_subdir='',
                remote_path='data/',
                files=(
                    File(name='EC-Pfam_calculated_associations.csv', key='ecdm_main', loader=load_ecdm_main_file),
                    File(name='OneToOne_ECPFamAssociations_march2016.csv', key='ecdm_1_to_1', loader=make_csv_loader())
                )
            ),
        )
    ),

    # TCDB:
    #
    # main reference: https://tcdb.org
    #
    # TCDB substrate file contains TC-ID of a transport protein and it's substrates' ChEBIs on each line:
    #
    # 3.A.1.208.25\tCHEBI:5437;glutathione|CHEBI:4050;cysteine
    # 1.C.57.1.3\tCHEBI:7990;peptide|CHEBI:25367;molecule
    # 9.B.470.1.2\tCHEBI:27333;xenobiotic
    # 1.A.27.1.3\tCHEBI:22563;anion
    SourceGroup(
        database='tcdb',
        base_url='https://tcdb.org',
        file_groups=(
            FileGroup(
                local_subdir='',
                remote_path='cgi-bin/substrates/',
                files=(
                    File(name='getSubstrates.py', key='tcdb_substrates', loader=load_tcdb_substrates_file),
                )  # TC <-> substrate mapping
            ),
            FileGroup(
                local_subdir='',
                remote_path='public/',
                files=(
                    File(name='tcdb', key='tcdb_proteins', loader=None),
                )  # All proteins in TCDB (FASTA format), not used further for now, more of a backup
            )
        )
    ),

    # ChEBI:
    #
    # ChEBI .sdf file is processed using RDKit (without any molecular preprocessing, raw data only)
    # and stored as pandas.DataFrame with the following columns:
    #
    # [
    # 'ChEBI ID', 'ChEBI NAME', 'STAR', 'MASS', 'FORMULA',
    # 'MONOISOTOPIC_MASS', 'SMILES', 'SYNONYM',
    # 'PubChem Substance Database Links', 'KEGG COMPOUND Database Links',
    # 'MetaboLights Database Links', 'ID', 'DEFINITION', 'INCHI', 'INCHIKEY',
    # 'IUPAC_NAME', 'BRENDA Database Links', 'SABIO-RK Database Links',
    # 'Reaxys Registry Numbers', 'PubMed Database Links',
    # 'Beilstein Registry Numbers', 'PubChem Compound Database Links',
    # 'CAS Registry Numbers', 'NMRShiftDB Database Links',
    # 'Gmelin Registry Numbers', 'BKMS-react Database Links',
    # 'BRENDA Ligand Database Links', 'KNApSAcK Database Links',
    # 'MetaCyc Database Links', 'CompTox Database Links',
    # 'LIPID MAPS Database Links', 'UniProt Database Links',
    # 'Rhea Database Links', 'SureChEMBL Database Links',
    # 'ChEMBL Database Links', 'BindingDB Database Links',
    # 'GeneOntology Database Links', 'PDB Database Links',
    # 'ChemSpider Database Links', 'LINCS Database Links',
    # 'HMDB Database Links', 'Wikipedia Database Links',
    # 'BPDB Database Links', 'SECONDARY_ID', 'DrugCentral Database Links',
    # 'KEGG GLYCAN Database Links', 'PDBeChem Database Links', 'CHARGE',
    # 'IEDB Database Links', 'Patent Database Links',
    # 'BioModels Database Links', 'Reactome Database Links', 'WURCS',
    # 'PPDB Database Links', 'FooDB Database Links',
    # 'SwissLipids Database Links', 'Agricola Database Links',
    # 'Virtual Metabolic Human Database Links', 'Golm Database Links',
    # 'Alan Wood's Pesticides Database Links', 'MassBank Database Links',
    # 'DrugBank Database Links', 'KEGG DRUG Database Links',
    # 'PubMed Central Database Links', 'Chinese Abstracts Database Links',
    # 'Expression Atlas Database Links', 'UM-BBD Database Links',
    # 'YMDB Database Links', 'GlyGen Database Links',
    # 'GlyTouCan Database Links', 'Carotenoids Database Database Links',
    # 'VSDB Database Links', 'MolBase Database Links',
    # 'The Signaling Pathways Project Database Links', 'ECMDB Database Links',
    # 'PPR Database Links', 'IntAct Database Links',
    # 'LIPID MAPS class Database Links', 'RESID Database Links',
    # 'WebElements Database Links', 'FAO/WHO standards Database Links',
    # 'ChemIDplus Database Links', 'CiteXplore Database Links',
    # 'SMID Database Links', 'CAS Database Links',
    # 'DrugBank Metabolite Database Links'
    # ]
    #
    # For the additional information, please check:
    # https://www.ebi.ac.uk/training/online/courses/chebi-the-online-chemical-dictionary-for-small-molecules/what-is-chebi/
    SourceGroup(
        database='chebi',
        base_url='https://ftp.ebi.ac.uk/pub/databases/chebi',
        file_groups=(
            FileGroup(
                local_subdir='',
                remote_path='SDF/',
                files=(
                    File(name='chebi.sdf.gz', key='chebi_data', loader=load_chebi_file),
                )
            ),
        )
    )
)
