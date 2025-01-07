#!/usr/bin/env python

__date__ = '2020-04-03'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import csv
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import yaml
import random
import warnings
import numpy as np
import pandas as pd
import scanpy as sc




# Set seed for reproducibility
seed_value = 0
y_n_print=False
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

# sc verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.settings.verbosity = 3
# sc.logging.print_versions()
# sc.settings.set_figure_params(dpi=80)

# mito_gene_list = sc.queries.mitochondrial_genes() # another query
# adata.var['gene_group__mito_transcript'] = [
#     x.startswith('MT-') for x in adata.var['gene_symbols']
# ]
# use this if var_names='gene_symbols' in sc.read_h5ad_mtx
# adata.var['mito_gene'] = [
#     x.startswith('MT-') for x in adata.var_names
# ]

# Label mitochondrial encoded proteins
#
# The below gene list was downloaded from on 3 Aug 2020:
# https://www.genenames.org/data/genegroup/#!/group/1974
gene_group__mito_protein = [
    'MT-ATP6',
    'MT-ATP8',
    'MT-CO1',
    'MT-CO2',
    'MT-CO3',
    'MT-CYB',
    'MT-ND1',
    'MT-ND2',
    'MT-ND3',
    'MT-ND4',
    'MT-ND4L',
    'MT-ND5',
    'MT-ND6'
]
# Label ribosomal protein genes
# Ribosomal protein: A ribosomal protein is any of
# the proteins that, in conjunction with rRNA, make up the ribosomal
# subunits involved in the cellular process of translation. A large
# part of the knowledge about these organic molecules has come from the
# study of E. coli ribosomes. Most ribosomal proteins have been
# isolated and specific antibodies have been produced. These, together
# with electronic microscopy and the use of certain reactives, have
# allowed for the determination of the topography of the proteins in
# ribosome. E.coli, other bacteria and Archaea have a 30S small subunit
# and a 50S large subunit, whereas humans and yeasts have a 40S small
# subunit and a 60S large subunit. Equivalent subunits are frequently
# numbered differently between bacteria, Archaea, yeasts and humans.
#
# The below gene list was downloaded from on 3 Aug 2020:
# https://www.genenames.org/data/genegroup/#!/group/1054
gene_group__ribo_protein = set([
    'DAP3',
    'FAU',
    'MRPL1',
    'MRPL1',
    'MRPL10',
    'MRPL10',
    'MRPL11',
    'MRPL11',
    'MRPL12',
    'MRPL12',
    'MRPL13',
    'MRPL13',
    'MRPL14',
    'MRPL14',
    'MRPL15',
    'MRPL15',
    'MRPL16',
    'MRPL16',
    'MRPL17',
    'MRPL17',
    'MRPL18',
    'MRPL18',
    'MRPL19',
    'MRPL19',
    'MRPL2',
    'MRPL2',
    'MRPL20',
    'MRPL20',
    'MRPL21',
    'MRPL21',
    'MRPL22',
    'MRPL22',
    'MRPL23',
    'MRPL23',
    'MRPL24',
    'MRPL24',
    'MRPL27',
    'MRPL27',
    'MRPL28',
    'MRPL28',
    'MRPL3',
    'MRPL3',
    'MRPL30',
    'MRPL30',
    'MRPL32',
    'MRPL32',
    'MRPL33',
    'MRPL33',
    'MRPL34',
    'MRPL34',
    'MRPL35',
    'MRPL35',
    'MRPL36',
    'MRPL36',
    'MRPL37',
    'MRPL37',
    'MRPL38',
    'MRPL38',
    'MRPL39',
    'MRPL39',
    'MRPL4',
    'MRPL4',
    'MRPL40',
    'MRPL40',
    'MRPL41',
    'MRPL41',
    'MRPL42',
    'MRPL42',
    'MRPL43',
    'MRPL43',
    'MRPL44',
    'MRPL44',
    'MRPL45',
    'MRPL45',
    'MRPL46',
    'MRPL46',
    'MRPL47',
    'MRPL47',
    'MRPL48',
    'MRPL48',
    'MRPL49',
    'MRPL49',
    'MRPL50',
    'MRPL50',
    'MRPL51',
    'MRPL51',
    'MRPL52',
    'MRPL52',
    'MRPL53',
    'MRPL53',
    'MRPL54',
    'MRPL54',
    'MRPL55',
    'MRPL55',
    'MRPL57',
    'MRPL57',
    'MRPL58',
    'MRPL9',
    'MRPS10',
    'MRPS10',
    'MRPS11',
    'MRPS11',
    'MRPS12',
    'MRPS12',
    'MRPS14',
    'MRPS14',
    'MRPS15',
    'MRPS15',
    'MRPS16',
    'MRPS16',
    'MRPS17',
    'MRPS17',
    'MRPS18A',
    'MRPS18A',
    'MRPS18B',
    'MRPS18B',
    'MRPS18C',
    'MRPS18C',
    'MRPS2',
    'MRPS2',
    'MRPS21',
    'MRPS21',
    'MRPS22',
    'MRPS22',
    'MRPS23',
    'MRPS23',
    'MRPS24',
    'MRPS24',
    'MRPS25',
    'MRPS25',
    'MRPS26',
    'MRPS26',
    'MRPS27',
    'MRPS27',
    'MRPS28',
    'MRPS28',
    'MRPS30',
    'MRPS30',
    'MRPS31',
    'MRPS31',
    'MRPS33',
    'MRPS33',
    'MRPS34',
    'MRPS34',
    'MRPS35',
    'MRPS35',
    'MRPS36',
    'MRPS36',
    'MRPS5',
    'MRPS6',
    'MRPS7',
    'MRPS9',
    'RPL10',
    'RPL10A',
    'RPL10L',
    'RPL11',
    'RPL12',
    'RPL13A',
    'RPL14',
    'RPL15',
    'RPL17',
    'RPL18A',
    'RPL19',
    'RPL21',
    'RPL22',
    'RPL23',
    'RPL23A',
    'RPL24',
    'RPL26',
    'RPL26L1',
    'RPL27',
    'RPL27A',
    'RPL28',
    'RPL29',
    'RPL3',
    'RPL30',
    'RPL31',
    'RPL32',
    'RPL34',
    'RPL35',
    'RPL35A',
    'RPL36',
    'RPL36A',
    'RPL36AL',
    'RPL37',
    'RPL37A',
    'RPL38',
    'RPL39',
    'RPL39L',
    'RPL3L',
    'RPL4',
    'RPL41',
    'RPL5',
    'RPL6',
    'RPL7',
    'RPL7A',
    'RPL7L1',
    'RPL8',
    'RPL9',
    'RPLP0',
    'RPLP1',
    'RPLP2',
    'RPS10',
    'RPS11',
    'RPS12',
    'RPS13',
    'RPS14',
    'RPS15',
    'RPS15A',
    'RPS16',
    'RPS17',
    'RPS18',
    'RPS19',
    'RPS2',
    'RPS20',
    'RPS21',
    'RPS23',
    'RPS24',
    'RPS25',
    'RPS26',
    'RPS27',
    'RPS27A',
    'RPS27L',
    'RPS28',
    'RPS29',
    'RPS3',
    'RPS3A',
    'RPS4X',
    'RPS4Y1',
    'RPS4Y2',
    'RPS5',
    'RPS6',
    'RPS7',
    'RPS8',
    'RPS9',
    'UBA52'
])
# Label ribosomal RNA
#
# The below gene list was downloaded from on 3 Aug 2020:
# https://www.genenames.org/data/genegroup/#!/group/848
gene_group__ribo_rna = [
    'MT-RNR1',
    'MT-RNR2',
    'RNA18S1',
    'RNA18S2',
    'RNA18S3',
    'RNA18S4',
    'RNA18S5',
    'RNA18SN1',
    'RNA18SN2',
    'RNA18SN3',
    'RNA18SN4',
    'RNA18SN5',
    'RNA28S1',
    'RNA28S2',
    'RNA28S3',
    'RNA28S4',
    'RNA28S5',
    'RNA28SN1',
    'RNA28SN2',
    'RNA28SN3',
    'RNA28SN4',
    'RNA28SN5',
    'RNA45S1',
    'RNA45S2',
    'RNA45S3',
    'RNA45S4',
    'RNA45S5',
    'RNA45SN1',
    'RNA45SN2',
    'RNA45SN3',
    'RNA45SN4',
    'RNA45SN5',
    'RNA5-8S1',
    'RNA5-8S2',
    'RNA5-8S3',
    'RNA5-8S4',
    'RNA5-8S5',
    'RNA5-8SN1',
    'RNA5-8SN2',
    'RNA5-8SN3',
    'RNA5-8SN4',
    'RNA5-8SN5',
    'RNA5S1',
    'RNA5S10',
    'RNA5S11',
    'RNA5S12',
    'RNA5S13',
    'RNA5S14',
    'RNA5S15',
    'RNA5S16',
    'RNA5S17',
    'RNA5S2',
    'RNA5S3',
    'RNA5S4',
    'RNA5S5',
    'RNA5S6',
    'RNA5S7',
    'RNA5S8',
    'RNA5S9',
    'RNR1',
    'RNR2',
    'RNR3',
    'RNR4',
    'RNR5'
]
# Label mitochondrial encoded trancripts
# This includes:
# * Mitochondrially encoded protein coding genes
# * Mitochondrially encoded transcribed regions
# * Mitochondrially encoded RNAs
#
# The below gene list was downloaded from on 3 Aug 2020:
# https://www.genenames.org/data/genegroup/#!/group/1972
gene_group__mito_transcript = [
    'MT-7SDNA',
    'MT-ATP6',
    'MT-ATP8',
    'MT-ATT',
    'MT-CO1',
    'MT-CO2',
    'MT-CO3',
    'MT-CSB1',
    'MT-CSB2',
    'MT-CSB3',
    'MT-CYB',
    'MT-HPR',
    'MT-HSP1',
    'MT-HSP2',
    'MT-LIPCAR',
    'MT-LSP',
    'MT-ND1',
    'MT-ND2',
    'MT-ND3',
    'MT-ND4',
    'MT-ND4L',
    'MT-ND5',
    'MT-ND6',
    'MT-OHR',
    'MT-OLR',
    'MT-RNR1',
    'MT-RNR2',
    'MT-RNR3',
    'MT-TA',
    'MT-TAS',
    'MT-TC',
    'MT-TD',
    'MT-TE',
    'MT-TER',
    'MT-TF',
    'MT-TFH',
    'MT-TFL',
    'MT-TFX',
    'MT-TFY',
    'MT-TG',
    'MT-TH',
    'MT-TI',
    'MT-TK',
    'MT-TL1',
    'MT-TL2',
    'MT-TM',
    'MT-TN',
    'MT-TP',
    'MT-TQ',
    'MT-TR',
    'MT-TS1',
    'MT-TS2',
    'MT-TT',
    'MT-TV',
    'MT-TW',
    'MT-TY'
]

def check_adata(adata, adata_id):
    """Check adata."""
    # get var names (genes) that are not unique
    vals, counts = np.unique(adata.var_names, return_counts=True)
    if np.sum(counts > 1):
        print('[{}] fixing {} duplicate var_names:\t{}'.format(
            adata_id,
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))
        adata.var_names_make_unique()

    # get obs_names (cell barcodes) that are not unique
    vals, counts = np.unique(adata.obs_names, return_counts=True)
    if np.sum(counts > 1):
        raise Exception('[{}] error {} duplicate obs_names:\t{}'.format(
            adata_id,
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))

    return(adata)


def scanpy_merge(
    h5ad_data,
    metadata,
    metadata_key,
    output_file,
    cellmetadata_filepaths=None,
    anndata_compression_opts=4,
    extra_metadata=None
):
    """Merge h5ad data.

    Parameters
    ----------
    h5ad_data : pandas.DataFrame
        Description of parameter `h5ad_data`.
    metadata : pandas.DataFrame
        Description of parameter `metadata`.
    output_file : string
        Description of parameter `output_file`.
    metadata_key : string
        Column in metadata that matches the "experiment_id" column in
        h5ad_data.
    anndata_compression_opts : int
        Anndata gzip compression level.

    Returns
    -------
    output_file : string
        output_file
    """

    # check the h5ad_data
    # check for required columns
    h5ad_data_required_cols = set(['experiment_id', 'data_path_h5ad_format'])
    if not h5ad_data_required_cols.issubset(h5ad_data.columns):
        raise Exception('Invalid h5ad_data.')
    # check no duplicate sample ids
    vals, counts = np.unique(h5ad_data['experiment_id'], return_counts=True)
    if np.sum(counts > 1):
        raise Exception('Error {} duplicate experiment_ids:\t{}'.format(
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))
    # check all files exist
    filt = h5ad_data['data_path_h5ad_format'].apply(
        lambda x: os.path.exists(x)
    )
    if np.sum(filt > 1):
        raise Exception('Error {} data_path_h5ad_format missing:\t{}'.format(
            np.sum(filt > 1),
            np.array2string(h5ad_data['data_path_h5ad_format'][filt])
        ))

    # Check the cellmetadata_filepaths if we have it
    if cellmetadata_filepaths is not None:
        # check all files exist
        filt = cellmetadata_filepaths['data_path_cellmetadata'].apply(
            lambda x: os.path.exists(x)
        )
        if np.sum(filt > 1):
            raise Exception(
                'Error {} data_path_cellmetadata missing:\t{}'.format(
                    np.sum(filt > 1),
                    np.array2string(
                        cellmetadata_filepaths['data_path_cellmetadata'][filt]
                    )
                )
            )


    # Init a dictionary to record the original number of cells per sample and
    # the number of cells after each filter.
    n_cells_dict = {}

    # Iterate over samples and load data
    adatasets = []
    adatasets__experiment_ids = []
    n_adatasets = 1
    
    for idx, row in h5ad_data.iterrows():
        print(row)
        idx1 = row['experiment_id'].split('__')[0]
        # Load the data
        adata = sc.read_h5ad(
            filename=row['data_path_h5ad_format'] + '/h5ad.h5ad' #,
            # var_names='gene_symbols',
            # ivar_names='gene_ids',
            # make_unique=False
        )

        adata_orig_cols = list(adata.obs.columns)
        adata_orig_cols.append("donor")
        adata = check_adata(adata, row['experiment_id'])

        # Record the total number of cells for this experiment_id
        n_cells_dict[row['experiment_id']] = {}
        n_cells_dict[row['experiment_id']]['before_filters'] = adata.n_obs
        if (len(adata.obs)<3):
            # we skip the donors that have 1 or 2 cells as this will cause issues down the line
            continue

        adata.var['gene_group__mito_transcript'] = [
            x in gene_group__mito_transcript for x in adata.var['gene_symbols']
        ]

        adata.var['gene_group__mito_protein'] = [
            x in gene_group__mito_protein for x in adata.var['gene_symbols']
        ]

        adata.var['gene_group__ribo_protein'] = [
            x in gene_group__ribo_protein for x in adata.var['gene_symbols']
        ]

        adata.var['gene_group__ribo_rna'] = [
            x in gene_group__ribo_rna for x in adata.var['gene_symbols']
        ]

        # Add in sample metadata.
        # NOTE: it would be more memory efficient to stash this in
        #       unstructured dict-like annotation (adata.uns)
        metadata_smpl = metadata[
            metadata[metadata_key] == row['experiment_id']
        ]
        try:
            extra_sample_metadata = extra_metadata[extra_metadata[metadata_key]==idx1]
            extra_sample_metadata= extra_sample_metadata.drop_duplicates(subset=['experiment_id'], keep='last')
        except:
            extra_sample_metadata = pd.DataFrame()
        if (len(extra_sample_metadata)>0):
            for col in extra_sample_metadata.columns:
                # print(col)
                if col in list(adata_orig_cols):
                    print(f' {col} already exist')
                else:
                    print(col)
                    # d1 = extra_sample_metadata[col].values
                    d2 = extra_sample_metadata[col].values
                    adata.obs[col] = np.repeat(list(set(d2)), adata.n_obs)
        adata_orig_cols = list(adata.obs.columns)
        
        for col in metadata_smpl.columns:
            if col in list(adata_orig_cols):
                print(f' {col} already exist')
            if (col == 'experiment_id'):
                print(col)
                adata.obs[col] = np.repeat(metadata_smpl[col].values, adata.n_obs)             
            else:
                print(col)
                adata.obs[col] = np.repeat(metadata_smpl[col].values, adata.n_obs)

        # Ensure we have experiment_in the final dataframe.
        if 'experiment_id' not in adata.obs.columns:
            adata.obs['experiment_id'] = adata.obs[metadata_key]

        # Add in per cell metadata if we have it.
        if cellmetadata_filepaths is not None:
            if row['experiment_id'] in cellmetadata_filepaths.index:
                cellmetadata = pd.read_csv(
                    cellmetadata_filepaths.loc[
                        row['experiment_id'], 'data_path_cellmetadata'
                    ],
                    sep='\t',
                    index_col='cell_barcode'
                )

                for col in cellmetadata.columns:
                    adata.obs[col] = cellmetadata.loc[adata.obs.index, col]

        if adata.n_obs > 0:
            adatasets.append(adata)
            adatasets__experiment_ids.append(row['experiment_id'])
            n_adatasets += 1
        else:
            raise Exception(
                'Error invalid h5ad_data file. Missing coluns.'
            )

    adata_merged = adatasets[0].concatenate(
        *adatasets[1:], join='outer',
        batch_categories=adatasets__experiment_ids
    )
    try:
        print(adata_merged.obs['batch'])
    except:
        _='only one variable, hence batch was not added'
        adata_merged.obs['batch']=adatasets__experiment_ids[0]
    adata_merged = check_adata(adata_merged, 'adata_merged')
    adata_merged.obs['experiment_id'] = adata_merged.obs['experiment_id'].str.split('__donor').str[0]
    # Re-calculate basic qc metrics of var (genes) for the whole dataset.
    # NOTE: we are only changing adata.var
    # obs_prior = adata_merged.obs.copy()
    try:
        sc.pp.calculate_qc_metrics(
            adata_merged,
            qc_vars=[
                'gene_group__mito_transcript',
                'gene_group__mito_protein',
                'gene_group__ribo_protein',
                'gene_group__ribo_rna'
            ],
            inplace=True
        )
        adata_merged.obs['pct_counts_gene_group__mito_transcript'].fillna(0, inplace=True)
        adata_merged.obs['pct_counts_gene_group__mito_protein'].fillna(0, inplace=True)
        adata_merged.obs['pct_counts_gene_group__ribo_rna'].fillna(0, inplace=True)
        adata_merged.obs['pct_counts_gene_group__ribo_protein'].fillna(0, inplace=True)
    except:
        _='most likely different data format such as ATAC which doesnt have gene IDs'
    # adata_merged.obs = obs_prior
    adata_merged.obs['log10_ngenes_by_count'] = np.log10(adata_merged.obs['n_genes_by_counts']) / np.log10(adata_merged.obs['total_counts'])
    adata_merged.write(
        '{}.h5ad'.format(output_file),
        compression='gzip',
        compression_opts=anndata_compression_opts
    )



def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge h5ad data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-txd', '--h5addata_file',
        action='store',
        dest='txd',
        required=True,
        help='File with the following headers: experiment_id\
            data_path_h5ad_format.'
    )

    parser.add_argument(
        '-mf', '--sample_metadata_file',
        action='store',
        dest='mf',
        required=True,
        help='File with metadata on samples matching experiment_id column in\
            h5addata_file. The column that links these two files is specified\
            in metadata_key.'
    )

    parser.add_argument(
        '-mcd', '--sample_metadata_columns_delete',
        action='store',
        dest='mcd',
        default='sample_status,study,study_id',
        help='Comma seperated list of columns to delete in\
            sample_metadata_file. If "" then no columns are deleted.\
            Note: whitespace should be represented with an underscore (_).\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-mk', '--metadata_key',
        action='store',
        dest='mk',
        default='experiment_id',
        help='Key to link metadata to h5addata_file experiment_id column.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cmf', '--cell_metadata_file',
        action='store',
        dest='cmf',
        default='',
        help='File with the following headers: experiment_id\
            data_path_cellmetadata, where data_path_cellmetadata is the path\
            to a tsv file contaning a cell_barcode column followed by other\
            columns to add to the annotations of each cell of that experiment\
            (e.g., doublet scores)\
            (default: None)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=2,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )


    parser.add_argument(
        '-em', '--extra_metadata',
        action='store',
        dest='extra_metadata',
        default='',
        help='Provide extra sample metadata file to be merged with the input'
    )

    parser.add_argument(
        '-pyml', '--params_yaml',
        action='store',
        dest='pyml',
        default='',
        help='YAML file containing cell filtering and downsampling (e.g.,\
            of total number of cells or reads per cell) parameters.\
            If file is not provided, no filtering or downsampling is\
            performed.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5ad appended.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--anndata_compression_opts',
        action='store',
        dest='anndata_compression_opts',
        default=4,
        type=int,
        help='Compression level in anndata. A larger value decreases disk \
            space requirements at the cost of compression time. \
            (default: %(default)s)'
    )

    options = parser.parse_args()
    print('\n0025-scanpy_merge_from_h5ad.py options: ' + str(options) + '\n')

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs

    # Load a file of the samples to analyse
    h5ad_data = pd.read_csv(options.txd, sep='\t')
    # h5ad_data = h5ad_data.rename(columns={
    #     'sanger_sample_id': 'experiment_id',
    #     'file_path': 'data_path_h5ad_format'
    # })
    h5ad_data_check = [
        'experiment_id',
        'data_path_h5ad_format'
    ]
    if not all(k in h5ad_data.columns for k in h5ad_data_check):
        raise Exception(
            'Error invalid h5ad_data file. Missing coluns.'
        )

    # Load the metadata
    metadata = pd.read_csv(options.mf, sep='\t')
    metadata.columns = metadata.columns.str.strip(
        ).str.replace(' ', '_').str.lower()

    # 
    try:
        extra_metadata = pd.read_csv(options.extra_metadata, sep='\t')
    except:
        extra_metadata = pd.DataFrame()

    # Delete the metadata columns that we do not want.
    if options.mcd != '':
        for i in options.mcd.split(','):
            if i in metadata.columns:
                metadata = metadata.drop(i, axis=1)

    # Make sure the matching key exists.
    if options.mk not in metadata.columns:
        raise Exception(
            'Error cannot find metadata_key in metadata.'
        )

    # Load cell metadata
    cellmetadata_filepaths = None
    if options.cmf != '':
        cellmetadata_filepaths = pd.read_csv(
            options.cmf,
            sep='\t',
            index_col='experiment_id'
        )

    # Run the merge function.
    scanpy_merge(
        h5ad_data,
        metadata,
        metadata_key=options.mk,
        output_file=options.of,
        cellmetadata_filepaths=cellmetadata_filepaths,
        anndata_compression_opts=options.anndata_compression_opts,
        extra_metadata= extra_metadata
    )



if __name__ == '__main__':
    main()
