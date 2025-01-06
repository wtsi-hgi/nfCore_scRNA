#!/usr/bin/env python

__date__ = '2023-03-13'
__version__ = '0.0.2'

import argparse
import scanpy as sc
from typing import Dict
import numpy as np
import tables
from distutils.version import LooseVersion
from typing import Dict
import tables
import scipy
import scipy.io
import scipy.sparse
import gzip
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import click
import logging
import os
import re
compression_opts = 'gzip'
filter_0_count_cells=False

def cellbender_to_tenxmatrix(adata,out_file='',out_dir='tenx_from_adata',verbose=True):
    """Write 10x like data from 10x H5.

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    output_file : string
        Description of parameter `output_file`.
    out_dir : string
        Description of parameter `out_dir`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    # Set up out_file
    if out_file != '':
        out_file = '{}-'.format(out_file)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # First filter out any cells that have 0 total counts
    if(filter_0_count_cells):
        zero_count_cells = adata.obs_names[np.where(adata.X.sum(axis=1) == 0)[0]]
        if verbose:
            print("Filtering {}/{} cells with 0 counts.".format(
                len(zero_count_cells),
                adata.n_obs
            ))
        adata = adata[adata.obs_names.difference(zero_count_cells, sort=False)]

    # Save the barcodes.
    out_f = os.path.join(
        out_dir,
        '{}barcodes.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
    pd.DataFrame(adata.obs.index).to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the features.
    # CellBender output: gene_symbols as index, enseble as `gene_ids`
    # Type is stored in adata.var.type
    out_f = os.path.join(
        out_dir,
        '{}features.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))

    # Check if read in by gene_symbol or gene_id

    if 'gene_symbols' in adata.var.columns:
        gene_var = 'gene_symbols'
    elif 'gene_ids' in adata.var.columns:
        gene_var = 'gene_ids'
    # elif 'id' in adata.var.columns:
    #     gene_var = 'id'
    # elif 'name' in adata.var.columns:
    #     gene_var = 'name'
    else:
        raise Exception(
            'Could not find "gene_symbols" or "gene_ids" in adata.var'
        )
    df_features = pd.DataFrame(
        data=[
            adata.var.index.values,
            adata.var.loc[:, gene_var].values,
            adata.var.feature_types.values
        ]
    )
    if df_features.shape[0]<df_features.shape[1]:
        df_features=df_features.T
    
    df_features.to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the count-adjusted matrix
    out_mtx = adata.X.transpose()
    if not isinstance(out_mtx, scipy.sparse.csr.csr_matrix):
        out_mtx = scipy.sparse.csr_matrix(out_mtx)
    out_f = os.path.join(
        out_dir,
        '{}matrix.mtx.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
        # print(out_mtx.dtype)  # it looks like even counts stored as float
    with gzip.open(out_f, 'wb', compresslevel=9) as f:
        scipy.io.mmwrite(
            f,
            out_mtx,
            comment='metadata_json: {"format_version": 2}'
            # field='real'  # can be integer if counts otherwise real
        )

    if verbose:
        print(
            'Done.'
        )

    return 0

def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d

def anndata_from_h5(file: str,analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent
            variables to be loaded properly into adata.obs and adata.obsm,
            rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """
    d = dict_from_h5(file)
    X = scipy.sparse.csc_matrix(
        (d.pop('data'), d.pop('indices'), d.pop('indptr')),
        shape=d.pop('shape')
    ).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print(
                'Warning: analyzed_barcodes_only=True, but the key ',
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" ',
                'is missing from the h5 file. ',
                'Will output all barcodes, and proceed as if ',
                'analyzed_barcodes_only=False'
            )

    print(d.keys())

    # Construct the count matrix.
    if 'gene_names' in d.keys():
        gene_symbols = d.pop('gene_names').astype(str)
    else:
        gene_symbols = d.pop('name').astype(str)
    adata = anndata.AnnData(
        X=X,
        obs={'barcode': d.pop('barcodes').astype(str)},
        var={
            'gene_ids': d.pop('id').astype(str),
            'gene_symbols': gene_symbols,
            'feature_type':d.pop('feature_type').astype(str),
        }
    )
    # adata = adata[:,adata.var.query('feature_type=="Gene Expression"').index]
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_ids', inplace=True)
    
    
    
    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print(
                'Unable to load data into AnnData: ', key, value, type(value)
            )

    if analyzed_barcodes_only:
        cols = adata.obs.columns[
            adata.obs.columns.str.startswith('barcodes_analyzed')
            | adata.obs.columns.str.startswith('barcode_indices')
        ]
        for col in cols:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata

# We load the cellbender data and check if it contains any antibody data. If it does we split the antibody and gex data and pass only the gex data to the cellbender.
# The Antibody data is subjected to the dsb.

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read Cellranger mtx data and split the citeseq and gex data.
            """
    )

    parser.add_argument(
        '-rw', '--raw_data',
        action='store',
        dest='raw_data',
        required=True,
        help='H5 AnnData file.'
    )
    
    parser.add_argument(
        '-o', '--outname',
        action='store',
        dest='outname',
        required=True,
        help='outputname'
    )   
    
    parser.add_argument(
        '-ha', '--hastag_labels',
        action='store',
        dest='hastag_labels',
        required=False,default=None,
        help='outputname'
    )   
    
    
    options = parser.parse_args()

    adata_cellranger_filtered = sc.read_10x_mtx(
        options.raw_data, var_names='gene_symbols', make_unique=True,
        cache=False, cache_compression=compression_opts,gex_only=False)
    all_feature_types = set(adata_cellranger_filtered.var['feature_types'])
    hashtags = set(options.hastag_labels.split(","))
    hashtags = ['Hashtag_.*']
    escaped_hashtags = [re.escape(tag) for tag in hashtags]
    matches = set(adata_cellranger_filtered.var.index[adata_cellranger_filtered.var.index.str.contains('|'.join(escaped_hashtags), regex=True)])
    matches2 = set(adata_cellranger_filtered.var.index[adata_cellranger_filtered.var.index.str.contains('|'.join(hashtags), regex=True)])
    combo = matches.union(matches2)
    if len(combo)>0:
        multiplexing_capure = adata_cellranger_filtered[:,list(combo)]
        multiplexing_capure = pd.DataFrame(multiplexing_capure.X.toarray(), index=multiplexing_capure.obs_names, columns=multiplexing_capure.var_names)
        multiplexing_capure.to_csv(f'{options.outname}__Multiplexing_Capture.tsv',sep='\t')
        
    for modality1 in set(adata_cellranger_filtered.var.feature_types):
        # {'Gene Expression', 'Multiplexing Capture', 'Antibody Capture'}
        print(f'---- Spliting {modality1}-----')
        modality = modality1.replace(" ",'_')
        adata_antibody = adata_cellranger_filtered[:,adata_cellranger_filtered.var.query(f'feature_types=="{modality1}"').index]
        # zero_count_cells = adata_antibody.obs_names[np.where(adata_antibody.X.sum(axis=1) == 0)[0]]
        # adata2 = adata_antibody[adata_antibody.obs_names.difference(zero_count_cells, sort=False)]
        # if(adata2.shape[0]>0):
        #     # Here we have actually captured some of the reads in the antibody dataset.

        if (modality=='Gene_Expression'):
            adata_antibody.write(
                f'{modality}-{options.outname}.h5ad',
                compression='gzip'
            )
        elif (modality=='Multiplexing_Capture'):
            all_indexes_multiplexing = set(adata_antibody.var.index).union(set(multiplexing_capure.columns))
            adata_antibody = adata_cellranger_filtered[:,list(all_indexes_multiplexing)]
            df = pd.DataFrame(adata_antibody.X.toarray(), index=adata_antibody.obs_names, columns=adata_antibody.var_names)
            df.to_csv(f'{options.outname}__{modality}.tsv',sep='\t')
            
        else:
            df = pd.DataFrame(adata_antibody.X.toarray(), index=adata_antibody.obs_names, columns=adata_antibody.var_names)
            df.to_csv(f'{options.outname}__{modality}.tsv',sep='\t')
            adata_antibody.var.index

        if len(all_feature_types)>1:
            _ = cellbender_to_tenxmatrix(
                adata_antibody,
                out_file='',
                out_dir=f'{options.outname}__{modality}'
            )
        else:
            os.system(f"ln -s {options.raw_data} {options.outname}__{modality}")
    
    # adata_gex = adata_cellranger_filtered[:,adata_cellranger_filtered.var.query('feature_types=="Gene Expression"').index]
    # adata_cellbender = anndata_from_h5('/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/5d/6a30871e864ed7bc03e949ef846a1d/cellbender_FPR_0.1_filtered.h5',
    #                                         analyzed_barcodes_only=True)

    
    
if __name__ == '__main__':
    main()