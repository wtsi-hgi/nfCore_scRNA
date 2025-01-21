#!/usr/bin/env python

__date__ = '2020-05-26'
__version__ = '0.0.1'

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import scvi
import argparse
import glob
import pandas as pd
import os

COLORS_LARGE_PALLETE = [
    '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
    '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
    '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
    '#77783C', '#EF4E22', '#AF1F26'
]

def save_plot(
    adata,
    out_file_base,
    color_var,
    colors_quantitative=True,
    colors_large_palette=COLORS_LARGE_PALLETE,
):
    
        
    def panel_grid(hspace, wspace, ncols, num_panels):
        """Init plot."""
        n_panels_x = min(ncols, num_panels)
        n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
        if wspace is None:
            #  try to set a wspace that is not too large or too small given the
            #  current figure size
            wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
        # each panel will have the size of rcParams['figure.figsize']
        fig = plt.figure(
            figsize=(
                n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
                n_panels_y * rcParams['figure.figsize'][1],
            )
        )
        left = 0.2 / n_panels_x
        bottom = 0.13 / n_panels_y
        gs = gridspec.GridSpec(
            nrows=n_panels_y,
            ncols=n_panels_x,
            left=left,
            right=1 - (n_panels_x - 1) * left - 0.01 / n_panels_x,
            bottom=bottom,
            top=1 - (n_panels_y - 1) * bottom - 0.1 / n_panels_y,
            hspace=hspace,
            wspace=wspace
        )
        return fig, gs

    fig, grid = panel_grid(
        hspace=0.125*2,
        wspace=None,
        ncols=4,
        num_panels=2
    )
    ax = plt.subplot(grid[0])
    # sc.pl.umap(
    #     SLEmap,
    #     color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2"],
    #     # Setting a smaller point size to get prevent overlap
    #     size=0.7,save='plots.pdf'
    # )

    # fig.savefig(
    #     'umap.png',
    #     #dpi=300,
    #     bbox_inches='tight'
    # )    
    
    """Save a plot."""


    if colors_quantitative is False:
        # Cast to category - required for booleans.
        adata.obs[color_var] = adata.obs[color_var].astype('category')
        n_categories = len(adata.obs[color_var].cat.categories)
        color_palette = None
        if n_categories <= len(plt.get_cmap('Dark2').colors):
            color_palette = 'Dark2'
        elif n_categories <= len(colors_large_palette):
            color_palette = colors_large_palette
    else:
        # Make sure that the numeric value is actually numeric
        adata.obs[color_var] = adata.obs[color_var].astype('double')
    legend_loc = 'right margin'
    color_palette = 'viridis'
    sc.pl.umap(
        adata=adata,
        color=color_var,
        palette=color_palette,
        alpha=0.4,
        title=color_var,
        show=False,
        ax=ax 
    )
    try:
        os.mkdir(f"{out_file_base}")
    except:
        _='exists'
        
    fig.savefig(
        f'{out_file_base}/{color_var}-umap.png',
        dpi=300,
        bbox_inches='tight'
    )


# This code will gather the most important graphs and place them in a summary folder.
parser = argparse.ArgumentParser(
    description="""
        Collect inputs
        """
)

parser.add_argument(
        '-h5ad_file', '--h5ad_file',
        action='store',
        dest='h5ad_file',
        required=True,
        help='Input h5ad_file after normalisation'
)

parser.add_argument(
    '-cq', '--colors_quantitative',
    action='store',
    dest='cq',
    default='',
    help='Comma seperated list of quantitative variable names for colors.\
        (default: "")'
)

parser.add_argument(
    '-cc', '--colors_categorical',
    action='store',
    dest='cc',
    default='',
    help='Comma seperated list of categorical variable names for colors.\
        (default: "")'
)

parser.add_argument(
    '-reduction_columns_cells', '--reduction_columns_cells',
    action='store',
    dest='reduction_columns_cells',
    default='cell_passes_qc,cell_passes_hard_filters',
    help='reduction_columns_cells'
)

parser.add_argument(
    '-reduction_columns_genes', '--reduction_columns_genes',
    action='store',
    dest='reduction_columns_genes',
    default='highly_variable',
    help='reduction_columns_genes'
)

parser.add_argument(
    '-gene_list_to_keep', '--gene_list_to_keep',
    action='store',
    dest='gene_list_to_keep',
    default='fake_fileinput.tsv',
    help='gene_list_to_keep'
)

options = parser.parse_args()

sc.set_figure_params(figsize=(5,5), dpi=150)

#import single cell data and CITE-seq data
# SLEmap = sc.read('adata-normalized.h5ad')
SLEmap = sc.read(options.h5ad_file, backed ='r')

all_cite_files = glob.glob("./*/*.matrix.csv")
CITE = pd.DataFrame()
for f1 in all_cite_files:
    c1 = pd.read_csv(f1,index_col=0)
    CITE=pd.concat([CITE,c1])

# merge citeseq data to obs of single cell data 
SLEmap.obs["Barcode"] = SLEmap.obs.index.str.split('__').str[0]
SLEmap.obs = SLEmap.obs.merge(CITE, left_on=['Barcode'],right_index=True,how='left', indicator=True)


# make citeseq data to obsm single cell data : required for totalVI
CITE_2 = SLEmap.obs[CITE.columns].copy()
SLEmap.obsm['protein_expression'] = CITE_2

cell_filters = options.reduction_columns_cells.split(';')
number_of_cell_filters = len(cell_filters)

gene_filters = options.reduction_columns_genes.split(';')
number_of_gene_filters = len(gene_filters)

if f"{options.gene_list_to_keep}"=='fake_fileinput.tsv':
    gene_list_to_keep = None
else:
    gene_list_to_keep = pd.read_csv(f"{options.gene_list_to_keep}",sep='\t')


# Combine cell filters dynamically
cell_condition = SLEmap.obs[cell_filters[0]]
for filter_name in cell_filters[1:]:
    cell_condition &= SLEmap.obs[filter_name]

# Combine gene filters dynamically
gene_condition = SLEmap.var[gene_filters[0]]
for filter_name in gene_filters[1:]:
    gene_condition &= SLEmap.var[filter_name]

# Add gene list filter if provided
if gene_list_to_keep is not None:
    gene_condition &= SLEmap.var_names.isin(gene_list_to_keep)

# Subset the dataset
SLEmap = SLEmap[cell_condition, gene_condition]


#run totalVI
SLEmap = SLEmap.to_memory().copy()
scvi.model.TOTALVI.setup_anndata(SLEmap, protein_expression_obsm_key="protein_expression",batch_key="experiment_id")
model = scvi.model.TOTALVI(SLEmap,latent_distribution="normal",n_layers_decoder=2)

model.train()
model.save("./scvi_model",adata=SLEmap, overwrite=True)


# Get latent expression from model: used for UMAP calculations
SLEmap.obsm["X_totalVI"] = model.get_latent_representation()


# Get denoised protein values: can fail with low memory: use 400 minimum for n_samples=25. 
# If it keeps failing lower n_samples or don't use it and use dsb values instead
rna, protein = model.get_normalized_expression(n_samples=25,return_mean=True)
SLEmap.obsm["denoised_protein"] = protein

sc.pp.neighbors(SLEmap, use_rep="X_totalVI",n_neighbors=15)

sc.tl.umap(SLEmap)

sc.set_figure_params(figsize=(5,5), dpi=150)



colors_quantitative = options.cq.split(',')
colors_categorical = options.cc.split(',')
adata = SLEmap
for color_var in colors_quantitative:
    try:       
        save_plot(
            adata=adata,
            out_file_base="figures",
            color_var=color_var,
            colors_quantitative=True,
        )
    except:
        print(f'{color_var} doesnt')
        
for color_var in colors_categorical:
    
    try:
        save_plot(
            adata=adata,
            out_file_base="figures",
            color_var=color_var,
            colors_quantitative=False,
        )
    except:
        print(f'{color_var} doesnt')

SLEmap.write("./totalVI_integrated.h5ad")