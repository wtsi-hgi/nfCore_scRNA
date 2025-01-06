#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.2'

import argparse
from distutils.version import LooseVersion
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import csv
import random
import numpy as np
import pandas as pd
import scrublet as scr

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import plotnine as plt9
import skimage.filters as skif
# from statannot import add_stat_annotation

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
# os.environ['PYTHONHASHSEED']=str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

sns.set(style='whitegrid')

# Set the default dpi
plt9.options.dpi = 100


def custom_cmap(rgb_list):
    """Make a custom cmap."""
    rgb_list = np.array(rgb_list)
    cmap = plt.cm.Reds
    cmap = cmap.from_list(rgb_list.shape[0], rgb_list)
    return cmap


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def plot_scrub_hist(
    scrub,
    threshold,
    output_file='scrublet_histogram',
    scale_y_log10=False,
    density=False,
    zscores=False
):
    """Plot better histogram of multiplets."""
    # Make a pandas dataframe of the data
    if not zscores:
        df_tmp1 = pd.DataFrame(
            data=scrub.doublet_scores_obs_,
            columns=['scores']
        )
        df_tmp2 = pd.DataFrame(
            data=scrub.doublet_scores_sim_,
            columns=['scores']
        )
        x_axis_label = 'Multiplet score'
    else:
        # Calculate z scores as in the scrublet source code (see call_doublets
        # function).
        df_tmp1 = pd.DataFrame(
            data=scrub.z_scores_,
            columns=['scores']
        )
        zscore_sim = (
            scrub.doublet_scores_sim_ - threshold
        ) / scrub.doublet_errors_sim_
        df_tmp2 = pd.DataFrame(
            data=zscore_sim,
            columns=['scores']
        )
        x_axis_label = 'Multiplet zscore'

    df_tmp1['type'] = 'Observed'
    df_tmp2['type'] = 'Simulated'

    df_plt = pd.concat([df_tmp1, df_tmp2])

    # Plot the data usig plotnine
    gplt = plt9.ggplot(df_plt, plt9.aes(
        x='scores'
    ))
    gplt = gplt + plt9.theme_bw()
    if density:
        gplt = gplt + plt9.geom_density(alpha=0.8)
    else:
        gplt = gplt + plt9.geom_histogram(alpha=0.8)
    if not zscores:
        gplt = gplt + plt9.geom_vline(xintercept=threshold, linetype='solid')
    if scale_y_log10:
        gplt = gplt + plt9.scale_y_continuous(
            trans='log10',
            # labels=comma_labels,
            minor_breaks=0
        )
    gplt = gplt + plt9.labs(
        x=x_axis_label,
        # y='Counts',
        title=''
    )
    gplt = gplt + plt9.facet_wrap(
        '~ {}'.format('type'),
        scales='free_y'
    )

    # Sort out whitespace issue in plotnine:
    # https://github.com/has2k1/plotnine/issues/106
    gplt = gplt + plt9.theme(subplots_adjust={'wspace': 0.35})

    gplt.save(
        '{}.png'.format(output_file),
        # dpi=300,
        width=4,
        height=2,
        limitsize=False
    )


def run_scrublet(
    adata,
    out_file_base='scrublet',
    expected_multiplet_rate=0.1,
    n_simulated_multiplet=100000,
    multiplet_threshold_method='threshold_li',
    scale_log10=False,
    verbose=True
):
    """Run scrublet.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object. Assume adata.X is counts.
    expected_multiplet_rate : integer
        Expected multiplet rate.
    out_file_base : string
        Tag for output file.
    n_simulated_multiplet : integer
        Number of multiplets to simulate.
    multiplet_threshold_method : string
        Method used to call multiplets.
    scale_log10 : string
        Scale simulated doublet scores to log10 before deriving the threshold.
    verbose : boolean
        Write extra info to standard out.

    Returns
    -------
    NULL
    """
    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Set color pallette for UMAP.
    zissou_palette_hex = ['#3B9AB2', '#EBCC2A', '#F21A00']
    zissou_palette_rgb = [colors.to_rgba(i) for i in zissou_palette_hex]

    # If no expected multiplet rate, then estimate the multiplet rate using
    # the coefficients from a lm predicting multiplet rate from recovered cells
    # in data distributed by 10x for their 3.1 chemistry.
    # https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry
    cells_recovered = len(adata)
    if expected_multiplet_rate == 0.0:
        multiplet_rate = 0.0007589 * cells_recovered + 0.0527214
        multiplet_rate = multiplet_rate / 100.0
    else:
        multiplet_rate = expected_multiplet_rate
    # expected_n_multiplets = multiplet_rate * cells_recovered
    if verbose:
        print('cells_input_data:\t{}'.format(cells_recovered))
        print('multiplet_rate:\t{}'.format(multiplet_rate))

    # The authors note that the method is not that sensitive to this parameter
    # https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb
    # From Chromium Single Cell 3' Reagent Kits User Guide (v2 Chemistry):
    #   https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC
    #   the expected multiplet rate is 0.06.
    # From Chromium Single Cell 3' Reagent Kits User Guide (v3.1 Chemistry):
    #   https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry
    #   the expected multiplet rate is ~3.9% for ~8000 input cells.
    adata.X.data[adata.X.data < 0] = 0
    scrub = scr.Scrublet(
        counts_matrix=adata.X,
        sim_doublet_ratio=n_simulated_multiplet/len(adata),  # Default is 2.0
        expected_doublet_rate=multiplet_rate,  # Default is 0.1
        random_state=0
    )

    # Run the scrublet pipeline with default parameters.
    # This function performs:
    # * Multiplet simulation
    # * Normalization, gene filtering, rescaling, PCA
    # * Multiplet score calculation
    # * Multiplet score threshold detection and multiplet calling
    try:
        multiplet_scores, predicted_multiplets = scrub.scrub_doublets(
            verbose=True,n_prin_comps=5
        )
    except:
        # this is to resolve an issue sometimes encountered when small nr o cells provided in adata
        # ValueError: n_components=30 must be between 1 and min(n_samples, n_features)=27 with svd_solver='arpack'
        multiplet_scores, predicted_multiplets = scrub.scrub_doublets(
            verbose=True,n_prin_comps=len(adata)-1
        )

    # Calculate the threshold for calling multiplets
    # The default method for scrublet is `threshold_minimum`, but in our
    # hands it looks like `threshold_otsu` looks better.
    # threshold = skif.threshold_minimum(scrub.doublet_scores_sim_)
    sim_doublet_scores = scrub.doublet_scores_sim_
    if scale_log10:
        sim_doublet_scores = np.log10(sim_doublet_scores)
    if multiplet_threshold_method == 'threshold_li':
        threshold = skif.threshold_li(
            sim_doublet_scores,
            initial_guess=skif.threshold_otsu(sim_doublet_scores)
        )
    elif multiplet_threshold_method == 'threshold_minimum':
        threshold = skif.threshold_minimum(sim_doublet_scores)
    elif multiplet_threshold_method == 'threshold_triangle':
        threshold = skif.threshold_triangle(sim_doublet_scores)
    elif multiplet_threshold_method == 'threshold_yen':
        threshold = skif.threshold_yen(sim_doublet_scores)
    elif multiplet_threshold_method == 'threshold_otsu':
        threshold = skif.threshold_otsu(sim_doublet_scores)
    if scale_log10:
        threshold = 10**threshold

    # For debug
    # print('threshold_li:\t{}'.format(
    #     skif.threshold_li(
    #         sim_doublet_scores,
    #         initial_guess=skif.threshold_otsu(sim_doublet_scores)
    #     )
    # ))
    # print('threshold_minimum:\t{}'.format(
    #     skif.threshold_minimum(sim_doublet_scores)
    # ))
    # print('threshold_triangle:\t{}'.format(
    #     skif.threshold_triangle(sim_doublet_scores)
    # ))
    # print('threshold_yen:\t{}'.format(
    #     skif.threshold_yen(sim_doublet_scores)
    # ))
    # print('threshold_otsu:\t{}'.format(
    #     skif.threshold_otsu(sim_doublet_scores)
    # ))
    # print('threshold_used:\t{}'.format(threshold))

    # Call multiplets using the otsu threshold
    predicted_multiplets = scrub.call_doublets(
        threshold=threshold,
        verbose=verbose
    )
    print(
        'If automatic threshold is poor, adjust threshold with',
        'scrub.call_doublets(threshold=<my_custom_threshold>)'
    )
    adata.obs['scrublet__predicted_multiplet'] = scrub.predicted_doublets_
    adata.obs['scrublet__multiplet_scores'] = scrub.doublet_scores_obs_
    adata.obs['scrublet__multiplet_zscores'] = scrub.z_scores_

    # Get the number of multiplets one would expect if the 10x prediction were
    # spot on. Taken from
    # https://github.com/vib-singlecell-nf/scrublet/blob/master/bin/sc_doublet_detection.py
    #
    # Estimate the multiplet rate using coefficients from a lm predicting
    # multiplet rate from recovered cells in data distributed by 10x for
    # their 3.1 chemistry.
    # https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry
    # cells_recovered = len(adata)
    # multiplet_rate = 0.0007589 * cells_recovered + 0.0527214
    # multiplet_rate = multiplet_rate / 100.0
    # multiplet_cells = adata.obs['scrublet__multiplet_scores'].sort_values(
    #     ascending=False
    # ).head(
    #     n=expected_n_doublets
    # ).index
    # adata.obs['scrublet__predicted_multiplet_based_10x_rate'] = False
    # adata.obs.loc[
    #     multiplet_cells,
    #     'scrublet__predicted_multiplet_based_10x_rate'
    # ] = True

    # Save the results.
    cols_save = [
        'scrublet__multiplet_scores',
        'scrublet__predicted_multiplet',
        'scrublet__multiplet_zscores'
    ]
    adata.obs[cols_save].to_csv(
        '{}-scrublet.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        index_label='cell_barcode',
        na_rep='',
        compression=compression_opts
    )

    # Plot a histogram of multiplets.
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        scale_y_log10=False,
        output_file='{}-histogram_multiplet_scores'.format(out_file_base)
    )
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        scale_y_log10=True,
        output_file='{}-histogram_multiplet_scores_log'.format(out_file_base)
    )
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        zscores=True,
        output_file='{}-histogram_multiplet_zscores'.format(out_file_base)
    )
    # plot_scrub_hist(
    #     scrub=scrub,
    #     threshold=threshold,
    #     density=True,
    #     scale_y_log10=False,
    #     output_file='{}-density_multiplet_scores'.format(out_file_base)
    # )
    # fig, ax = scrub.plot_histogram(
    #     scale_hist_obs='linear',
    #     scale_hist_sim='linear'
    # )
    # fig.savefig(
    #     '{}-histogram_multiplet_scores_v2.pdf'.format(out_file_base),
    #     # dpi=300,
    #     bbox_inches='tight'
    # )
    # plt.close(fig)

    # Plot the average number of UMIs in multiplets vs singlets.
    if 'total_counts' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    fig, ax = plt.subplots(figsize=(3, 3.5))
    ax = sns.boxplot(
        data=adata.obs,
        x='scrublet__predicted_multiplet',
        y='total_counts'
        # hue='scrublet__predicted_multiplet'
    )
    ax.set_yscale('log', basey=10)
    ax.set(xlabel='Predicted multiplet', ylabel='Number of molecules')
    # NOTE: I could not get p-value annotation to work.
    # ax, test_results = add_stat_annotation(
    #     ax,
    #     data=adata.obs,
    #     x='scrublet__predicted_multiplet',
    #     y='total_counts',
    #     #hue='scrublet__predicted_multiplet',
    #     box_pairs=[('True', 'False')],
    #     test='Mann-Whitney',
    #     text_format='full',
    #     loc='inside'
    # )
    # plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
    fig.savefig(
        '{}-boxplot_total_umi_counts.png'.format(out_file_base),
        # dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)  # Close the figure.

    # NOTE: Removed UMAP embedding on 30/06/2020 because it would not
    # work with singularity.
    # Plot UMAP embedding.
    # scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_))
    # fig, ax = scrub.plot_embedding(
    #     'UMAP',
    #     marker_size=1.5,
    #     color_map=custom_cmap(zissou_palette_rgb),
    #     order_points=True
    # )
    # fig.savefig(
    #     '{}-umap_multiplet_scores.png'.format(out_file_base),
    #     # dpi=300,
    #     bbox_inches='tight'
    # )
    # plt.close(fig)
    #
    # fig, ax = scrub.plot_embedding(
    #     'UMAP',
    #     score='zscore',
    #     marker_size=1.5,
    #     color_map=custom_cmap(zissou_palette_rgb),
    #     order_points=True
    # )
    # fig.savefig(
    #     '{}-umap_multiplet_zscores.png'.format(out_file_base),
    #     # dpi=300,
    #     bbox_inches='tight'
    # )
    # plt.close(fig)

    # Plot tSNE embedding.
    # scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read 10x data or AnnData object and annotates multiplets.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-txd', '--tenxdata_dir',
        action='store',
        dest='txd',
        default='',
        required=False,
        help='Path to directory with data in 10x matrix format.'
    )

    parser.add_argument(
        '-txh5', '--tenxdata_h5',
        action='store',
        dest='txh5',
        default='',
        required=False,
        help='10x h5 format file.'
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        default='',
        required=False,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-emr', '--expected_multiplet_rate',
        action='store',
        dest='emr',
        default=0.1,
        type=float,
        help='Expected multiplet rate. If 0.0 then automatically predict\
            multiplet rate based on the number of recovered cells using\
            published multiplet rates from 10x v3.1 chemistry. In tests, the\
            scrublet default of 0.1 was more aggressive in calling multiplets\
            than the automatic prediction and was able to capture a greater\
            fraction of simulated multiplets (with the automatically\
            calculated multiplet threshold). Note that the final multiplet\
            scores are calculated with z-scores such that the\
            expected_multiplet_rate does not have a huge effect on the final\
            calls.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-nsm', '--n_simulated_multiplet',
        action='store',
        dest='nsm',
        default=100000,
        type=int,
        help='Number of simulated multiplets. A larger number provides more\
            accurate multiplet calls at the cost of execution time.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-mtm', '--multiplet_threshold_method',
        action='store',
        dest='mtm',
        default='threshold_li',
        type=str,
        help='Method used to call multiplets threshold from simulated doublets \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--scale_log10',
        action='store_true',
        dest='scale_log10',
        help='Scale simulated doublet scores to log10 before deriving the \
            threshold. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='scrublet',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Load the data.
    if options.txd != '':
        adata = sc.read_10x_mtx(
            path=options.txd,
            # var_names='gene_symbols',
            var_names='gene_ids',
            make_unique=False
        )
    elif options.txh5 != '':
        adata = sc.read_10x_h5(path=options.txd)
        adata.var['gene_symbols'] = adata.var.index
        adata.var.index = adata.var['gene_ids'].values
        del adata.var['gene_ids']
    elif options.h5 != '':
        adata = sc.read_h5ad(filename=options.h5)
        print(
            'Scrublet uses counts. Assuming adata.X are counts'
        )
    else:
        raise Exception(
            'Error invalid input.'
        )

    run_scrublet(
        adata=adata,
        out_file_base=options.of,
        expected_multiplet_rate=options.emr,
        n_simulated_multiplet=options.nsm,
        multiplet_threshold_method=options.mtm,
        scale_log10=options.scale_log10,
        verbose=verbose
    )


if __name__ == '__main__':
    main()
