#!/usr/bin/env python


__date__ = '2021-01-20'
__version__ = '0.0.1'
# sklearn version used = 0.24.2
import argparse
import random
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
from distutils.version import LooseVersion
import numpy as np
import scanpy as sc
import re
import pandas as pd
from sklearn.svm import OneClassSVM
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
import matplotlib
matplotlib.use('Agg')
y_n_print=False
# import matplotlib.pyplot as plt
# matplotlib.style.use('ggplot')
import seaborn as sns

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


class onesidemad():
    def __init__(self, thresh):
        self.thresh=thresh.reshape(1, -1)
        self.median=None
        self.mad=None

    def fit(self, X):
        if X.size==0:
            return
        self.median=np.median(X, axis=0).reshape(1, -1)
        self.mad=np.median(np.abs(X - self.median), axis=0).reshape(1, -1)
    
    def stats(self):
        if self.median is None or self.mad is None:
            raise Exception("need to run fit() before")
        
        return self.median.ravel(), self.mad.ravel(), self.thresh.ravel(), (self.median + self.thresh*self.mad).ravel()

    def predict(self, X):
        if self.median is None or self.mad is None:
            raise Exception("need to run fit() before")
        
        return (np.sign(self.thresh) * (X - self.median - self.thresh*self.mad) > 0).any(axis=1)

    def decision_function(self,cols):
        tresholds_set = (self.median + self.thresh*self.mad).ravel()
        thresh = self.thresh.ravel()
        decission_string = ""
        for i,tr1 in enumerate(thresh):
            if tr1<0:
                trs=f"<{tresholds_set[i]}"
            else:
                trs=f">{tresholds_set[i]}"
            decission_string = f"{decission_string}{tr1}({cols[i]}{trs})"
        return decission_string

    def fit_predict(self, X):
        self.fit(X)
        return self.predict(X)

# Fit the model
def perform_adaptiveQC_Filtering(clf,adata,method,metadata_columns):
    # We return labels of the Treue/False of failing QC
    if method == 'LocalOutlierFactor':
        # use fit_predict to compute the predicted labels of the training
        # samples (when LOF is used for outlier detection, the estimator has
        # no predict, decision_function and score_samples methods).
        f = clf.fit_predict(
            adata.obs[metadata_columns].values
        ) == 1
        predicted_scores = clf.negative_outlier_factor_
    elif method == 'IsolationForest':
        f = clf.fit_predict(
            adata2.obs['pct_counts_gene_group__mito_transcript'].values
        ) == 1
        predicted_scores = clf.decision_function(adata.obs[metadata_columns].values) 
    elif method == 'MAD':
        metadata_columns2 = metadata_columns.copy()
        metadata_columns2.pop()
        f = clf.fit_predict(
            adata.obs[metadata_columns2].values
        ) == 0
        predicted_scores = clf.decision_function(metadata_columns2)
        
    else:
        f = clf.fit_predict(
            adata.obs[metadata_columns].values
        ) == 1
        predicted_scores = clf.decision_function(adata.obs[metadata_columns].values) 

    return predicted_scores,f

def generate_plots(adata,cell_qc_column,metadata_columns,metadata_columns_original,of):
    # Plot the identified outliers
    
   
    # print("Making plot")

    sns_plot = sns.PairGrid(
        adata.obs[set(metadata_columns)],
        hue=cell_qc_column,
        height=2.5,
        diag_sharey=False
    )
    sns_plot.map_upper(
        sns.scatterplot,
        marker='+',
        alpha=0.05,
        s=8,
        edgecolor=None
    )
    sns_plot.add_legend()
    for lh in sns_plot._legend.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]
    sns_plot.map_diag(sns.kdeplot)
    sns_plot.map_lower(
        sns.kdeplot,
        levels=3
    )
    sns_plot.savefig(f'{of}-outlier_cells__{cell_qc_column}.png')

    # Plot the cell density
    # print("Making density plot")
    sns_plot = sns.PairGrid(
        adata.obs[set(metadata_columns_original)],
        height=2.5,
        diag_sharey=False
    )
    sns_plot.map_upper(
        sns.kdeplot,
        cmap="viridis",
        shade=True
    )
    sns_plot.add_legend()
    for lh in sns_plot._legend.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]
    sns_plot.map_diag(sns.kdeplot)
    sns_plot.map_lower(
        sns.kdeplot,
        cmap="viridis",
        shade=True
    )
    sns_plot.savefig(f'{of}-cell_desity__{cell_qc_column}.png')
    
            
def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Performs automatic outlier detection over an anndata dataset,
            plotting the outlier cells.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )
    
    parser.add_argument(
        '-gt', '--gt_match_file',
        action='store',
        dest='gt_match_file',
        default=None,
        required=False,
        help='H5 AnnData file.'
    )    
    
    parser.add_argument(
        '-pt', '--patterns_exclude',
        action='store',
        dest='patterns_exclude',
        default=None,
        required=False,
        help='H5 AnnData file.'
    )    
    
    parser.add_argument(
        '--outliers_fraction',
        action='store',
        dest='outliers_fraction',
        default=0.0,
        type=float,
        help='Anticipated fraction of outlier cells. If 0.0, then runs sklearn \
            methods with "auto" as the anticipated number of outlier cells.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--metadata_columns',
        action='store',
        dest='metadata_columns',
        default='pct_counts_gene_group__mito_transcript,log1p_total_counts,log1p_n_genes_by_counts',
        help='Columns to use for outliers.\
            (default: %(default)s)'
    )
    
    parser.add_argument(
        '--MAD_thresholds',
        action='store',
        dest='MAD_thresholds',
        type=str,
        default='-5,-5,5',
        help='Columns to use for outliers.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--cell_qc_column',
        action='store',
        dest='cell_qc_column',
        default='cell_passes_qc',
        help='If column exists, cells are first filtered for those that \
            evaluate to true in this column, then this column is updated \
            based on the QC filtering. If this column does not exist then \
            it is added. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--method',
        action='store',
        dest='method',
        default='IsolationForest',
        help='Method for outlier detection. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--max_samples',
        action='store',
        dest='max_samples',
        default=0.1,
        type=float,
        help='The fraction of cells to draw from X to train each estimator. \
            Only valid if method == IsolationForest. \
            (default: %(default)s)'
    )
    # TODO: add option where user can say "run per each experiment id"

    parser.add_argument(
        '--cell_filtered_per_experiment_file',
        action='store',
        dest='cell_filtered_per_experiment',
        default='None',
        help='File detailing samples filtered per experiment. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='outliers',
        help='Basename of output files. \
            (default: %(default)s)'
    )
    parser.add_argument(
        '-fs', '--filter_strategy',
        action='store',
        dest='outlier_filtering_strategys',
        default='all_together',
        help='By default we filter the outliers for all the cells provided together, \
            but we may want to run it per celltype'
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

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)
    adata2 = sc.read_h5ad(filename='/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego_2025/results_cb3/handover/merged_h5ad/1.pre_QC_adata.h5ad')
    adata.obs['cell_id'] = adata.obs.index
    
    # Here we add an adaptive QC per Column
    # Drop out previous QCed cells
    
    # if cell_qc_column in adata.obs.columns:
    #     n_cells_original = adata.shape[0]
    #     adata = adata[adata.obs[cell_qc_column], :]
    #     print('Filtered out {} previously flagged cells using {}'.format(
    #         n_cells_original - adata.shape[0],
    #         cell_qc_column
    #     ))
    # else:
    #     adata.obs[cell_qc_column] = True
    # Get ballpark number of outliers
    outliers_fraction = options.outliers_fraction
    n_cells = adata.shape[0]
    n_outliers = outliers_fraction * n_cells
    if n_outliers == 0:
        outliers_fraction = 'auto'

    # Get a list of the data to use to calculate outliers
    metadata_columns = options.metadata_columns.split(',')
    # metadata_columns = [
    #     'pct_counts_gene_group__mito_transcript',
    #     #'total_counts',
    #     'log1p_total_counts',
    #     #'n_genes_by_counts',
    #     'log1p_n_genes_by_counts'
    # ]

    method = options.method
    
    
    methods = re.split(';|,',options.method)
    outlier_filtering_strategys = options.outlier_filtering_strategys
    # outlier_filtering_strategys = 'Azimuth:L0_predicted.celltype.l2;all_together;all_together::exclude'
    outlier_filtering_strategys = re.split(';|,',outlier_filtering_strategys)   
    metadata_columns_original = metadata_columns.copy() 
    
    for method in methods:
        print(f'-------Running {method} outlier filtering strategy-------')
        nmethod = method
        
        if method == 'LocalOutlierFactor':
            cell_qc_column = nmethod+':'+options.cell_qc_column
            # fit the model for outlier detection (default)
            clf = LocalOutlierFactor(
                #n_neighbors=100,
                # contamination=outliers_fraction
            )
        elif method == 'IsolationForest':
            max_samples = options.max_samples
            # if max_samples == 0.0:
            #     if n_cells < 1000:
            #         max_samples = 250
            #     else:
            #         max_samples = 0.1
            cell_qc_column = options.cell_qc_column
            print("Using max_samples of:\t{}".format(max_samples))
            clf = IsolationForest(
                #n_estimators=500,
                max_samples=max_samples,
                warm_start=False,
                contamination=outliers_fraction,
                random_state=0,
                bootstrap=True
            )
        elif method == 'EllipticEnvelope':
            cell_qc_column = nmethod+':'+options.cell_qc_column
            if outliers_fraction == 'auto':
                outliers_fraction = 0.1
            clf = EllipticEnvelope(
                contamination=outliers_fraction
            )
        elif method == 'OneClassSVM':
            cell_qc_column = nmethod+':'+options.cell_qc_column
            clf = OneClassSVM(
                # nu=n_outliers,
                # kernel="rbf",
                # gamma=0.1
            )
        elif method == 'MAD':
            mad_thresh = np.array([ int(x) for x in re.split(';|,',options.MAD_thresholds)])
            nmethod = method+options.MAD_thresholds
            cell_qc_column = nmethod+':'+options.cell_qc_column
            clf = onesidemad(mad_thresh)
        else:
            raise ValueError('ERROR: invalid method.')
        print(cell_qc_column)
        metadata_columns=metadata_columns_original.copy()
        metadata_columns.append(cell_qc_column)
        # We perform the adaptive qc either for all data together or per user defined column of unique values.

        matching = [s for s in outlier_filtering_strategys if "::" in s]
        for m1 in matching:
            # duplicate the strategy if we have flagged it as an extra step by adding ::no_exclude flag
            m2 = m1.split('::')[0]
            if m1 =='all_together::exclude':
                adata.obs[m1] = 'all_cells'
            else:
                adata.obs[m1] = m2
            
        # We load the GT match file and determine the celline match that needs to be subjected to adaptive qc independently
        if(options.patterns_exclude):
            ######## The folowing bit of code takes the GT match outputs and utilises this to run adaptive qc on cellines independently - only if the celline is expected.
            if(options.gt_match_file!='fake_file.fq'):
                GT_match_file = pd.read_csv(options.gt_match_file,sep='\t')
                GT_patterns = options.patterns_exclude.split(';')
                for outlier_filtering_strategy in outlier_filtering_strategys:
                    strategy = outlier_filtering_strategy
                    if strategy!='all_together':
                        for pattern in GT_patterns:
                            Matches = GT_match_file[GT_match_file['donor_gt'].str.contains(pattern)]
                            if len(Matches)>0:
                                for i,match in Matches.iterrows():
                                    donor=match['donor_query']
                                    pool=match['pool']
                                    expected = match['Match Expected']
                                    if expected:
                                        filter_query = f"convoluted_samplename=='{pool}' and Donor=='{donor}'"
                                        idx1 = adata.obs.query(filter_query).index
                                        try:
                                            adata.obs[strategy] = adata.obs[strategy].cat.add_categories(pattern)
                                        except:
                                            _='category exists already'
                                        adata.obs.loc[idx1,strategy]=pattern
            ########                 
        # all_index = pd.DataFrame(adata.obs.index,columns=['col'])
        # all_together = all_indexes.str[0]+'-'+all_indexes.str[1]+'-'+all_indexes.str[2]
        
        for outlier_filtering_strategy in outlier_filtering_strategys:
            metadata_columns = metadata_columns_original.copy()
            if method == 'IsolationForest' and outlier_filtering_strategy == 'all_together':
                if outlier_filtering_strategy == 'all_together':
                    cell_qc_column = options.cell_qc_column
                    cell_qc_column_score = options.cell_qc_column+':score'
                else:
                    cell_qc_column = f'{options.cell_qc_column}-{nmethod}-per:{outlier_filtering_strategy}'
                    cell_qc_column_score = f'{options.cell_qc_column}-{nmethod}-per:{outlier_filtering_strategy}:score'                      
            else:
                cell_qc_column = f'{options.cell_qc_column}-{nmethod}-per:{outlier_filtering_strategy}'
                cell_qc_column_score = f'{options.cell_qc_column}-{nmethod}-per:{outlier_filtering_strategy}:score' 
            metadata_columns.append(cell_qc_column)               
            if (outlier_filtering_strategy == 'all_together'):

                adata.obs[cell_qc_column] = True
                adata.obs[cell_qc_column_score] = None
                prediction_score, fail_pass = perform_adaptiveQC_Filtering(clf,adata,method,metadata_columns)
                adata.obs[cell_qc_column] = fail_pass
                adata.obs[cell_qc_column_score] = prediction_score
            else:
                # cell_qc_column = f'{cell_qc_column}-per:{outlier_filtering_strategy}'
                # cell_qc_column_score = f"{cell_qc_column}:score"
                

                try:
                    os.mkdir(f'per_celltype_outliers__{outlier_filtering_strategy}')
                except:
                    _='exists already'
                try:
                    outlier_strategy_cols = set(adata.obs[outlier_filtering_strategy])
                except:
                    print(f'{outlier_filtering_strategy} - user provided col doesnt exist')
                    continue
                
                adata.obs[cell_qc_column] = True
                adata.obs[cell_qc_column_score] = ''
                
                for subset_id_for_ad_qc in outlier_strategy_cols:
                    subset_ad = adata[adata.obs[outlier_filtering_strategy]==subset_id_for_ad_qc]
                    print(f'filtering:{subset_id_for_ad_qc} {outlier_filtering_strategy}')    
                    if(len(subset_ad)>100):
                        # We only perform adaptive qc when there is at least 100 cells, otherwise we assume that all pass
                        prediction_score, fail_pass = perform_adaptiveQC_Filtering(clf,subset_ad,method,metadata_columns)
                        adata.obs.loc[subset_ad.obs.index,cell_qc_column]=fail_pass
                        subset_ad.obs.loc[subset_ad.obs.index,cell_qc_column]=fail_pass
                        
                        adata.obs.loc[subset_ad.obs.index,cell_qc_column_score]=prediction_score
                        subset_ad.obs.loc[subset_ad.obs.index,cell_qc_column_score]=prediction_score
                        
                    else:
                        print(f'For a category {subset_id_for_ad_qc} we have only {len(subset_ad)} cells and as its not sufficient ammount to estimate distributions we assuma all pass QC')
                    # subset_ad.uns['cell_outlier_estimator'] = method
                    of = f'per_celltype_outliers__{outlier_filtering_strategy}/{subset_id_for_ad_qc}---{options.of}'
                    generate_plots(subset_ad,cell_qc_column,metadata_columns,metadata_columns_original,of)
                    del subset_ad
            adata.uns['cell_outlier_estimator'] = method  
            adata.obs[['cell_id', cell_qc_column,cell_qc_column_score]].to_csv(
                f'{options.of}-outliers_filtered__{cell_qc_column}.tsv',
                sep='\t',
                index=False,
                header=True
            )

            # Calculate cell_filtered_per_experiment
            filter_columns = [
                'experiment_id',
                'filter_type',
                'n_cells_left_in_adata'
            ]
            if options.cell_filtered_per_experiment == 'None':
                df_cell_filt_per_exp = adata.obs['experiment_id'].value_counts()
                df_cell_filt_per_exp = df_cell_filt_per_exp.rename_axis(
                    'experiment_id'
                ).reset_index(name='n_cells_left_in_adata')
                df_cell_filt_per_exp['filter_type'] = 'before_filters'
                df_cell_filt_per_exp = df_cell_filt_per_exp[filter_columns]
            else:
                # Load the samples filtered per experiment file:
                df_cell_filt_per_exp = pd.read_csv(
                    options.cell_filtered_per_experiment,
                    sep="\t"
                )
                filt = df_cell_filt_per_exp['filter_type'] != 'after_filters'
                df_cell_filt_per_exp = df_cell_filt_per_exp.loc[filt, :]
                
            # Now calculate the n cells left after all filters
            adata_after_filters = adata.obs.loc[adata.obs[cell_qc_column], :]
            df_cells_filtered = adata_after_filters['experiment_id'].value_counts()
            df_cells_filtered = df_cells_filtered.rename_axis(
                'experiment_id'
            ).reset_index(name='n_cells_left_in_adata')
            df_cells_filtered['filter_type'] = '{} {} outlier_{}'.format(
                'filter__all_samples',
                'after_outlier_filter',
                method
            )
            df_cells_filtered = df_cells_filtered[filter_columns]
            df_cell_filt_per_exp = df_cell_filt_per_exp.append(
                df_cells_filtered,
                ignore_index=True
            )
            df_cells_filtered['filter_type'] = 'after_filters'
            df_cell_filt_per_exp = df_cell_filt_per_exp.append(
                df_cells_filtered,
                ignore_index=True
            )
            # Save the final dataframe
            df_cell_filt_per_exp.to_csv(
                f'{options.of}-cell_filtered_per_experiment__{cell_qc_column}.tsv',
                sep='\t',
                index=False,
                header=True
            )
            
            generate_plots(adata,cell_qc_column,metadata_columns,metadata_columns_original,options.of)

    
    # adata.write(
    #     '{}.h5ad'.format(options.of),
    #     compression='gzip',
    #     compression_opts=options.anndata_compression_opts
    # )
    
    

    
    # Write adata matrix for each of the 
    # adata.obs['convoluted_samplename']

if __name__ == '__main__':
    main()
