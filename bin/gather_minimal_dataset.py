#!/usr/bin/env python3

__date__ = '2022-21-03'
__version__ = '0.0.2'
# DEBUG = True

DEFAULT_SAMPLE_SIZE_PERC = 10
SAMPLE_SIZE_MIN = 20
import logging
import sys
import argparse
import random
import h5py
import pandas
import anndata
import re
import scanpy as sc
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import glob
import pandas as pd
import numpy as np
import statistics
import h5py
import scanpy
from datetime import date
write_h5=False
today = date.today()
date_now = today.strftime("%Y-%m-%d")
date_of_transfer = today.strftime("%Y-%m-")
month = today.strftime('%m')
day='30'

if (month=='03'):
    day='31'
elif(month=='01'):
    day='30'
elif(month=='02'):
    day='29'
date_of_transfer+=day

ANNDATA_FILE_QC = "adata.h5ad"
DATA_DIR_AZIMUTH = "azimuth"
AZIMUTH_ASSIGNMENTS_FNSUFFIX = '_predicted_celltype_l2.tsv'
SCRUBLET_ASSIGNMENTS_FNSUFFIX = 'scrublet.tsv'

COLUMNS_AZIMUTH = {
    
    'Azimuth:L0_predicted.celltype.l2': 'l2.mapped.azimuth.celltyp.l0',
    'Azimuth:L1_predicted.celltype.l2': 'l2.mapped.azimuth.celltyp.l1',
    'Azimuth:mapping.score.celltype.l1': 'azimuth.pred.score.l1',
    'Azimuth:mapping.score.celltype.l2': 'azimuth.pred.score.l2',
    'Azimuth:mapping.score.celltype.l3': 'azimuth.pred.score.l3',
    'Azimuth:predicted.celltype.l1': 'azimuth.celltyp.l1',
    'Azimuth:predicted.celltype.l1.score': 'azimuth.pred.score.l1',
    'Azimuth:predicted.celltype.l2': 'azimuth.celltyp.l2',
    'Azimuth:predicted.celltype.l2.score': 'azimuth.pred.score.l2',
    'Azimuth:predicted.celltype.l3': 'azimuth.celltyp.l3',
    'Azimuth:predicted.celltype.l3.score': 'azimuth.pred.score.l3',
    'Azimuth:predicted.celltype.l2': 'azimuth.celltyp.l2',
    'Azimuth:predicted.celltype.l3': 'azimuth.celltyp.l3',
    'Azimuth:L1_predicted.celltype.l2': 'azimuth.celltyp.l1',
    'Azimuth:predicted.celltype.l2.score': 'azimuth.pred.score.l2',
    'Azimuth:mapping.score.celltype.l2': 'azimuth.map.score',
    'Celltypist*':'Celltypist*'
    }

COLUMNS_DECONV = {
    'donor_id': 'vireo.donor.id',
    'prob_max': 'vireo.prob.max',
    'prob_doublet': 'vireo.prob.doublet'
    }
COLUMNS_QC = {
    'cell_passes_hard_filters': 'cell_passes_hard_filters', 'State':'State','Donor id':'Donor id','Vacutainer ID':'Vacutainer ID', 'live_cell_count':'live_cell_count', 'viability':'viability', 'SITE':'SITE', 'GEM_BATCH':'GEM_BATCH', 'Gender':'Gender', 'DRAW_DATE':'DRAW_DATE', 'customer_measured_volume':'customer_measured_volume',
    'cell_passes_qc': 'qc.filter.pass',"id_pool_lims":'lims.pool.id','chromium_channel':'chromium.run.id',
    'cell_passes_qc:score':'qc.filter.pass:score',
    'cell_passes_qc-per:all_together::exclude':'qc.filter.pass.spikein_exclude',
    'cell_passes_qc-per:all_together::exclude:score':'qc.filter.pass.spikein_exclude:score',
    'cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2':'qc.filter.pass.AZ:L0',
    'cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2:score':'qc.filter.pass.AZ:L0:score',
    'total_counts': 'qc.umi.count.total',
    'total_counts_gene_group__mito_transcript': 'qc.umi.count.mt',
    'pct_counts_gene_group__mito_transcript': 'qc.umi.perc.mt',
    'pct_counts_in_top_500_genes':'pct_counts_in_top_500_genes',
    'n_genes_by_counts': 'qc.genes.detected.count',
    'Azimuth:L0_Azimuth:predicted.celltype.l2':'azimuth.celltyp.l0.mapped_fromL2',
    'Azimuth:L1_Azimuth:predicted.celltype.l2':'azimuth.celltyp.l1.mapped_fromL2',
    'total_counts_gene_group__mito_transcript':'total_counts_gene_group__mito_transcript',
    'pct_counts_gene_group__mito_transcript':'pct_counts_gene_group__mito_transcript',
    'total_counts_gene_group__mito_protein':'total_counts_gene_group__mito_protein',
    'pct_counts_gene_group__mito_protein':'pct_counts_gene_group__mito_protein',
    'Celltypist:Immune_All_High:predicted_labels':'celltypist.celltyp.Immune_All_High',
    'Celltypist:Immune_All_Low:predicted_labels':'celltypist.celltyp.Immune_All_Low',
    'cellbender_latent_probability':'cellbender.latent.probability',
    'scrublet__multiplet_scores':'scrublet__multiplet_scores',
    'scrublet__predicted_multiplet':'scrublet__predicted_multiplet',
    'scrublet__multiplet_zscores':'scrublet__multiplet_zscores',
    'log10_ngenes_by_count':'log10_ngenes_by_count',
    'pct_counts_in_top_50_genes':'pct_counts_in_top_50_genes',
    'pct_counts_in_top_100_genes':'pct_counts_in_top_100_genes',
    'pct_counts_in_top_200_genes':'pct_counts_in_top_200_genes',
    'pct_counts_in_top_500_genes':'pct_counts_in_top_500_genes',
    'S_score':'S_score', 'G2M_score':'G2M_score',
    'phase':'phase'
    }
COLUMNS_CELLBENDER = {'cellbender_latent_probability': 'cellbender.latent.probability'}
COLUMNS_DATASET = {
    'experiment_id': 'experiment.id',
    'tranche.id':'tranche.id',
    'chromium_run_id': 'chromium.run.id',
    'chromium_lane': 'chromium.lane',
    'instrument':'instrument'
    }

COLUMNS_SCRUBLET = {
    'scrublet__multiplet_scores': 'scrublet.scores',
    'scrublet__predicted_multiplet': 'scrublet.multiplet',
    'scrublet__multiplet_zscores': 'scrublet.zscores',
    'scds_DropletType':'scds.multiplet','scds_score':'scds.score',
    'scDblFinder_DropletType':'scDblFinder.multiplet','scDblFinder_Score':'scDblFinder.score',
    'DoubletDecon_DropletType':'DoubletDecon.multiplet',
    'DoubletFinder_DropletType':'DoubletFinder.multiplet','DoubletFinder_score':'DoubletFinder.score'
    }

COLUMNS_OUTPUT = \
    {**COLUMNS_DATASET, **COLUMNS_CELLBENDER, **COLUMNS_DECONV, **COLUMNS_QC, **COLUMNS_AZIMUTH}
COLUMNS_OUTPUT_WITH_SCRUBLET = \
    {**COLUMNS_DATASET, **COLUMNS_CELLBENDER, **COLUMNS_DECONV, **COLUMNS_SCRUBLET, **COLUMNS_QC, **COLUMNS_AZIMUTH}

def get_df_from_mangled_index(df, expid):
    idx = df.index.str.split(pat='-{}__'.format(expid))
    try:
        xf = pandas.DataFrame.from_records(idx, columns = ('barcode', 'donor'), index = df.index)
        if xf.shape[0] != df.shape[0]:
            sys.exit("ERROR: when untangling mangled index.")
    except:
        df['barcode_idx']= df.index+'-'+expid+'__'+df['donor_id'].astype(str)
        df2 = df.reset_index() 
        df2 = df2.rename({'index':'barcode','donor_id':'donor'},axis=1)
        xf = df2[['barcode','donor']]
        xf.index=df2['barcode']
    return xf

def load_deconv_file_table(fnam):
    df = pandas.read_table(fnam)
    df = df.rename(columns={'experiment_id':'experiment_mangled_id'})
    ix = df.experiment_mangled_id.str.split(pat='__')
    jx = df.h5ad_filepath.transform(lambda a: os.path.basename(a))
    xf = pandas.DataFrame.from_records(ix, columns = ('exp_id', 'donor_id'), index = df.index)
    df.insert(1, 'experiment_id', xf['exp_id'])
    df.insert(2, 'donor_id', xf['donor_id'])
    df.insert(3, 'file_name_h5ad', jx)
    return df

def load_raw_file_table(fnam):
    df = pandas.read_table(fnam, index_col = 'experiment_id')
    df = df[['data_path_raw_h5']]
    x = df['data_path_raw_h5'].transform(lambda a: os.path.basename(a))
    df.insert(1, 'file_name_raw_h5', x)
    return df


def anndata_from_h5(
    file: str,
    analyzed_barcodes_only: bool = True
 ) -> 'anndata.AnnData':
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
            'gene_symbols': gene_symbols
        }
    )
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


def load_scrublet_assignments(expid, datadir_scrublet):
    filpath = None
    fnam = '{}{}'.format(expid, SCRUBLET_ASSIGNMENTS_FNSUFFIX)
    fnam2 = '{}{}'.format(expid, SCRUBLET_ASSIGNMENTS_FNSUFFIX.replace('-',''))
    fnam3 = '{}{}'.format(expid, 'scrublet.tsv.gz')
    for fn in os.listdir(datadir_scrublet):
        if fn == fnam or fn == fnam2 or fn == fnam3:
            filpath = os.path.join(datadir_scrublet, fn)
            break
 
    if not filpath:
        sys.exit("ERROR: could not find filename '{}' in direcotry {}\n"
            .format(fnam, datadir_scrublet))
    sys.stderr.write("loading scrublet annotation from file {} ...\n".format(filpath))
    scb = pandas.read_table(filpath).set_index('cell_barcode', drop = True)
    return scb

def fetch_cellbender_annotation(dirpath, expid,Resolution):
    
    try:
        h5_path = f"{args.results_dir}/{os.path.dirname(dirpath)}/cellbender_FPR_{Resolution}_filtered.h5"
        f = h5py.File(h5_path, 'r')
    except:
        try:
            h5_path = glob.glob(f"{os.path.dirname(dirpath)}/cellbender*{Resolution}_filtered.h5")[0]
        except:
            Resolution=Resolution.replace('pt','.')
            h5_path = glob.glob(f"{os.path.dirname(dirpath)}/cellbender*{Resolution}_filtered.h5")[0]
        f = h5py.File(h5_path, 'r')
    # ad = scanpy.read_10x_h5(h5_path, genome='background_removed')
    # interesting data is in /matrix/barcodes and matrix/latent_cell_probability
   
    df = pandas.DataFrame({
        "barcodes":f['/matrix/barcodes'],
        "cellbender_latent_probability":f['/matrix/latent_cell_probability']
        })
    bc = df['barcodes'].transform(lambda a: a.decode("ascii"))
    df["barcodes"] = bc
    f.close()
    return df.set_index("barcodes", drop = True)

def get_lane_and_runid_from_experiment_id(df, insert_pos = 0):
    x = df['experiment_id'].transform(lambda a: a.split('_lane_'))
    xf = pandas.DataFrame.from_records(x, columns = ('chromium_run_id', 'chromium_lane'), index = df.index)
    df.insert(insert_pos, 'chromium_run_id', xf['chromium_run_id'])
    df.insert(insert_pos + 1, 'chromium_lane', xf['chromium_lane'])
    return df

def write_subsample(adat, oufnam, sample_size_perc = DEFAULT_SAMPLE_SIZE_PERC):
    n_cells = adat.shape[0]
    n_cells_subsample = int(n_cells * sample_size_perc / 100);
    if n_cells_subsample < SAMPLE_SIZE_MIN:
        n_cells_subsample = SAMPLE_SIZE_MIN
    if n_cells_subsample >= n_cells:
        sys.stderr.write("# write_subsample: unchanged ...\n")
        ad = adat
    else:
        sys.stderr.write("# write_subsample sample {:d} cells ({:d}% of {:d})\n".
            format(n_cells_subsample, sample_size_perc, n_cells))
        sv = random.sample(range(n_cells), n_cells_subsample)
        sv.sort()
        ad = adat[sv,:]

    sys.stderr.write("# writing subsample file {} ...\n".format(oufnam))
    ad.obs.index.name = 'barcode'
    if write_h5:
        ad.write(oufnam + '.h5ad')
    ad.obs.to_csv(oufnam + '.tsv', sep = "\t", na_rep = "N/A")
    return

def gather_donor(donor_id, ad, ad_lane_raw, qc_obs, columns_output = COLUMNS_OUTPUT,outdir = os.curdir,oufh = sys.stdout,lane_id=1):
    
    oufnam = "{}.{}".format(expid, donor_id)
    sys.stderr.write("processing {} {} ...\n".format(expid, donor_id))
    if write_h5:
        oufh.write("{}\t{}\t{}.h5ad\t{}.tsv\n".format(expid, donor_id, oufnam, oufnam))

    # loading deconvoluted dataset
    try:
        ad=ad.to_memory()
    except:
        _='already loaded'
    ad.var.index.name = "ensembl_id"
    idx_all = list(set(ad.obs.index).intersection(set(ad_lane_raw.obs.index)))
    ad = ad[idx_all]
    ad.raw = ad_lane_raw[idx_all]
    if donor_id != "unassigned" and donor_id != "doublet":

        dt = qc_obs[qc_obs.donor == donor_id]
        try:
            dt = dt.set_index(dt['barcode'])
        except:
            _='barcode is already the index'
        if len(set(ad.obs.index).intersection(dt.index))<1:
            dt.index = dt.index+'-'+dt['experiment_id'].astype(str)+'__'+dt['donor'].astype(str)
        df = pandas.concat([ad.obs, dt], axis = 1, join = 'inner')
        df =df.loc[:,~df.columns.duplicated(keep='last')]
        df['tranche.id']=args.experiment_name
        colnams = list(columns_output.keys())
        colnams_overlap = sorted(set(colnams).intersection(set(df.columns)))
        ad.obs = df[colnams_overlap].rename(columns = columns_output)
        dt=df[colnams_overlap].rename(columns = columns_output)
        # Stats
        print('Performing the stats analysis')
        experiment_id = args.experiment_name
        pool_id = list(set(df.experiment_id))[0]
        try:
            chromium_channel_number = list(set(df.chromium_channel))[0]
        except:
            chromium_channel_number=' '
        donor_id = list(set(df.donor_id))[0]
    else:
        dt = ad.obs
        experiment_id=' ';pool_id=' ';chromium_channel_number=' ';donor_id=' '

    dt.index.name = 'barcode'
    ad.obs.index.name = 'barcode'
    ad.obs = ad.obs.loc[:,~ad.obs.columns.duplicated()]
    dt = dt.loc[:,~dt.columns.duplicated()].copy()
    dt[set(dt.columns)].to_csv(os.path.join(outdir, oufnam + '.tsv'), sep = "\t", na_rep = "N/A")
    sys.stderr.write("writing file {} ...\n".format(oufnam))
    
    if write_h5:
        path1=os.path.join(outdir, oufnam + '.h5ad')
        try:
            ad.obs['qc.filter.pass.AZ:L0'] = ad.obs['qc.filter.pass.AZ:L0'].astype('bool')
        except:
            pass
        ad.obs['cell_passes_hard_filters'] = ad.obs['cell_passes_hard_filters'].astype('bool')
        ad.obs['qc.filter.pass'] = ad.obs['qc.filter.pass'].astype('bool')
        ad.write(path1,compression='gzip')

    return {
        'Experiment ID':experiment_id,
        'Pool ID':pool_id,
        'Chromium channel number':chromium_channel_number,
        'Donor id':donor_id
    }

def gather_pool(expid, args, df_raw, df_cellbender, adqc, oufh = sys.stdout,lane_id=1,Resolution='0pt5'):
    
    # Get the merged in metadata
    outdir = f'{args.outdir}/{expid}'
    try:
        os.mkdir(outdir)
    except:
        print('dir exists')

    ######################
    #Cellranger datasets
    ######################
    columns_output = {**COLUMNS_DATASET, **COLUMNS_DECONV, **COLUMNS_QC}
    #Reading unfiltered raw cellranger datasets
    try:
        adata_cellranger_raw = scanpy.read_10x_mtx(f"{df_raw.loc[expid, 'data_path_10x_format']}/raw_feature_bc_matrix")
        try:
            if write_h5:
                os.symlink(f"{df_raw.loc[expid, 'data_path_10x_format']}/raw_feature_bc_matrix.h5", f"./{outdir}/Cellranger_raw_feature_bc_matrix__{expid}.h5")
        except:
            print('File already linked')
    except:
        adata_cellranger_raw = scanpy.read_10x_mtx(f"{df_raw.loc[expid, 'data_path_10x_format']}/raw_feature_bc_matrix")
        try:
            if write_h5:
                os.symlink(f"{df_raw.loc[expid, 'data_path_10x_format']}/raw_feature_bc_matrix.h5", f"./{outdir}/Cellranger_raw_feature_bc_matrix__{expid}.h5")
        except:
            print('File already linked')

     
    # Reading filtered cellranger files

    adata_cellranger_filtered = scanpy.read_10x_mtx(f"{df_raw.loc[expid, 'data_path_10x_format']}/filtered_feature_bc_matrix")

    try:
        if write_h5:
            os.symlink(f"{df_raw.loc[expid, 'data_path_10x_format']}/filtered_feature_bc_matrix.h5",f"{outdir}/Cellranger_filtered_feature_bc_matrix__{expid}.h5")  
    except:
        print('File already linked')

               
    zero_count_cells_cellranger_raw = adata_cellranger_raw.obs_names[np.where(adata_cellranger_raw.X.sum(axis=1) == 0)[0]]
    ad_lane_raw = adata_cellranger_raw[adata_cellranger_raw.obs_names.difference(zero_count_cells_cellranger_raw, sort=False)]
    scanpy.pp.calculate_qc_metrics(adata_cellranger_raw, inplace=True)
    df_total_counts = pd.DataFrame(data= adata_cellranger_raw.obs.sort_values(by=['total_counts'], ascending=False).total_counts)
    df_total_counts['barcode_row_number'] = df_total_counts.reset_index().index + 1 
    df_total_counts['barcodes'] = df_total_counts.index
    df_total_counts_cellranger_raw = df_total_counts
    df_total_counts_cellranger_raw['dataset']='Cellranger Raw'

    if df_cellbender is not None and (len(df_cellbender)!=0):
        df_cellbender = df_cellbender.reset_index()
        df_cellbender = df_cellbender.drop_duplicates(subset=['experiment_id'])
        df_cellbender = df_cellbender.set_index('experiment_id')
        df2 = glob.glob(f'{args.results_dir}/preprocessing/data_modalities_split/preprocess/{expid}/Gene_Expression-{expid}.h5ad')[0]

        f=df_cellbender.loc[expid, 'data_path_10x_format']
        if (type(f) == str):
            f=[f]
        for id in f:
            print(id)
            cell_bender_path = id
        try:
            ad_lane_filtered = anndata.read_h5ad(df2)
            if write_h5:
                try:
                    path2 = os.path.realpath(df2)
                    os.symlink(path2, f"./{outdir}/Cellbender_filtered_{Resolution}__{expid}.h5ad")
                except:
                    print('File already linked')
        except:
            ad_lane_filtered = scanpy.read_10x_mtx(cell_bender_path)  
            if write_h5:
                ad_lane_filtered.write(
                f"./{outdir}/Cellbender_filtered_{Resolution}__{expid}.h5ad" ,
                compression='gzip'
                )      

        columns_output = {**columns_output, **COLUMNS_CELLBENDER}
    else:
        ad_lane_filtered = adata_cellranger_filtered
        df_cellbender=None
        cell_bender_path=None

    # Removing Zero count cells from the matices
    zero_count_cells = ad_lane_filtered.obs_names[np.where(ad_lane_filtered.X.sum(axis=1) == 0)[0]]
    ad_lane_filtered = ad_lane_filtered[ad_lane_filtered.obs_names.difference(zero_count_cells, sort=False)]
    zero_count_cells = adata_cellranger_filtered.obs_names[np.where(adata_cellranger_filtered.X.sum(axis=1) == 0)[0]]
    adata_cellranger_filtered = adata_cellranger_filtered[adata_cellranger_filtered.obs_names.difference(zero_count_cells, sort=False)]

    #############
    #Cellranger Metrics Datasheet
    #############

    metrics = pd.read_csv(df_raw.loc[expid, 'data_path_10x_format']+'/metrics_summary.csv')

    
    try:
        Donor_Cohort_Assignments = pd.read_csv(f'{args.results_dir}/gtmatch/{expid}/{expid}_gt_donor_assignments.csv')
        All_assignments = pd.read_csv(f'{args.results_dir}/gtmatch/assignments_all_pools.tsv',sep='\t')
        All_assignments = All_assignments.set_index(All_assignments['pool']+'__'+All_assignments['donor_gt'].astype(str).str.replace('^0*', '', regex=True).str.replace('.*THP1.*', 'THP1', regex=True).str.replace('.*U937.*', 'U937', regex=True))
        All_assignments['tp2'] = All_assignments.index
        all_dubs =  list(set(All_assignments[All_assignments['tp2'].duplicated()].index))
        all_dubs2 =  All_assignments.loc[all_dubs]
        All_assignments = All_assignments.drop_duplicates(subset=['tp2'],keep=False)
    except:
        Donor_Cohort_Assignments = pd.DataFrame(columns=['panel'])
        All_assignments= pd.DataFrame(columns=['panel'])

    if (args.extra_meta):
        Metadata = pd.read_csv(f"{args.extra_meta}",sep='\t')
        Metadata = Metadata.set_index(Metadata['Pool ID']+'__'+Metadata['donor'].astype(str).str.replace('^0*', '', regex=True).str.replace('.*THP1.*', 'THP1', regex=True).str.replace('.*U937.*', 'U937', regex=True))
        Metadata['donor2']=All_assignments['donor_query']
        Metadata = Metadata.set_index(Metadata['Pool ID']+'__'+Metadata['donor2'])
        
        if len(all_dubs2)>0:
            for i,n in all_dubs2.iterrows():
                print(i)
                new =Metadata[Metadata['experiment_id']==i]
                new['donor2']=n['donor_query']
                Metadata = Metadata.append(new, ignore_index=True)
            Metadata = Metadata.drop_duplicates()
        Metadata = Metadata.set_index(Metadata['Pool ID']+'__'+Metadata['donor2'])
        Metadata['tp2'] = Metadata.index
        Metadata = Metadata.drop_duplicates(subset=['tp2'],keep=False)
        Metadata = Metadata.drop('tp2',axis=1)
        Metadata = Metadata.drop('donor2',axis=1)
        Metadata = Metadata.drop('donor',axis=1)
        Metadata = Metadata.drop('experiment_id',axis=1)
    else:
        Metadata = pd.DataFrame()

    #############
    #Cell-type assignments
    #############

    azt = pd.read_csv(f'{args.results_dir}/celltype_assignemt/All_Celltype_Assignments.tsv',sep='\t',index_col=0)
    azt_cols_to_add = azt.columns[azt.columns.str.contains('Azimuth')]
    ct_cols_to_add = azt.columns[azt.columns.str.contains('Celltypist')]
    sc_cols_to_add = azt.columns[azt.columns.str.contains('scpred_prediction')]
    for i3 in set(azt_cols_to_add) - set(columns_output.keys()):
        columns_output = {**columns_output,  **{i3:i3}}
    for i3 in set(sc_cols_to_add) - set(columns_output.keys()):
        columns_output = {**columns_output,  **{i3:i3}}
    for i3 in set(ct_cols_to_add) - set(columns_output.keys()):
        columns_output = {**columns_output,  **{i3:i3}}
        
        
    cols = pd.DataFrame(adqc.obs.columns)
    cols =cols[cols[0].str.contains('cell_passes_qc')]
    for i3 in set(cols[0]) - set(columns_output.keys()):
        columns_output = {**columns_output,  **{i3:i3}}   
    # scpred_to_add = azt.columns[azt.columns.str.contains('Scpred')]
    ##########################
    # Scrublet
    #########################
    doublet_data = glob.glob(f'{args.results_dir}/doublet_detection/doublet_results_combined/*.tsv')
    doublet_data_combined = pd.DataFrame()
    for f1 in doublet_data:
        print(f1)
        pool_name = f1.split('__')[0].split('/')[-1]
        d2 = pd.read_csv(f1,sep='\t')
        d2['Exp']=pool_name
        doublet_data_combined = pd.concat([doublet_data_combined,d2])
    doublet_data_combined = doublet_data_combined.drop_duplicates(subset='barcodes')
    scb = doublet_data_combined.set_index('barcodes')
    

    columns_output = {**columns_output,  **COLUMNS_SCRUBLET}    
    columns_output = {k: v for k, v in columns_output.items() if 'idx1' not in k and 'idx' not in v}

    ############################################################
    # Loading deconvoluted data including unassigned and doublets
    ###########################################

    s = adqc.obs['convoluted_samplename'] == expid
    ad = adqc[s]
    df = get_df_from_mangled_index(ad.obs, expid)
    #df.insert(0, 'barcode', df.index.values)
    df_pre = pandas.concat([df,ad.obs], axis = 1)
    df_pre['mengled_index'] = df_pre.index
    df = df_pre.set_index("barcode", drop = True)

    if cell_bender_path is not None:
        # cellbender removes the barcodes - 
        dfcb = fetch_cellbender_annotation(cell_bender_path, expid,Resolution)
        dc = pandas.concat([df, dfcb], axis = 1, join = 'inner')
        if dc.shape[0] != df.shape[0]:
            sys.exit("ERROR: barcodes missing in cellbender file.")
        df = dc.copy()
        dc=dc.set_index('mengled_index')
        
    obsqc = df
    all_QC_lane = ad

    donor_tables=pd.DataFrame([])
    for d1 in set(obsqc['donor']):
        donor_table={}
        donor_table['experiment_id']=expid
        donor_table['donor_id']=d1
        donor_table['file_path_h5ad']='all_QC_lane'
        df_donors=pd.DataFrame([donor_table])
        try:
            donor_tables=pd.concat([donor_tables,df_donors])
        except:
            donor_tables=donor_table
    df_donors = donor_tables
    df_donors = df_donors.reset_index(drop=True)
    if scb is not None:
        obsqc = pandas.concat([obsqc,scb], axis = 1, join = 'inner')

    if len(Metadata)>0:
        obsqc['barcode']= obsqc.index
        obsqc = obsqc.set_index(obsqc['convoluted_samplename'].astype(str)+'__'+obsqc['donor_id'].astype(str))
        obsqc['chromium_channel']=Metadata['chromium_channel']
        obsqc.update(Metadata)
        for col in Metadata.columns:
            if col not in obsqc.columns:
                obsqc[col] = Metadata[col]
    if len(All_assignments)>0:
        All_assignments = All_assignments.set_index(All_assignments['pool']+'__'+All_assignments['donor_query'])
        obsqc['Donor id']=All_assignments['donor_gt original']
        obsqc['Vacutainer ID']=All_assignments['donor_gt']
    fctr = 0

    #####################
    #Performing Calculations and gathering data
    #####################

    try:
        
        Count_of_UKB = int(all_QC_lane.obs['nr_ukbb_samples'].astype(str)[0])
        Count_of_ELGH = int(all_QC_lane.obs['nr_elgh_samples'].astype(str)[0])
        Count_of_Spikeins = int(all_QC_lane.obs['nr_spikeins'].astype(str)[0])

    except:
        Count_of_UKB = 0    
        Count_of_ELGH = 0    
        Count_of_Spikeins = 0

    try:
        date_of_sequencing = all_QC_lane.obs['last_updated'][0]
    except:
        date_of_sequencing = 'Sequencing date not vailable'

    try:
        try:
            Machine_id = all_QC_lane.obs['instrument'][0]
        except:
            Machine_id = all_QC_lane.obs['instrument_name'][0]
    except:
        Machine_id = 'Machine_id not vailable'

    try:
        Run_ID = str(all_QC_lane.obs['id_run'][0]) 
    except:
        Run_ID = 'Run_ID not vailable'
    try:
        chromium_channel = str(all_QC_lane.obs['chromium_channel'][0]) 
    except:
        chromium_channel = 'Run_ID not vailable'
        
    Donors = list(df_donors.donor_id)

    try:
        Donors.remove('doublet')
    except:
        _ = 'no doublets detected'
    try:
        Donors.remove('unassigned')
    except:
        _ = 'no unassigned detected'

    Donors_in_pool = len(Donors)
    try:
        Number_of_Reads = int(metrics['Number of Reads'].values[0].replace(',','')) # NOT SURE HOW THIS IS CALCULATED.
    except:
        try:
             Number_of_Reads = int(metrics['Number of reads'].values[0].replace(',','')) 
        except:
            Number_of_Reads = None
    try:
        Fraction_Reads_in_Cells = metrics['Fraction Reads in Cells'].values[0]
    except:
        try:
            Fraction_Reads_in_Cells = metrics['Sequencing saturation'].values[0]
        except:
            Fraction_Reads_in_Cells = None
    try:
        Mean_reads_per_cell = int(metrics['Mean Reads per Cell'].values[0].replace(',',''))
    except:
        try:
            Mean_reads_per_cell = int(metrics['Mean reads per cell'].values[0].replace(',',''))
        except:
            Mean_reads_per_cell= None
    
    f = pd.DataFrame(adata_cellranger_filtered.X.sum(axis=1))
    Median_UMI_Counts_per_cellranger= statistics.median(f[f>0][0])
    Mean_Reads = statistics.mean(f[f>0][0])

    f = pd.DataFrame(ad_lane_raw.X.sum(axis=1))
    Median_UMI_Counts_per_before_filter= statistics.median(f[f>0][0])

    f = pd.DataFrame(all_QC_lane.X.sum(axis=1))
    Median_UMI_Counts_per_Cell_after_all_filter= statistics.median(f[f>0][0])

    f = pd.DataFrame(ad_lane_filtered.X.sum(axis=1))
    Median_UMI_Counts_per_Cell_after_cellbender= statistics.median(f[f>0][0])

    f = pd.DataFrame(ad_lane_filtered.X.sum(axis=0)).T
    Median_UMI_Counts_per_Gene = statistics.median(f[f[0]>0][0])
    try:
        Valid_Droplet_percentage = metrics['Valid Barcodes'].values[0]
    except:
        _='metric not available'
    df1 = ad_lane_filtered.to_df()
    Number_of_cells = len(set(df1.index))
    Total_UMIs_before_10x_filter = np.sum(ad_lane_raw.X) #this may be after the normalisation

    Total_UMIs_after_cellbender_filter = np.sum(ad_lane_filtered.X) #This is more 27840
    Total_UMIs_after_cellbender = sum(all_QC_lane.obs['total_counts']) #This is less 22817

    Droplets_removed_by_filtering = len(set(ad_lane_raw.obs.index)-set(ad_lane_filtered.obs.index))
    Total_Drroplets_before_10x_filtering = len(set(pd.DataFrame(ad_lane_raw.obs).index))
    Doublets_donor = 0
    Unassigned_donor = 0
    try:
        Cells_before_QC_filters=len(all_QC_lane.obs['cell_passes_qc'])
    except:
        # Cell QC column was not defined, sellecting the first entry as the main.
        cell_passes_qc = [item for item in all_QC_lane.obs.columns if 'cell_passes_qc' in item and 'score' not in item][0]
        all_QC_lane.obs['cell_passes_qc']=all_QC_lane.obs[cell_passes_qc]
        Cells_before_QC_filters=len(all_QC_lane.obs['cell_passes_qc'])
        adqc.obs['cell_passes_qc']=adqc.obs[cell_passes_qc]
    
    Cells_passing_QC=len(all_QC_lane.obs[all_QC_lane.obs['cell_passes_qc']])
    Cells_failing_QC=len(all_QC_lane.obs[all_QC_lane.obs['cell_passes_qc']==False])

    UMIS_mapped_to_mitochondrial_genes = sum(all_QC_lane.obs['total_counts_gene_group__mito_transcript'])
    UMIS_mapped_to_ribo_genes = sum(all_QC_lane.obs['total_counts_gene_group__ribo_protein'])
    UMIS_mapped_to_ribo_rna = sum(all_QC_lane.obs['total_counts_gene_group__ribo_rna'])
    UMIs_mapped_to_genes = sum(all_QC_lane.obs['total_counts'])

   
    Total_donors_deconvoluted_in_pool = len(Donor_Cohort_Assignments[Donor_Cohort_Assignments['panel']!='NONE'].index)
    ELGH_Donors_Deconvoluted_in_pool = len(Donor_Cohort_Assignments[Donor_Cohort_Assignments['panel']=='GT_ELGH'].index)
    UKB_Donors_Deconvoluted_in_pool = len(Donor_Cohort_Assignments[Donor_Cohort_Assignments['panel']=='GT_UKBB'].index)
    Spikeins_Deconvoluted_in_pool = len(Donor_Cohort_Assignments[Donor_Cohort_Assignments['panel']=='GT_cell_lines'].index)
    if (Total_donors_deconvoluted_in_pool==0):
        Total_donors_deconvoluted_in_pool=1
    data_donor =[]

    data_donor_for_stats ={ 'cells before QC filters':[],
                            'cells failing QC':[],
                            'cells passing QC':[],
    }
    all_probs = pd.DataFrame()
    all_probs = pd.DataFrame()

    Tranche_Pass_Fail='PASS'
    Tranche_Failure_Reason =' '
    try:
        if (float(Fraction_Reads_in_Cells.strip('%'))<=70):
            Tranche_Pass_Fail='FAIL'
            Tranche_Failure_Reason +='Fraction of reads in cells for pool<=70%; '
    except:
        _='cant validate reasoning'
    try:
        if (Mean_reads_per_cell<=25000):
            Tranche_Pass_Fail='FAIL'
            Tranche_Failure_Reason +='Mean reads per cell for all cells in pool <=25000; '
    except:
        _='cant validate reasoning'        
    try:
        Summary_check = pd.read_csv(f'{args.results_dir}/deconvolution/vireo_raw/{expid}/vireo_{expid}/donor_ids.tsv',sep='\t')
    except:
        try:
            Summary_check = pd.read_csv(f'{args.results_dir}/deconvolution/vireo_raw/vireo_{expid}/donor_ids.tsv',sep='\t')
        except:
            try:
                Summary_check = pd.read_csv(f'{args.results_dir}/deconvolution/vireo_raw/existing_vireo/{expid}/donor_ids.tsv',sep='\t')
            except:
                Summary_check = pd.read_csv(f'{args.results_dir}/deconvolution/vireo_raw/existing_vireo/{expid}/vireo_{expid}/donor_ids.tsv',sep='\t')
            
    try:
        Doublets_donor = len(Summary_check[Summary_check['donor_id']=='doublet'])
    except:
        Doublets_donor = 0
    try:
        Unassigned_donor = len(Summary_check[Summary_check['donor_id']=='unassigned'])
    except:
        Unassigned_donor = 0
    for i in df_donors.index:
        print(i)
        # if df_donors.loc[i]["donor_id"] != 'donor2':
        #     continue
        # else:
        #     _=''
        Donor_Stats=[]
        row = df_donors.loc[i]
        print(row)
        path1 = row['file_path_h5ad']
        if(path1=='all_QC_lane'):
            tp1 = pd.DataFrame(all_QC_lane.obs)
            tp1 = tp1[tp1['donor_id']==row['donor_id']].index
            Deconvoluted_Donor_Data = adqc[tp1]
            Donor_barcodes = Deconvoluted_Donor_Data.obs.index.str.split('-').str[:2]
            Donor_barcodes = Donor_barcodes.str[0]+'-'+Donor_barcodes.str[1]
            Deconvoluted_Donor_Data.obs.index = Donor_barcodes
        else:
            path1 = re.sub('.*/results/', 'results/', path1)
            Deconvoluted_Donor_Data = anndata.read_h5ad(path1)
            Donor_barcodes = Deconvoluted_Donor_Data.obs.index

        #################
        #Deconvolution data
        #################
        
        # Issue with the all_QC_lane is that they are filtered and the unassigned cells are removed - we can merge them back together and 
        if (row["donor_id"] == 'unassigned'):
            Unassigned_donor = len(Deconvoluted_Donor_Data.obs)
        elif (row["donor_id"] == 'doublet'):
            Doublets_donor = len(Deconvoluted_Donor_Data.obs)
        else:
            data = all_QC_lane.obs
            s = data['donor_id']
            s=s.reset_index()
            s = s.set_index('donor_id')
            s.columns=['index']
            try:
                Mengled_barcodes_donor = list(s.loc[row["donor_id"]]['index'])
            except:
                continue
            donor_number =row["donor_id"].replace('donor','') #Todo - will need to change upon Vireo runs with genotype, can just pick it as an fctr
            donor_id = row["donor_id"]
            intersect_set = set(Donor_barcodes).intersection(set(ad_lane_filtered.obs.index))

            all_probs = pd.concat([all_probs,pd.DataFrame(Deconvoluted_Donor_Data.obs['prob_doublet'])])
            Donor_qc_files = Deconvoluted_Donor_Data
            UMIs = np.sum(Donor_qc_files.X)
            Donor_cells_for_donor=len(Deconvoluted_Donor_Data.obs)
            Donor_cells_passes_qc = len(Deconvoluted_Donor_Data.obs[Deconvoluted_Donor_Data.obs['cell_passes_qc']])
            Donor_cells_fails_qc = len(Deconvoluted_Donor_Data.obs[Deconvoluted_Donor_Data.obs['cell_passes_qc']==False])
            
            data_donor_for_stats['cells before QC filters'].append(Donor_cells_for_donor)
            data_donor_for_stats['cells failing QC'].append(Donor_cells_fails_qc)
            data_donor_for_stats['cells passing QC'].append(Donor_cells_passes_qc)

            Donor_cell_assignments = azt.loc[Mengled_barcodes_donor] #for this have to figure out when the cell type is unasigned.
            Cell_types_detected = len(set(Donor_cell_assignments['Azimuth:predicted.celltype.l2']))
            Donor_UMIS_mapped_to_mitochondrial_genes = sum(Donor_qc_files.obs['total_counts_gene_group__mito_transcript'])
            Donor_UMIS_mapped_to_ribo_genes = sum(Donor_qc_files.obs['total_counts_gene_group__ribo_protein'])
            Donor_UMIS_mapped_to_ribo_rna = sum(Donor_qc_files.obs['total_counts_gene_group__ribo_rna'])

            Donor_matrix = Donor_qc_files.to_df()
            Donor_matrix_gene_sums = Donor_matrix.sum()
            genes_detected_with_counts_greater_than_0 = len(Donor_matrix_gene_sums[Donor_matrix_gene_sums>0])
            genes_with_UMI_count_larger_than_3 =  len(Donor_matrix_gene_sums[Donor_matrix_gene_sums>=3])
            Median_UMIs_per_gene= statistics.median(pd.DataFrame(Donor_qc_files.X.sum(axis=0))[0])
            Median_UMIs_per_cell= statistics.median(pd.DataFrame(Donor_qc_files.X.sum(axis=1))[0])

            Cell_numbers = ''
            for type1 in set(Donor_cell_assignments['Azimuth:predicted.celltype.l2']):
                nr_cells_of_this_type = len(Donor_cell_assignments[Donor_cell_assignments['Azimuth:predicted.celltype.l2']==type1])
                Cell_numbers+=f"{type1}:{nr_cells_of_this_type} ; "


            Donor_Stats = gather_donor(
                donor_id=row["donor_id"],
                ad= Donor_qc_files,
                ad_lane_raw=ad_lane_raw,
                qc_obs = obsqc,
                columns_output = columns_output,
                outdir = outdir,
                oufh = oufh,
                lane_id=lane_id
            )
            Donor_Stats['Experiment ID']=args.experiment_name
            if Donor_Stats['Donor id']!='':
                # Only generate donor stats for the donors excluding unasigned and doublets. 
                Pass_Fail='PASS'
                Failure_Reason =' '

                if (Median_UMIs_per_cell<=400):
                    Pass_Fail='FAIL'
                    Failure_Reason +='Median_UMIs_per_cell<=400; '
                if (Donor_UMIS_mapped_to_mitochondrial_genes/UMIs>=0.5):
                    Pass_Fail='FAIL'
                    Failure_Reason +='Donor_UMIS_mapped_to_mitochondrial_genes/UMIs>=0.5; '
                if (Donor_cells_passes_qc<=200):
                    Pass_Fail='FAIL'
                    Failure_Reason +='Donor_cells_passes_qc<=200; '
                if (Donor_cells_for_donor<=400):
                    Pass_Fail='FAIL'
                    Failure_Reason +='Donor_cells_for_donor<=400; '
            
                try:
                    Date_sample_received = ' '
                except:
                    Date_sample_received = 'No sample info available'   

                try:
                    Date_of_sample_sequencing = ' '
                except:
                    Date_of_sample_sequencing = 'No sample info available'   

            
                Donor_Stats_extra = {
                    'Donor id':donor_id,
                    'Donors in pool':Donors_in_pool,
                    'Median UMIs per cell':Median_UMIs_per_cell,
                    'Nr UMIS mapped to mitochondrial genes':Donor_UMIS_mapped_to_mitochondrial_genes,
                    'Nr cells passes qc':Donor_cells_passes_qc,
                    'Total Nr cells for donor':Donor_cells_for_donor, 
                    'Date sample received':Date_sample_received, #Do this
                    'Date of sample sequencing':Date_of_sample_sequencing, 
                    'Donor_number': donor_number,
                    'Nr UMIs':UMIs,
                    'Median UMIs per gene':Median_UMIs_per_gene,
                    'Nr UMIS mapped to ribo genes':Donor_UMIS_mapped_to_ribo_genes,
                    'Nr UMIS mapped to ribo rna':Donor_UMIS_mapped_to_ribo_rna,
                    'Total Nr cells fails qc':Donor_cells_fails_qc,
                    'Genes detected with counts > 0':genes_detected_with_counts_greater_than_0,
                    'Genes with UMI count >= 3':genes_with_UMI_count_larger_than_3,
                    'Overall Pass Fail':Pass_Fail,
                    'Failure reason':Failure_Reason,
                    'Date of data transfer':date_of_transfer, # make this the last day of the month
                    'QC Report end date':date_now,
                    'Cell types detected':Cell_types_detected,
                    'Cell type numbers':Cell_numbers,
                    
                }

                Donor_Stats.update(Donor_Stats_extra)
                data_donor.append(Donor_Stats)

        fctr += 1
    all_probs = all_probs[~all_probs.index.duplicated(keep='first')]
    azt['prob_doublet']=all_probs['prob_doublet']
    Donor_df = pd.DataFrame(data_donor)

    try:
        Median_cells_passes_qc=statistics.median(data_donor_for_stats['cells passing QC'])
        Median_cells_fails_qc=statistics.median(data_donor_for_stats['cells failing QC'])
        Median_Nr_cells_for_donor=statistics.median(data_donor_for_stats['cells before QC filters'])
        Stdev_cells_passes_qc=statistics.stdev(data_donor_for_stats['cells passing QC'])
        Stdev_cells_fails_qc=statistics.stdev(data_donor_for_stats['cells failing QC'])
        Stdev_Nr_cells_for_donor=statistics.stdev(data_donor_for_stats['cells before QC filters'])
    except:
        Median_cells_passes_qc=data_donor_for_stats['cells passing QC'][0]
        Median_cells_fails_qc=data_donor_for_stats['cells failing QC'][0]
        Median_Nr_cells_for_donor=data_donor_for_stats['cells before QC filters'][0]
        Stdev_cells_passes_qc=0
        Stdev_cells_fails_qc=0
        Stdev_Nr_cells_for_donor=0       
    Percentage_of_unassigned_cells = Unassigned_donor/(Donors_in_pool+Doublets_donor+Unassigned_donor)*100  

    data_tranche = {
        'Experiment id':args.experiment_name,
        'Pool id':list(set(Donor_df['Pool ID']))[0],
        'Machine id':Machine_id, #Change this - it will come from the extra metadata file if available
        'Run id':Run_ID, #Generate - feed in from extra metadate if available
        'Date of sample sequencing':date_of_sequencing,#Change this
        'Donors in pool':Donors_in_pool,
        'Total donors deconvoluted in pool':Total_donors_deconvoluted_in_pool,
        'UKB donors expected in pool':Count_of_UKB, 
        'UKB donors deconvoluted in pool':UKB_Donors_Deconvoluted_in_pool, 
        'ELGH donors expected in the pool':Count_of_ELGH, 
        'ELGH donors deconvoluted in the pool':ELGH_Donors_Deconvoluted_in_pool, 
        'Spikeins expected in the pool':Count_of_Spikeins, 
        'Spikeins deconvoluted in the pool':Spikeins_Deconvoluted_in_pool, 
        'Chromium channel number':chromium_channel,
        'Number of Reads':Number_of_Reads,
        'Fraction Reads in Cells':Fraction_Reads_in_Cells,
        'Mean Reads per Cell':Mean_reads_per_cell,

        'Median UMI Counts per Droplet before filter':Median_UMI_Counts_per_before_filter,
        'Median UMI Counts per Droplet after Cellranger filter':Median_UMI_Counts_per_cellranger,
        'Median UMI Counts per Droplet after Cellbender filter':Median_UMI_Counts_per_Cell_after_cellbender,
        'Median UMI Counts per Cell after Cellbender filter; doublet removal; unassigned removal':Median_UMI_Counts_per_Cell_after_all_filter,
        'Total UMIs before filter':Total_UMIs_before_10x_filter,
        'Total UMIs after Cellbender filter':Total_UMIs_after_cellbender_filter,
        'Total UMIs after cellbender filter; doublet removal; unassigned removal':Total_UMIs_after_cellbender,
        
        'UMIS mapped to mitochondrial genes':UMIS_mapped_to_mitochondrial_genes,
        'UMIs mapped to genes':UMIs_mapped_to_genes,
        'UMIS mapped to ribo genes':UMIS_mapped_to_ribo_genes,
        'UMIS mapped to ribo rna':UMIS_mapped_to_ribo_rna,
        'Percentage of unassigned cells':Percentage_of_unassigned_cells,

        'Droplets before filtering':Total_Drroplets_before_10x_filtering,
        'Empty droplets - removed by filtering':Droplets_removed_by_filtering,
        'Droplets identified as doublet':Doublets_donor,
        'Total Droplets with a single cell':Cells_before_QC_filters+Unassigned_donor,
        'Droplets with donor unassigned':Unassigned_donor+Doublets_donor,
        'Number of unassigned cells':Unassigned_donor,
        'Cells before QC filters':Cells_before_QC_filters,
        'Total Cells failing QC':Cells_failing_QC,
        'Total Cells passing QC':Cells_passing_QC,
        'Total Droplets with donor assignment':Cells_before_QC_filters,

        'Median cells passes qc':Median_cells_passes_qc,
        'Median cells fails qc':Median_cells_fails_qc,
        'Median Nr cells for donor':Median_Nr_cells_for_donor,
        'Stdev cells passes qc':Stdev_cells_passes_qc,
        'Stdev cells fails qc':Stdev_cells_fails_qc,
        'Stdev Nr cells for donor':Stdev_Nr_cells_for_donor,
        'QC Report end date':date_now,
        'Tranche Pass/Fail':Tranche_Pass_Fail,
        'Tranche Failure Reason':Tranche_Failure_Reason
    }
    return fctr, data_tranche, data_donor

def set_argument_parser():

    parser = argparse.ArgumentParser(description="gather minimal dataset for scRNAseq handover")
    parser.add_argument("--output-dir", "-o",
        default=os.curdir,
        help="output directory",
        dest="outdir")
    parser.add_argument("--results_dir", required = True,
                        help="Folder where all the results are stored",
                        dest='results_dir')
    parser.add_argument("--input_table", required = True,
                    help="The input folder used",
                    dest='input_table')
    parser.add_argument("--write_h5", required = True,
                help="Should we write the h5ad files for each donor?",
                dest='write_h5')
    parser.add_argument("--extra_metadata", required = False,
                    help="The input folder used",
                    dest='extra_metadata')
    parser.add_argument("--cellbender", required = True,
                help="Cellbender paths if not cellranger filtering used",
                dest='cellbender')
    parser.add_argument("--resolution", required = True,
                help="Cellbender resolution used",
                dest='resolution')
    parser.add_argument("--experiment_name", required = False,default='default',
            help="Cellbender resolution used",
            dest='experiment_name')
    parser.add_argument("--extra_meta", required = False,default=None,
            help="Extra Meta To merge in final files",
            dest='extra_meta')
                
    return parser.parse_args()

if __name__ == '__main__':
    args = set_argument_parser()

    if args.outdir != os.curdir and not os.access(args.outdir, os.F_OK):
        os.mkdir(args.outdir)
        os.mkdir(f"{args.outdir}_summary")
    if (args.write_h5=='false'):
        write_h5=False
    else:
       write_h5=True 

    oufh = open(os.path.join(args.outdir, "files.tsv"), 'w')
    oufh.write("experiment_id\tdonor_id\tfilename_h5ad\tfilename_annotation_tsv\n")
    df_raw = pandas.read_table(args.input_table, index_col = 'experiment_id')
    if (args.cellbender)=='cellranger':
        # Here we do not use cellbender and go with default cellranger
        df_cellbender = None
    else:
        # Here we have run the cellbender as par of pipeline. 
        # cellbender/*/cellbender-epochs_*/cellbender-FPR_0pt01-filtered_10x_mtx
        file_path = glob.glob(f'{args.results_dir}/nf-preprocessing/cellbender/*/cellbender-epochs_*/*{args.resolution}*10x_mtx*')
        file_path2 = glob.glob(f'{args.results_dir}/nf-preprocessing/cellbender/*/*{args.resolution}*10x_mtx*')
        joined_file_paths = file_path+file_path2
        df_cellbender = pd.DataFrame(joined_file_paths,columns=['data_path_10x_format'])
        df_cellbender['experiment_id']=df_cellbender['data_path_10x_format'].str.split('/').str[-3]
        df_cellbender= df_cellbender.set_index('experiment_id')
    
    Resolution = args.resolution
    
    # Load the final QCd dataset
    a1 = glob.glob(f'{args.results_dir}/*/*/*outlier_filtered_adata.h5ad')
    a2 =glob.glob(f'{args.results_dir}/*/*outlier_filtered_adata.h5ad')
    a3 =glob.glob(f'{args.results_dir}/*outlier_filtered_adata.h5ad')
    all_inter = list(set(a3).union(set(a2)).union(set(a1)))[0]
    adqc = anndata.read_h5ad(all_inter, backed='r')
    
    if 'log1p_cp10k' not in adqc.layers:
        adqc = adqc.to_memory()
        normalized_counts = sc.pp.normalize_total(
            adqc,
            target_sum=1e4,
            exclude_highly_expressed=False,
            key_added='normalization_factor',  # add to adata.obs
            inplace=False
        )['X']
        log1p_cp10k = np.log1p(normalized_counts)
        adqc.layers['log1p_cp10k'] = log1p_cp10k
        adqc.uns['log1p_cp10k'] = {'transformation': 'ln(CP10k+1)'}
        del normalized_counts
        del log1p_cp10k
    
    try:
        a1 = glob.glob(f'{args.results_dir}/*/*/*adata-normalized_*.h5ad')
        a2 =glob.glob(f'{args.results_dir}/*/*adata-normalized_*.h5ad')
        a3 =glob.glob(f'{args.results_dir}/*adata-normalized_*.h5ad')
        all_inter = list(set(a3).union(set(a2)).union(set(a1)))[0]
        adqc_norm = anndata.read_h5ad(all_inter, backed='r')
        # Here we want to add any metadata columns that may be missed for donors. 
        adqc.obs['phase'] = adqc_norm.obs['phase']
        adqc.obs['S_score'] = adqc_norm.obs['S_score']
        adqc.obs['G2M_score'] = adqc_norm.obs['G2M_score']
    except:
        _='not normalised'
        
    fctr = 0
    data_tranche_all=[]
    data_donor_all=[]
    count = 1
    All_probs_and_celltypes = pd.DataFrame()
    Sample_metadata = pd.DataFrame()
    
    # SETTING TRANCHE NAME
    # try:
    #     adqc.obs['cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2'] = adqc.obs['cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2:score'].astype(float,errors='ignore')
    # except:
    #     _='no values associated'
    # try:
    #     adqc.obs['cell_passes_qc-per:all_together::exclude:score']= adqc.obs['cell_passes_qc-per:all_together::exclude:score'].astype(float,errors='ignore')
    # except:
    #     _='no values associated'
    
    # Now we calculate all the statistics for each of the pools.
    for expid in df_raw.index:
        print(expid)
        # if expid != 'CRD_CMB13637311':
        #     continue
        s = adqc.obs['convoluted_samplename'] == expid
        ad = adqc[s]
        if ad.n_obs == 0:
            continue #Here no cells has passed the qc thresholds.
        nf, data_tranche, data_donor = gather_pool(expid, args, df_raw, df_cellbender, ad, oufh = oufh, lane_id=count,Resolution=Resolution)
        data_tranche_all.append(data_tranche)
        data_donor_all= data_donor_all+data_donor
        count += 1
        fctr += nf

    Donor_Report = pd.DataFrame(data_donor_all)
    Tranche_Report = pd.DataFrame(data_tranche_all)

    Donor_Report['Pool_ID.Donor_Id']=Donor_Report['Pool ID']+'_'+Donor_Report['Donor id']
    Donor_Report=Donor_Report.set_index('Pool_ID.Donor_Id')
    Donor_Report.to_csv(f'{args.outdir}_summary/{args.experiment_name}_Donor_Report.tsv',sep='\t')
    Tranche_Report.to_csv(f'{args.outdir}_summary/{args.experiment_name}_Tranche_Report.tsv',sep='\t',index=False)
    oufh.close()
    sys.stderr.write("# wrote {:d} files to directory {:s}\n".format(fctr, args.outdir))
    exit(0)