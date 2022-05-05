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
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import glob
import pandas as pd
import numpy as np
import tables
from plot_cellranger_vs_cellbender import anndata_from_h5
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
AZIMUTH_ASSIGNMENTS_FNSUFFIX = '_predicted_celltype_l2.tsv.gz'
SCRUBLET_ASSIGNMENTS_FNSUFFIX = '-scrublet.tsv.gz'

COLUMNS_AZIMUTH = {
    'predicted.celltype.l2': 'azimuth.celltyp.l2',
    'predicted.celltype.l2.score': 'azimuth.pred.score.l2',
    'mapping.score': 'azimuth.map.score'
    }
COLUMNS_DECONV = {
    'donor_id': 'vireo.donor.id',
    'prob_max': 'vireo.prob.max',
    'prob_doublet': 'vireo.prob.doublet'
    }
COLUMNS_QC = {
    'cell_passes_qc': 'qc.filter.pass',
    'total_counts': 'qc.umi.count.total',
    'total_counts_gene_group__mito_transcript': 'qc.umi.count.mt',
    'pct_counts_gene_group__mito_transcript': 'qc.umi.perc.mt',
    'n_genes_by_counts': 'qc.genes.detected.count'
    }
COLUMNS_CELLBENDER = {'cellbender_latent_probability': 'cellbender.latent.probability'}
COLUMNS_DATASET = {
    'experiment_id': 'experiment.id',
    'chromium_run_id': 'chromium.run.id',
    'chromium_lane': 'chromium.lane'
    }
COLUMNS_SCRUBLET = {
    'scrublet__multiplet_scores': 'scrublet.scores',
    'scrublet__predicted_multiplet': 'scrublet.multiplet',
    'scrublet__multiplet_zscores': 'scrublet.zscores'
    }
COLUMNS_OUTPUT = \
    {**COLUMNS_DATASET, **COLUMNS_CELLBENDER, **COLUMNS_DECONV, **COLUMNS_QC, **COLUMNS_AZIMUTH}
COLUMNS_OUTPUT_WITH_SCRUBLET = \
    {**COLUMNS_DATASET, **COLUMNS_CELLBENDER, **COLUMNS_DECONV, **COLUMNS_SCRUBLET, **COLUMNS_QC, **COLUMNS_AZIMUTH}

def get_df_from_mangled_index(df, expid):
    idx = df.index.str.split(pat='-{}__'.format(expid))
    xf = pandas.DataFrame.from_records(idx, columns = ('barcode', 'donor'), index = df.index)
    if xf.shape[0] != df.shape[0]:
        sys.exit("ERROR: when untangling mangled index.")
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

def gather_azimuth_annotation(expid, datadir_azimuth, index_label = None):
    # e.g. A4C06803ACD34DFB-adata_franke_Pilot_3_lane_3_predicted_celltype.tsv.gz
    filpath = None
    expid2=expid
    fnsfx = '_{}{}'.format(expid, AZIMUTH_ASSIGNMENTS_FNSUFFIX)
    for fn in os.listdir(datadir_azimuth):
        if fn.endswith(fnsfx):
            filpath = os.path.join(datadir_azimuth, fn)
            break
    if not filpath:
        expid ='full'
        fnsfx = '_{}{}'.format(expid, AZIMUTH_ASSIGNMENTS_FNSUFFIX)
        for fn in os.listdir(datadir_azimuth):
            if fn.endswith(fnsfx):
                filpath = os.path.join(datadir_azimuth, fn)
                break
            
    if not filpath:
        sys.exit("ERROR: could not find filename suffix '{}' in direcotry {}\n"
            .format(fnsfx, datadir_azimuth))
    azt = pandas.read_table(filpath)
    if (expid=='full'):
        expid=expid2
        azt = azt[azt.index.str.contains(expid2)]
        # filter to only the experiment id inputs.
    df = get_df_from_mangled_index(azt, expid)
    azt.insert(0, "barcode", df["barcode"])
    azt.insert(1, "donor", df["donor"])
    azt.insert(2, "experiment_id", expid)
    if index_label is not None and index_label == "barcode":
        azt.insert(0, "mangled_cell_id", df.index)
        azt = azt.set_index("barcode", drop = True)
    return azt

def load_scrublet_assignments(expid, datadir_scrublet):
    filpath = None
    fnam = '{}{}'.format(expid, SCRUBLET_ASSIGNMENTS_FNSUFFIX)
    for fn in os.listdir(datadir_scrublet):
        if fn == fnam:
            filpath = os.path.join(datadir_scrublet, fn)
            break
 
    if not filpath:
        sys.exit("ERROR: could not find filename '{}' in direcotry {}\n"
            .format(fnam, datadir_scrublet))
    sys.stderr.write("loading scrublet annotation from file {} ...\n".format(filpath))
    scb = pandas.read_table(filpath).set_index('cell_barcode', drop = True)
    return scb

def fetch_qc_obs_from_anndata(adqc, expid, df_cellbender = None):

    s = adqc.obs['convoluted_samplename'] == expid
    if s.shape[0] < 1:
        sys.exit("ERROR: No QC data found for experiment_id = '{:s}'"
            .format(expid))
    ad = adqc[s]
    df = get_df_from_mangled_index(ad.obs, expid)
    #df.insert(0, 'barcode', df.index.values)
    df_pre = pandas.concat([df,ad.obs], axis = 1)
    df_pre['mengled_index'] = df_pre.index
    df = df_pre[['barcode', 'donor',
        'n_genes_by_counts', 'log1p_n_genes_by_counts',
        'total_counts', 'log1p_total_counts',
        'total_counts_gene_group__mito_transcript', 'pct_counts_gene_group__mito_transcript',
        'cell_passes_qc','mengled_index'
    ]].set_index("barcode", drop = True)

    if df_cellbender is not None:
        # cellbender removes the barcodes - 
        dfcb = fetch_cellbender_annotation(df_cellbender, expid)
        dc = pandas.concat([df, dfcb], axis = 1, join = 'inner')
        if dc.shape[0] != df.shape[0]:
            sys.exit("ERROR: barcodes missing in cellbender file.")
        df = dc.copy()
        dc=dc.set_index('mengled_index')
        
        ad.obs['cellbender_latent_probability']=dc['cellbender_latent_probability']
    return df,ad

def fetch_cellbender_annotation(df_cellbender, expid):
    dirpath = df_cellbender.loc[expid, 'data_path_10x_format']
    try:
        h5_path = f"{args.results_dir}/{os.path.dirname(dirpath)}/cellbender_FPR_0pt05_filtered.h5"
        f = h5py.File(h5_path, 'r')
    except:
        h5_path = f"{os.path.dirname(dirpath)}/cellbender_FPR_0pt05_filtered.h5"
        f = h5py.File(h5_path, 'r')
    # ad = scanpy.read_10x_h5(h5_path, genome='background_removed')
    # interesting data is in /matrix/barcodes and matrix/latent_cell_probability
    f = h5py.File(h5_path, 'r')
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

def gather_donor(donor_id, ad, ad_lane_raw, azimuth_annot, qc_obs, columns_output = COLUMNS_OUTPUT,outdir = os.curdir,oufh = sys.stdout,lane_id=1):
    
    oufnam = "{}.{}".format(expid, donor_id)
    sys.stderr.write("processing {} {} ...\n".format(expid, donor_id))
    if write_h5:
        oufh.write("{}\t{}\t{}.h5ad\t{}.tsv\n".format(expid, donor_id, oufnam, oufnam))

    # loading deconvoluted dataset
    

    ad.var.index.name = "ensembl_id"
    ad.raw = ad_lane_raw[ad.obs.index, :]
    if donor_id != "unassigned" and donor_id != "doublet":
        # add annotation from QC
        df = pandas.concat([ad.obs, azimuth_annot.loc[azimuth_annot.donor == donor_id]], axis = 1, join = 'outer')
        df = df[['experiment_id'] + list(COLUMNS_DECONV.keys()) + list(COLUMNS_AZIMUTH.keys())]
        try:
            df = get_lane_and_runid_from_experiment_id(df, insert_pos = 1)
        except:
            # here we do not know the lane ID.
            df['chromium_lane']=lane_id
            df['chromium_run_id']=df['experiment_id'][0]



        dfqc = qc_obs[qc_obs.donor == donor_id]
        dt = pandas.concat([df,dfqc], axis = 1, join = 'inner')
        if dt.shape[0] != ad.obs.shape[0]:
            sys.exit("ERROR: Number of cells in file {:s} changed from {:d} to {:d}\n"
                .format(oufnam, ad.obs.shape[0], df.shape[0]))

        colnams = list(columns_output.keys())
        ad.obs = dt[colnams].rename(columns = columns_output)
        dt = pandas.concat([df, dfqc], axis = 1, join = 'outer')[colnams]
        dt.rename(columns = columns_output, inplace = True)


            # Stats
        print('Performing the stats analysis')
        experiment_id = list(set(df.experiment_id))[0]
        pool_id = list(set(df.chromium_run_id))[0] #???
        chromium_channel_number = list(set(df.chromium_run_id))[0]
        donor_id = list(set(df.donor_id))[0]
    else:
        dt = ad.obs
        experiment_id='';pool_id='';chromium_channel_number='';donor_id=''

    dt.index.name = 'barcode'
    ad.obs.index.name = 'barcode'
    dt.to_csv(os.path.join(outdir, oufnam + '.tsv'), sep = "\t", na_rep = "N/A")
    sys.stderr.write("writing file {} ...\n".format(oufnam))
    if write_h5:
        ad.write(os.path.join(outdir, oufnam + '.h5ad'))

    return {
        'Experiment ID':experiment_id,
        'Pool ID':pool_id,
        'Chromium channel number':chromium_channel_number,
        'Donor id':donor_id
    }

def gather_pool(expid, args, df_raw, df_cellbender, adqc, oufh = sys.stdout,lane_id=1):
    
    ######################
    #Cellranger datasets
    ######################
    columns_output = {**COLUMNS_DATASET, **COLUMNS_DECONV, **COLUMNS_QC}
    #Unfiltered
    compression_opts = 'gzip'
    adata_cellranger_raw = scanpy.read_10x_mtx(f"{df_raw.loc[expid, 'data_path_10x_format']}/raw_feature_bc_matrix")
    zero_count_cells_cellranger_raw = adata_cellranger_raw.obs_names[np.where(adata_cellranger_raw.X.sum(axis=1) == 0)[0]]
    ad_lane_raw = adata_cellranger_raw[adata_cellranger_raw.obs_names.difference(zero_count_cells_cellranger_raw, sort=False)]
    scanpy.pp.calculate_qc_metrics(adata_cellranger_raw, inplace=True)
    df_total_counts = pd.DataFrame(data= adata_cellranger_raw.obs.sort_values(by=['total_counts'], ascending=False).total_counts)
    df_total_counts['barcode_row_number'] = df_total_counts.reset_index().index + 1 
    df_total_counts['barcodes'] = df_total_counts.index
    df_total_counts_cellranger_raw = df_total_counts
    df_total_counts_cellranger_raw['dataset']='Cellranger Raw'

    if df_cellbender is not None:
        try:
            # depends whether the absolute or relative path was recorded.
            ad_lane_filtered = scanpy.read_10x_mtx(f"{df_cellbender.loc[expid, 'data_path_10x_format']}")
        except:
            ad_lane_filtered = scanpy.read_10x_mtx(f"{args.results_dir}/{df_cellbender.loc[expid, 'data_path_10x_format']}")
        
        dfcb = fetch_cellbender_annotation(df_cellbender, expid)
        columns_output = {**columns_output, **COLUMNS_CELLBENDER}
    else:
        ad_lane_filtered = scanpy.read_10x_mtx(f"{df_raw.loc[expid, 'data_path_10x_format']}/filtered_feature_bc_matrix")
        df_cellbender=None
    
    zero_count_cells_cellranger_filtered = ad_lane_filtered.obs_names[np.where(ad_lane_filtered.X.sum(axis=1) == 0)[0]]
    ad_lane_filtered = ad_lane_filtered[ad_lane_filtered.obs_names.difference(zero_count_cells_cellranger_filtered, sort=False)]
    adata_cellranger_filtered=ad_lane_filtered
    scanpy.pp.calculate_qc_metrics(adata_cellranger_filtered, inplace=True)
    df_total_counts = pd.DataFrame(data= adata_cellranger_filtered.obs.sort_values(by=['total_counts'], ascending=False).total_counts)
    df_total_counts['barcodes'] = df_total_counts.index
    df_total_counts['barcode_row_number'] = df_total_counts.reset_index().index + 1 
    df_total_counts_cellranger_filtered = df_total_counts
    df_total_counts_cellranger_filtered['dataset'] = 'Cellranger Filtered'

    #############
    #Cellranger Metrics Datasheet
    #############
    metrics = pd.read_csv(df_raw.loc[expid, 'data_path_10x_format']+'/metrics_summary.csv')
    
    #############
    #Azimuth cell-type assignments
    #############
    datadir_azimuth = f'{args.results_dir}/celltype/azimuth' 
    if os.path.isdir(datadir_azimuth):
        try:
            azt = gather_azimuth_annotation(
                expid, datadir_azimuth=datadir_azimuth,
                index_label = 'barcode')
            columns_output = {**columns_output, **COLUMNS_AZIMUTH}
        except:
            azt = None
    else:
        azt = None

    ##########################
    # Scrublet
    #########################
    datadir_scrublet=f'{args.results_dir}/multiplet.method=scrublet'
    if os.path.isdir(datadir_scrublet):
        # Scrublet loading QC
        scb = load_scrublet_assignments(
            expid,
            datadir_scrublet=datadir_scrublet
        )
        columns_output = {**columns_output,  **COLUMNS_SCRUBLET}
    else:
        scb = None
        
    
    ############################################################
    # Loading deconvoluted data including unassigned and doublets
    ###########################################
    datadir_deconv=f'{args.results_dir}/deconvolution/split_donor_h5ad'
    donor_table = os.path.join(datadir_deconv, expid, "{}.donors.h5ad.tsv".format(expid))
    df_donors = pandas.read_table(donor_table, header=None, names=("experiment_id", "donor_id", "file_path_h5ad"))
    obsqc,all_QC_lane = fetch_qc_obs_from_anndata(adqc, expid, df_cellbender = df_cellbender)

    if scb is not None:
        obsqc = pandas.concat([obsqc,scb], axis = 1, join = 'outer')
    
    fctr = 0

    #####################
    #Performing Calculations and gathering data
    #####################

    Raw_counts_data_per_lane = ad_lane_raw
    Per_lane_QC_File_data = obsqc
    Azimuth_Cell_Assignments_data = azt
    Deconvoluted_Donor_Data_sheet = df_donors
    All_AnnData_QC_lane = adqc
    Adata_counts = adqc.to_df()
    
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
    Number_of_Reads = int(metrics['Number of Reads'].values[0].replace(',','')) # NOT SURE HOW THIS IS CALCULATED.
    Fraction_Reads_in_Cells = metrics['Fraction Reads in Cells'].values[0]
    Median_reads_per_cell = int(metrics['Mean Reads per Cell'].values[0].replace(',',''))

    f = pd.DataFrame(ad_lane_filtered.X.sum(axis=1))
    Median_UMI_Counts_per_Cell_10x_filter= statistics.median(f[f>0][0])
    f = pd.DataFrame(ad_lane_filtered.X.sum(axis=1))
    Median_UMI_Counts_per_Cell= statistics.median(f[f>0][0])
    f = pd.DataFrame(ad_lane_filtered.X.sum(axis=0)).T
    Median_UMI_Counts_per_Gene = statistics.median(f[f[0]>0][0])
    Valid_Droplet_percentage = metrics['Valid Barcodes'].values[0]
    df1 = ad_lane_filtered.to_df()
    Number_of_cells = len(set(df1.index))
    Total_UMIs_before_10x_filter = np.sum(ad_lane_raw.X) #this may be after the normalisation
    Total_UMIs_after_cellbender_filter = np.sum(ad_lane_filtered.X)
    Droplets_removed_by_filtering = len(set(ad_lane_raw.obs.index)-set(ad_lane_filtered.obs.index))
    Total_Drroplets_before_10x_filtering = len(set(pd.DataFrame(ad_lane_raw.obs).index))
    Doublets_donor = 0
    Unassigned_donor = 0
    Cells_before_QC_filters=len(all_QC_lane.obs['cell_passes_qc'])
    Cells_passing_QC=len(all_QC_lane.obs[all_QC_lane.obs['cell_passes_qc']])
    Cells_failing_QC=len(all_QC_lane.obs[all_QC_lane.obs['cell_passes_qc']==False])
    try:
        Azimuth_Cell_Assignments_data=Azimuth_Cell_Assignments_data.set_index('mangled_cell_id')
        all_QC_lane.obs['predicted celltype']=Azimuth_Cell_Assignments_data['predicted.celltype.l2']
    
    except:
        print('skipped az')
    UMIS_mapped_to_mitochondrial_genes = sum(all_QC_lane.obs['total_counts_gene_group__mito_transcript'])
    UMIS_mapped_to_ribo_genes = sum(all_QC_lane.obs['total_counts_gene_group__ribo_protein'])
    UMIS_mapped_to_ribo_rna = sum(all_QC_lane.obs['total_counts_gene_group__ribo_rna'])
    UMIs_mapped_to_genes = sum(all_QC_lane.obs['total_counts'])
    Total_UMIs_after_cellbender = sum(all_QC_lane.obs['total_counts'])

    data_donor =[]

    data_donor_for_stats ={ 'cells before QC filters':[],
                            'cells failing QC':[],
                            'cells passing QC':[],
    }
    all_probs = pd.DataFrame()
    for i in df_donors.index:
        # feeds in the individual assignments here.
        row = df_donors.loc[i]
        path1 = row['file_path_h5ad']
        #################
        #Deconvolution data
        #################
        Deconvoluted_Donor_Data = anndata.read_h5ad(path1)
        Donor_barcodes = Deconvoluted_Donor_Data.obs.index
        
        # issue with the all_QC_lane is that they are filtered and the unassigned cells are removed - ve can merge them back together and 
        if (row["donor_id"] == 'unassigned'):
            Unassigned_donor = len(Deconvoluted_Donor_Data.obs)
        elif (row["donor_id"] == 'doublet'):
            Doublets_donor = len(Deconvoluted_Donor_Data.obs)
        else:
            data = all_QC_lane.obs
            s = data['donor_id']
            s=s.reset_index()
            s = s.set_index('donor_id')
            Mengled_barcodes_donor = list(s.loc[row["donor_id"]]['index'])
            donor_number =row["donor_id"].replace('donor','') #Todo - will need to change upon Vireo runs with genotype, can just pick it as an fctr
            donor_id = row["donor_id"]
            intersect_set = set(Donor_barcodes).intersection(set(ad_lane_filtered.obs.index))

            all_probs = pd.concat([all_probs,pd.DataFrame(Deconvoluted_Donor_Data.obs['prob_doublet'])])
            Donor_qc_files = all_QC_lane[Mengled_barcodes_donor]
            UMIs = np.sum(Donor_qc_files.X)
            Donor_cells_for_donor=len(all_QC_lane[Mengled_barcodes_donor].obs)
            Donor_cells_passes_qc = len(all_QC_lane[Mengled_barcodes_donor].obs[all_QC_lane.obs['cell_passes_qc']])
            Donor_cells_fails_qc = len(all_QC_lane[Mengled_barcodes_donor].obs[all_QC_lane.obs['cell_passes_qc']==False])
            
            data_donor_for_stats['cells before QC filters'].append(Donor_cells_for_donor)
            data_donor_for_stats['cells failing QC'].append(Donor_cells_passes_qc)
            data_donor_for_stats['cells passing QC'].append(Donor_cells_fails_qc)

            Donor_cell_assignments = Azimuth_Cell_Assignments_data.loc[Mengled_barcodes_donor] #for this have to figure out when the cell type is unasigned.
            Cell_types_detected = len(set(Donor_cell_assignments['predicted.celltype.l2']))
            Donor_UMIS_mapped_to_mitochondrial_genes = sum(Donor_qc_files.obs['total_counts_gene_group__mito_transcript'])
            Donor_UMIS_mapped_to_ribo_genes = sum(Donor_qc_files.obs['total_counts_gene_group__ribo_protein'])
            Donor_UMIS_mapped_to_ribo_rna = sum(Donor_qc_files.obs['total_counts_gene_group__ribo_rna'])

            Donor_matrix = Donor_qc_files.to_df()
            Donor_matrix_gene_sums = Donor_matrix.sum()
            genes_detected_with_counts_greater_than_0 = len(Donor_matrix_gene_sums[Donor_matrix_gene_sums>0])
            genes_with_UMI_count_larger_than_3 =  len(Donor_matrix_gene_sums[Donor_matrix_gene_sums>=3])
            Median_UMIs_per_gene= statistics.median(pd.DataFrame(Donor_qc_files.X.sum(axis=1))[0])
            Median_UMIs_per_cell= statistics.median(pd.DataFrame(Donor_qc_files.X.sum(axis=1))[0])

            Cell_numbers = ''
            for type in set(Donor_cell_assignments['predicted.celltype.l2']):
                nr_cells_of_this_type = len(Donor_cell_assignments[Donor_cell_assignments['predicted.celltype.l2']==type])
                Cell_numbers+=f"{type}:{nr_cells_of_this_type} ; "


        Donor_Stats = gather_donor(
            row["donor_id"],
            Deconvoluted_Donor_Data,
            ad_lane_raw,
            azimuth_annot = azt,
            qc_obs = obsqc,
            columns_output = columns_output,
            outdir = args.outdir,
            oufh = oufh,
            lane_id=lane_id
        )

        if Donor_Stats['Donor id']!='':
            # Only generate donor stats for the donors excluding unasigned and doublets. 
            Pass_Fail='PASS'
            Failure_Reason =''

            if (Median_UMIs_per_gene<=400):
                Pass_Fail='FAIL'
                Failure_Reason +='Median_UMIs_per_gene<=400; '
            if (Donor_UMIS_mapped_to_mitochondrial_genes/UMIs>=0.5):
                Pass_Fail='FAIL'
                Failure_Reason +='Donor_UMIS_mapped_to_mitochondrial_genes/UMIs>=0.5; '
            if (Donor_cells_passes_qc<=200):
                Pass_Fail='FAIL'
                Failure_Reason +='Donor_cells_passes_qc<=200; '
            if (Donor_cells_for_donor<=400):
                Pass_Fail='FAIL'
                Failure_Reason +='Donor_cells_for_donor<=400; '
        
            Donor_Stats_extra = {
                'Donor id':donor_id,
                'Donors in pool':Donors_in_pool,
                'Median UMIs per cell':Median_UMIs_per_cell,
                'Nr UMIS mapped to mitochondrial genes':Donor_UMIS_mapped_to_mitochondrial_genes,
                'Nr cells passes qc':Donor_cells_passes_qc,
                'Total Nr cells for donor':Donor_cells_for_donor, 
                'Date sample received':'comes from extra metadata', #Do this
                'Date of sample sequencing':'comes from extra metadata', 
                'Donor_number': donor_number,
                'Nr UMIs':UMIs,
                'Median UMIs per gene':Median_UMIs_per_gene,
                'Nr UMIS mapped to ribo genes':Donor_UMIS_mapped_to_ribo_genes,
                'Nr UMIS mapped to ribo rna':Donor_UMIS_mapped_to_ribo_rna,
                'Total Nr cells fails qc':Donor_cells_fails_qc,
                'Genes detected with counts > 0':genes_detected_with_counts_greater_than_0,
                'Genes with UMI count >= 3':genes_with_UMI_count_larger_than_3,
                'Cell type numbers':Cell_numbers,
                'Cell types detected':Cell_types_detected,
                'Overall Pass Fail':Pass_Fail,
                'Failure reason':Failure_Reason,
                'Date of data transfer':date_of_transfer, # make this the last day of the month
                'QC Report end date':date_now,
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
        'Experiment id':expid,
        'Pool id':list(set(Donor_df['Pool ID']))[0],
        'Machine id':'XXXX', #Change this - it will come from the extra metadata file if available
        'Run id':'Comes from the fetch', #Generate - feed in from extra metadate if available
        'Date of sample sequencing':'Comes from the LIMs or Sequencescape',#Change this
        'Total Droplets with donor assignment':Cells_before_QC_filters,
        'Droplets identified as doublet':Doublets_donor,
        'Droplets with donor unassigned':Unassigned_donor,
        'UKB donors deconvoluted in pool':0, #change this - would be better if this was an automated input.  ie - if 2 vcfs fed in then 2 matches
        'ELGH donors in the pool':0, #change this - would be better if this was an automated input - ie - if 2 vcfs fed in then 2 matches
        'Chromium channel number':list(set(Donor_df['Chromium channel number']))[0],
        'Donors in pool':Donors_in_pool,
        'Number of Reads':Number_of_Reads,
        'Fraction Reads in Cells':Fraction_Reads_in_Cells,
        'Mean Reads per Cell':Median_reads_per_cell,
        'Median UMI Counts per Cell after 10x filter':Median_UMI_Counts_per_Cell_10x_filter,
        'Median UMI Counts per Cell after Cellbender':Median_UMI_Counts_per_Cell,
        'Total UMIs before filter':Total_UMIs_before_10x_filter,
        'Total UMIs after Cellbender filter':Total_UMIs_after_cellbender_filter,
        'UMIS mapped to mitochondrial genes':UMIS_mapped_to_mitochondrial_genes,
        'UMIs mapped to genes':UMIs_mapped_to_genes,
        'UMIS mapped to ribo genes':UMIS_mapped_to_ribo_genes,
        'UMIS mapped to ribo rna':UMIS_mapped_to_ribo_rna,
        'Percentage of unassigned cells':Percentage_of_unassigned_cells,
        'Droplets before filtering':Total_Drroplets_before_10x_filtering,
        'Empty droplets - removed by filtering':Droplets_removed_by_filtering,
        'Total Droplets with a single cell':Number_of_cells,
        'Cells before QC filters':Cells_before_QC_filters,
        'Total Cells failing QC':Cells_failing_QC,
        'Total Cells passing QC':Cells_passing_QC,
        'Median cells passes qc':Median_cells_passes_qc,
        'Median cells fails qc':Median_cells_fails_qc,
        'Median Nr cells for donor':Median_Nr_cells_for_donor,
        'Stdev cells passes qc':Stdev_cells_passes_qc,
        'Stdev cells fails qc':Stdev_cells_fails_qc,
        'Stdev Nr cells for donor':Stdev_Nr_cells_for_donor,
        'Total UMIs after cellbender':Total_UMIs_after_cellbender,
        'QC Report end date':date_now
    }
    return fctr, data_tranche, data_donor,azt

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
    parser.add_argument("--extra_metadata", required = False,
                    help="The input folder used",
                    dest='extra_metadata')
    parser.add_argument("--cellbender", required = True,
                help="Cellbender paths if not cellranger filtering used",
                dest='cellbender')
    parser.add_argument("--resolution", required = True,
                help="Cellbender resolution used",
                dest='resolution')
    return parser.parse_args()

if __name__ == '__main__':
    args = set_argument_parser()

    if args.outdir != os.curdir and not os.access(args.outdir, os.F_OK):
        os.mkdir(args.outdir)
        os.mkdir(f"{args.outdir}_summary")

    oufh = open(os.path.join(args.outdir, "files.tsv"), 'w')
    oufh.write("experiment_id\tdonor_id\tfilename_h5ad\tfilename_annotation_tsv\n")
    df_raw = pandas.read_table(args.input_table, index_col = 'experiment_id')
    if (args.cellbender)=='cellranger':
        # here we do not use cellbender and go with default cellranger
        df_cellbender = None
    elif (args.cellbender)=='cellbender':
        # here we have run the cellbender as par of pipeline. 
        file_path = glob.glob(f'{args.results_dir}/nf-preprocessing/cellbender/qc_cluster_input_files/*{args.resolution}*')[0]
        df_cellbender = pandas.read_table(file_path, index_col = 'experiment_id')
    else:
        # this is existing_cellbender, hence using this input
        df_cellbender = pandas.read_table(f'{args.cellbender}', index_col = 'experiment_id')

    adqc = anndata.read_h5ad(f'{args.results_dir}/adata.h5ad')
    fctr = 0
    data_tranche_all=[]
    data_donor_all=[]
    count = 1
    All_probs_and_celltypes = pd.DataFrame()


    args.extra_metadata = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/ELGH_fech/results/yascp_inputs/Extra_Metadata.tsv'   
    Sample_metadata = pd.DataFrame()
    try:
        if args.extra_metadata:
            Extra_Metadata = pd.read_csv('/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/ELGH_fech/results/yascp_inputs/Extra_Metadata.tsv',sep='\t')
            # Extra_Metadata = pd.read_csv(args.extra_metadata,sep='\t')
            Extra_Metadata=Extra_Metadata.set_index('sanger_sample_id')
            Mappings_between_sanger_sampe_and_NS = Extra_Metadata['public_name']
            Library_IDs = pd.read_csv('/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/ELGH_fech/results/yascp_inputs/Library_IDs.csv',sep='\t')
            Sample_Manifest = pd.read_csv('/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/ELGH_fech/results/yascp_inputs/Sample_Manifest.csv',sep='\t')
            Sample_Manifest=Sample_Manifest.set_index('Barcode (scan)')
            Library_IDs = Library_IDs.set_index('S2-046 ID')
            
            for Lid in Mappings_between_sanger_sampe_and_NS.iteritems():
                print(Lid)
                Samples = Library_IDs[Library_IDs['library ID:'].str.contains(Lid[1])]
                Samples['Sanger_id']=Lid[0]
                Samples['NS_ID']=Lid[1]

                intersect = set(Sample_Manifest.index).intersection(set(Samples.index))
                missing = set(Samples.index)- set(Sample_Manifest.index)

                # now combine the missing and existing values.
                All_Vals = pd.DataFrame()
                if len(intersect)>0:
                    Recieved_sample_info = Sample_Manifest.loc[intersect]
                    Recieved_sample_info = Recieved_sample_info[['Volume (ul)','Issue (optional)','Date received']]
                    All_Vals = pd.concat([All_Vals,Recieved_sample_info],axis=0)
                    
                if len(missing)>0:
                    # samples not recorded in sampe reception.
                    missing = intersect
                    Recieved_sample_info_missing = pd.DataFrame(columns={'Volume (ul)','Issue (optional)','Date received'},index=[list(missing)])
                    All_Vals = pd.concat([All_Vals,Recieved_sample_info_missing],axis=0)
                # 
                S1 = Samples.join(Recieved_sample_info)
                # Sample_metadata = Sample_metadata.merge(Samples)
                Sample_metadata=pd.concat([Sample_metadata,S1])
                print('Done')
                # set(Samples.index)
            # Read in sample manifest
            # Read in Library IDs
            print('Done')
            Sample_metadata2= Sample_metadata[['NS_ID','Volume (ul)','Issue (optional)','Date received','Time blood samples taken','Date and time of last meal:','Sanger_id']]
            Sample_metadata2.Sanger_id in df_raw.index
            Sample_metadata2 = Sample_metadata2.reset_index()
            S2 = Sample_metadata2.set_index('Sanger_id')
            S3 = S2.loc[set(df_raw.index)]
            df_raw.index
            S3.to_csv('Sample_data3.tsv',sep='\t')

            
    except:
        print('not working')

    


    for expid in df_raw.index:
        nf, data_tranche, data_donor, azt = gather_pool(expid, args, df_raw, df_cellbender, adqc, oufh = oufh, lane_id=count)
        # add the stuff to the adata.
        azt=azt.set_index('mangled_cell_id')
        All_probs_and_celltypes=pd.concat([All_probs_and_celltypes,azt])
        data_tranche_all.append(data_tranche)
        data_donor_all= data_donor_all+data_donor
        count += 1
        fctr += nf
    
    Donor_Report = pd.DataFrame(data_donor_all)
    Tranche_Report = pd.DataFrame(data_tranche_all)




    Donor_Report['Donor id_Experiment ID']=Donor_Report['Donor id']+'_'+Donor_Report['Experiment ID']
    Donor_Report=Donor_Report.set_index('Donor id_Experiment ID')
    Donor_Report.to_csv(f'{args.outdir}_summary/Donor_Report.tsv',sep='\t')
    Tranche_Report.to_csv(f'{args.outdir}_summary/Tranche_Report.tsv',sep='\t',index=False)
    oufh.close()
    sys.stderr.write("# wrote {:d} files to directory {:s}\n".format(fctr, args.outdir))
    exit(0)
