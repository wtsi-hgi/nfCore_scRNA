
include {AZIMUTH} from '../azimuth/main'
include {CELLTYPIST} from '../celltypist/main'
include {SPLIT_BATCH_H5AD} from '../split_batch_h5ad/main'
include {CELLTYPE_FILE_MERGE} from './functions'

workflow CELL_TYPE_ASSIGNEMT{
    
    take:
        file__anndata_merged
        file__cells_filtered
        
    main:
        // if (params.split_ad_per_bach){
            log.info '---Splitting the assignment for each batch---'
            SPLIT_BATCH_H5AD(file__anndata_merged,params.split_ad_per_bach)
            SPLIT_BATCH_H5AD.out.sample_file.view()
            // Here we may want to not split it and just pass in an entire h5ad file for annotations.
            SPLIT_BATCH_H5AD.out.sample_file
                .splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{ch_experiment_filth5}
            SPLIT_BATCH_H5AD.out.files_anndata_batch.flatMap().set{ch_batch_files}
        // }else{
        //     log.info '---Running assignment for full ad---'
        //     PREPERE_H5AD_FOR_CELLTYPE(file__anndata_merged)
        //     ch_batch_files = file__anndata_merged
        //     ch_experiment_filth5 = Channel.from("dummy").map { dummy -> tuple("full_h5ad") }
        //     ch_experiment_filth5=ch_experiment_filth5.combine(file__anndata_merged)
        //     ch_experiment_filth5.view()
        //     ch_batch_files.view()
        // }

        Channel.fromList(params.celltypist.models)
            .set{ch_celltypist_models}
        ch_experiment_filth5.combine(ch_celltypist_models)

        
        AZIMUTH(params.output_dir,ch_batch_files)
        CELLTYPIST(ch_experiment_filth5.combine(ch_celltypist_models))
        CELLTYPE_FILE_MERGE(AZIMUTH.out.predicted_celltype_labels.collect(),CELLTYPIST.out.sample_predicted_labels_csv.collect(),file__anndata_merged)
        file__anndata_merged2=CELLTYPE_FILE_MERGE.out.file__anndata_merged2

    emit:
        file__anndata_merged2
}