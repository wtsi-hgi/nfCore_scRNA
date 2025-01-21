
include {AZIMUTH} from "$projectDir/modules/nf-core/modules/azimuth/main"
include {CELLTYPIST} from "$projectDir/modules/nf-core/modules/celltypist/main"
include {SPLIT_BATCH_H5AD} from "$projectDir/modules/nf-core/modules/split_batch_h5ad/main"
include {KERAS_CELLTYPE} from "$projectDir/modules/nf-core/modules/keras_celltype/main"
include {PSEUDOBULK_AGREGATION_PER_CELLTYPE} from "$projectDir/modules/nf-core/modules/pseudobulk_aggregation_per_celltype/main"
include {CELLTYPE_FILE_MERGE} from './functions'

workflow CELL_TYPE_ASSIGNEMT{
    
    take:
        file__anndata_merged
        file__cells_filtered
        
    main:
        // if (params.split_ad_per_bach){
        log.info '---Splitting the assignment for each batch---'
        file__anndata_merged3 = file__anndata_merged
        SPLIT_BATCH_H5AD(file__anndata_merged,params.split_ad_per_bach)
        // SPLIT_BATCH_H5AD.out.sample_file.view()
        // Here we may want to not split it and just pass in an entire h5ad file for annotations.
        SPLIT_BATCH_H5AD.out.sample_file
            .splitCsv(header: true, sep: "\t", by: 1)
            .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{ch_experiment_filth5}
        SPLIT_BATCH_H5AD.out.az_sample_file
            .splitCsv(header: true, sep: "\t", by: 1)
            .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{az_ch_experiment_filth5}

        SPLIT_BATCH_H5AD.out.files_anndata_batch.flatMap().set{ch_batch_files}
        if (params.celltype_assignment.run_keras){
            KERAS_CELLTYPE(ch_experiment_filth5) 
            all_extra_fields = KERAS_CELLTYPE.out.predicted_celltype_labels.collect()
        }else{
            all_extra_fields = Channel.of()
        }
        

        if (params.celltype_assignment.run_azimuth){
            AZIMUTH(params.outdir,ch_batch_files)
            az_out = AZIMUTH.out.predicted_celltype_labels.collect()
        }else{
            az_out = Channel.of("$projectDir/assets/fake_file2.fq")
        }

        if (params.celltype_assignment.run_celltypist){
            Channel.fromList(params.celltypist.models)
                .set{ch_celltypist_models}
            CELLTYPIST(az_ch_experiment_filth5.combine(ch_celltypist_models))
            ct_out = CELLTYPIST.out.sample_predicted_labels_csv.collect()
        }else{
            ct_out = Channel.of()
        }
        
        CELLTYPE_FILE_MERGE(az_out,ct_out,all_extra_fields,file__anndata_merged)
        
        file__anndata_merged2=CELLTYPE_FILE_MERGE.out.file__anndata_merged2
        PSEUDOBULK_AGREGATION_PER_CELLTYPE(file__anndata_merged2,CELLTYPE_FILE_MERGE.out.celltype_assignments)

    emit:
        file__anndata_merged2
}