process TOTAL_VI_INTEGRATION{
    
    if (params.utilise_gpu){
        label 'process_low'
    }else{
        label 'process_medium'
    }
    memory { 
            sizeInGB = adata.size() / 1e9 * 2 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }

    publishDir  path: "${outdir_prev}/totalVi",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.total_vi_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(adata)
        path(citedata)
        val(outdir_prev)
        val(colors_quantitative)
        val(colors_categorical)
        val(reduction_columns_cells)
        val(reduction_columns_genes)
        path(gene_list_to_keep)
    output:
        path("./figures"), emit: figs, optional: true
        path("./scvi_model"), emit: scvi_model, optional: true
        path("./totalVI_integrated.h5ad"), emit: totalVI_integrated, optional: true
        
    script:
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative '${colors_quantitative}'"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical '${colors_categorical}'"
        }

        """
            totalVI.py -h5ad_file ${adata} ${cmd__colors_quant} ${cmd__colors_cat} --reduction_columns_cells '${reduction_columns_cells}' --reduction_columns_genes '${reduction_columns_genes}' --gene_list_to_keep ${gene_list_to_keep}
        """
}
