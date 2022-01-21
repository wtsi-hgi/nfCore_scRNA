def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process OUTLIER_FILTER {
    // Takes annData object, plots predicted outlier cells
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    tag "${samplename}"
    
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    
    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__cells_filtered)
        val(metadata_columns)
        val(method)
        val(outliers_fraction)
        val(max_samples)
        val(anndata_compression_opts)

    output:
        path("${runid}-adata.h5ad", emit: anndata)
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Append run_id to output file.
        outfile = "${runid}-adata"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "filter_outlier_cells: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        0026-filter_outlier_cells.py \
            --h5_anndata ${file__anndata} \
            --cell_filtered_per_experiment_file ${file__cells_filtered} \
            --outliers_fraction 0 \
            --metadata_columns ${metadata_columns} \
            --cell_qc_column cell_passes_qc \
            --method ${method} \
            --outliers_fraction ${outliers_fraction} \
            --max_samples ${max_samples} \
            --output_file ${outfile} \
            --anndata_compression_opts ${anndata_compression_opts}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
