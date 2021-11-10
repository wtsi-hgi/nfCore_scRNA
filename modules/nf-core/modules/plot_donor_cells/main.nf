
process PLOT_DONOR_CELLS {

    tag "${sample_donor_summary_tsv}"
    
    label 'process_low'
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    publishDir "${params.outdir}/plots/", mode: "${params.plot_donor_ncells.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.indexOf(".pdf") > 0 ? filename.replaceFirst("outputs/","") : "$filename"}
    
    when: 
    params.plot_donor_ncells.run

    input: 
    path(sample_donor_summary_tsv)

    output: 
    path("outputs/*.pdf"), emit: sample_pdf

    script:
    """
        python plot_donor_ncells.py \\
        --output_dir \$PWD/outputs \\
        --sample_donor_summary_tsv ${sample_donor_summary_tsv} \\
        --plotnine_dpi ${params.plot_donor_ncells.plotnine_dpi}
    """
}
