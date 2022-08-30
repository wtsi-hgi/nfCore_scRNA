process CREATE_ARTIFICIAL_BAM_CHANNEL_process{
    
    // This will be a csv file that will be later converted to the channel as per data emmited from main_deconvolution channel
    scratch false      // use tmp directory
    label 'process_low'
    input:
        path(outdir)
        path(fech_folder)

    output:
        val(outdir, emit: ch_experiment_bam_bai_barcodes)

    script:
        """
            echo 'lets do this'
        """
        
}



workflow CREATE_ARTIFICIAL_BAM_CHANNEL{
    take:
        input_channel
    main:
        // CREATE_ARTIFICIAL_BAM_CHANNEL_process(params.outdir,input_channel)
        // ch_experiment_bam_bai_barcodes=CREATE_ARTIFICIAL_BAM_CHANNEL_process.out.ch_experiment_bam_bai_barcodes

        input_channel
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam","${params.outdir}/deconvolution/vireo/${row.experiment_id}/donor_ids.tsv")}
            .set{ch_experiment_bam_bai_barcodes}
        ch_experiment_bam_bai_barcodes.view()
    // ch_experiment_bam_bai_barcodes
    //               .map { samplename, bam_file, bai_file, barcodes_tsv_gz -> tuple(samplename, file(bam_file)) }
    //               .combine(vireo_out_sample_donor_ids, by: 0 )
    //               .set { ch_experiment_bam_vireo_donor_ids }
    emit:
        ch_experiment_bam_bai_barcodes

}