process CAPTURE_VIREO{
  label 'process_tiny'
  input:
    path(vireo_location)
   
  output:
    // tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool
    path("output_vireo.csv"),emit:vireo_loc
  script:
  """
    for OUTPUT in \$(ls -d ${vireo_location}/*/)
    do
    samplename1=\$(echo \$OUTPUT | sed 's/${vireo_location}\\///g')
    samplename1=\${samplename1:0:-1}
    echo "\$samplename1 \$PWD/${vireo_location}/\$samplename1/\${samplename1}_headfix_vireo.vcf.gz \$PWD/${vireo_location}/\$samplename1/\${samplename1}_headfix_vireo.vcf.gz.tbi" >> output_vireo.csv
    done
  """    
}



process VIREO {
    tag "${samplename}"
    label 'process_high'
    publishDir "${params.outdir}/deconvolution/vireo/${samplename}/",  mode: "${params.vireo.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("vireo_${samplename}/","") }



    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

     when:
      params.vireo.run

    input:
      tuple val(samplename), path(cell_data), val(n_pooled), path(donors_gt_vcf), path(donor_gt_csi)

    output:
      tuple val(samplename), path("vireo_${samplename}/*"), emit: output_dir
      tuple val(samplename), path("vireo_${samplename}/donor_ids.tsv"), emit: sample_donor_ids
      tuple val(samplename), path("vireo_${samplename}/GT_donors.vireo.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
      tuple val(samplename), path("vireo_${samplename}/GT_donors.vireo.vcf.gz"), emit: infered_vcf
      path("vireo_${samplename}/${samplename}.sample_summary.txt"), emit: sample_summary_tsv
      path("vireo_${samplename}/${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv
      tuple  val(samplename), path("vireo_${samplename}/GT_donors.vireo.vcf.gz"), path("vireo_${samplename}/${samplename}.sample_summary.txt"),path("vireo_${samplename}/${samplename}__exp.sample_summary.txt"),path("vireo_${samplename}/donor_ids.tsv"),path(vcf_file),path(donor_gt_csi), emit: all_required_data

    script:
      vcf_file = ""
      if (params.genotype_input.vireo_with_gt){
        vcf = " -d ${donors_gt_vcf} --forceLearnGT"
        vcf_file = donors_gt_vcf
        com2 = "cd vireo_${samplename} && ln -s ../${donors_gt_vcf} GT_donors.vireo.vcf.gz"
        com2 = ""
      }else{
         vcf = ""
         vcf_file = donors_gt_vcf
         com2 = ""
      }

    """

      umask 2 # make files group_writable

      
      vireo -c $cell_data -N $n_pooled -o vireo_${samplename} ${vcf} -t GT --randSeed 1 -p $task.cpus --nInit 200
      # add samplename to summary.tsv,
      # to then have Nextflow concat summary.tsv of all samples into a single file:
      gzip vireo_${samplename}/GT_donors.vireo.vcf || echo 'vireo_${samplename}/GT_donors.vireo.vcf already gzip'
      cat vireo_${samplename}/summary.tsv | \\
        tail -n +2 | \\
        sed s\"/^/${samplename}\\t/\"g > vireo_${samplename}/${samplename}.sample_summary.txt

      cat vireo_${samplename}/summary.tsv | \\
        tail -n +2 | \\
        sed s\"/^/${samplename}__/\"g > vireo_${samplename}/${samplename}__exp.sample_summary.txt
      ${com2}
    """
}
