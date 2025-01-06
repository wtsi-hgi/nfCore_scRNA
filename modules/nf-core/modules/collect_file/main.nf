process collect_file{
  label 'process_tiny'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.nf_scrna_qc_sif_container}"
  } else {
    container "wtsihgi/nf_scrna_qc:6bb6af5"
  }
    // In nf there is a function collectFile - however if you provide a symlinked file directory tusing nf function will overwrite the source instead of replacing the file
    // This snipped is a replication of the function but as a nf module and hence the problem is avoided.

  publishDir  "${outpath2}/",
              saveAs: {filename ->
                  if ("${outpath}" == "0") {
                      null
                  }else {
                      filename
                  }
              },
              mode: "${params.copy_mode}",
              overwrite: "true"

  input:
    path(files_to_concentrate)
    val(name)
    val(outpath)
    val(header_size)
    val(seed)
   
  output:
    path("${name}"),emit:output_collection

  
  script:

    if ("${outpath}" == "0") {
        outpath2="${params.outdir}"
    }else {
        outpath2=outpath
    }
    file1 = files_to_concentrate[0]
    files= files_to_concentrate.join(" ")
    if (seed!='' ){
        seed_inset = "-seed '${seed}'"
    }else{
        seed_inset = ""
    }
    if (header_size>0){
        header_insert = "head -n ${header_size} ${file1} > ${name}"
        
    }else{
        header_insert = ""
    }

    header_insert

    """
        combine.py -d '${files}' -o ${name} ${seed_inset}
    """    
}