#!/usr/bin/env nextflow

// NOTE: label 'big_mem' may be useful at some point


def random_hex(n) {
  Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
  echo_mode = true
}


// process AMBIENTNESS_QUANTIFICATION{
//     label 'process_medium'
  
//     if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//       container "${params.nf_cellbender_container}"
//       maxRetries = 1
//     } else {
//       container "wtsihgi/nf_cellbender_container:3cc9983"
//     }

// }


process cellbender__rb__get_input_cells {

  label 'process_low'
  
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.nf_cellbender_container}"
    //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/wtsihgi_nf_cellbender_v1.2.img"
    maxRetries = 1
    // workdir /tmp

    
  } else {
    container "wtsihgi/nf_cellbender_container:3cc9983"
  }

  // Calculates thresholds for input cells of cellbender__remove_background
  // ------------------------------------------------------------------------
 // use tmp directory

  publishDir  "${outdir}",
      mode: "${params.cellsnp.copy_mode}",
      overwrite: "true"

  input:
    val(outdir_prev)
    tuple(
    val(experiment_id),
    path(file_10x_barcodes),
    path(file_10x_features),
    path(file_10x_matrix),val(ncells)
    )
    val(estimate_params_umis)
    

  output:
    val(outdir, emit: outdir)
    tuple(
    val(experiment_id),
    path(file_10x_barcodes),
    path(file_10x_features),
    path(file_10x_matrix),
    path("${outfile}-expected_cells.txt"),
    path("${outfile}-total_droplets_included.txt"),
    emit: cb_input
    )
    path(
    "${outfile}-expected_cells.txt",
    emit: expected_cells
    )
    path(
    "${outfile}-total_droplets_included.txt",
    emit: total_droplets_include
    )
    path("${outfile}-cell_estimate_cutoff.tsv.gz")
    path("${outfile}-total_droplets_cutoff.tsv.gz")
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

  script:

    outdir = "${outdir_prev}/${experiment_id}"
    outdir = "${outdir}/cellbender-estimate_ncells_nemptydroplets"
    outfile = "umi_count_estimates"
    cmd__expected_ncells = ""

    cmd__droplets_include = ""

    if ("${outdir}" == "${outdir_prev}/${experiment_id}/cellbender-params_setbyuser") {
    outdir = "${outdir_prev}/${experiment_id}/cellbender-estimated_ncells_nemptydroplets"
    }

    if (params.cellbender_rb.estimate_params_umis.value.method_estimate_ncells=='expected'){
        ncells2=ncells.toInteger()-1000
        cell_numbers = "--expected_ncells ${ncells2}"
    }else{
        cell_numbers = ""
    }

    """

    rm -fr plots
    mkdir txd_input
    ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

    get_estimates_from_umi_counts.py ${cell_numbers} --tenxdata_path txd_input --output_file ${outfile} --expected_nemptydroplets_umi_cutoff   ${estimate_params_umis.expected_nemptydroplets_umi_cutoff} --method_estimate_ncells   ${estimate_params_umis.method_estimate_ncells} --lower_bound_umis_estimate_ncells   ${estimate_params_umis.lower_bound_umis_estimate_ncells} --method_estimate_nemptydroplets   ${estimate_params_umis.method_estimate_nemptydroplets} --lower_bound_umis_estimate_nemptydroplets   ${estimate_params_umis.lower_bound_umis_estimate_nemptydroplets} --upper_bound_umis_estimate_nemptydroplets   ${estimate_params_umis.upper_bound_umis_estimate_nemptydroplets} --estimate_nemptydroplets_add_umifactor   ${estimate_params_umis.estimate_nemptydroplets_umi_add_factor} --estimate_nemptydroplets_subtract_dropletfactor   ${estimate_params_umis.estimate_nemptydroplets_subtract_cell_factor} --estimate_nemptydroplets_min_nemptydroplets ${estimate_params_umis.estimate_nemptydroplets_min_drop}    ${cmd__expected_ncells} ${cmd__droplets_include}

    mkdir plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true
    """
}


process cellbender__preprocess_output{
    label 'process_low'
    tag "${experiment_id}_cb"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.nf_cellbender_container}"
      //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/wtsihgi_nf_cellbender_v1.2.img"
      maxRetries = 1
      // memory = 250.GB
      cpus = 1
    } else {
      container "wtsihgi/nf_cellbender_container:3cc9983"
    }
    publishDir  path: "${outdir}",
        saveAs: {filename ->
          if (filename == "barcodes.tsv.gz") {
          null
          } else if(filename.equalsIgnoreCase("features.tsv.gz")){
            null
          }
          else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
          null
          } else {
          filename.replaceAll("-", "")
          }
        },
        mode: "${params.cellsnp.copy_mode}",
        overwrite: "true" 
 
  input:

    tuple(
      val(experiment_id),
      val(fpr),
      val(outdir),
      path(filtered),
      val(cb_params)
    )

    tuple(
      val(outdir),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix)
      
    )

    tuple(
      val(experiment_id),
      val(outdir),
      path(expected_cells),
      path(total_droplets_include)
    )



  output:
    tuple(val(experiment_id),path("cellbender-FPR_${params.cellbender_resolution_to_use}-filtered_10x_mtx"), emit: alternative_input)
    tuple(val(experiment_id),path("cellbender-FPR_${params.cellbender_resolution_to_use}-unfiltered_10x_mtx"), emit: alternative_input_raw)
    path("*filtered_10x_mtx/barcodes.tsv.gz", emit: tenx_barcodes)
    path("*filtered_10x_mtx/features.tsv.gz", emit: tenx_features)
    path("*filtered_10x_mtx/matrix.mtx.gz", emit: tenx_matrix)
    path("Warnings.log") optional true
    path(
      "${outfile}-filtered_10x_mtx-file_list.tsv",
       emit: results_list
    )  
    tuple(
      val(experiment_id),
      val(outdir),
      path("*_unfiltered.h5"),
      path(expected_cells),
      path(total_droplets_include),
      emit: experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
    )
    tuple(
      val(outdir),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      path("*_filtered.h5"),
      emit: cb_plot_input
    )

    tuple(
      val(experiment_id),
      val(outdir),
      emit: out_paths
    )

  script:
    // 032-clean_cellbender_results.py --nf_outdir_tag results/nf-preprocessing/cellbender/Card_Val11650914/cellbender-epochs_250__learnrt_1Eneg7__zdim_100__zlayer_500__lowcount_10 --cb_outfile_tag cellbender --experiment_id Card_Val11650914 --fpr '0.01 0.05 0.1' --cb_params cellbender_params-epochs_250__learnrt_1Eneg7__zdim_100__zlayer_500__lowcount_10
    // cp cellbender-filtered_10x_mtx-file_list.tsv AA868D1A8EB8249-cellbender-filtered_10x_mtx-file_list.tsv

    outfile = "cellbender"
    """
      export LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so
      for i in \$(ls *.h5); do
        echo \$i
        out_file=\$(echo \$i | sed s/"\\."/"pt"/g | sed s/"pth5"/"\\.h5"/)
        # If outfile does not exist, move i to out_file
        [ ! -f \$out_file ] && mv \$i \$out_file
        echo \$out_file
      done
      032-clean_cellbender_results.py --nf_outdir_tag ${outdir} --cb_outfile_tag ${outfile} --experiment_id ${experiment_id} --fpr '${fpr}' --cb_params ${cb_params}
      cp ${outfile}-filtered_10x_mtx-file_list.tsv ${outfile}-filtered_10x_mtx-file_list.tsv || echo 'same file'
    """

}

process cellbender__remove_background {
  // Remove ambient RNA
  // ------------------------------------------------------------------------
  //cache false    // cache results from run
  //maxForks 2   // hard to control memory usage. limit to 3 concurrent
// cb_plot_input,out_paths,results_list,experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude

  tag "${experiment_id}_cb"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    if (params.cellbender_v == '0.3.1'){
      container "${params.nf_cellbender_container_032}"
    }else{
      container "${params.nf_cellbender_container}"
    }
    
    maxRetries = 1
	  cpus = 1

  } else {
    if (params.cellbender_v == '0.3.1'){
      container "us.gcr.io/broad-dsde-methods/cellbender:0.3.1"
    }else{
      container "wtsihgi/nf_cellbender_container:3cc9983"
    }
    
  }
    
  // set LD_PRELOAD to fix mkl / anaconda conflict
  // cf. https://stackoverflow.com/questions/36659453/intel-mkl-fatal-error-cannot-load-libmkl-avx2-so-or-libmkl-def-so
  // containerOptions '--env LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so'

  //     // use GPU
  if (params.utilise_gpu){
    label 'gpu'

	// only one label here, otherwise bsub -R -M will be doubled..
	// label 'process_high_memory'
    gpu_text_info = '--cuda'
  }else{
    label 'process_medium'
  }

  // scratch false    // use tmp directory
  // echo false   // echo output from script

  publishDir  path: "${outdir}",
      saveAs: {filename ->
        if (filename == "barcodes.tsv.gz") {
        null
        } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
        null
        } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
        null
        } else {
        filename.replaceAll("-", "")
        }
      },
      mode: "${params.cellsnp.copy_mode}",
      overwrite: "true"

  input:
    val(outdir_prev)
    tuple(
      val(experiment_id),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      path(expected_cells),
      path(total_droplets_include),
      val(low_count_threshold),
      val(epochs),
      val(learning_rate),
      val(zdims),
      val(zlayers),
    )
    // each epochs
    // each learning_rate
    // each zdims
    // each zlayers
    // each low_count_threshold
    each fpr

  output:
    val(outdir, emit: outdir)
    tuple(
      val(experiment_id),
      val(fpr),
      val(outdir),
      path("*.h5"),val(cb_params),
      emit: cleanup_input
    )
    tuple(val(experiment_id),path("cellbender_FPR_${params.cellbender_resolution_to_use.replaceAll('pt', '.')}_filtered.h5"),emit: cb_to_use_downstream)

    tuple(
      val(outdir),
      val(experiment_id),
      val(epochs),
      val(learning_rate),
      val(zdims),
      val(zlayers),
      val(low_count_threshold),
      val(fpr),
      emit: cb_results_details
    )

    tuple(
      val(outdir),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      emit: cb_plot_input
    )

    path("${outfile}_cell_barcodes.csv", emit: barcodes)
    
    path("*.h5", emit: filtered_h5s)
    path(file_10x_barcodes, emit: raw_tenx_barcodes)
    path(file_10x_features, emit: raw_tenx_features)
    path(file_10x_matrix, emit: raw_tenx_matrix)
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

    tuple(
      val(experiment_id),
      val(outdir),
      
      path(expected_cells),
      path(total_droplets_include),
      emit: experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
    )

    tuple(
      val(experiment_id),
      val(outdir),
      emit: out_paths
    )

  script:

    if (params.utilise_gpu){
      gpu_text_info = '--cuda'
    }else{
      gpu_text_info = "--cpu-threads ${task.cpus}"
    }

    if (params.cellbender_v == '0.3.1'){
      option1='--checkpoint-mins 100'
    }else{
      option1=''
    }






    outdir = "${outdir_prev}/${experiment_id}"
    lr_string = "${learning_rate}".replaceAll("\\.", "pt")
    lr_string = "${lr_string}".replaceAll("-", "neg")
    fpr_string = "${fpr}".replaceAll("\\.", "pt").replaceAll(" ", "_")
    cb_params = "cellbender_params"
    cb_params = "${cb_params}-epochs_${epochs}"
    cb_params = "${cb_params}__learnrt_${lr_string}"
    cb_params = "${cb_params}__zdim_${zdims}"
    cb_params = "${cb_params}__zlayer_${zlayers}"
    cb_params = "${cb_params}__lowcount_${low_count_threshold}"
    outdir = "${outdir}/${cb_params}".replaceAll("cellbender_params","cellbender")
    outfile = "cellbender"
    
    """
    echo ${experiment_id}
    # LD_PRELOAD to fix mkl/anaconda python error
    # cf. https://stackoverflow.com/questions/36659453/intel-mkl-fatal-error-cannot-load-libmkl-avx2-so-or-libmkl-def-so
    export LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so
    #// nvidia-smi
    rm -fr plots
    mkdir txd_input
    ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz
    export TMPDIR=\$PWD
    cellbender remove-background --input txd_input ${gpu_text_info} ${option1} --output ${outfile} --expected-cells \$(cat ${expected_cells}) --total-droplets-included \$(cat ${total_droplets_include}) --model full --z-dim ${zdims} --z-layers ${zlayers} --low-count-threshold ${low_count_threshold} --epochs ${epochs} --learning-rate ${learning_rate} --fpr ${fpr}
    # If outfile does not have h5 appended to it, move it.
    [ -f ${outfile} ] && mv ${outfile} ${outfile}.h5

    mkdir plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true
    """
}

process cellbender__remove_background__qc_plots {
  label 'process_low'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.nf_scrna_qc_sif_container}"
    //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/wtsihgi_nf_cellbender_v1.2.img"
  } else {
    container "wtsihgi/nf_scrna_qc:6bb6af5"
  }
  

  // QC plots from cellbdender
  // ------------------------------------------------------------------------
  //tag { output_dir }
  //cache false    // cache results from run
  //maxForks 2   // hard to control memory usage. limit to 3 concurrent
  // scratch false    // use tmp directory
  // echo echo_mode   // echo output from script

  publishDir  path: "${outdir}",
      saveAs: {filename ->
        if (filename.equalsIgnoreCase("barcodes.tsv.gz")) {
        null
        } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
        null
        } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
        null
        } else {
        filename.replaceAll("-", "")
        }
      },
      mode: "${params.copy_mode}",
      overwrite: "true"

  input:
    tuple(
    val(outdir_prev),
    path(file_10x_barcodes),
    path(file_10x_features),
    path(file_10x_matrix),
    file(h5_filtered_cellbender)
    )

  output:
    val(outdir, emit: outdir)
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

  script:
    outdir = "${outdir_prev}" // /${experiment_id}"
    h5_filtered_cellbender = h5_filtered_cellbender.join(",")
    """
    echo "outdir: ${outdir}"
    mkdir -p txd_input
    ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz
    # Make a file with list of our files
    echo "${h5_filtered_cellbender}" | sed s/","/"\\n"/g > files.txt
    for i in \$(cat files.txt); do
    echo \$i
    out_file=\$(echo \$i | sed s/".h5"//)
    035-analyse_cellbender_results.py   --tenxdata_path txd_input   --h5_cellbender \$i   --output_file cellbender_results-\$out_file   --number_cpu ${task.cpus}
    done
    rm files.txt
    mkdir -p plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true
    """
}

process capture_cellbender_files{

  publishDir  path: "${outdir}/cellbender",
  saveAs: {filename ->
        if (filename.contains("captured")) {
          null
        } else if (filename.contains(".h5")){
          null
        }
        else {
          filename.replaceAll("tmp1234/cellbender/", "")
        }
      },
        mode: "${params.copy_mode}",
    overwrite: "true"
  label 'process_tiny'



  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.nf_scrna_qc_sif_container}"
    //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/wtsihgi_nf_cellbender_v1.2.img"
  } else {
    container "wtsihgi/nf_scrna_qc:6bb6af5"
  }

  // cache false
  input:
    path(cellbender_location) //cant use path here - needs a full path
    val(outdir)
    path(input_file)
  output:
    path("tmp1234/cellbender/*") optional true
    path("captured/*/*FPR_${params.cellbender_resolution_to_use}*"),emit:alt_input optional true
    path("cellbender/*/*/cellbender_FPR_${params.cellbender_resolution_to_use.replaceAll('pt', '.')}_filtered.h5"),emit:cb_to_use_downstream optional true
    // cellbender/*/*/cellbender_FPR_0pt1_filtered.h5
    // tuple(val(experiment_id),path("cellbender_FPR_${params.cellbender_rb_resolution_to_use.replaceAll('pt', '.')}_filtered.h5"),emit: cb_to_use_downstream)

    
  script:
  """
  
  
    mkdir tmp1234
    mkdir tmp1234/cellbender
    cd tmp1234/cellbender
    for d in ../../${cellbender_location}/* ; do
        echo \$d
        basename1=\$(basename \$d)
        if grep -q \$basename1 ../../${input_file}
        then
          ln -s \$d .
        else
            f='no'
        fi        
    done
    cd ../..
    capture_res_files_cb.py -res ${params.cellbender_resolution_to_use}
  """    

}

process cellbender__remove_background__qc_plots_2 {

  label 'process_low'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_sif_container}"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
  }

  // Second set of QC plots from cellbdender.
  // This task compare the cellbender output with both the cellranger filtered and cellragner raw outputs
  // ------------------------------------------------------------------------
  tag { "$experiment_id" }
  publishDir "${outdir2}/$experiment_id/compare_cellranger/", pattern: "fpr_${fpr}/${experiment_id}/*.png", 
    saveAs: {filename ->
    filename.replaceAll("fpr_${fpr}/${experiment_id}/", "fpr_${fpr}/")
  },
    mode: "${params.cellsnp.copy_mode}",
    overwrite: "true"
  
  //cache false    // cache results from run
  //maxForks 2   // hard to control memory usage. limit to 3 concurrent
  // scratch false    // use tmp directory
  // echo echo_mode   // echo output from script
  input:
    tuple val(experiment_id), val(outdir), path(cellbender_unfiltered_h5s), path(expectedcells), path(totaldropletsinclude), path(raw_cellranger_mtx), path(filtered_cellranger_mtx), val(fpr)
    val(outdir2)
  output:
    val(outdir, emit: outdir)
    path("fpr_${fpr}/${experiment_id}/*.png"), emit: plots_png 

  script:
  """
  fprid=\$(echo $fpr | sed s'/\\./pt/'g)
  cellbender_unfiltered_h5=cellbender_FPR_\${fprid}_unfiltered.h5

  n_expected_cells=\$(cat $expectedcells)
  n_total_droplets_included=\$(cat $totaldropletsinclude)

  037-plot_cellranger_vs_cellbender.py \\
    --samplename \"${experiment_id}\" \\
    --raw_cellranger_mtx \"${raw_cellranger_mtx}\" \\
    --filtered_cellranger_mtx \"${filtered_cellranger_mtx}\" \\
    --cellbender_unfiltered_h5 \"\$cellbender_unfiltered_h5\" \\
    --fpr \"${fpr}\" \\
    --n_expected_cells \"\${n_expected_cells}\" \\
    --n_total_droplets_included \"\${n_total_droplets_included}\" \\
    --out_dir \$PWD
  """
}


process cellbender__gather_qc_input {


  label 'process_low'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.nf_scrna_qc_sif_container}"
    //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
  } else {
    container "wtsihgi/nf_cellbender_container:3cc9983"
  }

  // Prepare cell bender output for qc_cluster pipeline. For each epoch and
  // learning rate, and fpr, gather 10x matrix files into format for
  // qc_cluster.
  // ------------------------------------------------------------------------

  publishDir  path: "${outdir}",
      saveAs: {filename -> filename.replaceAll("-", "")},
      mode: "${params.cellsnp.copy_mode}",
      overwrite: "true"

  input:
    val(outdir_prev)
    file(cb_results_tsvs)

  output:
    val(outdir, emit: outdir)
    path("*.tsv", emit: qc_input_files)
    path("file_paths_10x-*${params.cellbender_resolution_to_use}.tsv", emit: celbender_path)
  script:
    
    outdir = "${outdir_prev}/qc_cluster_input_files"
    cb_results_tsvs = cb_results_tsvs.join(",")

    """
        045-prepare_nf_qc_cluster_input.py --cb_results_tsvs ${cb_results_tsvs}
    """
}
