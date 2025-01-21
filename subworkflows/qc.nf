
// Load base.config by default for all pipelines - typically included in the nextflow config.
// Modules to include.
include {OUTLIER_FILTER; MERGE_OUTLIER_FILES} from "$projectDir/modules/nf-core/modules/outlier_filter/main"
include {PLOT_STATS} from "$projectDir/modules/nf-core/modules/plot_stats/main"
include {ESTIMATE_PCA_ELBOW} from "$projectDir/modules/nf-core/modules/estimate_pca_elbow/main"
include {SUBSET_PCS} from "$projectDir/modules/nf-core/modules/subset_pcs/main"
include {NORMALISE_AND_PCA; PCA} from "$projectDir/modules/nf-core/modules/normalise_and_pca/main"
include {HARMONY} from "$projectDir/modules/nf-core/modules/harmony/main"
include {BBKNN} from "$projectDir/modules/nf-core/modules/bbknn/main"
include {ADD_EXTRA_METADATA_TO_H5AD} from "$projectDir/modules/nf-core/modules/adata_manipulations/main"
include {LISI} from "$projectDir/modules/nf-core/modules/lisi/main"
include {UMAP; UMAP as UMAP_HARMONY; UMAP as UMAP_BBKNN;} from "$projectDir/modules/nf-core/modules/umap/main"
include {CLUSTERING; CLUSTERING as CLUSTERING_HARMONY; CLUSTERING as CLUSTERING_BBKNN;} from "$projectDir/modules/nf-core/modules/clustering/main"
include {CELL_HARD_FILTERS} from "$projectDir/modules/nf-core/modules/cell_hard_filters/main"
include {DONT_INTEGRATE} from "$projectDir/modules/nf-core/modules/reduce_dims/main"
include {TOTAL_VI_INTEGRATION} from "$projectDir/modules/nf-core/modules/totalVi/main"
include { DSB_PROCESS; PREPROCESS_PROCESS; DSB_INTEGRATE; MULTIMODAL_INTEGRATION; VDJ_INTEGRATION } from '../modules/nf-core/modules/citeseq/main'

workflow qc {
    take:
        file__anndata_merged
        file__cells_filtered
        gt_outlier_input
        channel_dsb
        vireo_paths
        assignments_all_pools
        matched_donors
        chanel_cr_outs
    main:
        log.info "--- Running QC metrics --- "
        // if(params.extra_metadata!=''){
        //     log.info '''--- Adding extra metadata to h5ad---'''
        //     ADD_EXTRA_METADATA_TO_H5AD(file__anndata_merged,params.extra_metadata)
        //     file__anndata_merged = ADD_EXTRA_METADATA_TO_H5AD.out.file__anndata
        // }else{
        //     log.info '''--- No extra metadata to add to h5ad ---'''
        // }

        if (params.cell_hard_filters){
            if(params.sample_qc.cell_filters.experiment.value != '' | params.sample_qc.cell_filters.all_samples.value != '' | params.sample_qc.downsample_cells_fraction.value != '' | params.sample_qc.downsample_cells_n.value != '' | params.sample_qc.downsample_feature_counts.value != ''){
                log.info """---Flagging/filtering hard filters.----"""
                CELL_HARD_FILTERS(file__anndata_merged,params.hard_filters_drop)
                file__anndata_merged = CELL_HARD_FILTERS.out.anndata
            }
        }

        
        //FILTERING OUTLIER CELLS
        if (params.filter_outliers) {
            log.info """---Running automatic outlier cell filtering.----"""
            OUTLIER_FILTER(
                params.outdir,
                file__anndata_merged,
                file__cells_filtered,
                params.sample_qc.cell_filters.filter_outliers.outliers_fraction,
                params.sample_qc.cell_filters.filter_outliers.max_samples,
                params.anndata_compression_opts,
                gt_outlier_input,
                params.sample_qc.gt_match_based_adaptive_qc_exclusion_pattern,
                params.sample_qc.cell_filters.filter_outliers.methods_thresholds
            )
            OUTLIER_FILTER.out.cells_filtered.collect().set{test2}

            MERGE_OUTLIER_FILES(file__anndata_merged,test2 )
            file__anndata_merged = MERGE_OUTLIER_FILES.out.anndata
            // file__cells_filtered = OUTLIER_FILTER.out.cells_filtered
                
        }
        
        
        if (params.normalise_andata){
            log.info """---Normalising data For data clustering and integration.----"""
            NORMALISE_AND_PCA(
                file__anndata_merged,
                params.normalise.mode,
                params.normalise.layer,
                params.normalise.minimum_number_of_cells_for_donor,
                params.normalise.gene_filters.variable_genes_exclude,
                params.normalise.genes_score,
                params.normalise.gene_filters.genes_exclude,
                params.normalise.gene_filters.genes_at_least_in_nr_cells,
                params.reduced_dims.vars_to_regress.value
            )



            if (params.citeseq){
                log.info """---Integrating data using Seurat integration method----"""
                NORMALISE_AND_PCA.out.sample_QCd_adata.flatten().map{sample -> tuple("${sample}".replaceFirst(/___sample_QCd_adata.h5ad/,"").replaceFirst(/.*\//,""),sample)}.set{alt_input}
                channel_dsb2 = channel_dsb.combine(alt_input, by: 0)
                DSB_PROCESS(channel_dsb2)


                if(params.totalVi.run_process){
                    TOTAL_VI_INTEGRATION(NORMALISE_AND_PCA.out.anndata,
                    DSB_PROCESS.out.citeseq_rsd.collect(),
                    NORMALISE_AND_PCA.out.outdir, 
                    params.umap.colors_quantitative.value, 
                    params.umap.colors_categorical.value, 
                    params.totalVi.reduction_columns_cells,
                    params.totalVi.reduction_columns_genes,
                    params.totalVi.gene_list_to_keep)
                }


                DSB_PROCESS.out.ch_for_norm.subscribe { println "1:: DSB_PROCESS.out.ch_for_norm: $it" }
                
                // PREPROCESS_PROCESS()
                // DSB_PROCESS.out.citeseq_rsd.subscribe { println "1:: DSB_PROCESS.out.citeseq_rsd: $it" }
                vireo_paths_map = vireo_paths.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), row)}
                // vireo_paths_map.subscribe { println "1:: vireo_paths_map $it" }
                // DSB_PROCESS.out.ch_for_norm.subscribe { println "1:: vireo_paths_map $it" }
                vireo_paths_map.combine(DSB_PROCESS.out.ch_for_norm, by: 0).set{norm_chanel}
                norm_chanel.combine(matched_donors).set{inp4}
                PREPROCESS_PROCESS(inp4,params.reduced_dims.vars_to_regress.value)

                if(params.seurat_integration.run_process){
                    DSB_INTEGRATE(
                        PREPROCESS_PROCESS.out.tmp_rsd.collect(),
                        params.reduced_dims.vars_to_regress.value,
                        params.reduced_dims.seurat_integration.k_anchor,
                        params.reduced_dims.seurat_integration.dims,
                        params.reduced_dims.seurat_integration.ndim_sct,
                        params.reduced_dims.seurat_integration.ndim_citeBgRemoved,
                        params.reduced_dims.seurat_integration.ndim_cite_integrated
                        )

                    MULTIMODAL_INTEGRATION(
                        DSB_INTEGRATE.out.tmp_rds_file,
                    )

                    VDJ_INTEGRATION(
                        chanel_cr_outs.collect(),
                        MULTIMODAL_INTEGRATION.out.wnn_integrated_file
                    )
                }
            }    


            andata = NORMALISE_AND_PCA.out.anndata
            outdir = NORMALISE_AND_PCA.out.outdir

        }else{
            andata = file__anndata_merged
            outdir = "${params.outdir}"
            LI4 = Channel.of([1, 'dummy_lisi'])
        }


        log.info """---Estimating PCA elbow.----"""
        PCA(andata,outdir,params.normalise.layer)
        ESTIMATE_PCA_ELBOW(
            PCA.out.outdir,
            PCA.out.anndata,
            params.reduced_dims.n_dims.add_n_to_estimate
        )

        if (params.reduced_dims.n_dims.auto_estimate) {
            log.info "n_pcs = automatically estimated."
            n_pcs = ESTIMATE_PCA_ELBOW.out.auto_elbow
        } else {
            log.info "n_pcs = Channel.from(params.reduced_dims.n_dims.value)"
            n_pcs = Channel.from(params.reduced_dims.n_dims.value)
        }

        SUBSET_PCS(
            PCA.out.outdir,
            PCA.out.anndata,
            PCA.out.metadata,
            PCA.out.pcs,
            PCA.out.param_details,
            n_pcs
        )
        
        PLOT_STATS(file__anndata_merged,
                    file__cells_filtered,
                    SUBSET_PCS.out.outdir,
                    SUBSET_PCS.out.anndata,
                    n_pcs)
                    
        file__anndata_merged = PCA.out.anndata
        
        LI4 = PLOT_STATS.out.LI

        if (params.cluster.known_markers.run_process) {
            channel__cluster__known_markers = Channel
                .fromList(params.cluster.known_markers.value)
                .map{row -> tuple(row.file_id, file(row.file))}
        } else {
            channel__cluster__known_markers = tuple('', '')
        }




        // "Correct" PCs using Harmony or BBKNN
        if (params.harmony.run_process) {
            HARMONY(
                PCA.out.outdir,
                PCA.out.anndata,
                PCA.out.metadata,
                PCA.out.pcs,
                PCA.out.param_details,
                n_pcs,
                Channel.fromList( params.harmony.variables_and_thetas.value)
            )

            UMAP_HARMONY(
                HARMONY.out.outdir,
                HARMONY.out.anndata,
                HARMONY.out.metadata,
                HARMONY.out.pcs,
                HARMONY.out.reduced_dims,
                "False",
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value,
                'harmony'
            )

            cluster_harmony__outdir = UMAP_HARMONY.out.outdir
            cluster_harmony__anndata = UMAP_HARMONY.out.anndata
            anndata =  UMAP_HARMONY.out.anndata
            outdir = UMAP_HARMONY.out.outdir
            cluster_harmony__metadata = UMAP_HARMONY.out.metadata
            cluster_harmony__pcs = UMAP_HARMONY.out.pcs
            cluster_harmony__reduced_dims = UMAP_HARMONY.out.reduced_dims
            

            // cluster_harmony__outdir.subscribe { println "cluster_harmony__outdir input: $it" }
            // cluster_harmony__anndata.subscribe { println "cluster_harmony__anndata input: $it" }
            // cluster_harmony__reduced_dims.subscribe { println "cluster_harmony__reduced_dims input: $it" }
            // cluster_harmony__metadata.subscribe { println "cluster_harmony__metadata input: $it" }
            // cluster_harmony__pcs.subscribe { println "cluster_harmony__pcs input: $it" }

            CLUSTERING_HARMONY(
                cluster_harmony__outdir,
                cluster_harmony__anndata,
                cluster_harmony__metadata,
                cluster_harmony__pcs,
                cluster_harmony__reduced_dims,
                "False",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy         
            )
            
            lisi_input2 = HARMONY.out.reduced_dims_params.collect()
            LI14 = CLUSTERING_HARMONY.out.dummy_output.collect()
            LI1 = LI14.mix(cluster_harmony__pcs)
                
        }else{
            lisi_input2 = Channel.of([1, 'dummy'])
            LI1 = Channel.of([1, 'dummy_harmony'])
        }

        if (params.bbknn.run_process) {

            if (params.bbknn.estimate_neighbors){
                // Here we use BBKNN dims for the UMAPS
                neighbors =  ["-1"]
                neighbors_umap = ["-1"]
                
            }else{
                neighbors = params.cluster.number_neighbors.value
                neighbors_umap = params.umap.n_neighbors.value
            }

            BBKNN(
                PCA.out.outdir,
                PCA.out.anndata,
                PCA.out.metadata,
                PCA.out.pcs,
                PCA.out.param_details,
                n_pcs,
                params.bbknn.batch_variable.value
            )

            UMAP_BBKNN(
                BBKNN.out.outdir,
                BBKNN.out.anndata,
                BBKNN.out.metadata,
                BBKNN.out.pcs,
                BBKNN.out.reduced_dims,
                "True",  // Don't look at the reduced_dims parameter
                neighbors,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value,
                'bbknn'
            )

            cluster_bbknn__outdir = UMAP_BBKNN.out.outdir
            cluster_bbknn__anndata = UMAP_BBKNN.out.anndata
            outdir = UMAP_BBKNN.out.outdir
            anndata = UMAP_BBKNN.out.anndata
            cluster_bbknn__metadata = UMAP_BBKNN.out.metadata
            cluster_bbknn__pcs = UMAP_BBKNN.out.pcs
            cluster_bbknn__reduced_dims = UMAP_BBKNN.out.reduced_dims

            CLUSTERING_BBKNN(
                cluster_bbknn__outdir,
                cluster_bbknn__anndata,
                cluster_bbknn__metadata,
                cluster_bbknn__pcs,
                cluster_bbknn__reduced_dims,
                "True",  // use_pcs_as_reduced_dims
                neighbors,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                neighbors_umap,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy
            )
            
            lisi_input3 = BBKNN.out.reduced_dims_params.collect()
            LI2 = CLUSTERING_BBKNN.out.dummy_output.collect()
                
        }else{
            lisi_input3 = Channel.of([1, 'dummy'])
            LI2 = Channel.of([1, 'dummy_bbknn'])
        }


        if (params.dont_integrate_just_cluster){
            DONT_INTEGRATE(
                PCA.out.outdir,
                PCA.out.anndata,
                PCA.out.metadata,
                PCA.out.pcs,
                PCA.out.param_details,
                n_pcs
            )

            UMAP(
                DONT_INTEGRATE.out.outdir,
                DONT_INTEGRATE.out.anndata,
                PCA.out.metadata,
                PCA.out.pcs,
                DONT_INTEGRATE.out.reduced_dims,
                "True",  // Don't look at the reduced_dims parameter
                params.cluster.number_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value,
                'bbknn'
            )
            cluster_outdir = UMAP.out.outdir
            cluster_anndata = UMAP.out.anndata
            cluster_metadata = UMAP.out.metadata
            cluster_pcs = UMAP.out.pcs
            cluster_reduced_dims = UMAP.out.reduced_dims
            
            CLUSTERING(
                cluster_outdir,
                cluster_anndata,
                cluster_metadata,
                cluster_pcs,
                cluster_reduced_dims,
                "True",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy
            )
        }

        if (params.lisi.run_process) {
            lisi_input = SUBSET_PCS.out.reduced_dims_params.collect()
            lisi_input_first = lisi_input.mix(lisi_input2)
            lisi_input_second = lisi_input_first.mix(lisi_input3)

            LISI(
                PCA.out.outdir,
                PCA.out.metadata,
                params.lisi.variables.value,
                lisi_input_second.collect()
            )
            
            LI3 = LISI.out.outdir.collect()
        }else{
            LI3 = Channel.of([1, 'dummy_lisi'])
        }
        LI5=LI1.combine(LI2)
        LI6=LI5.combine(LI3)
        LI=LI6.combine(LI4)
    emit:
        LI
        file__anndata_merged
        
}
