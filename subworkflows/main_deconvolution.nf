nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP;capture_cellsnp_files;DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION; ASSESS_CALL_RATE } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE } from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include { VIREO;POSTPROCESS_SUMMARY;VIREO as VIREO_SINGLE_DONOR; REMOVE_DUPLICATED_DONORS_FROM_GT;VIREO_SUBSAMPLING;VIREO_SUBSAMPLING_PROCESSING; GENOTYPE_MATCHER; CAPTURE_VIREO } from "$projectDir/modules/nf-core/modules/vireo/main"
include { SUBSET_BAM_PER_BARCODES_AND_VARIANTS } from "$projectDir/modules/nf-core/modules/subset_bam_per_barcodes_and_variants/main"
include { FREEBAYES } from "$projectDir/modules/nf-core/modules/freebayes/main"
include { GUZIP_VCF } from "$projectDir/modules/nf-core/modules/guzip_vcf/main"
include { PREPROCESS_GENOTYPES } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { SOUPORCELL } from "$projectDir/modules/nf-core/modules/souporcell/main"
include { SPLIT_DONOR_H5AD; PREP_ASSIGNMENTS_FILE } from "$projectDir/modules/nf-core/modules/split_donor_h5ad/main"
include { PLOT_DONOR_CELLS } from "$projectDir/modules/nf-core/modules/plot_donor_cells/main"
include {SOUPORCELL_VS_VIREO} from "$projectDir/modules/nf-core/modules/plot_souporcell_vs_vireo/main"
include { REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE; ENHANCE_VIREO_METADATA_WITH_DONOR } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { match_genotypes } from './match_genotypes'
include {ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/nf-core/modules/genotypes/main"
include {SUBSET_WORKF} from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include {REPLACE_GT_DONOR_ID2; REPLACE_GT_DONOR_ID2 as REPLACE_GT_DONOR_ID_SUBS } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { STAGE_FILE } from "$projectDir/subworkflows/local/retrieve_recourses"
include {GT_MATCH_POOL_IBD } from "$projectDir/modules/nf-core/modules/genotypes/main"

include {VIREO_GT_FIX_HEADER; VIREO_GT_FIX_HEADER as VIREO_GT_FIX_HEADER_SUBS; MERGE_GENOTYPES_IN_ONE_VCF; VIREO_ADD_SAMPLE_PREFIX; MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN; MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_SUBSET} from "$projectDir/modules/nf-core/modules/genotypes/main"
include {collect_file as collect_file1;
        collect_file as collect_file2;
        collect_file as collect_file3;
        collect_file as collect_file4;
        collect_file as collect_file5;
        collect_file as collect_file6;
        collect_file as collect_file7;
        collect_file as collect_file8;
        collect_file as collect_file9;
        collect_file as collect_file10} from "$projectDir/modules/nf-core/modules/collect_file/main"

workflow  main_deconvolution {

    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        scrublet_paths
        ch_ref_vcf2
        genome

    main:
		log.info "#### running DECONVOLUTION workflow #####"


        vcf_candidate_snps = STAGE_FILE(params.cellsnp.vcf_candidate_snps)
        // genome.subscribe { println "genome: $it" }

        if (params.genotype_input.run_with_genotype_input) {
            // We have to produce a single vcf file for each individual pool.
            // Therefore we create 2 channels:
            // 1) All the expected vcf ids listed in the donor table
            log.info "# selecting SNPs and subsetting GT #"
            Channel.fromPath(params.input_data_table,      
                followLinks: true,
                checkIfExists: true
                ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
                .set { donors_in_pools }
            
            // This will subsequently result in a joint vcf file per cohort per donors listed for each of the pools that can be used in VIREO and/or GT matching algorythm.
            ch_ref_vcf = PREPROCESS_GENOTYPES(ch_ref_vcf2)
            SUBSET_WORKF(ch_ref_vcf,donors_in_pools,'AllExpectedGT',genome)
            merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
            merged_expected_genotypes2 = merged_expected_genotypes.combine(vcf_candidate_snps)
            GT_MATCH_POOL_IBD(SUBSET_WORKF.out.samplename_subsetvcf_ibd,'Withing_expected','Expected')

            DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION(params.add_snps_to_pile_up_based_on_genotypes_provided,merged_expected_genotypes2)
            cellsnp_panels = DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.cellsnp_pool_panel
            
            informative_uninformative_sites = DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.informative_uninformative_sites

        }else{
            ch_ref_vcf = Channel.of()
        }

        ch_poolid_donor_assignment = Channel.empty()
        ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorstbi,donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}
        ch_experiment_bam_bai_barcodes.combine(ch_experiment_npooled, by: 0).set{cellsnp_with_npooled}

        if (params.genotype_input.run_with_genotype_input) {
            // Here we provide sites to pile up the snps within the pool. 
            // As a starting point instead of the default cellsnp panel users are encouraged to:
            if (params.provide_within_pool_donor_specific_sites_for_pilup){
                cellsnp_with_npooled_pre = cellsnp_with_npooled.join(cellsnp_panels, remainder: true)
            }else{
                cellsnp_with_npooled_pre = cellsnp_with_npooled.combine(vcf_candidate_snps)
            }
            cellsnp_with_npooled_pre2 = cellsnp_with_npooled_pre.combine(vcf_candidate_snps)
            cellsnp_with_npooled2 = cellsnp_with_npooled_pre2.map { experiment, bam, bai, barcodes, nrPool, panel, default_cellsnp ->
                                                                    [
                                                                        experiment,
                                                                        bam,
                                                                        bai,
                                                                        barcodes,
                                                                        nrPool,
                                                                        panel == null ? default_cellsnp : panel
                                                                    ]
                                                                }
        }else{
            cellsnp_with_npooled2 = cellsnp_with_npooled.combine(vcf_candidate_snps)
        }

        log.info('Capturing some of the existing CELLSNP files')
        capture_cellsnp_files(params.existing_cellsnp)
        capture_cellsnp_files.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], "${row[1]}")}
            .set{cellsnp_output_dir1}
        cellsnp_output_dir1.join(cellsnp_with_npooled2, remainder: true).set{filter_channel}
        // filter_channel.subscribe { println "cellsnp_output_dir1 filter_channel: $it" }
        filter_channel.filter{ it[1] == null }.map{row -> tuple(row[0], row[2],row[3],row[4],row[5],row[6])}.set{cellsnp_with_npooled}
        
        capture_cellsnp_files.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], file("${row[1]}/cellSNP.cells.vcf.gz"))}
            .set{cellsnp_cell_vcfs}           
        capture_cellsnp_files.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0])}
            .set{cellsnp_cap_ids} 

        CELLSNP(cellsnp_with_npooled)
        if (params.genotype_input.run_with_genotype_input) {
            // Here we assess how many informative sites has been called on. 
            CELLSNP.out.cell_vcfs.combine(DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.informative_uninformative_sites, by: 0).set{assess_call_rate_input}
            ASSESS_CALL_RATE(assess_call_rate_input)
            collect_file9(ASSESS_CALL_RATE.out.variants_description.collect(),"all_variants_description.tsv",params.outdir+'/deconvolution/concordances',1,'')
        }
        for_bam_pileups = CELLSNP.out.for_bam_pileups
        
        cellsnp_cap_ids.combine(ch_experiment_bam_bai_barcodes, by: 0).combine(vcf_candidate_snps).set{cellsnp_cap_ids2}
        cellsnp_cap_ids2.map{row -> tuple(row[0], row[4],row[1])}.set{cellsnp_cap_ids2}
        for_bam_pileups2 = cellsnp_cap_ids2.mix(for_bam_pileups)

        cellsnp_output_dir2 = CELLSNP.out.cellsnp_output_dir
        cellsnp_output_dir=cellsnp_output_dir1.concat(cellsnp_output_dir2)
        cellsnp_cell_vcfs2=cellsnp_cell_vcfs.concat(CELLSNP.out.cell_vcfs)

        cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf2}
        full_vcf2.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf2}
         // Here we run Vireo software to perform the donor deconvolution. Note that we have coded the pipeline to be capable in using
         // the full genotypes as an input and also subset to the individuals provided as an input in the donor_vcf_ids column.
        if (params.genotype_input.vireo_with_gt && params.genotype_input.run_with_genotype_input) {
            log.info "---running Vireo with genotype input----"
            // for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
            // Here we subset the genotypes. This happens if the input.nf contains subset_genotypes = true
            // log.info "---We are using subset genotypes running Vireo----"
            // We need to make sure that the expected genotypes dont contain repeated genotypes - donors sequenced twice.

            if (params.genotype_phenotype_mapping_file!=''){
                REMOVE_DUPLICATED_DONORS_FROM_GT(merged_expected_genotypes,params.genotype_phenotype_mapping_file,params.input_data_table)
                merged_expected_genotypes2= REMOVE_DUPLICATED_DONORS_FROM_GT.out.merged_expected_genotypes
            }else{
                merged_expected_genotypes2 = merged_expected_genotypes
            }

            cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                .combine(merged_expected_genotypes2, by: 0).set{full_vcf_pre}

            full_vcf_pre.join(full_vcf2, remainder: true).set{filter_channel_vcf}
            // just take the null chanel, and add that to full vcf 
            filter_channel_vcf.filter{ it[1] == null }.map{row -> tuple(row[0], row[2],row[3],row[4],row[5])}.set{full_vcf4}
            full_vcf_pre.mix(full_vcf4).set{full_vcf}
        } else {
            // Vireo can also be run without the genotypes, and the performance is equally good then running with, however have to be careful in what sites file is provided.
            // Here we run it without the genotypes and the donors are labeled as donor0, donor1 etc, dependant on the number of donors set in the input file n_pooled column
            log.info "-----running Vireo without genotype input----"
            full_vcf = full_vcf2
        }


        CAPTURE_VIREO(params.existing_vireo)
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/donor_ids.tsv")}
            .set{vireo_out_sample_donor_ids_cap}  
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}")}
            .set{vireo_paths_0} 
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/GT_donors.vireo.vcf.gz")}
            .set{vireo_out_sample_donor_vcf_cap}  
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/summary.tsv")}
            .set{vireo_out_sample_summary_tsv_cap} 
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/GT_donors.vireo.vcf.gz", "${row}/donor_ids.tsv")}
            .set{all_required_data_cap_pre} 


        // When all the channels are prpeared accordingly we exacute the vireo with the prpeared channel.
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
        full_vcf2.join(vireo_out_sample_donor_ids_cap, remainder: true).set{filter_channel_vireo}

        all_required_data_cap_pre.join(full_vcf2, remainder: false).set{all_required_data_cap_pre}
        all_required_data_cap_pre.map{row -> tuple(row[0], row[1],row[2],row[5],row[6])}.set{all_required_data_cap}

        filter_channel_vireo.filter{ it[5] == null }.map{row -> tuple(row[0], row[1],row[2],row[3],row[4])}.set{full_vcf3}
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}
        
        not_deconvoluted.join(for_bam_pileups2, remainder: false).set{for_bam_pileups_pre}
        for_bam_pileups_pre.map{row -> tuple(row[0], row[5],row[6],row[7])}.set{for_bam_pileups_3}

        Channel.of(1..params.vireo.subsample_times).set{itterations}
        full_vcf2.combine(itterations).set{vireo_extra_repeats}

        VIREO(full_vcf3)
        VIREO.out.summary.mix(vireo_out_sample_summary_tsv_cap).set{summary_files}
        POSTPROCESS_SUMMARY(summary_files)
        


        //But to make it consistent we still randomly assign donor IDs per pool for each of the names.
        vireo_with_gt = Channel.of(params.genotype_input.vireo_with_gt)
        VIREO.out.all_required_data.set{replacement_input3}
        replacement_input3.mix(all_required_data_cap).set{replacement_input2}

        
        replacement_input2.combine(POSTPROCESS_SUMMARY.out.summary_tsvs, by: 0).set{replacement_input}

        replacement_input.combine(vireo_with_gt).set{vir_repl_input}

        REPLACE_GT_DONOR_ID2(vir_repl_input)
        // combine the outputs of the capture_vireo and Vireo
        vireo_paths_0.concat(REPLACE_GT_DONOR_ID2.out.output_dir).set{vireo_paths}

        VIREO_GT_FIX_HEADER(REPLACE_GT_DONOR_ID2.out.infered_vcf,genome)
        VIREO_ADD_SAMPLE_PREFIX(VIREO_GT_FIX_HEADER.out.infered_vcf)
        

        vireo_out_sample_donor_vcf = VIREO_GT_FIX_HEADER.out.infered_vcf
        vireo_out_sample_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample_summary_tsv
        vireo_out_sample__exp_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample__exp_summary_tsv
        vireo_out_sample_donor_ids = REPLACE_GT_DONOR_ID2.out.sample_donor_ids

        

        // if (params.souporcell.run){
        //     // YASCP pipeline is also capable in running SOUPORCELL instead of VIREO. If activated SOUPORCELL will be used.
        //     // yascp currently doesnt have an option to take souporcell assignments as downstream instead of vireo but this will be added shortly.
        //     // This runs the Souporcell Preprocessing

        //     if (params.souporcell.use_raw_barcodes) {
        //         // read raw cellranger barcodes per pool for souporcell
        //         log.info """Here use the raw barcodes"""
        //     } else if (! params.souporcell.use_raw_barcodes) {
        //         ch_experiment_bam_bai_barcodes
        //             .combine(ch_experiment_npooled, by: 0)
        //             .set{ch_experiment_bam_bai_barcodes_npooled}
        //     }

        //     // This runs the Souporcell
        //     // Similarly to the VIREO Soupocell can be run with and without genotypes and folowing prpeares the inputs accordingly to each option.
        //     // Soupocell cant digest a gz file gence we extract the data.
        //     if (params.genotype_input.run_with_genotype_input) {

        //         if (params.genotype_input.subset_genotypes){
        //             // this will run the soupocell with the subset genotypes. This happens if the input.nf contains subset_genotypes = true
        //             log.info "---We are using subset genotypes running Suporcell----"
        //             GUZIP_VCF(SUBSET_GENOTYPE.out.samplename_subsetvcf)
        //             // ch_experiment_bam_bai_barcodes_npooled.combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
        //             ch_experiment_bam_bai_barcodes_npooled.combine(GUZIP_VCF.out.souporcell_vcf, by: 0).set{full_vcf}

        //         }else{
        //             log.info "---We are using a full genotype input for Suporcell----"
        //             // this however currently doesnt work and the individuals have to be provided.
        //             // here just add the full vcf path to each of the ch_experiment_bam_bai_barcodes_npooled
        //             GUZIP_VCF(tuple('full_vcf', file(params.genotype_input.full_vcf_file)))
        //             GUZIP_VCF.out.souporcell_vcf.map { sample, vcf -> vcf }.set{vcf_file}
        //             ch_experiment_bam_bai_barcodes_npooled.map {
        //                 samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
        //                 tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters)
        //                 }.set{full_vcf}
        //             full_vcf.combine(vcf_file).set{full_vcf}

        //         }
        //         full_vcf.combine(donors_in_lane, by: 0).set{full_vcf}
        //     }
        //     else{
        //         // Soupocell can also run without the genotypes. This will prpeare channels to run it withoutt.
        //         log.info "-----running Suporcell without genotype input----"
        //         // here make add an empty entry [] to the ch_experiment_bam_bai_barcodes_npooled
        //         ch_experiment_bam_bai_barcodes_npooled.map {
        //             samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
        //             tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters,[],[])
        //             }.set{full_vcf}
        //     }

        //     // When all the channels are prepeared then we can run the soupocell accordingly.
        //     SOUPORCELL(full_vcf,
        //         Channel.fromPath(params.souporcell.reference_fasta).collect())
        //     // Regardless if the Soupocell is run with or without genotypes we still need to match the donor ids with the cluster ids since this does not happen automatically is Soupocell.
        // }

        // if (params.souporcell.run && params.vireo.run) {
        //     // This is an orringinal comparison of Vireo and soupocell, this needs to be cleaned up.
        //     SOUPORCELL_VS_VIREO(
        //         vireo_out_sample_donor_ids
        //         .combine(SOUPORCELL.out.souporcell_output_files.map {a,b,c,d -> tuple(a,b)},
        //             by: 0))
        // }


        ch_experiment_bam_bai_barcodes
            .map { samplename, bam_file, bai_file, barcodes_tsv_gz -> tuple(samplename, file(bam_file)) }
            .combine(vireo_out_sample_donor_ids, by: 0 )
            .set { ch_experiment_bam_vireo_donor_ids }

        // If sample is not deconvoluted we will use scrublet to detect the doublets and remove them.
        not_deconvoluted.map{ experiment, donorsvcf, npooled,t,t2 -> tuple(experiment, 'None')}.set{not_deconvoluted2}
        split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
        split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)
        // combining these 2 channels in one
        split_channel3 = split_channel.mix(split_channel2)
        // adding the scrublet paths to the channel.
        split_channel4 = split_channel3.combine(scrublet_paths, by: 0)
        split_channel5 = split_channel4.map{
            val_sample, val_donor_ids_tsv, val_filtered_matrix_h5, path_scrublet ->
            [  val_sample,file(val_donor_ids_tsv),file(val_filtered_matrix_h5),path_scrublet,params.outdir]
        }

        SPLIT_DONOR_H5AD(split_channel5)
        PREP_ASSIGNMENTS_FILE(split_channel5)
        cell_assignments = PREP_ASSIGNMENTS_FILE.out.cell_assignments
        if (params.bam_pileup_per_donor){
            bam_subset_chanel = SPLIT_DONOR_H5AD.out.sample_donor_level_barcodes.transpose().combine(for_bam_pileups_3, by: 0)
            SUBSET_BAM_PER_BARCODES_AND_VARIANTS(bam_subset_chanel)  // This process subsets the ba

            chromosomes =  Channel.of(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)
            freebayes_pre = SUBSET_BAM_PER_BARCODES_AND_VARIANTS.out.freebayes_input

            freebayes_in = freebayes_pre.combine(chromosomes)

            FREEBAYES(freebayes_in,genome)
            vireo_out_sample_donor_vcf = FREEBAYES.out.freebayes_vcf

            FREEBAYES.out.gt_pool.groupTuple(by:0).set{fbb1}
            fbb1.map{row->tuple(row[0], row[1], row[2], row[3][0])}.set{fbb1}
            MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN(fbb1,'freebayes')
            MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN.out.gt_pool.groupTuple(by:0).set{fbb2}
            MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES(fbb2,'freebayes')
            gt_math_pool_against_panel_input2 = MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES.out.gt_pool.combine(ch_ref_vcf)
            vir_inp =  MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES.out.vir_input

        }else{
            gt_math_pool_against_panel_input2 = Channel.of()
            vir_inp =  Channel.of()
        }
        // vir_inp.subscribe { println "vir_inp: $it" }
        vir_inp2 = vir_inp.collect()
        vir_inp3 = vireo_paths.collect()
        vireo_paths2 = gt_matcher_inp = vir_inp3.mix(vir_inp2)
              
        GENOTYPE_MATCHER(gt_matcher_inp)

        matched_donors = GENOTYPE_MATCHER.out.matched_donors

        if (params.genotype_input.run_with_genotype_input) {
            if (params.do_vireo_subsampling){
                VIREO_SUBSAMPLING(vireo_extra_repeats)
                VIREO_SUBSAMPLING.out.output_dir.concat(VIREO.out.output_dir_subsampling).set{tuple_1}
                tuple_1.groupTuple(by:0).set{vspp0}
                VIREO_SUBSAMPLING_PROCESSING(vspp0)
                subsampling_donor_swap = VIREO_SUBSAMPLING_PROCESSING.out.subsampling_donor_swap
            }else{
                subsampling_donor_swap = Channel.from("$projectDir/assets/fake_file.fq")
            }

            

            VIREO_GT_FIX_HEADER.out.gt_pool
                .combine(ch_ref_vcf)
                .set { gt_math_pool_against_panel_input3 }
            gt_math_pool_against_panel_input = gt_math_pool_against_panel_input3.mix(gt_math_pool_against_panel_input2)

            match_genotypes(vireo_out_sample_donor_vcf,merged_expected_genotypes,VIREO_GT_FIX_HEADER.out.gt_pool,gt_math_pool_against_panel_input,genome,ch_ref_vcf,cellsnp_cell_vcfs2,cell_assignments,subsampling_donor_swap,informative_uninformative_sites)
            gt_matches = match_genotypes.out.donor_match_table.collect()

            ENHANCE_STATS_GT_MATCH(match_genotypes.out.donor_match_table_enhanced,params.input_data_table)
            collect_file1(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/deconvolution/vireo_processed',1,'')
            collect_file10(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/deconvolution/gtmatch',1,'')
            assignments_all_pools = collect_file1.out.output_collection
            gt_matches = Channel.from("$projectDir/assets/fake_file.fq")
            assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
        }else{
            gt_matches = Channel.from("$projectDir/assets/fake_file.fq")
            assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
        }

        // collect file paths to h5ad files in tsv tables:
        // header_seed="experiment_id\tdonor\th5ad_filepath"
        // collect_file2(SPLIT_DONOR_H5AD.out.donors_h5ad_tsv.collect(),"donors_h5ad.tsv",0,0,header_seed)

        // header_seed="experiment_id\th5ad_filepath"
        // collect_file3(SPLIT_DONOR_H5AD.out.exp__donors_h5ad_tsv.collect(),"exp__donors_h5ad.tsv",0,0,header_seed)

        // header_seed="experiment_id\tdonor\th5ad_filepath"
        // collect_file4(SPLIT_DONOR_H5AD.out.donors_h5ad_assigned_tsv.collect(),"donors_h5ad_assigned.tsv",0,0,header_seed)

        header_seed="experiment_id\th5ad_filepath"
        out_h5ad = collect_file5(SPLIT_DONOR_H5AD.out.exp__donors_h5ad_assigned_tsv.collect(),"exp__donors_h5ad_assigned.tsv",0,0,header_seed)

        // header_seed="experiment_id\th5ad_filepath"
        // collect_file6(SPLIT_DONOR_H5AD.out.h5ad_tsv.collect(),"cellranger_as_h5ad.tsv",0,0,header_seed)

        header_seed="experiment_id\tdonor\tn_cells"
        ch_vireo_donor_n_cells_tsv = collect_file7(REPLACE_GT_DONOR_ID2.out.sample_summary_tsv.collect(),"vireo_donor_n_cells.tsv",0,0,header_seed)

        // paste experiment_id and donor ID columns with __ separator
        header_seed="experiment_id\tn_cells"
        vireo_out_sample__exp_summary_tsv = collect_file8(SPLIT_DONOR_H5AD.out.donor_n_cells.collect(),"vireo_exp__donor_n_cells.tsv",0,0,header_seed)

        if (params.genotype_input.run_with_genotype_input & params.genotype_input.posterior_assignment) {
            if (params.extra_sample_metadata!='' & params.add_donor_metadata){
                // Here we have sample level metadata but we have chosen to keep the donor ids.
                // in this scenario we enhance the donor level vireo metadata file to add the donor metadata to the h5ads and eventually to the donor and teanche report
                ENHANCE_VIREO_METADATA_WITH_DONOR(params.extra_sample_metadata,vireo_out_sample__exp_summary_tsv,REPLACE_GT_DONOR_ID.out.assignments.collect())
                vireo_out_sample__exp_summary_tsv = ENHANCE_VIREO_METADATA_WITH_DONOR.out.replaced_vireo_exp__donor_n_cells_out
            }
        }
        PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)

    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf
        ch_poolid_csv_donor_assignments = ch_poolid_donor_assignment
        sample_possorted_bam_vireo_donor_ids = ch_experiment_bam_vireo_donor_ids
        assignments_all_pools
        vireo_paths2
        matched_donors
}
