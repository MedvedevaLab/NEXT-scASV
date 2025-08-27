#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import workflows from separate files
include { SPLIT_WORKFLOW } from './split'
include { REMAP_WORKFLOW } from './remap'
include { CALL_WORKFLOW } from './call'
include { COUNT_READS_WORKFLOW } from './count_reads'
include { FIND_ASE_WORKFLOW } from './find_ase'

// Process to create temporary directory
process CREATE_TEMP_AND_OUT_DIR {
    output:
    path "temp_dir", emit: temp_dir
    
    script:
    """
    mkdir -p ${params.outdir}
    mkdir -p ${params.tmp_dir}
    ln -sf ${params.tmp_dir} temp_dir
    
    echo "Temporary directory created: ${params.tmp_dir}"
    echo "Output directory created: ${params.outdir}"
    """
}

workflow {
    // Define output directories for each step
    def split_output = "${params.outdir}/split"
    def remap_output = "${params.outdir}/align"
    def call_output = "${params.outdir}"
    def count_reads_output = "${params.outdir}/count_reads"
    def final_vcf = "${params.outdir}/vcf_by_sample"
    def find_ase_output = "${params.outdir}/find_ase"
    
    // Create temporary folder before all steps
    CREATE_TEMP_AND_OUT_DIR()
    temp_dir_channel = CREATE_TEMP_AND_OUT_DIR.out.temp_dir
    
    // Step 1: Split FASTQ files
    if (params.run_split) {
        log.info "Running SPLIT step"
        SPLIT_WORKFLOW(
            params.input,
            split_output
        )
        split_result = SPLIT_WORKFLOW.out
    } else {
        // Create empty channel for split result when not running split
        split_result = Channel.empty()
    }
    
    // Step 2: Remap reads
    if (params.run_remap) {
        log.info "Running REMAP step"
        
        // Use split result channel if split was run, otherwise use directory path
        remap_input = params.run_split ? split_result : split_output
        
        REMAP_WORKFLOW(
            remap_input,
            remap_output
        )
        remap_result = REMAP_WORKFLOW.out
    } else {
        // Create empty channel for remap result when not running remap
        remap_result = Channel.empty()
    }
    
    // Step 3: Call variants
    if (params.run_call) {
        log.info "Running CALL step"
        
        // Use remap result channel if remap was run, otherwise use directory path
        call_input = params.run_remap ? remap_result : remap_output
        
        CALL_WORKFLOW(
            call_input,
            call_output
        )
        call_result = CALL_WORKFLOW.out.filtered_results
        merged_variants = CALL_WORKFLOW.out.merged_variants
    } else {
        // Create empty channel for call result when not running call
        call_result = Channel.empty()
    }
    
    // Step 4: Count reads
    if (params.run_count_reads) {
        log.info "Running COUNT_READS step"
        
        // Use result channels if previous steps were run, otherwise use directory paths
        count_bams_input = params.run_remap ? remap_result : remap_output
        count_calls_input = params.run_call ? call_result : call_output
        count_merged_input = params.run_call ? merged_variants : call_output
        
        COUNT_READS_WORKFLOW(
            count_bams_input,
            count_calls_input,
            count_merged_input,
            count_reads_output
        )
        count_reads_result = COUNT_READS_WORKFLOW.out.sample_vcfs
    } else {
        // Create empty channel for count reads result when not running count reads
        count_reads_result = Channel.empty()
    }
    
    // Step 5: Find ASE (includes annotation)
    if (params.run_find_ase) {
        log.info "Running FIND_ASE step (includes annotation)"
        
        // Use count reads result channel if count reads was run, otherwise use final_vcf directory path
        ase_input = params.run_count_reads ? count_reads_result : final_vcf
        
        FIND_ASE_WORKFLOW(
            ase_input,
            find_ase_output
        )
        find_ase_result = FIND_ASE_WORKFLOW.out.fit_results
        annotated_tables = FIND_ASE_WORKFLOW.out.annotated_tables
        annotation_plots = FIND_ASE_WORKFLOW.out.plots
    } else {
        // Create empty channel for find ASE result when not running find ASE
        find_ase_result = Channel.empty()
    }
}