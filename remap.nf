#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { REMAP_FASTQ } from './modules/remap_fastq'
include { PROCESS_BAM } from './modules/process_bam'
include { BAM_TO_FASTQ } from './modules/bam_to_fastq'

// Validate required parameters
def checkRequiredParams() {
    def missing = []
    if (params.meta_json == null) missing << "meta_json"
    if (params.hisat2_index == null) missing << "hisat2_index"
    if (params.genome == null) missing << "genome"
    
    if (missing.size() > 0) {
        error "Check the config file. Missing required parameter(s): ${missing.join(', ')}"
    }
}

workflow REMAP_WORKFLOW {
    take:
        input_source // Can be input_dir (string) or channel from SPLIT_WORKFLOW
        output_dir

    main:
        log.info "Output directory: ${output_dir}"
        // Check required parameters
        checkRequiredParams()
        
        // Determine if input is a directory path or a channel from previous workflow
        if (input_source instanceof String || input_source instanceof groovy.lang.GString) {
            // Independent mode: parse metadata and find FASTQ files
            log.info "Running REMAP_WORKFLOW in independent mode with input directory: ${input_source}"
            
            // Create channel from metadata entries
            Channel
                .fromPath(params.meta_json)
                .map { meta_file ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_entries_list ${meta_file}"
                    def output = cmd.execute().text
                    return output
                }
                .splitCsv(sep: '\t', header: false)
                .map { sample_id, pat_id, group ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_sample_name '${sample_id}' '${group}'"
                    def sample_name = cmd.execute().text.trim().replaceAll("'", '')
                    def r1_file = file("${input_source}/${sample_name}.R1.fastq")
                    def r2_file = file("${input_source}/${sample_name}.R2.fastq")
                    def remap_out_path = "${output_dir}"
                    
                    // Check if files exist
                    if (!r1_file.exists() || !r2_file.exists()) {
                        log.warn "Missing FASTQ file(s) ${input_source}/${sample_name}.R1.fastq"
                        return null
                    }
                    
                    // Use pat_id downstream as the key (CALL/COUNT_READS expect pat_id)
                    [sample_name, pat_id, r1_file, r2_file, remap_out_path]
                }
                .filter { it != null }
                .set { remap_ch }
                
        } else {
            // Chained mode: use channel from SPLIT_WORKFLOW
            log.info "Running REMAP_WORKFLOW in chained mode with input channel"
            
            input_source
                .map { sample_name, sample_id, group, r1_fastq, r2_fastq ->
                    def remap_out_path = "${output_dir}"
                    // Chained mode must pass pat_id as sample_id (2nd element) to stay consistent
                    [sample_name, sample_id, r1_fastq, r2_fastq, remap_out_path]
                }
                .set { remap_ch }
        }

        // Remap each split fastq file
        // REMAP_FASTQ output: [sample_name, sample_id, bam, bai]
        remapped_ch = REMAP_FASTQ(remap_ch)
        
        // Create channel from BAM entries (for existing BAM files)
        Channel
            .fromPath(params.meta_json)
            .map { meta_file ->
                def cmd = "python3 ${projectDir}/bin/meta_parser.py create_bam_entries ${meta_file}"
                def output = cmd.execute().text
                return output
            }
            .splitCsv(sep: '\t', header: false)
            .map { sample_name, sample_id, bam_path, bai_path ->
                def bam_file = file(bam_path)
                def bai_file = file(bai_path)
                
                // Check if files exist
                if (!bam_file.exists() || !bai_file.exists()) {
                    log.warn "Missing BAM file(s) for ${sample_name}"
                    return null
                }
                
                [sample_name, sample_id, bam_file, bai_file]
            }
            .filter { it != null }
            .set { bam_ch }
        
        // Prepare channels for PROCESS_BAM
        // From remapped BAMs
        remapped_ch
            .map { sample_name, sample_id, bam, bai ->
                [sample_name, sample_id, bam]
            }
            .set { remapped_process_ch }
        
        // From existing BAMs
        bam_ch
            .map { sample_name, sample_id, bam, bai ->
                [sample_name, sample_id, bam]
            }
            .set { existing_bam_process_ch }
        
        // Combine both channels
        process_bam_ch = remapped_process_ch.mix(existing_bam_process_ch)
        
        // Process the BAMs
        processed_bams = PROCESS_BAM(process_bam_ch)
    
    emit:
        processed_bams
}
