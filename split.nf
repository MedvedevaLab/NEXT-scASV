#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { SPLIT_FASTQ } from './modules/split_fastq'

// Validate required parameters
def checkRequiredParams() {
    def missing = []
    if (params.split_table == null) missing << "split_table"
    if (params.meta_json == null) missing << "meta_json" 
    if (params.input == null) missing << "input"
    if (params.outdir == null) missing << "outdir"
    
    if (missing.size() > 0) {
        error "Check the config file. Missing required parameter(s): ${missing.join(', ')}"
    }
}

workflow SPLIT_WORKFLOW {
    take:
        input_dir
        output_dir
    
    main:
        // Check required parameters
        checkRequiredParams()
        
        // Parse split_table to get sample_id/group combinations
        def sampleGroups = []
        def splitTableContent = file(params.split_table).text
        
        // Skip header and process each line
        splitTableContent.tokenize('\n')
                        .findAll { it.trim() && !it.trim().startsWith('barcode') } // Skip header and empty lines
                        .collect { line ->
                            def fields = line.trim().split('\t')
                            if (fields.size() >= 4) {
                                // Handle tables with a leading index column
                                [fields[2], fields[3]] // [sample_id, group]
                            } else if (fields.size() >= 3) {
                                [fields[1], fields[2]] // [sample_id, group]
                            }
                        }
                        .unique()
                        .each { sample_id, group ->
                            sampleGroups << [sample_id, group, file(params.split_table)]
                        }
        
        // Create channel from sample/group combinations
        sample_group_ch = Channel.fromList(sampleGroups)
        
        // Run the splitting process for each sample/group combination
        // SPLIT_FASTQ output: [sample_id, group, r1_fastq, r2_fastq]
        split_results = SPLIT_FASTQ(
            sample_group_ch,
            input_dir,
            output_dir
        )
        
        // Transform output to structured format for downstream workflows
        structured_output = split_results.map { sample_id, group, r1_fastq, r2_fastq ->
            // Avoid parsing sample_id/group from filenames; sample_id may itself contain '-'
            def sample_name = "${sample_id}-${group}"
            [sample_name, sample_id, group, r1_fastq, r2_fastq]
        }
    
    emit:
        structured_output
}