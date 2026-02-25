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

        // Build sample -> sample_id mapping from meta.json
        def metaObj = new groovy.json.JsonSlurper().parseText(file(params.meta_json).text)
        def sampleToSampleId = [:]
        metaObj.each { k, v ->
            def sample = k?.toString()
            def sample_id = v?.pat_id?.toString()
            if (sample && sample_id) sampleToSampleId[sample] = sample_id
        }
        
        // Parse split_table to get sample/group combinations
        def sampleGroups = []
        def splitTableContent = file(params.split_table).text
        
        // Skip header and process each line
        splitTableContent.tokenize('\n')
                        .findAll { it.trim() && !it.trim().startsWith('barcode') } // Skip header and empty lines
                        .collect { line ->
                            def fields = line.split('\t')
                            if (fields.size() >= 3) {
                                [fields[1], fields[2]] // [sample, group]
                            }
                        }
                        .unique()
                        .each { sample, group ->
                            sampleGroups << [sample, group, file(params.split_table)]
                        }
        
        // Create channel from sample/group combinations
        sample_group_ch = Channel.fromList(sampleGroups)
        
        // Run the splitting process for each sample/group combination
        // SPLIT_FASTQ output: [sample, group, r1_fastq, r2_fastq]
        split_results = SPLIT_FASTQ(
            sample_group_ch,
            input_dir,
            output_dir
        )
        
        // Transform output to structured format for downstream workflows
        structured_output = split_results.map { sample, group, r1_fastq, r2_fastq ->
            // Avoid parsing sample/group from filenames; sample may itself contain '-'
            def sample_name = "${sample}-${group}"
            def sample_id = sampleToSampleId[sample?.toString()]
            if (!sample_id) {
                error "Missing sample_id for sample='${sample}' in meta_json='${params.meta_json}'"
            }
            // Emit: [sample_name, sample, sample_id, group, r1_fastq, r2_fastq]
            [sample_name, sample, sample_id, group, r1_fastq, r2_fastq]
        }
    
    emit:
        structured_output
}