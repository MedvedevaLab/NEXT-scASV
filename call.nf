#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { CALL_VARIANTS } from './modules/call_variants'
include { MERGE_VARIANTS } from './modules/merge_variants'
include { FILTER_VARIANTS } from './modules/filter_variants'

// Validate required parameters
def checkRequiredParams() {
    def missing = []
    if (params.meta_json == null) missing << "meta_json"
    if (params.genome == null) missing << "genome"
    if (params.ref_genome == null) missing << "ref_genome"
    if (params.dbsnp == null) missing << "dbsnp"
    if (params.autosomes == null) missing << "autosomes"
    if (params.min_GQ == null) missing << "min_GQ"
    if (params.min_DP == null) missing << "min_DP"
    if (params.min_AD == null) missing << "min_AD"

    if (missing.size() > 0) {
        error "Check the config file. Missing required parameter(s): ${missing.join(', ')}"
    }
}


workflow CALL_WORKFLOW {
    take:
        input_source // Can be input_dir (string) or channel from REMAP_WORKFLOW
        output_dir
    
    main:
        // Check required parameters
        checkRequiredParams()

        // Parse metadata using Python script
        def meta = Channel.fromPath(params.meta_json)
            .map { meta_file ->
                def cmd = "python3 ${projectDir}/bin/meta_parser.py parse_meta_json ${meta_file}"
                def output = cmd.execute().text
                def parsed = new groovy.json.JsonSlurper().parseText(output)
                println "Parsed metadata: ${parsed}"
                return parsed
            }
            .first()
        
        // Determine if input is a directory path or a channel from previous workflow
        if (input_source instanceof String || input_source instanceof groovy.lang.GString) {
            // Independent mode: parse metadata and find BAM files
            log.info "Running CALL_WORKFLOW in independent mode with input directory: ${input_source}"
            
            // Create channel from metadata entries
            Channel
                .fromPath(params.meta_json)
                .map { meta_file ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_entries_list ${meta_file}"
                    def output = cmd.execute().text
                    return output
                }
                .splitCsv(sep: '\t', header: false)
                .map { sample_id, group ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_sample_name '${sample_id}' '${group}'"
                    def sample_name = cmd.execute().text.trim().replaceAll("'", '')
                    def bam_file = file("${input_source}/${sample_name}_marked_filtered.bam")
                    def bai_file = file("${input_source}/${sample_name}_marked_filtered.bam.bai")
                    
                    // Check if files exist
                    if (!bam_file.exists() || !bai_file.exists()) {
                        log.warn "Missing BAM file(s) for ${sample_name} (${bam_file})"
                        return null
                    }
                    
                    [sample_name, sample_id, bam_file, bai_file]
                }
                .filter { it != null }
                .set { filtered_bams }
                
        } else {
            // Chained mode: use channel from REMAP_WORKFLOW
            log.info "Running CALL_WORKFLOW in chained mode with input channel"
            
            input_source
                .map { sample_name, bam_file, bai_file ->
                    // Extract sample_id from sample_name (assuming format like "sample_id+group")
                    def sample_id = sample_name.split('\\-')[0]
                    [sample_name, sample_id, bam_file, bai_file]
                }
                .set { filtered_bams }
        }

        // Get chromosomes for variant calling
        chrom_ch = Channel.fromPath(params.autosomes)
            .splitCsv(header:true, sep:'\t')
            .map { row -> row.chroms }

        // Collect all BAM files for variant calling
        all_bams = filtered_bams
            .map { sample_name, sample_id, bam, bai -> bam }
            .collect()
        
        // Call variants for each chromosome
        variants = CALL_VARIANTS(chrom_ch, all_bams)
        
        // Collect variants for merging
        variants_to_merge = variants.collect()
        
        // Merge all variant calls
        merged_variants = MERGE_VARIANTS(variants_to_merge)

        // Prepare input for FILTER_VARIANTS
        sample_ids = filtered_bams
            .map { sample_name, sample_id, bam, bai -> sample_id }
            .unique()

        filter_input = sample_ids
            .combine(merged_variants)
            .map { sample_id, snps_vcf, snps_csi -> 
                tuple(sample_id, snps_vcf, snps_csi)
            }
        
        // Filter variants for each sample
        filtered_results = FILTER_VARIANTS(filter_input)
    
    emit:
        filtered_results
        merged_variants
}
