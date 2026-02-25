#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { CALL_VARIANTS } from './modules/call_variants'
include { CALL_VARIANTS_CELLSNP } from './modules/call_variants_cellsnp'
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
                .map { sample_id, pat_id, group ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_sample_name '${sample_id}' '${group}'"
                    def sample_name = cmd.execute().text.trim().replaceAll("'", '')
                    def bam_file = file("${input_source}/${sample_name}_marked_filtered.bam")
                    def bai_file = file("${input_source}/${sample_name}_marked_filtered.bam.bai")
                    
                    // Check if files exist
                    if (!bam_file.exists() || !bai_file.exists()) {
                        log.warn "Missing BAM file(s) for ${sample_name} (${bam_file})"
                        return null
                    }
                    log.info "✓ File exists: ${sample_name}"
                    // Use pat_id downstream
                    tuple(sample_name, pat_id.toString(), bam_file, bai_file)
                }
                .filter { it != null }
                .set { filtered_bams }
                
        } else {
            // Chained mode: use channel from REMAP_WORKFLOW
            log.info "Running CALL_WORKFLOW in chained mode with input channel"
            
            input_source
                .map { sample_name, sample_id, bam_file, bai_file ->
                    // Chained mode must provide pat_id as the 2nd element
                    tuple(sample_name, sample_id.toString(), bam_file, bai_file)
                }
                .set { filtered_bams }
        }

        // Get chromosomes for variant calling
        //chrom_ch = Channel.fromPath(params.autosomes)
        //    .splitCsv(header:true, sep:'\t')
        //    .map { row -> row.chroms }
	chrom_ch = Channel.fromPath(params.autosomes)
	    .splitText()
	    .map { line -> line.trim() }
            .filter { it && !it.startsWith('#') }
            .view { chrom -> "INPUT_CHROM: ${chrom}" }

        // Collect BAM paths (bcftools caller only needs BAMs)
        all_bams = filtered_bams
            .map { sample_name, pat_id, bam, bai -> bam }
            .collect()

        // Early sanity check: make failures actionable (e.g. if something like 'R' sneaks in)
        all_bams_checked = all_bams.map { bams ->
            if( !(bams instanceof List) )
                error "CALL_WORKFLOW expected a List of BAM paths, but got ${bams?.getClass()?.name}: ${bams}"
            def bad = bams.findAll { x ->
                !(x instanceof java.nio.file.Path) && !(x instanceof File)
            }
            if( bad && bad.size() > 0 )
                error "CALL_WORKFLOW got non-path values in BAM list (showing up to 10): ${bad.take(10)}"
            return bams
        }
        
        // Call variants for each chromosome
        if (params.variant_caller == 'cellsnp') {
            log.info "Calling variants with cellsnp-lite (Mode 2b) per chromosome"

            // For cellsnp-lite (-S/-i), we must keep BAMs and sample IDs in the exact same order.
            bam_id_pairs = filtered_bams
                .map { sample_name, pat_id, bam, bai -> tuple(pat_id.toString(), bam) }
                .collect()

            sorted_bam_id_pairs = bam_id_pairs
                .map { pairs ->
                    return pairs.toList().sort { a, b -> a[0].toString() <=> b[0].toString() }
                }

            all_bams_cellsnp = sorted_bam_id_pairs
                .map { pairs -> pairs.collect { it[1] } }

            all_sample_ids = sorted_bam_id_pairs
                .map { pairs -> pairs.collect { it[0].toString() } }

            // Also validate the derived list
            all_bams_cellsnp_checked = all_bams_cellsnp.map { bams ->
                if( !(bams instanceof List) )
                    error "CALL_WORKFLOW expected a List of BAM paths for cellsnp-lite, but got ${bams?.getClass()?.name}: ${bams}"
                def bad = bams.findAll { x ->
                    !(x instanceof java.nio.file.Path) && !(x instanceof File)
                }
                if( bad && bad.size() > 0 )
                    error "CALL_WORKFLOW got non-path values in cellsnp BAM list (showing up to 10): ${bad.take(10)}"
                return bams
            }

            variants = CALL_VARIANTS_CELLSNP(chrom_ch, all_bams_cellsnp_checked, all_sample_ids)
        } else {
            log.info "Calling variants with bcftools mpileup/call per chromosome"
            variants = CALL_VARIANTS(chrom_ch, all_bams_checked)
        }
        
        // Collect variants for merging
        variants_to_merge = variants.collect()
        
        // Merge all variant calls
        merged_variants = MERGE_VARIANTS(variants_to_merge)

        // Prepare input for FILTER_VARIANTS
        sample_ids = filtered_bams
            .map { sample_name, pat_id, bam, bai -> pat_id.toString() }
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
