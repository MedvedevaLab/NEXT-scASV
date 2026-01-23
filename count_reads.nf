#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import meta parser
include { GENERATE_H5_TABLES } from './modules/mapp_bias_and_allele_count/generate_h5_tables'
include { REMAP_BAMFILES } from './modules/mapp_bias_and_allele_count/remap_bam'
include { COUNT_READS } from './modules/mapp_bias_and_allele_count/count_reads'
include { RECODE_VCF } from './modules/mapp_bias_and_allele_count/recode_vcf'
include { CONCAT_AND_SPLIT_VCF } from './modules/mapp_bias_and_allele_count/concat_and_split_vcf'

// Validate required parameters
def checkRequiredParams() {
    def missing = []
    if (params.meta_json == null) missing << "meta_json"
    if (params.workdir == null) missing << "workdir"
    if (params.nuclear_chroms == null) missing << "nuclear_chroms"
    if (params.genome_chrom_sizes == null) missing << "genome_chrom_sizes"
    if (params.hisat2_index == null) missing << "hisat2_index"

    if (missing.size() > 0) {
        error "Check the config file. Missing required parameter(s): ${missing.join(', ')}"
    }
}

workflow COUNT_READS_WORKFLOW {
    take:
        bams_source // Can be directory (string) or channel from REMAP_WORKFLOW
        called_source // Can be directory (string) or channel from CALL_WORKFLOW
        merged_variants // Channel from CALL_WORKFLOW with merged VCF files
        output_dir
    
    main:
        // Check required parameters
        checkRequiredParams()
        
        // Determine input types
        def bams_is_string = bams_source instanceof String || bams_source instanceof groovy.lang.GString
        def called_is_string = called_source instanceof String || called_source instanceof groovy.lang.GString
        def merged_is_string = merged_variants instanceof String || merged_variants instanceof groovy.lang.GString
        
        // Handle all possible combinations
        if (bams_is_string && called_is_string) {
            // Both independent: use directory paths
            log.info "Running COUNT_READS_WORKFLOW with both inputs as directory paths"
            processed_bams_dir = bams_source
            called_dir = called_source
            merged_dir = "${called_dir}/Called_SNPs"

            // Create BAM file paths from metadata
            sample_ch = Channel.fromPath(params.meta_json)
                .flatMap { meta_file ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_bam_bed_entries ${meta_file} ${processed_bams_dir} ${called_dir}/filtered_SNPs"
                    def output = cmd.execute().text
                    def entries = new groovy.json.JsonSlurper().parseText(output)
                    return entries
                }
                .filter { it != null && it[0] != null }  // Remove any null entries
                .map { tuple ->
                    return tuple
                }
                
        } else if (bams_is_string && !called_is_string) {
            // BAMs from directory, called from channel
            log.info "Running COUNT_READS_WORKFLOW with BAMs from directory and called results from channel"
            processed_bams_dir = bams_source
            
            // Create BAM file paths from metadata
            sample_ch = Channel.fromPath(params.meta_json)
                .flatMap { meta_file ->
                    def cmd = "python3 ${projectDir}/bin/meta_parser.py create_bam_bed_entries ${meta_file} ${processed_bams_dir} ${params.outdir}/filtered_SNPs"
                    def output = cmd.execute().text
                    def entries = new groovy.json.JsonSlurper().parseText(output)
                    return entries
                }
                .filter { it != null && it[0] != null }  // Remove any null entries
                .map { tuple ->
                    return tuple
                }
                
                
        } else {
            // Both from channels (fully chained mode)
            log.info "Running COUNT_READS_WORKFLOW with both inputs from channels"
            
            // Combine processed BAMs with called results
            // bams_source: [sample_name, sample_id, bam_file, bai_file]
            // called_source: [sample_id, bed_file, tbi_file] from FILTER_VARIANTS
            
            // Transform BAM channel to have sample_id as key for joining
            bam_keyed = bams_source
                .map { sample_name, sample_id, bam_file, bai_file ->
                    [sample_id, sample_name, bam_file, bai_file]
                }
            
            // Transform called_source to have sample_id as key
            bed_keyed = called_source
                .map { sample_id, bed_file, tbi_file ->
                    [sample_id, bed_file, tbi_file]
                }
            
            // Join BAM and BED channels by sample_id
            sample_ch = bam_keyed
                .combine(bed_keyed)
                .filter { bam_sample_id, sample_name, bam_file, bai_file, bed_sample_id, bed_file, tbi_file ->
                    bam_sample_id == bed_sample_id  // Only keep matching sample_ids
                }
                .map { bam_sample_id, sample_name, bam_file, bai_file, bed_sample_id, bed_file, tbi_file ->
                    [sample_name, bam_sample_id, bam_file, bai_file, bed_file, tbi_file]
                }
        }

        // Handle merged variants - create genotype_channel based on input type
        if (merged_is_string) {
            // merged_variants is a string (directory path)
            merged_dir = "${merged_variants}/Called_SNPs"
            genotype_channel = Channel.fromPath("${merged_dir}/all.filtered.annotated.snps.vcf.gz")
                .map { vcf_file -> 
                    def csi_file = file("${vcf_file}.csi")
                    [vcf_file, csi_file]
                }
        } else {
            // merged_variants is a channel from CALL_WORKFLOW
            genotype_channel = merged_variants
        }

        // Handle H5 files - either use existing ones or generate new ones
        h5_tables = params.h5_files ? 
            Channel.fromPath("${params.h5_files}/*.h5").collect() :
            GENERATE_H5_TABLES(
                genotype_channel,
                file(params.genome_chrom_sizes)
            ).collect()

        // Remap BAM files
        remapped = REMAP_BAMFILES(
            sample_ch,
            file(params.nuclear_chroms),
            h5_tables
        )

        // Combine remapped BAMs with filtered variants for counting
        if (called_is_string) {
            // Called from directory: use bed files from specified directory
            count_input = remapped.map { sample_name, sample_id, bam, bai, passing_bam, passing_bai ->
                def bed = file("${called_source}/filtered_SNPs/${sample_id}.bed.gz")
                def tbi = file("${bed}.tbi")
                [sample_name, bam, bai, passing_bam, passing_bai, bed, tbi]
            }
        } else {
            // Called from channel: construct bed file paths from output directory
            count_input = remapped.map { sample_name, sample_id, bam, bai, passing_bam, passing_bai ->
                def bed = file("${params.outdir}/filtered_SNPs/${sample_id}.bed.gz")
                def tbi = file("${bed}.tbi")
                [sample_name, bam, bai, passing_bam, passing_bai, bed, tbi]
            }
        }

        // Count reads for each sample
        counted = COUNT_READS(count_input)

        // Collect all count files
        read_files = counted
            .map { sample_name, counts_bed, counts_tbi -> 
                [counts_bed, counts_tbi]
            }
            .flatten()
            .collect()


        chrom_ch = Channel.fromPath(params.autosomes)
            .splitCsv(header:true, sep:'\t')
            .map { row -> row.chroms }
        
        // Recode VCF - combine VCF with each chromosome to create 22 processes
        vcf_chrom_combinations = genotype_channel
            .combine(chrom_ch)
            .map { vcf_file, csi_file, chrom ->
                [vcf_file, csi_file, chrom]
            }
        
        recoded_vcfs = RECODE_VCF(
            vcf_chrom_combinations.map { vcf_file, csi_file, chrom -> vcf_file },
            read_files,
            vcf_chrom_combinations.map { vcf_file, csi_file, chrom -> csi_file },
            file(params.meta_json),
            vcf_chrom_combinations.map { vcf_file, csi_file, chrom -> chrom }
        )
        
        // Concatenate chromosome VCFs and split by sample
        sample_vcfs = CONCAT_AND_SPLIT_VCF(
            recoded_vcfs.collect(),
            file(params.meta_json)
        )
    
    emit:
        sample_vcfs
} 