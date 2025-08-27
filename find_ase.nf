#!/usr/bin/env nextflow

include { FIT_ASE_MODEL } from './modules/find_ase'
include { ANNOTATE_PVALUE_TABLES } from './modules/annotation'

workflow FIND_ASE_WORKFLOW {
    take:
        vcf_source // Can be vcf_dir (string) or channel from COUNT_READS_WORKFLOW
        output_dir
    
    main:
        // Determine if input is a directory path or a channel from previous workflow
        if (vcf_source instanceof String || vcf_source instanceof groovy.lang.GString) {
            // Independent mode: use directory path to get sample VCF files
            log.info "Running FIND_ASE_WORKFLOW in independent mode with VCF directory: ${vcf_source}"
            
            // Get sample VCF files from the final_vcf directory
            sample_vcfs = Channel.fromPath("${vcf_source}/*.vcf.gz")
                
        } else {
            // Chained mode: use channel from COUNT_READS_WORKFLOW
            log.info "Running FIND_ASE_WORKFLOW in chained mode with input channel"
            
            // The vcf_source is a channel of sample VCF files from CONCAT_AND_SPLIT_VCF
            sample_vcfs = vcf_source
        }
        
        // Use sample VCF files directly for ASE analysis
        project_name = "rna_project_${params.model}"
        
        // Collect all VCF files into a single list before passing to FIT_ASE_MODEL
        sample_vcfs
            .collect()
            .set { all_vcfs }
        
        // FIT_ASE_MODEL now handles everything including group-specific combining
        fit_results = FIT_ASE_MODEL(all_vcfs, project_name, params.model)
        
        // Extract project directory from fit_results for annotation
        project_dir = fit_results
            .map { proj_name, model_type, proj_dir, fit_file, init_file, json_file, test_file, comb_file ->
                return proj_dir  // This is the project directory path
            }
        
        // Run annotation on the project directory
        ANNOTATE_PVALUE_TABLES(project_dir)
        
    emit:
        fit_results = fit_results
        annotated_tables = ANNOTATE_PVALUE_TABLES.out.annotated_tables
        plots = ANNOTATE_PVALUE_TABLES.out.plots
}