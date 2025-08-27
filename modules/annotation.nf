#!/usr/bin/env nextflow

process ANNOTATE_PVALUE_TABLES {
    tag "Annotating p-value tables in ${project_dir}"
    publishDir "${params.outdir}/annotated", mode: 'copy'
    debug true
    
    input:
    path project_dir
    
    output:
    path "*_annotated.tsv", emit: annotated_tables
    path "*_annotation_plot.png", emit: plots
    
    script:
    """
    echo "Processing all TSV files in ${project_dir}/pvalues/"

    python3 ${projectDir}/bin/annotation/annotation_pipeline.py \
        --pvalues_dir "${project_dir}/pvalues" \
        --annotation_meta_dir ${params.annotation_meta_dir ?: "${projectDir}/bin/annotation"}
    """
}