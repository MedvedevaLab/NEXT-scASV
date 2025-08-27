process RECODE_VCF {
    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    // Split by chromosome for parallelization
    input:
        path vcf_file
        path read_files
        path csi_file
        path meta_json
        val chrom

    output:
        path "allele_counts.${chrom}.vcf.gz*"

    script:
    """
    
    
    # Create sample map file using meta_parser
    python3 ${projectDir}/bin/meta_parser.py create_sample_map_entries ${meta_json} > sample_map.tsv

    echo "Processing chromosome: ${chrom}"
    ${projectDir}/bin/recode_vcf.py \
        ${vcf_file} \
        sample_map.tsv \
        allele_counts.${chrom}.vcf \
        --chrom ${chrom}

    bgzip -c allele_counts.${chrom}.vcf > allele_counts.${chrom}.vcf.gz
    bcftools index allele_counts.${chrom}.vcf.gz
    """
}