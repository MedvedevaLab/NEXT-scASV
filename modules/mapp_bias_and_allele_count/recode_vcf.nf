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
    
    # Keep only entries that have counts files staged in this work dir
    : > sample_map.filtered.tsv
    while IFS=\$'\\t' read -r sample genotype counts; do
        if [[ -f "\${counts}" ]]; then
            printf "%s\\t%s\\t%s\\n" "\${sample}" "\${genotype}" "\${counts}" >> sample_map.filtered.tsv
        else
            echo "WARNING: counts file not found for \${sample} (\${counts})" >&2
        fi
    done < sample_map.tsv
    mv sample_map.filtered.tsv sample_map.tsv

    if [[ ! -s sample_map.tsv ]]; then
        echo "ERROR: sample_map.tsv is empty after filtering missing counts files." >&2
        exit 1
    fi

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
