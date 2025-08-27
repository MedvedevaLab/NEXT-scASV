process CONCAT_AND_SPLIT_VCF {
    publishDir "${params.outdir}/vcf_by_sample", mode: 'copy'
    
    input:
        path vcf_files
        path meta_json

    output:
        path "*.vcf.gz"

    script:
    """
    # Sort VCF files by chromosome order for proper concatenation
    ls allele_counts.*.vcf.gz | sort -V > vcf_list.txt
    
    # Concatenate all chromosome VCFs
    bcftools concat -f vcf_list.txt -O z -o concatenated_allele_counts.vcf.gz
    bcftools index concatenated_allele_counts.vcf.gz
    
    # Extract sample names from VCF header
    bcftools query -l concatenated_allele_counts.vcf.gz > sample_names.txt
    
    # Split VCF by sample
    while read sample; do
        echo "Processing sample: \$sample"
        bcftools view -c1 --output-type v -s \$sample concatenated_allele_counts.vcf.gz \
            | bcftools filter  -i 'AD[0:1] >= 5 && AD[0:0] >= 5'  --output-type z -o \${sample}.vcf.gz
        
    done < sample_names.txt
    
    # Clean up temporary files
    rm concatenated_allele_counts.vcf.gz*
    rm vcf_list.txt
    rm sample_names.txt
    """
} 