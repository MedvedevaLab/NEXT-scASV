process MERGE_VARIANTS {
    publishDir "${params.outdir}/Called_SNPs", mode: 'copy'
    
    input:
        path '*'

    output:
        tuple path("all.filtered.annotated.snps.vcf.gz"), 
              path("all.filtered.annotated.snps.vcf.gz.csi")
        // tuple path("all.filtered.annotated.vcf.gz"), path("all.filtered.annotated.vcf.gz.csi"), path("all.filtered.annotated.snps.vcf.gz"), path("all.filtered.annotated.snps.vcf.gz.csi")

    script:
    """
    # Create sorted list of VCF files
    ls *.filtered.annotated.vcf.gz > vcf_files.txt
    sort -k1,1 -V vcf_files.txt > vcf_files_sorted.txt

    # Merge variants
    bcftools concat \
        --output-type z \
        -f vcf_files_sorted.txt \
        > all.filtered.annotated.vcf.gz

    bcftools index all.filtered.annotated.vcf.gz

    # Filter for biallelic SNPs
    bcftools view \
        -m2 -M2 -v snps \
        --output-type z \
        all.filtered.annotated.vcf.gz \
        > all.filtered.annotated.snps.vcf.gz

    bcftools index all.filtered.annotated.snps.vcf.gz
    """
} 