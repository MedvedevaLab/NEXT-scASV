process CALL_VARIANTS_CELLSNP {
    // Cellsnp-lite Mode 2b: well-based single cells / bulk without given SNPs (-R/-b not used)
    // Emits a per-chromosome, dbSNP-annotated VCF with per-sample genotypes (GT/AD/DP, etc.)
    debug true

    input:
        val chrom
        path bams
        val sample_ids

    output:
        path("${chrom}.filtered.annotated.vcf.gz*")

    script:
    """
    set -euo pipefail

    # Create input lists in matching order
    : > bam_list.txt
    for bam in ${bams}; do
        realpath "\$bam" >> bam_list.txt
    done

    : > sample_list.txt
    for sid in ${sample_ids}; do
        echo "\$sid" >> sample_list.txt
    done

    # Run cellsnp-lite (Mode 2b): no -R, no -b, use sample list + BAM list
    # Use -f so REF is fetched from the reference genome; infer ALT from data.
    mkdir -p cellsnp_out

    cellsnp-lite \\
        -S bam_list.txt \\
        -i sample_list.txt \\
        -O cellsnp_out \\
        -p ${task.cpus} \\
        --minMAF ${params.cellsnp_minMAF ?: 0.1} \\
        --minCOUNT ${params.cellsnp_minCOUNT ?: 100} \\
        --chrom ${chrom} \\
        --cellTAG None \\
        --UMItag None \\
        --gzip \\
        --genotype \\
        -f ${params.ref_genome}

    # Convert/normalize into the same naming convention as the bcftools caller so downstream steps stay unchanged
    # Keep only biallelic SNPs and annotate dbSNP IDs.
    bcftools view \\
        -m2 -M2 -v snps \\
        --output-type z \\
        -o ${chrom}.filtered.vcf.gz \\
        cellsnp_out/cellSNP.cells.vcf.gz

    bcftools index ${chrom}.filtered.vcf.gz

    bcftools annotate -r ${chrom} \\
        -a ${params.dbsnp} \\
        --columns ID \\
        --output-type z \\
        ${chrom}.filtered.vcf.gz \\
    > ${chrom}.filtered.annotated.vcf.gz

    bcftools index ${chrom}.filtered.annotated.vcf.gz

    # Cleanup
    rm -rf cellsnp_out ${chrom}.filtered.vcf.gz ${chrom}.filtered.vcf.gz.csi bam_list.txt sample_list.txt
    """
}


