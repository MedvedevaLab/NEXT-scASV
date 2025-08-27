process CALL_VARIANTS {
    // publishDir "${params.outdir}/Called_SNPs", mode: 'copy'
    debug true  // Add debug output
    
    input:
        val chrom
        path bams

    output:
        path("${chrom}.filtered.annotated.vcf.gz*")

    script:
    """
    # Create BAM list file
    for bam in ${bams}; do
        realpath \$bam >> bam_list.txt
    done

    # Call variants
    bcftools mpileup -r ${chrom} \
        --fasta-ref ${params.ref_genome} \
        --redo-BAQ \
        --adjust-MQ 50 \
        --gap-frac 0.05 \
        --max-depth 1000 \
        --max-idepth 2000 \
        --annotate FORMAT/DP,FORMAT/AD \
        --bam-list bam_list.txt \
        --output-type u \
        --threads ${task.cpus} \
    | bcftools call \
        --keep-alts \
        --multiallelic-caller \
        --format-fields GQ \
        --output-type v \
        --threads ${task.cpus} \
    | bcftools filter -i"INFO/DP>=10" \
        --output-type z - \
    | bcftools norm \
        --check-ref x -m - \
        --fasta-ref ${params.ref_genome} \
    | bcftools annotate -x ^INFO/DP \
    | bcftools +fill-tags -- \
    > ${chrom}.filtered.vcf

    # Compress and index
    bcftools view ${chrom}.filtered.vcf -Oz -o ${chrom}.filtered.vcf.gz
    bcftools index ${chrom}.filtered.vcf.gz

    # Annotate with dbSNP
    bcftools annotate -r ${chrom} \
        -a ${params.dbsnp} \
        --columns ID \
        --output-type z \
        ${chrom}.filtered.vcf.gz \
    > ${chrom}.filtered.annotated.vcf.gz

    bcftools index ${chrom}.filtered.annotated.vcf.gz

    # Cleanup
    rm ${chrom}.filtered.vcf*
    """
} 