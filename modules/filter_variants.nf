process FILTER_VARIANTS {
    publishDir "${params.outdir}/filtered_SNPs", mode: 'copy'
    tag "${sample_id}"
    debug true  // Add debug output
    
    input:
        tuple val(sample_id), path(genotype_file), path(genotype_file_csi)

    output:
        tuple val(sample_id), path("${sample_id}.bed.gz"), path("${sample_id}.bed.gz.tbi")

    script:
    """
    bcftools query \
        -s ${sample_id} \
        -i'GT="1/1" || GT="1/0" || GT="0/1"' \
        -f'%CHROM\\t%POS0\\t%POS\\t%CHROM:%POS:%REF:%ALT\\t%REF\\t%ALT\\t[%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}]\\n' \
        ${genotype_file} \
    | awk -v OFS="\\t" \
        -v min_GQ=${params.min_GQ} -v min_AD=${params.min_AD} -v min_DP=${params.min_DP} \
        '\$8<min_GQ { next; } \$9<min_DP { next; } \
        (\$7=="0/1" || \$7=="1/0" || \$7=="0|1" || \$7=="1|0") && (\$10<min_AD || \$11<min_AD) { next; } \
        { print; }' \
    | sort-bed - \
    | grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn | grep -v KI | grep -v GL \
    | bgzip -c > ${sample_id}.bed.gz

    tabix -p bed ${sample_id}.bed.gz
    """
}