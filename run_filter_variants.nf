#!/usr/bin/env nextflow

// Parameters
params.sample_name = "KR_SGI_B007_L001_5GEX_H086"
params.bam = ""
params.bai = ""
params.genotype_file = "/mnt/flashgpu/lab2/change_meta/ase_pipeline/tests/real_var3/results/Called_SNPs/all.filtered.annotated.snps.vcf.gz"
params.min_DP = 0
params.min_AD = 0
params.min_GQ = 0
params.outdir = "results"


process FILTER_VARIANTS {
    publishDir "filter_check", mode: 'move'
    tag "${sample_name}"
    debug true  // Add debug output
    
    input:
        tuple val(sample_name), val(bam), val(bai), path(genotype_file), val(min_DP), val(min_AD), val(min_GQ)

    output:
        tuple val(sample_name), path("${sample_name}.bed.gz"), path("${sample_name}.bed.gz.tbi")

    script:
    """
    bcftools query \
        -s ${sample_name} \
        -i'GT="1/1" || GT="1/0" || GT="0/1"' \
        -f'%CHROM\\t%POS0\\t%POS\\t%CHROM:%POS:%REF:%ALT\\t%REF\\t%ALT\\t[%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}]\\n' \
        ${genotype_file} \
    | awk -v OFS="\\t" \
        -v min_GQ=${min_GQ} -v min_AD=${min_AD} -v min_DP=${min_DP} \
        '\$8<min_GQ { next; } \$9<min_DP { next; } \
        (\$7=="0/1" || \$7=="1/0" || \$7=="0|1" || \$7=="1|0") && (\$10<min_AD || \$11<min_AD) { next; } \
        { print; }' \
    | sort-bed - \
    | grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn | grep -v KI | grep -v GL \
    | bgzip -c > ${sample_name}.bed.gz

    tabix -p bed ${sample_name}.bed.gz
    """
}

workflow {
    // Create input channel
    input_ch = Channel.of([
        "KR_SGI_B007_L001_5GEX_H086",
        "",
        "",
        file("/mnt/flashgpu/lab2/change_meta/ase_pipeline/tests/real_var3/results/Called_SNPs/all.filtered.annotated.snps.vcf.gz"),
        0,
        0,
        0
    ])
    
    // Run the process
    FILTER_VARIANTS(input_ch)
    
    // Print completion message
    FILTER_VARIANTS.out.view { sample, bed, tbi ->
        "Completed filtering for sample: ${sample}"
    }
} 