process CONCAT_VCF {
    cpus { params.concat_vcf_cpus ?: 2 }

    input:
        path vcf_files

    output:
        path "concatenated_allele_counts.vcf.gz"

    script:
    """
    set -euo pipefail

    # Sort VCF files by chromosome order for proper concatenation
    ls allele_counts.*.vcf.gz | sort -V > vcf_list.txt

    # Concatenate all chromosome VCFs
    bcftools concat -f vcf_list.txt -O z -o concatenated_allele_counts.vcf.gz
    bcftools index -t concatenated_allele_counts.vcf.gz
    """
}

process VCF_SAMPLE_LIST {
    input:
        path concat_vcf

    output:
        stdout

    script:
    """
    set -euo pipefail
    bcftools query -l ${concat_vcf} \
        | tr -d '\r' \
        | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*\$//'
    """
}

process SPLIT_VCF {
    publishDir "${params.outdir}/vcf_by_sample", mode: 'copy'
    tag { sample }

    input:
        tuple path(concat_vcf), val(sample)

    output:
        path "${sample}.vcf.gz"

    script:
    """
    set -euo pipefail
    echo "Processing sample: ${sample}"
    bcftools view -c1 --output-type v -s "${sample}" ${concat_vcf} \
        | bcftools filter -i 'FORMAT/AD[0:1] >= 5 && FORMAT/AD[0:0] >= 5' --output-type z -o "${sample}.vcf.gz"
    """
}