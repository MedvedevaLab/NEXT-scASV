process COUNT_READS {
    publishDir "${params.outdir}/count_reads", mode: 'copy'
    tag "${sample_name}"

    input:
        tuple val(sample_name), path(bam), path(bai), path(passing_bam), path(passing_bai), path(bed), path(tbi)

    output:
        tuple val(sample_name), path("${sample_name}.counts.bed.gz"), path("${sample_name}.counts.bed.gz.tbi")

    script:
    """
    ${projectDir}/bin/count_tags_pileup.py \
        ${bed} ${bam} ${passing_bam} \
    | sort-bed - | bgzip -c > ${sample_name}.counts.bed.gz

    tabix -p bed ${sample_name}.counts.bed.gz
    """
}