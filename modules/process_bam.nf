process PROCESS_BAM {
    publishDir "${params.outdir}/align", mode: 'copy'
    tag "${sample_name}"
    debug true  // Add debug output

    input:
        tuple val(sample_name), path(bam)

    output:
        tuple val(sample_name), path("${sample_name}_marked_filtered.bam"), path("${sample_name}_marked_filtered.bam.bai")

    script:
    """
    samtools sort -n -@ ${task.cpus} -o ${sample_name}_sorted.bam -O bam ${bam}

    # Filter reads
    python3 ${projectDir}/bin/bwa_script/filter_reads.py \
        ${sample_name}_sorted.bam \
        ${sample_name}_marked.bam \
        ${params.genome}

    samtools sort -@ ${task.cpus} -l0 ${sample_name}_marked.bam \
        | samtools view -b -F 512 - \
        > ${sample_name}_marked_filtered.bam

    samtools index ${sample_name}_marked_filtered.bam
    """
} 