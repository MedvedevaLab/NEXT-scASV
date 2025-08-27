process BAM_TO_FASTQ {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    tag "${sample_name}"

    input:
        tuple val(sample_name), path(bam)

    output:
        // tuple val(sample_name), path("${sample_name}_R1.fastq.gz"), path("${sample_name}_R2.fastq.gz")
        path("multiqc_report.html"), emit: multiqc_report

    script:
    """
    # Convert BAM to FASTQ
    samtools fastq -1 ${sample_name}_R1.fastq.gz -2 ${sample_name}_R2.fastq.gz ${bam}

    # Run FastQC on the FASTQ files
    fastqc ${sample_name}_R1.fastq.gz ${sample_name}_R2.fastq.gz -o fastqc_results

    # Run MultiQC
    multiqc fastqc_results -o .
    """
} 