process REMAP_FASTQ {
    // publishDir "${out_path}", mode: 'move'
    tag "${sample_name}"
    debug true  // Add debug output

    input:
        tuple val(sample_name), val(sample_id), path(r1_fastq), path(r2_fastq), val(out_path)

    output:
        tuple val(sample_name), val(sample_id), path("${sample_name}_dedup_with_tag.bam"), path("${sample_name}_dedup_with_tag.bam.bai")

    script:
    // def (base_name, group_name) = sample_name.split('\\+')
    """
    # Trim adapters add after debag: --trimmed-only
    cutadapt -g TTTCTTATATGGG \
        -o ${sample_name}_R1_trimm.fastq.gz \
        -p ${sample_name}_R2_trimm.fastq.gz -j ${task.cpus} \
        ${r1_fastq} ${r2_fastq}

    # Map with HISAT2
    /mnt/flashgpu/lab2/Allele-specifity/hisat2/hisat2 --very-sensitive -p ${task.cpus} --seed 13 \
        --rg-id "ID:${sample_name}" --rg "SM:${sample_id}" \
        -x ${params.hisat2_index} \
        -1 ${sample_name}_R1_trimm.fastq.gz \
        -2 ${sample_name}_R2_trimm.fastq.gz \
        | samtools view -b \
        > ${sample_name}_remaped.bam

    # Sort and filter
    samtools sort -@ ${task.cpus} -l0 ${sample_name}_remaped.bam \
        | samtools view -@ ${task.cpus} -b -F 512 - \
        > ${sample_name}_remaped_sorted.bam
    
    samtools index -@ ${task.cpus} ${sample_name}_remaped_sorted.bam

    # UMI deduplication
    umi_tools dedup \
        --chimeric-pairs=discard \
        --unpaired-reads=discard \
        --paired \
        --per-cell \
        --log=${sample_name}_dedup.log \
        --temp-dir ${params.tmp_dir} \
        --stdin ${sample_name}_remaped_sorted.bam \
        > ${sample_name}_dedup.bam

    # Add read group tags
    samtools addreplacerg \
        -r "ID:${sample_name}" \
        -r "SM:${sample_id}" \
        -o ${sample_name}_dedup_with_tag.bam \
        ${sample_name}_dedup.bam

    samtools index ${sample_name}_dedup_with_tag.bam

    # Cleanup intermediate files
    rm ${sample_name}_R*_trimm.fastq.gz ${sample_name}_remaped*.bam ${sample_name}_dedup.bam
    """
}