process REMAP_BAMFILES {
    publishDir "${params.outdir}/reremapped", mode: 'copy'
    tag "${sample_name}"

    input:
        tuple val(sample_name), val(sample_id), path(bam), path(bai), path(bed), path(tbi)
        path nuclear_chroms
        path h5_files

    output:
        tuple val(sample_name), val(sample_id), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), 
              path("${sample_name}.passing.bam"), path("${sample_name}.passing.bam.bai")

    script:
    """
    samtools sort -@ ${task.cpus} -o sorted.bam -O bam ${bam}

    python3 ${projectDir}/bin/WASP/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir \${PWD} \
        --snp_tab snp_tab.h5 \
        --snp_index snp_index.h5 \
        --haplotype haplotypes.h5 \
        --samples ${sample_id} \
        sorted.bam

    hisat2 --very-sensitive -p ${task.cpus} --seed 13 \
        -x ${params.hisat2_index} \
        -1 sorted.remap.fq1.gz \
        -2 sorted.remap.fq2.gz \
        | samtools sort -n -@ ${task.cpus} \
        > remapped.bam

    python3 ${projectDir}/bin/bwa_script/filter_reads.py \
        remapped.bam \
        remapped.marked.bam \
        ${nuclear_chroms}

    samtools sort -@ ${task.cpus} -l0 remapped.marked.bam \
        | samtools view -b -F 512 - \
        > remapped.marked.filtered.bam

    python3 ${projectDir}/bin/WASP/mapping/filter_remapped_reads.py \
        sorted.to.remap.bam \
        remapped.marked.filtered.bam \
        remapped.keep.bam

    samtools merge -f reads.before.bam sorted.bam
    samtools sort -@${task.cpus} -o reads.before.sorted.bam reads.before.bam

    samtools merge -f reads.passing.bam remapped.keep.bam sorted.keep.bam
    samtools sort -@${task.cpus} -o reads.passing.sorted.bam reads.passing.bam

    mv reads.before.sorted.bam ${sample_name}.bam
    samtools index ${sample_name}.bam
    mv reads.passing.sorted.bam ${sample_name}.passing.bam
    samtools index ${sample_name}.passing.bam

    rm sorted.bam
    rm remapped.bam
    rm remapped.marked.bam
    rm remapped.marked.filtered.bam
    rm remapped.keep.bam
    rm reads.before.bam
    rm reads.passing.bam
    rm sorted.to.remap.bam
    """
}