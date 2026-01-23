process SPLIT_FASTQ {
    publishDir "${output_dir}", mode: 'copy'
    tag "${sample_id}:${group}"
    // errorStrategy 'ignore'  // Continue workflow even if this process fails
    debug true  // Add debug output
    
    input:
        tuple val(sample_id), val(group), path(split_table)
        val(input_dir)
        val(output_dir)

    output:
        tuple val(sample_id), val(group), path("${sample_id}-${group}*.R1.fastq"), path("${sample_id}-${group}*.R2.fastq")

    script:
    """
    # Get file path from meta.json
    FILE_PATH=\$(python3 -c "import json; f=open('${params.meta_json}'); meta=json.load(f); print(meta['${sample_id}']['path'])")
    
    # Run with error logging
    python3 ${projectDir}/bin/split_fastq.py \
        --sample_id ${sample_id} \
        --group ${group} \
        --split_table ${split_table} \
        --input_dir ${input_dir} \
        --file_path \${FILE_PATH} \
        --outdir \${PWD}
    
    """
} 