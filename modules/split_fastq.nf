process SPLIT_FASTQ {
    publishDir "${output_dir}", mode: 'copy'
    tag "${sample}:${group}"
    // errorStrategy 'ignore'  // Continue workflow even if this process fails
    debug true  // Add debug output
    
    input:
        tuple val(sample), val(group), path(split_table)
        val(input_dir)
        val(output_dir)

    output:
        tuple val(sample), val(group), path("${sample}-${group}*.R1.fastq"), path("${sample}-${group}*.R2.fastq")

    script:
    """
    # Get file path from meta.json
    FILE_PATH=\$(python3 -c "import json; f=open('${params.meta_json}'); meta=json.load(f); print(meta['${sample}']['path'])")
    
    # Run with error logging
    python3 ${projectDir}/bin/split_fastq.py \
        --sample ${sample} \
        --group ${group} \
        --split_table ${split_table} \
        --input_dir ${input_dir} \
        --file_path \${FILE_PATH} \
        --outdir \${PWD}
    
    """
} 