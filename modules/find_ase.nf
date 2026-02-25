process FIT_ASE_MODEL {
    tag "Fitting ASE model for ${project_name}"
    publishDir "${params.outdir}/find_ase", mode: 'copy'
    debug true
    
    input:
    path input_files
    val project_name
    val model_type
    
    output:
    tuple val(project_name), val(model_type),
        path("rna_project_${model_type}", type: 'dir'),
        path("rna_project_${model_type}.fit.lzma"),
        path("rna_project_${model_type}.init.lzma"),
        path("rna_project_${model_type}.json"),
        path("rna_project_${model_type}.test.lzma"),
        path("rna_project_${model_type}.comb.lzma")
    
    script:
    """    
    # Create a file list of all VCF files (escape \\n to avoid Groovy newline interpolation)
    printf '%s\\n' ${input_files} > tmp_input_files.tsv
    if [[ ! -s tmp_input_files.tsv ]]; then
        echo "ERROR: tmp_input_files.tsv is empty; no VCF files were provided to FIT_ASE_MODEL" >&2
        exit 1
    fi
    echo "all vcf model"
    cat tmp_input_files.tsv

    # Keep only non-empty VCF files; mixalime fit fails on empty datasets.
    : > tmp_input_nonempty.tsv
    while IFS= read -r vcf_file; do
        if zcat "\${vcf_file}" | awk 'BEGIN{ok=1} /^#/ {next} {ok=0; exit 0} END{exit ok}'; then
            echo "\${vcf_file}" >> tmp_input_nonempty.tsv
        fi
    done < tmp_input_files.tsv
    mv tmp_input_nonempty.tsv tmp_input_files.tsv
    if [[ ! -s tmp_input_files.tsv ]]; then
        echo "ERROR: all input VCF files are empty after per-sample filtering." >&2
        exit 1
    fi
    mapfile -t INPUT_VCFS < tmp_input_files.tsv
    
    mixalime create "${project_name}" "tmp_input_files.tsv"
    echo "Project created"
    mixalime fit "${project_name}" ${model_type}
    echo "model fitted"
    mixalime test "${project_name}"
    echo "model tested"

    mixalime combine --subname all_groups -g tmp_input_files.tsv "${project_name}"

    mixalime export all "${project_name}" "${project_name}"
    mixalime plot all "${project_name}" "${project_name}"
    echo "exported and plotted"
    
    # Get all group IDs from metadata and run combine for each group
    python3 ${projectDir}/bin/meta_parser.py get_group_ids ${params.meta_json} > group_ids.txt
    
    echo "Processing groups:"
    cat group_ids.txt
    
    while IFS= read -r group_id; do
        echo "Processing group: \$group_id"
        
        # Filter VCF files for this group
        python3 ${projectDir}/bin/meta_parser.py filter_vcf_files_by_group ${params.meta_json} \$group_id "\${INPUT_VCFS[@]}" > \${group_id}_tmp_input.tsv
        
        echo "Filtered VCF files for group \$group_id:"
        cat \${group_id}_tmp_input.tsv
        
        # Run combine only when group has at least one non-empty VCF
        if [[ -s \${group_id}_tmp_input.tsv ]]; then
            mixalime combine --subname \$group_id -g \${group_id}_tmp_input.tsv "${project_name}"
            echo "Group \$group_id combined"
        else
            echo "Skipping group \$group_id: no matching non-empty VCF files"
        fi
        
        # Clean up group-specific temp file
        rm \${group_id}_tmp_input.tsv
    done < group_ids.txt
    
    # Export all results (including group-specific ones)
    mixalime export all "${project_name}" "${project_name}"
    echo "All groups exported"
    
    # Clean up
    rm tmp_input_files.tsv group_ids.txt
    """
}