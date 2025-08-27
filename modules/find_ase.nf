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
    # Create a file list of all VCF files
    printf '%s\n' ${input_files} > tmp_input_files.tsv
    echo "all vcf model"
    cat tmp_input_files.tsv
    
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
        python3 ${projectDir}/bin/meta_parser.py filter_vcf_files_by_group ${params.meta_json} \$group_id ${input_files} > \${group_id}_tmp_input.tsv
        
        echo "Filtered VCF files for group \$group_id:"
        cat \${group_id}_tmp_input.tsv
        
        # Run combine for this group
        mixalime combine --subname \$group_id -g \${group_id}_tmp_input.tsv "${project_name}"
        echo "Group \$group_id combined"
        
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