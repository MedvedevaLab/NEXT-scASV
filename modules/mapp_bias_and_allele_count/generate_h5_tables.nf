process GENERATE_H5_TABLES {

    publishDir "${params.outdir}/h5_tables", mode: 'copy'

    input:
        tuple path(vcf_file), path(vcf_csi)
        path chrom_sizes

    output:
        path '*.h5'

    script:
    """
    tabix -p vcf ${vcf_file}
    chroms=("\$(tabix -l ${vcf_file})")
    for chrom in \${chroms[@]}; do
        echo \${chrom}
        bcftools view -r \${chrom} -Oz ${vcf_file} > \${chrom}.chrom.vcf.gz
        bcftools index \${chrom}.chrom.vcf.gz
    done

    gzip -c ${chrom_sizes} > chrom_sizes.txt.gz

    ${projectDir}/bin/WASP/snp2h5/snp2h5 --chrom chrom_sizes.txt.gz \
        --format vcf \
        --haplotype haplotypes.h5 \
        --snp_index snp_index.h5 \
        --snp_tab snp_tab.h5 \
        *.chrom.vcf.gz
    """
}