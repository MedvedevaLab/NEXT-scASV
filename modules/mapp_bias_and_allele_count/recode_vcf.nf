process RECODE_VCF {
    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    // Split by chromosome for parallelization
    input:
        path vcf_file
        path read_files
        path csi_file
        path meta_json
        val chrom

    output:
        path "allele_counts.${chrom}.vcf.gz*"

    script:
    """
    
    
    # Build both mappings and auto-select the one matching samples present in VCF header.
    python3 ${projectDir}/bin/meta_parser.py create_sample_map_entries ${meta_json} "sample_name" > sample_map.sample_name.tsv
    python3 ${projectDir}/bin/meta_parser.py create_sample_map_entries ${meta_json} "sample_id" > sample_map.sample_id.tsv

    python3 - "${vcf_file}" <<'PY'
import pathlib
import sys
import pysam

vcf_path = sys.argv[1]
vcf_samples = set(pysam.VariantFile(vcf_path, mode='r').header.samples)

def norm_sample(x: str) -> str:
    # Some historical VCFs have malformed sample labels like "[SAMPLE,".
    # Normalize both sides for matching, then write back the raw VCF sample label.
    return x.strip().strip(",").strip("[]").strip()

norm_to_raw = {}
for raw in vcf_samples:
    norm_to_raw.setdefault(norm_sample(raw), raw)

candidates = [
    ("sample_name", "sample_map.sample_name.tsv"),
    ("sample_id", "sample_map.sample_id.tsv"),
]

best_key = None
best_rows = []

for key, map_path in candidates:
    rows = []
    with open(map_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            sample, genotype_id, counts = line.split("\\t")
            if not pathlib.Path(counts).exists():
                continue
            raw_genotype = norm_to_raw.get(norm_sample(genotype_id))
            if raw_genotype is None:
                continue
            rows.append((sample, raw_genotype, counts))
    if len(rows) > len(best_rows):
        best_key = key
        best_rows = rows

if not best_rows:
    raise SystemExit(
        "ERROR: Could not build non-empty sample_map.tsv that matches VCF samples and existing counts files."
    )

with open("sample_map.tsv", "w", encoding="utf-8") as out:
    for sample, genotype_id, counts in best_rows:
        out.write(f"{sample}\\t{genotype_id}\\t{counts}\\n")

print(f"Using genotype_key={best_key}; matched_rows={len(best_rows)}", file=sys.stderr)
PY

    if [[ ! -s sample_map.tsv ]]; then
        echo "ERROR: sample_map.tsv is empty after filtering missing counts files." >&2
        exit 1
    fi

    echo "Processing chromosome: ${chrom}"
    ${projectDir}/bin/recode_vcf.py \
        ${vcf_file} \
        sample_map.tsv \
        allele_counts.${chrom}.vcf \
        --chrom ${chrom}

    bgzip -c allele_counts.${chrom}.vcf > allele_counts.${chrom}.vcf.gz
    bcftools index allele_counts.${chrom}.vcf.gz
    """
}