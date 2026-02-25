# NEXT-scASV: step-by-step tutorial

This tutorial is aligned with `nextflow_full_local.config` and `KR_SGI_B007_L002_5GEX_H076_meta.json`.
Goal: user downloads the Kaggle dataset and can run without manual path edits.

Dataset:
`https://www.kaggle.com/datasets/andreylalaley/next-scasv-tutorial-dataset`

## 0) Prepare folders

```bash
mkdir -p test_out/split test_tmp ref/hisat2_grch38
```

## 1.0) Download data with `curl`

Download:

```bash
curl -L -o next-scasv-tutorial-dataset.zip  https://www.kaggle.com/api/v1/datasets/download/andreylalaley/next-scasv-tutorial-dataset
```

## 1.1) Unpack dataset and place split FASTQ files

```bash
unzip -o next-scasv-tutorial-dataset.zip
tar -xzf splitted_test.tar.gz -C test_out/split
```

Expected files in project root after unzip:

- `splitted_test.tar.gz`
- `00-common_all.vcf` (or `00-common_all.vcf.gz`)
- `00-common_all.vcf.gz.tbi`
- `GRCh38.primary_assembly.genome.all.txt`
- `GRCh38.primary_assembly.genome.autosomes.txt`
- `GRCh38.primary_assembly.genome.chrom_sizes`
- `GRCh38.primary_assembly.genome.nuclear.txt`
- `KR_SGI_B007_L002_5GEX_H076_meta.json`
- `KR_SGI_B007_L002_5GEX_H076_split.tsv`
- `nextflow_full_local.config`

## 1.2) Normalize dbSNP file name (if needed)

`nextflow_full_local.config` uses `00-common_all.vcf.gz`.
If Kaggle provided only `00-common_all.vcf`, compress and index it:

```bash
if [ -f "00-common_all.vcf" ] && [ ! -f "00-common_all.vcf.gz" ]; then
  bgzip -c 00-common_all.vcf > 00-common_all.vcf.gz
  bcftools index -t 00-common_all.vcf.gz
fi
```

## 1.3) HISAT2 reference index

If `ref/hisat2_grch38/GRCh38.*.ht2` already exists, skip this step.
Otherwise:

```bash
wget -O ref/hisat2_grch38/grch38_genome.tar.gz \
  https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzf ref/hisat2_grch38/grch38_genome.tar.gz -C ref/hisat2_grch38
```

`nextflow_full_local.config` expects:
`params.hisat2_index = "${projectDir}/ref/hisat2_grch38/GRCh38"`

## 1.4) Reference FASTA for `genome` and `ref_genome`

Both parameters should point to the same GRCh38 FASTA file:

- `params.genome = "${projectDir}/ref/GRCh38.primary_assembly.genome.fa"`
- `params.ref_genome = "${projectDir}/ref/GRCh38.primary_assembly.genome.fa"`

```bash
wget -O ref/GRCh38.primary_assembly.genome.fa.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip -f ref/GRCh38.primary_assembly.genome.fa.gz
samtools faidx ref/GRCh38.primary_assembly.genome.fa
```

## 2) Set up runtime

### 2A) Conda

```bash
conda env create -f as_env.yml
conda activate as_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
```

Install Nextflow into `bin/`:

```bash
curl -s https://get.nextflow.io > nextflow
mv nextflow bin/nextflow
chmod +x bin/nextflow
```

### 2B) Docker (alternative)

```bash
docker build -t next-scasv:latest .
```

## 3) Run

### 3.1) Use preconfigured file

No copy required. This config is already aligned:
`nextflow_full_local.config`

### 3.2) Start with Conda

```bash
./bin/nextflow run main.nf -c nextflow_full_local.config
```

### 3.3) Start in Docker container (alternative)

```bash
docker run --rm -it \
  -v /mnt/flashgpu/lab2/Allele-specifity/rebultle_version/NEXT-scASV:/app \
  next-scasv:latest \
  nextflow run main.nf -c nextflow_full_local.config
```

## Result check

Outputs are written to:
- `test_out`
- `test_out/work`
- reports under `test_out/reports` (if enabled)
