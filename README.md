## NEXT-scASV

### Abstract
NEXT-scASV (Nextflow pipeline for Allele-Specific Variant calling), a scalable Nextflow pipeline for calling allele-specific variants (ASVs) from 5′ single-cell RNA sequencing data. NEXT-scASV integrates established tools for alignment (HISAT2), reads deduplication  (umi_tools), variant calling (bcftools mpileup), reference bias correction (WASP), and advanced statistical modeling (MIXALIME). The pipeline consists of five main subflows: (i) Data Splitting, (ii) Alignment and Filtering, (iii) Variant Calling, (iv) Allelic Read Counting, and (v) ASV Calling. 

<div align="center">
<img src="scheme.png" alt="Pipeline scheme" width="50%" />
</div>

### 1) Installation

- **Create Conda env (recommended)**
  - Using the provided env file:
    ```bash
    conda env create -f as_env.yml
    conda activate as_env
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
    ```
  - If `conda activate` is unavailable in your shell:
    ```bash
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate as_env
    ```
  - Create manually (older setup):
    ```bash
    conda create -n as_env python==3.10 -y
    conda activate as_env
    conda install -y bcftools samtools hdf5 pytables numpy bioconda::sinto bioconda::bedops bioconda::hisat2 -c bioconda -c conda-forge
    python -m pip install --no-cache-dir pysam mixalime
    conda install -y -c conda-forge openjdk=17
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
    ```
  - Install WASP under `bin/`:
    ```bash
    git clone https://github.com/bmvdgeijn/WASP bin/WASP
    sed -i 's/np.int/int/g' bin/WASP/mapping/snptable.py
    ```
  - Nextflow requires Java 17+ (the env installs `openjdk=17`).

- **Add Nextflow to `bin/`**
  - Download Nextflow and place it in this repo's `bin` folder:
    ```bash
    curl -s https://get.nextflow.io > nextflow
    mv nextflow bin/nextflow
    chmod +x bin/nextflow
    ```

- **Docker**
  - A `Dockerfile` is provided. It installs Nextflow, HISAT2, and WASP.
  - Build the image from the project root:
    ```bash
    docker build -t asv-pipeline:latest .
    ```
  - Quick sanity check:
    ```bash
    docker run --rm asv-pipeline:latest nextflow -version
    ```
  - You can then run Nextflow with Docker enabled via `-with-docker` (the default `nextflow_dev.config` already enables Docker; adjust as needed).
  - The image sets `JAVA_HOME`, `PATH`, and `NXF_HOME` for non-root runs.

### 2) Pipeline steps (overview)

- **Split** (`split.nf`): barcode-based split of 5′ scRNA FASTQ into per-group files using `bin/split_fastq.py`.
- **Remap** (`remap.nf`): adapter trimming, HISAT2 alignment, UMI deduplication, and tagged BAM output.
- **Call** (`call.nf`): variant calling per chromosome (`bcftools` or `cellsnp-lite`) and merge/filtering.
- **Count reads** (`count_reads.nf`): remapping bias correction (WASP) and allele count tables.
- **Find ASE** (`find_ase.nf`): MIXALIME modeling and report/plot generation.

### 3) Supplementary files

- Load supplementary annotation/model files from the archive (to be provided later). Place annotation files under `bin/annotation/`.

### 4) Configure `nextflow.config`

Modify `nextflow.config` in the project root. Key parameters (see `nextflow_dev.config` for full list and defaults):

- **I/O and temp**
  - `params.input`: Input directory. 
  - `params.outdir`: Output directory.
  - `params.tmp_dir`: Temporary directory (writable, large enough).

- **Metadata**
  - `params.split_table`: Path to split table TSV (check `tests/split_table.tsv` for expample).
  - `params.meta_json`: Path to sample metadata JSON (check `tests/meta.json` for example).

- **References**
  - `params.genome`, `params.ref_genome`: reference FASTA paths (as required by your workflows).
  - `params.hisat2_index`: Prefix to HISAT2 index.
  - `params.dbsnp`: dbSNP VCF (bgzipped and indexed).
  - `params.autosomes`, `params.nuclear_chroms`, `params.genome_chrom_sizes`: text tables files with chrom lists and sizes (each chrom on the separate line).

- **Workflow toggles**
  - `params.run_split`, `params.run_remap`, `params.run_call`, `params.run_count_reads`, `params.run_find_ase`: set `true/false` to enable steps.

- **Calling/Filtering**
  - `params.min_GQ`, `params.min_DP`, `params.min_AD`: genotype quality/depth thresholds.
  - `params.variant_caller`: variant calling backend for the CALL step (`bcftools` or `cellsnp`).
    - For `cellsnp`, the pipeline runs **cellsnp-lite Mode 2b** per chromosome (no `-R`, no `-b`, `--cellTAG None --UMItag None`) and uses `params.6` / `params.cellsnp_minCOUNT`.

- **Count/ASE**
  - `params.h5_files_b`: optional H5 directory (set `null` to generate from scratch).
  - `params.model`: ASE model, e.g. `BetaNB`.
  - `params.annotation_meta_dir`: directory with annotation pickles; `null` to use defaults under `bin/annotation/`.

- **Process/Cluster**
    ```
    process {
        executor = 'slurm'
        queue = '<name>'
        memory = '64 GB'
        cpus = 8
    }
    ```

### 5) Run

6
```bash
./bin/nextflow run main.nf -c nextflow.config
```

- If running with Docker explicitly:
```bash
./bin/nextflow run main.nf -c nextflow.config -with-docker
```

Outputs will be created under `params.outdir` with subfolders for each step (e.g., `split/`, `align/`, `count_reads/`, `find_ase/`) and reports in `reports/` if enabled.

### 6) Test data

- Test FASTQ data are in `test_data/`.
- Use `KR_SGI_B007_L002_5GEX_H076_meta.json` and `KR_SGI_B007_L002_5GEX_H076_split.tsv`.
- The split table may include a leading numeric index column (as in the provided TSV).
- A minimal split-only config is available in `tests/nextflow_test.config`.
  ```bash
  ./bin/nextflow run main.nf -c tests/nextflow_test.config
  ```
  - Outputs will land in `test_out/` by default.

### 7) Full local run (after split)

- Fill in required reference paths in `tests/nextflow_full_local.config`.
- The following reference files are already included under `test_data/ref/`:
  - `GRCh38.primary_assembly.genome.chrom_sizes`
  - `GRCh38.primary_assembly.genome.all.txt` (used as `autosomes`)
  - `GRCh38.primary_assembly.genome.nuclear.txt`
- You still need to download:
  - GRCh38 primary assembly FASTA (use as `genome` and `ref_genome`)
  - HISAT2 GRCh38 index (set `hisat2_index` to the index prefix)
  - dbSNP VCF for GRCh38 (bgzipped + indexed; set `dbsnp`)
- Example download commands:
  ```bash
  # GRCh38 FASTA
  curl -L -o GRCh38.primary_assembly.genome.fa.gz https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip -f GRCh38.primary_assembly.genome.fa.gz

  # HISAT2 GRCh38 index (prebuilt)
  mkdir -p ref/hisat2_grch38
  curl -L -o ref/hisat2_grch38/grch38_hisat2.tar.gz \
    https://genome-idx.s3.amazonaws.com/hisat/grch38_hisat2.tar.gz
  tar -xzf ref/hisat2_grch38/grch38_hisat2.tar.gz -C ref/hisat2_grch38

  # Or build locally from the FASTA (takes longer)
  hisat2-build ref/GRCh38.primary_assembly.genome.fa ref/hisat2_grch38/GRCh38

  # dbSNP VCF (GRCh38) + index
  curl -L -o dbsnp_155_GRCh38.vcf.gz \
    https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz
  tabix -f -p vcf dbsnp_155_GRCh38.vcf.gz

  # If you already have 00-common_all.vcf.gz, use it for `dbsnp` instead
  tabix -f -p vcf 00-common_all.vcf.gz
  ```
- Run with:
  ```bash
  ./bin/nextflow run main.nf -c tests/nextflow_full_local.config
  ```
