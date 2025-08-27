## NEXT-ASV Pipeline

End-to-end Nextflow pipeline for ASV detection

### 1) Installation

- **Create Conda env (recommended)**
  - Using the provided env file:
    ```bash
    conda env create -f as_env.yml
    conda activate as_env
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
    ```
  - Create manually (older setup):
    ```bash
    conda create -n as_env python==3.10 -y
    conda activate as_env
    conda install -y bcftools samtools hdf5 pytables numpy bioconda::sinto bioconda::bedops -c bioconda -c conda-forge
    python -m pip install --no-cache-dir pysam mixalime
    conda install -y -c conda-forge openjdk=11.0.20
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
    ```
  - Follow `https://github.com/bmvdgeijn/WASP` installation.
  - In `WASP/mapping/snptable.py` change `np.int` → `int` (NumPy ≥ 1.20 deprecates `np.int`).

- **Add Nextflow to `bin/`**
  - Download Nextflow and place it in this repo's `bin` folder:
    ```bash
    curl -s https://get.nextflow.io > nextflow
    mv nextflow bin/nextflow
    chmod +x bin/nextflow
    ```

- **Docker**
  - A `Dockerfile` is provided. Build the image from the project root:
    ```bash
    docker build -t asv-pipeline:latest .
    ```
  - You can then run Nextflow with Docker enabled via `-with-docker` (the default `nextflow_dev.config` already enables Docker; adjust as needed).

### 2) Supplementary files

- Load supplementary annotation/model files from the archive (to be provided later). Place annotation files under `bin/annotation/`.

### 3) Configure `nextflow.config`

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

### 4) Run

From the project root:
```bash
./bin/nextflow run main.nf -c nextflow.config
```

- If running with Docker explicitly:
```bash
./bin/nextflow run main.nf -c nextflow.config -with-docker
```

Outputs will be created under `params.outdir` with subfolders for each step (e.g., `split/`, `align/`, `count_reads/`, `find_ase/`) and reports in `reports/` if enabled.