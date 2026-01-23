FROM continuumio/miniconda3

# Copy project files
COPY . /app/

# Create conda environment
RUN conda env create -f /app/as_env.yml

# Set working directory
WORKDIR /app

# Make Nextflow executable
RUN chmod +x /app/bin/nextflow

# Activate the conda environment
SHELL ["conda", "run", "-n", "as_env", "/bin/bash", "-c"]

# Default command
CMD ["./bin/nextflow", "run", "main.nf", "-c", "nextflow.config", "--help"]