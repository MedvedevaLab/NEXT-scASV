FROM continuumio/miniconda3

# Copy project files
COPY . /app/

# Create conda environment
RUN conda env create -f /app/as_env.yml

# Set working directory
WORKDIR /app

# Make scripts executable
RUN chmod +x /app/run_workflow.sh && chmod +x /app/run.sh && chmod +x /app/bin/nextflow

# Activate the conda environment
SHELL ["conda", "run", "-n", "as_env", "/bin/bash", "-c"]

# Default command
CMD ["./run.sh", "--help"] 