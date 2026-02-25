FROM continuumio/miniconda3

# Copy project files
COPY . /app/

# Install curl/git for Nextflow + WASP
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl git \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
RUN conda env create -f /app/as_env.yml

# Ensure tools from the env and repo are on PATH
ENV PATH="/opt/conda/envs/as_env/bin:/app/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/envs/as_env/lib:${LD_LIBRARY_PATH}"

# Set working directory
WORKDIR /app

# Make Nextflow executable and apply WASP compatibility patch if needed
RUN chmod +x /app/bin/nextflow \
    && if [ -f /app/bin/WASP/mapping/snptable.py ] && grep -q 'np.int' /app/bin/WASP/mapping/snptable.py; then \
         sed -i 's/np.int/int/g' /app/bin/WASP/mapping/snptable.py; \
       fi

# Default command
CMD ["nextflow", "run", "main.nf", "-c", "nextflow.config", "--help"]