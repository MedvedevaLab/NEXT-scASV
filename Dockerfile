FROM continuumio/miniconda3

# Copy project files
COPY . /app/

# Install curl/git for Nextflow + WASP
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl git \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
RUN conda env create -f /app/as_env.yml \
    && conda clean -afy

# Set working directory
WORKDIR /app

# Install WASP under /app/bin/WASP and apply numpy fix
RUN git clone https://github.com/bmvdgeijn/WASP /app/bin/WASP \
    && sed -i 's/np.int/int/g' /app/bin/WASP/mapping/snptable.py

# Install Nextflow into /app/bin
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /app/bin/nextflow \
    && chmod +x /app/bin/nextflow

# Make conda env available for runtime commands
ENV JAVA_HOME="/opt/conda/envs/as_env"
ENV PATH="/opt/conda/envs/as_env/bin:/app/bin:${PATH}"
ENV NXF_HOME="/app/.nextflow"

RUN mkdir -p /app/.nextflow

# Default command
CMD ["nextflow", "run", "main.nf", "-c", "nextflow.config", "--help"]