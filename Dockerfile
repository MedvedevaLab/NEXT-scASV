FROM continuumio/miniconda3

# Copy project files
COPY . /app/

# Install curl/git for Nextflow + WASP
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl git \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
RUN conda env create -f /app/as_env.yml

# Ensure env works and is used by default (even in bash -lc)
RUN /opt/conda/envs/as_env/bin/python -c "import numpy, pandas, scipy, pysam; print('as_env imports: OK')"
ENV CONDA_AUTO_ACTIVATE_BASE=false
ENV PATH="/opt/conda/envs/as_env/bin:/app/bin:/opt/conda/bin:${PATH}"
RUN printf '%s\n' \
      'source /opt/conda/etc/profile.d/conda.sh' \
      'conda activate as_env' \
    > /etc/profile.d/activate-as_env.sh \
    && ln -sf /opt/conda/envs/as_env/bin/python /usr/local/bin/python \
    && ln -sf /opt/conda/envs/as_env/bin/pip /usr/local/bin/pip

# Download Nextflow into /app/bin (optionally pin version at build time)
ARG NEXTFLOW_VERSION=""
RUN mkdir -p /app/bin \
    && if [ -n "${NEXTFLOW_VERSION}" ]; then export NXF_VER="${NEXTFLOW_VERSION}"; fi \
    && LD_LIBRARY_PATH="" /usr/bin/curl -fsSL https://get.nextflow.io | /bin/bash \
    && mv nextflow /app/bin/nextflow \
    && chmod +x /app/bin/nextflow \
    && ln -sf /app/bin/nextflow /usr/local/bin/nextflow

# Download WASP into /app/bin/WASP (pin default tag; override at build time)
ARG WASP_VERSION="v0.3.4"
RUN if [ ! -d /app/bin/WASP ]; then \
      if [ -n "${WASP_VERSION}" ]; then \
        git clone --depth 1 --branch "${WASP_VERSION}" https://github.com/bmvdgeijn/WASP.git /app/bin/WASP; \
      else \
        git clone --depth 1 https://github.com/bmvdgeijn/WASP.git /app/bin/WASP; \
      fi; \
    fi

# Set working directory
WORKDIR /app

# Apply WASP compatibility patch if needed
RUN if [ -f /app/bin/WASP/mapping/snptable.py ] && grep -q 'np.int' /app/bin/WASP/mapping/snptable.py; then \
         sed -i 's/np.int/int/g' /app/bin/WASP/mapping/snptable.py; \
       fi

# Default command
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "as_env"]
CMD ["nextflow", "run", "main.nf", "-c", "nextflow.config", "--help"]