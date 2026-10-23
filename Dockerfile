# Container for the Nextflow workflow, built from the same pixi.lock as local development,
# so the container and `pixi run hcv-cluster` use identical tool versions.
#
#   docker build -t hcv-cluster .
#   nextflow run main.nf -profile docker --input sequences.fasta
FROM ghcr.io/prefix-dev/pixi:0.67.0

# procps provides `ps`, which Nextflow uses to collect task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hcv-cluster

# Install the locked environment first so code changes don't invalidate this layer.
COPY pixi.toml pixi.lock ./
RUN pixi install --locked --environment pipeline && pixi clean cache --yes

COPY *.py ./
COPY reference_data ./reference_data

ENV PATH=/opt/hcv-cluster/.pixi/envs/pipeline/bin:$PATH
