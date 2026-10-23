# Container for the Nextflow workflow, built from the same pixi.lock as local development,
# so the container and `pixi run hcv-cluster` use identical tool versions.
#
#   docker build --build-arg GIT_COMMIT=$(git rev-parse HEAD) -t hcv-cluster:<version> .
#   nextflow run main.nf -profile docker --input sequences.fasta   (uses hcv-cluster:<manifest version>)
FROM ghcr.io/prefix-dev/pixi:0.67.0

# procps provides `ps`, which Nextflow uses to collect task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hcv-cluster

# Install the locked environment first so code changes don't invalidate this layer.
COPY pixi.toml pixi.lock ./
RUN pixi install --locked --environment pipeline && pixi clean cache --yes

# The package carries its own reference data (hcv_cluster/data).
COPY hcv_cluster ./hcv_cluster

# Recorded in run_info.json (the image has no .git directory).
ARG GIT_COMMIT=unknown

ENV PATH=/opt/hcv-cluster/.pixi/envs/pipeline/bin:$PATH \
    PYTHONPATH=/opt/hcv-cluster \
    HCV_CLUSTER_COMMIT=$GIT_COMMIT
