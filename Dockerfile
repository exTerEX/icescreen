FROM condaforge/mambaforge:latest

ARG ICESCREEN_VERSION=1.3.3

LABEL maintainer="Andreas Sagen"
LABEL description="ICEscreen environment for the advena Nextflow pipeline"
LABEL org.opencontainers.image.source="https://github.com/exterex/icescreen"
LABEL org.opencontainers.image.licenses="AGPL-3.0"

COPY environment.yaml /tmp/environment.yaml

RUN mamba env create -f /tmp/environment.yaml \
    && mamba clean -afy && rm /tmp/environment.yaml

ENV ICESCREEN_VERSION=${ICESCREEN_VERSION}
# TODO: Dynamically add correct conda environment path to PATH variable
ENV PATH=/opt/conda/envs/icescreen_env_1-3-3/bin:$PATH

WORKDIR /data
