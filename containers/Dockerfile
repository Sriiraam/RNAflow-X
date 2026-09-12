FROM condaforge/miniforge3@sha256:23c769cbfc96d261aa6bd9dbbcb8c0a6ca6a66887410141076fa79c334e82c81

LABEL org.opencontainers.image.title="RNAFlowX"
LABEL org.opencontainers.image.description="Reproducible bulk RNA-seq analysis environment"
LABEL org.opencontainers.image.version="1.1.0"
LABEL org.opencontainers.image.authors="Sriram B"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash=5.2.21-2ubuntu4 \
        ca-certificates=20260601~24.04.1 \
        procps=2:4.0.4-4ubuntu3.3 \
        && rm -rf /var/lib/apt/lists/*

# Conda is used only inside the Docker image to install the
# reproducible bioinformatics runtime. Local development uses Python venv.
COPY containers/locks/conda-linux-64.lock /tmp/conda-linux-64.lock

RUN conda install -y --file /tmp/conda-linux-64.lock && \
    conda clean -afy && \
    rm /tmp/conda-linux-64.lock

RUN fastqc --version && \
    fastp --version && \
    salmon --version && \
    multiqc --version && \
    Rscript -e 'library(tximport); library(DESeq2); library(optparse); library(clusterProfiler); library(enrichplot); library(org.Hs.eg.db); library(ggplot2); cat("R dependencies OK\n")'

WORKDIR /workspace

CMD ["/bin/bash"]
