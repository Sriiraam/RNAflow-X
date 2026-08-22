FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.title="RNAFlowX"
LABEL org.opencontainers.image.description="Reproducible bulk RNA-seq analysis environment"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.authors="Sriram B"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        procps \
        && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
    -c conda-forge \
    -c bioconda \
    fastqc \
    fastp \
    salmon \
    multiqc \
    r-base \
    r-optparse \
    bioconductor-tximport \
    bioconductor-deseq2 \
    bioconductor-clusterprofiler \
    bioconductor-enrichplot \
    bioconductor-org.hs.eg.db \
    r-ggplot2 \
    && conda clean -afy

RUN fastqc --version && \
    fastp --version && \
    salmon --version && \
    multiqc --version && \
    Rscript -e 'library(tximport); library(DESeq2); library(optparse); library(clusterProfiler); library(enrichplot); library(org.Hs.eg.db); library(ggplot2); cat("R dependencies OK\n")'

WORKDIR /workspace

CMD ["/bin/bash"]
