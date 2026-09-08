# RNAFlowX Project Plan and Implementation Status

## 1. Project Overview

RNAFlowX is a reproducible and containerized bulk RNA-seq analysis and workflow-engineering platform built with Nextflow DSL2.

The project was designed to demonstrate an end-to-end production-oriented bioinformatics workflow under resource-constrained local execution while remaining portable to additional execution environments.

## 2. Core Analysis Workflow

The implemented workflow supports:

1. Validated paired-end FASTQ input through a samplesheet.
2. Raw-read quality control with FastQC.
3. Adapter trimming and quality filtering with fastp.
4. Post-trimming quality control.
5. Transcript-level quantification with Salmon.
6. Transcript-to-gene aggregation with tximport.
7. Differential expression analysis with DESeq2.
8. Functional enrichment analysis.
9. GO, KEGG, and GSEA downstream analysis.
10. Analytical visualization and scientific reporting.

## 3. Workflow Engineering

RNAFlowX includes:

- Nextflow DSL2 workflow orchestration
- Modular workflow organization
- Explicit configuration hierarchy
- Docker containerization
- Local execution profile
- SLURM execution profile
- Azure execution configuration
- Git/GitHub version control
- Automated CI/testing
- pytest project validation
- Nextflow execution reports
- Timeline, trace, and DAG generation
- Benchmarking utilities
- Software and parameter provenance
- Reproducibility documentation

## 4. Data and Reporting Layer

The completed project additionally includes:

- SQLite-based structured result storage
- Interactive Streamlit dashboard
- Differential-expression visualization
- Functional-enrichment visualization
- GSEA visualization
- Downloadable analytical outputs
- R Markdown reporting
- Quarto publication-style reporting
- CITATION.cff project citation metadata

## 5. Resource-Constrained Design

The primary development environment was a local Linux/WSL system with approximately:

- 8 GB total RAM
- approximately 6 GB normally available RAM
- approximately 1 TB storage
- no dependency on paid cloud compute

To support this environment, the workflow uses:

- Salmon rather than genome-scale STAR alignment
- Conservative process resource limits
- Limited process concurrency
- A small four-sample paired-end dataset
- Avoidance of unnecessary large intermediate files

The primary dataset contains four paired-end RNA-seq samples with approximately 308.36 MB of compressed sequencing data.

## 6. Frozen Technical Decisions

| Component | Implementation |
|---|---|
| Organism | Homo sapiens |
| Genome assembly | GRCh38.p14 |
| Annotation | GENCODE v50 |
| Transcript reference | GENCODE v50 transcript FASTA |
| Gene annotation | GENCODE v50 basic annotation GFF3 |
| Quantification | Salmon |
| Gene aggregation | tximport |
| Differential expression | DESeq2 |
| Functional analysis | GO / KEGG / GSEA |
| Workflow engine | Nextflow DSL2 |
| Containerization | Docker |
| Primary execution | Local Linux/WSL |
| Additional execution configs | SLURM / Azure |
| Samples | 4 |
| Experimental design | 2 control + 2 treatment |
| Input size | ~308.36 MB compressed |
| Reporting | R Markdown / Quarto / HTML |
| Interactive presentation | Streamlit |
| Structured data storage | SQLite |
| Testing | pytest + CI |
| Benchmarking | Implemented |
| Observability | Nextflow report / timeline / trace / DAG |
| Provenance | Implemented |

## 7. Design Principles

RNAFlowX prioritizes:

- Reproducibility
- Modularity
- Traceability
- Resource efficiency
- Portability
- Version control
- Containerized execution
- Explicit configuration
- Automated testing
- Scientific reporting
- Execution observability
- Measurable performance

## 8. Completion Status

The planned RNAFlowX workflow has been implemented.

Core analytical objectives were completed, including QC, preprocessing, quantification, gene-level aggregation, differential expression, downstream functional analysis, visualization, and reporting.


RNAFlowX therefore represents both a bulk RNA-seq analysis workflow and a portfolio demonstration of reproducible bioinformatics workflow engineering.
