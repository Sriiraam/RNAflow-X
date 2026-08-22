# RNAFlowX Methodology

## 1. Workflow Overview

RNAFlowX implements a modular bulk RNA-seq workflow using Nextflow DSL2.

```text
Paired-end FASTQ
       │
       ▼
   FastQC (Raw)
       │
       ▼
     FastP
       │
       ▼
 FastQC (Trimmed)
       │
       ▼
     Salmon
       │
       ▼
    tximport
       │
       ▼
     DESeq2
       │
       ├── PCA
       ├── MA Plot
       ├── Volcano Plot
       └── Differential Expression Tables
       │
       ├────────► MultiQC
       │
       ▼
 GO / KEGG / GSEA
       │
       ▼
Scientific Reports
       │
       ▼
Streamlit Dashboard
```

## 2. Quality Control

FastQC evaluates sequencing-read quality before and after FastP preprocessing.

MultiQC consolidates quality-control information into a unified report.

## 3. Read Preprocessing

FastP performs adapter removal and quality filtering on paired-end sequencing reads.

## 4. Transcript Quantification

Salmon performs alignment-free transcript-level abundance estimation using the configured transcriptome reference index.

## 5. Gene-Level Aggregation

tximport imports Salmon transcript estimates and aggregates them into gene-level counts using the transcript-to-gene mapping.

## 6. Differential Expression

DESeq2 performs differential-expression analysis using the gene-level count matrix and experimental metadata.

Major outputs include:

- differential-expression results
- significant-gene results
- normalized counts
- PCA
- MA plot
- volcano plot

## 7. Functional Analysis

The downstream R analysis layer performs:

- GO Biological Process enrichment
- GO Molecular Function enrichment
- GO Cellular Component enrichment
- KEGG pathway enrichment
- Gene Set Enrichment Analysis (GSEA)

GO, KEGG and GSEA are currently downstream analytical steps and are not native Nextflow processes.

## 8. Reporting and Visualization

RNAFlowX provides:

- MultiQC quality reporting
- Quarto scientific reporting
- R Markdown reporting
- SQLite analytical storage
- Interactive Streamlit dashboard

## 9. Execution Architecture

RNAFlowX supports:

- Local execution
- Docker execution
- SLURM/HPC-ready configuration
- Azure Batch-ready configuration

Container portability has additionally been demonstrated using a local Kubernetes/kind environment.

## 10. Reproducibility

RNAFlowX supports reproducibility through:

- Nextflow DSL2 workflow orchestration
- Docker containerization
- version-controlled configuration
- GitHub Actions CI/CD
- automated pytest validation
- execution reports and timelines
- task-level tracing
- workflow DAG generation
- environment provenance capture
- reproducible benchmarking
