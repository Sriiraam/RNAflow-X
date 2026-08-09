# RNAflow-X Project Plan

## 1. Project Overview

RNAflow-X is a reproducible, containerized Bulk RNA-seq analysis workflow designed to demonstrate production-oriented bioinformatics workflow engineering using a resource-constrained execution environment.

The project focuses on building a complete end-to-end RNA-seq workflow rather than demonstrating only individual bioinformatics tools.

The workflow will integrate:

- Nextflow DSL2
- Salmon
- FastQC
- FastP
- R
- DESeq2
- Python
- Docker
- Git/GitHub
- Structured execution logging
- Benchmarking
- Metadata and provenance tracking
- Automated reporting
- Reproducible configuration

## 2. Project Objective

The primary objective is to develop a professional Bulk RNA-seq workflow capable of:

1. Accepting paired-end FASTQ input through a validated samplesheet.
2. Performing raw-read quality control.
3. Performing adapter and quality trimming.
4. Performing transcript-level quantification using Salmon.
5. Generating gene-level abundance estimates.
6. Performing differential expression analysis using DESeq2.
7. Producing publication-style analytical visualizations.
8. Generating execution and performance metrics.
9. Capturing software, parameter, input and reference provenance.
10. Producing a reproducible final analysis report.

## 3. Resource Constraints

The workflow is intentionally designed to run on a local machine with approximately:

- 8 GB total RAM
- approximately 6 GB available RAM during normal operation
- approximately 1 TB storage
- no dependency on paid cloud infrastructure

The project therefore deliberately avoids:

- Full GRCh38 genome indexing
- STAR genome indexing
- Large cloud compute instances
- Large-scale datasets
- Unnecessary intermediate files

The target dataset is limited to four paired-end RNA-seq samples with approximately 308.36 MB of compressed sequencing data.

## 4. Frozen Technical Decisions

| Component | Decision |
|---|---|
| Organism | Homo sapiens |
| Genome assembly | GRCh38.p14 |
| Annotation release | GENCODE v50 |
| Transcript reference | GENCODE v50 transcript FASTA |
| Gene annotation | GENCODE v50 basic annotation GFF3 |
| Quantification | Salmon |
| Differential expression | DESeq2 |
| Workflow engine | Nextflow DSL2 |
| Containerization | Docker |
| Primary execution environment | Local Linux/WSL |
| Cloud execution | Not required |
| Samples | 4 |
| Experimental design | 2 control + 2 treatment |
| Input size | ~308.36 MB compressed |
| Reporting | HTML + graphical outputs |
| Benchmarking | Required |
| Observability | Required |
| Provenance | Required |

## 5. Design Philosophy

The project prioritizes:

- Reproducibility
- Modularity
- Traceability
- Resource efficiency
- Version control
- Containerized execution
- Explicit configuration
- Automated reporting
- Measurable performance

The workflow should remain understandable to another bioinformatician without requiring knowledge of the original development environment.

## 6. Success Criteria

RNAflow-X will be considered complete when:

- The pipeline executes successfully from a clean environment.
- All samples pass through the automated workflow.
- QC outputs are generated.
- Salmon quantification completes successfully.
- Differential expression analysis completes successfully.
- Benchmark metrics are captured.
- Nextflow execution reports are generated.
- Software and parameter provenance is recorded.
- A final analytical report is generated.
- The repository can be reproduced from the documented instructions.
- The workflow is version-controlled through Git.