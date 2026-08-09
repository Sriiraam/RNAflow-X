# RNAflow-X — Frozen Technical Decisions

**Project:** RNAflow-X  
**Pipeline type:** Bulk RNA-seq  
**Organism:** Homo sapiens  
**Status:** FROZEN  
**Purpose:** Reproducible, production-style bulk RNA-seq workflow for portfolio and technical demonstration.

---

## 1. Project Scope

RNAflow-X is a reproducible bulk RNA-seq analysis pipeline designed to demonstrate industry-relevant bioinformatics workflow engineering practices.

The project focuses on:

- Nextflow DSL2 workflow development
- Containerized execution with Docker
- Reproducible configuration
- RNA-seq quality control
- Read preprocessing
- Transcript quantification
- Gene-level quantification
- Differential expression analysis
- MultiQC reporting
- Reproducibility and provenance
- Execution observability
- Benchmarking
- Git/GitHub-based project governance

The pipeline is intentionally designed to operate on a small-to-medium benchmark dataset so that the complete workflow can be executed locally without paid cloud infrastructure.

---

# 2. Dataset Decision

## GEO Study

**GEO accession:** `GSE342612`

**Study:** Transcriptomic and H3K27ac chromatin responses of human microglia to acute PFOS exposure and recovery.

Only the **RNA-seq** component is used by RNAflow-X.

The H3K27ac CUT&Tag data are outside the scope of this project.

## Biological system

- Organism: Homo sapiens
- Cell line: HMC3
- Cell type: microglial cell
- Platform: Illumina NextSeq 550
- Library layout: Paired-end
- Read length: approximately 72 bp
- Assay: RNA-Seq

---

# 3. Frozen Experimental Design

The project uses four biological samples:

| Sample | GEO accession | SRA run | Condition | Role |
|---|---|---|---|---|
| GSM9934646 | GSM9934646 | SRR40038350 | Vehicle control, 24 h | Control |
| GSM9934647 | GSM9934647 | SRR40038349 | Vehicle control, 24 h | Control |
| GSM9934653 | GSM9934653 | SRR40038343 | 50 µM PFOS, 24 h | Treatment |
| GSM9934652 | GSM9934652 | SRR40038344 | 50 µM PFOS, 24 h | Treatment |

### Experimental comparison

```text
Control
  ├── Vehicle replicate 1
  └── Vehicle replicate 2

Treatment
  ├── 50 µM PFOS replicate 1
  └── 50 µM PFOS replicate 2
