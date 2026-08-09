# Dataset Specification

## 1. Dataset Source

RNAflow-X uses publicly available RNA-seq data from the NCBI Gene Expression Omnibus (GEO) and Sequence Read Archive (SRA).

### GEO Series

**GSE342612**

### BioProject

**PRJNA1508658**

### Organism

Homo sapiens

### Cell model

HMC3 human microglial cells

### Instrument

Illumina NextSeq 550

### Library layout

Paired-end

### Assay

RNA-seq

## 2. Biological Question

The dataset investigates transcriptional responses of human HMC3 microglial cells following acute PFOS exposure.

For RNAflow-X, a simplified two-group comparison is used:

- Vehicle control
- 50 µM PFOS treatment

This provides a clear binary differential-expression design suitable for validating the workflow.

## 3. Selected Samples

Four biological samples are used.

### Control

| Sample | Run | GEO Sample | Condition | Approx. compressed size |
|---|---|---|---|---:|
| control_rep1 | SRR40038349 | GSM9934647 | Vehicle control, 24 h | 136.01 MB |
| control_rep2 | SRR40038350 | GSM9934646 | Vehicle control, 24 h | 49.12 MB |

### Treatment

| Sample | Run | GEO Sample | Condition | Approx. compressed size |
|---|---|---|---|---:|
| pfos50_rep1 | SRR40038343 | GSM9934653 | 50 µM PFOS, 24 h | 61.96 MB |
| pfos50_rep2 | SRR40038344 | GSM9934652 | 50 µM PFOS, 24 h | 61.27 MB |

## 4. Dataset Size

Total compressed sequencing data:

**~308.36 MB**

This size was deliberately selected to allow complete end-to-end execution on local hardware without paid cloud infrastructure.

## 5. Experimental Design

The simplified design is:

```text
Control
├── control_rep1
└── control_rep2

Treatment
├── pfos50_rep1
└── pfos50_rep2