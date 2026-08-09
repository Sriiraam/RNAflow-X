# RNAflow-X Methodology

## 1. Workflow Overview

RNAflow-X follows the following analytical workflow:

```text
Raw FASTQ
   |
   v
FastQC
   |
   v
FastP
   |
   v
Salmon Quantification
   |
   v
Transcript Abundance
   |
   v
tximport
   |
   v
Gene-level Count Matrix
   |
   v
DESeq2
   |
   +----> PCA
   |
   +----> MA Plot
   |
   +----> Volcano Plot
   |
   +----> Heatmap
   |
   +----> Differential Expression Tables
   |
   v
Final Report