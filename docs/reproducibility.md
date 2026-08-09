# Reproducibility Specification

## 1. Objective

RNAflow-X is designed so that another user can reproduce the analysis using the documented inputs, reference files, software versions and pipeline configuration.

## 2. Version Control

The project is maintained using Git and GitHub.

The Git commit associated with each analysis run will be recorded.

## 3. Software Reproducibility

Software versions will be explicitly recorded.

Docker images will use pinned versions rather than the `latest` tag.

## 4. Reference Reproducibility

The following reference information will be recorded:

- Genome assembly
- GENCODE release
- Reference filename
- Reference checksum
- Annotation filename
- Annotation checksum

## 5. Input Provenance

Every FASTQ input will be associated with:

- Sample ID
- SRA run accession
- GEO accession
- Biological condition
- Replicate
- File name
- File checksum

## 6. Parameter Provenance

The final execution metadata will record the relevant pipeline parameters.

## 7. Execution Metadata

RNAflow-X will generate structured metadata containing:

```text
Pipeline name
Pipeline version
Git commit
Execution date
Nextflow version
Container version
Input files
Reference files
Reference versions
Pipeline parameters
Execution environment