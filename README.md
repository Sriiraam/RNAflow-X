# RNAFlowX

## Reproducible Bulk RNA-seq Analysis & Engineering Platform

RNAFlowX is a modular bulk RNA-sequencing analysis platform built with **Nextflow DSL2**. It combines reproducible workflow orchestration, sequencing quality control, transcript quantification, differential-expression analysis, functional enrichment, scientific reporting, and an interactive Streamlit dashboard.

The project demonstrates both **bioinformatics analysis** and **production-oriented workflow engineering**, rather than functioning as a single-use RNA-seq script.

---

## Overview

RNAFlowX processes paired-end RNA-seq data through quality control, preprocessing, transcript quantification, gene-level aggregation, differential-expression analysis, reporting, downstream pathway analysis, and interactive visualization.

The current implementation contains two connected analytical layers.

### Automated Nextflow Workflow

```text
Paired-end FASTQ
       │
       ▼
    FastQC
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
       ├──────────────► MultiQC
       │
       ▼
Differential Expression Results
```

### Downstream Analysis & Presentation

```text
DESeq2 Results
      │
      ├──► GO Enrichment
      ├──► KEGG Enrichment
      └──► GSEA
              │
              ▼
       Scientific Reports
              │
              ▼
       Streamlit Dashboard
```

> **Current architecture note:** GO, KEGG, and GSEA are implemented through the downstream R analysis layer. They are not yet orchestrated as native Nextflow processes.

---

## Key Features

- Modular **Nextflow DSL2** workflow architecture
- Raw and post-trimming **FastQC**
- Adapter trimming and quality filtering with **FastP**
- Alignment-free transcript quantification with **Salmon**
- Transcript-to-gene aggregation using **tximport**
- Differential-expression analysis with **DESeq2**
- PCA, MA, and volcano visualizations
- GO Biological Process, Molecular Function, and Cellular Component enrichment
- KEGG pathway enrichment
- Gene Set Enrichment Analysis (**GSEA**)
- Unified **MultiQC** quality reporting
- Interactive multi-page **Streamlit dashboard**
- Quarto and R Markdown scientific reporting
- Local execution profile
- Docker configuration layer
- Version-controlled workflow configuration
- Reference-data checksum tracking
- Structured engineering documentation

---

## Experimental Design

RNAFlowX currently uses a deliberately constrained public human bulk RNA-seq dataset so the complete workflow can be executed on modest local hardware.

| Property | Value |
|---|---|
| GEO Series | **GSE342612** |
| BioProject | **PRJNA1508658** |
| Organism | *Homo sapiens* |
| Cell Model | HMC3 human microglial cells |
| Sequencing Platform | Illumina NextSeq 550 |
| Library Layout | Paired-end |
| Assay | RNA-seq |
| Comparison | Vehicle control vs 50 µM PFOS |
| Exposure | 24 hours |
| Biological Samples | 4 |

### Selected Samples

| Group | Sample | SRA Run |
|---|---|---|
| Control | `control_rep1` | SRR40038349 |
| Control | `control_rep2` | SRR40038350 |
| PFOS 50 µM | `pfos50_rep1` | SRR40038343 |
| PFOS 50 µM | `pfos50_rep2` | SRR40038344 |

Total compressed sequencing input is approximately **308 MB**.

The dataset size was intentionally constrained to support complete local execution without depending on paid cloud infrastructure.

Detailed dataset documentation is available in [`docs/dataset.md`](docs/dataset.md).

---

## Technology Stack

| Layer | Technologies |
|---|---|
| Workflow Orchestration | Nextflow DSL2 |
| Quality Control | FastQC, MultiQC |
| Preprocessing | FastP |
| Quantification | Salmon |
| Count Aggregation | tximport |
| Differential Expression | DESeq2 |
| Functional Analysis | clusterProfiler, GO, KEGG, GSEA |
| Statistical Programming | R |
| Dashboard | Streamlit, Pandas, Plotly |
| Scientific Reporting | Quarto, R Markdown |
| Version Control | Git / GitHub |
| Container Configuration | Docker |

---

## Repository Architecture

```text
RNAFlowX/
│
├── main.nf
├── nextflow.config
│
├── workflows/
│   ├── rnaseq.nf
│   ├── qc.nf
│   ├── quantification.nf
│   └── counting.nf
│
├── modules/
│   ├── qc/
│   │   └── fastqc.nf
│   ├── preprocessing/
│   │   └── fastp.nf
│   ├── quantification/
│   │   └── salmon.nf
│   ├── counting/
│   │   └── tximport.nf
│   ├── differential_expression/
│   │   └── deseq2.nf
│   └── reporting/
│       └── multiqc.nf
│
├── bin/
│   ├── run_tximport.R
│   ├── run_deseq2.R
│   └── run_enrichment.R
│
├── conf/
│   ├── params.config
│   ├── base.config
│   ├── local.config
│   └── docker.config
│
├── assets/
│   ├── samplesheet.csv
│   ├── logo.svg
│   └── hero_dna.svg
│
├── data/
│   ├── metadata.csv
│   └── reference/
│
├── streamlit_app/
│   ├── app.py
│   ├── pages/
│   │   ├── 1_Quality_Control.py
│   │   ├── 2_Differential_Expression.py
│   │   ├── 3_Functional_Enrichment.py
│   │   ├── 4_GSEA.py
│   │   └── 5_Downloads.py
│   └── assets/
│
├── reports/
│   ├── RNAFlowX_Final_Report.Rmd
│   ├── RNAFlowX_Publication.qmd
│   └── RNAFlowX_Report.qmd
│
├── docs/
│   ├── dataset.md
│   ├── methodology.md
│   ├── reference.md
│   ├── reproducibility.md
│   ├── benchmarking.md
│   ├── frozen_decisions.md
│   └── project_plan.md
│
├── tests/
├── LICENSE
└── CITATION.cff
```

The separation between orchestration, analytical processes, statistical scripts, configuration, reporting, and visualization is intentional. This keeps the workflow modular and allows individual components to evolve independently.

---

## Workflow Architecture

### 1. Input

RNAFlowX consumes a sample manifest describing paired-end FASTQ files.

Default manifest:

```text
assets/samplesheet.csv
```

Experimental metadata are maintained separately:

```text
data/metadata.csv
```

This separates sequencing-file discovery from experimental-design information.

### 2. Quality Control

Raw sequencing reads are evaluated using **FastQC**.

After preprocessing, FastQC is executed again on trimmed reads to provide before/after quality assessment.

### 3. Preprocessing

**FastP** performs adapter removal and read-quality filtering.

FastP HTML and JSON reports are collected for downstream reporting.

### 4. Transcript Quantification

Filtered reads are quantified against a pre-built transcriptome index using **Salmon**.

This provides alignment-free transcript abundance estimation.

### 5. Gene-Level Aggregation

**tximport** imports Salmon transcript abundances and aggregates them to gene-level counts using the configured transcript-to-gene mapping.

### 6. Differential Expression

**DESeq2** performs differential-expression analysis using the sample metadata and generated count matrix.

Outputs include:

- complete differential-expression results
- significant-gene results
- normalized counts
- PCA
- MA plot
- volcano plot
- serialized DESeq2 object
- analysis summary

### 7. Unified QC Reporting

FastQC, FastP, and Salmon reporting artifacts are collected by **MultiQC** into a consolidated HTML quality report.

### 8. Functional Analysis

DESeq2 results are subsequently processed by the R enrichment layer to generate:

- GO Biological Process
- GO Molecular Function
- GO Cellular Component
- KEGG pathways
- GSEA

---

## Configuration

RNAFlowX separates configuration from workflow implementation.

Primary Nextflow configuration:

```text
nextflow.config
```

Parameter configuration:

```text
conf/params.config
```

Current default parameters include:

```text
samplesheet  = assets/samplesheet.csv
salmon_index = data/reference/salmon_index
tx2gene      = data/reference/tx2gene.tsv
metadata     = data/metadata.csv
outdir       = results
```

Execution-specific configuration is separated into:

```text
conf/local.config
conf/docker.config
```

---

## Running RNAFlowX

### Requirements

Core runtime requirements include:

- Linux or WSL2
- Java
- Nextflow >= 24.10.0
- FastQC
- FastP
- Salmon
- MultiQC
- R
- DESeq2
- tximport

Additional R packages are required for GO, KEGG, and GSEA analysis.

### Local Execution

From the repository root:

```bash
nextflow run main.nf -profile local
```

Nextflow loads the project parameters and local execution configuration automatically.

### Resume an Interrupted Run

Nextflow caching allows previously completed processes to be reused:

```bash
nextflow run main.nf -profile local -resume
```

### Docker Profile

A Docker execution profile is defined through:

```text
conf/docker.config
```

The Docker/container layer remains part of the project's ongoing reproducibility validation and should not yet be interpreted as a fully validated production container deployment.

---

## Pipeline Outputs

Pipeline outputs are organized under:

```text
results/
```

Major output categories include:

```text
results/
├── counting/
├── differential_expression/
├── enrichment/
├── multiqc/
└── quantification/
```

### Quantification

The quantification layer produces Salmon abundance estimates for each biological sample.

### Count Matrix

tximport generates the gene-level count matrix used by DESeq2.

### Differential Expression

The DESeq2 output directory contains:

```text
differential_expression.csv
significant_genes.csv
normalized_counts.csv
PCA.png
MA_plot.png
volcano_plot.png
dds.rds
deseq2_summary.txt
```

### Functional Enrichment

Functional analysis produces separate results for:

```text
GO/
├── BP/
├── CC/
└── MF/

KEGG/

GSEA/
```

GO and KEGG results are additionally separated into:

- all significant genes
- upregulated genes
- downregulated genes

---

## Interactive Streamlit Dashboard

RNAFlowX includes an interactive multi-page Streamlit interface built directly on pipeline-generated outputs.

### Dashboard Pages

1. **Project Overview** — workflow, experimental design, architecture, and engineering principles.
2. **Quality Control** — sample-level sequencing quality metrics and FastQC/MultiQC status.
3. **Differential Expression** — PCA, MA plot, volcano plot, statistics, and DEG tables.
4. **Functional Enrichment** — GO and KEGG enrichment results.
5. **GSEA** — ranked gene-set enrichment results and enrichment statistics.
6. **Downloads & Reports** — access to major analytical outputs and reports.

### Run the Dashboard

From the repository root:

```bash
streamlit run streamlit_app/app.py
```

The dashboard is a presentation layer over pipeline-generated outputs rather than a separate analytical workflow.

---

## Scientific Reporting

RNAFlowX includes scientific reporting implementations under:

```text
reports/
```

The project currently contains both **Quarto** and **R Markdown** reporting approaches.

The reporting layer connects:

```text
Methods
   │
   ▼
Quality Control
   │
   ▼
Differential Expression
   │
   ▼
Functional Enrichment
   │
   ▼
Biological Interpretation
```

MultiQC separately provides the unified sequencing-quality report.

---

## Reproducibility Strategy

RNAFlowX treats reproducibility as an engineering requirement rather than only a documentation concern.

Current reproducibility controls include:

- fixed dataset selection
- documented sample provenance
- frozen technical decisions
- explicit sample manifest
- separate experimental metadata
- version-controlled Nextflow modules
- version-controlled R scripts
- configuration profiles
- fixed reference resources
- reference checksum tracking
- deterministic output organization
- software-version reporting
- Nextflow execution caching

Detailed documentation:

- [Frozen Technical Decisions](docs/frozen_decisions.md)
- [Dataset Specification](docs/dataset.md)
- [Methodology](docs/methodology.md)
- [Reference Resources](docs/reference.md)
- [Reproducibility](docs/reproducibility.md)
- [Benchmarking](docs/benchmarking.md)

---

## Engineering Principles

### Modularity

Each major analytical stage is isolated into a reusable Nextflow DSL2 module.

### Separation of Concerns

Workflow orchestration, statistical analysis, configuration, scientific reporting, and dashboard presentation remain separate layers.

### Reproducibility

Dataset selection, reference resources, metadata, parameters, software execution, and output organization are explicitly controlled.

### Local-First Development

The benchmark dataset is deliberately constrained so development and validation can occur on modest hardware without requiring paid cloud infrastructure.

### Extensibility

Execution profiles and modular processes provide a foundation for future HPC, container, and infrastructure integration.

---

## Implementation Status

### Implemented

- [x] Nextflow DSL2 architecture
- [x] Modular workflow organization
- [x] Raw FastQC
- [x] FastP preprocessing
- [x] Post-trimming FastQC
- [x] Salmon transcript quantification
- [x] tximport gene-level aggregation
- [x] DESeq2 differential expression
- [x] PCA visualization
- [x] MA visualization
- [x] Volcano visualization
- [x] MultiQC reporting
- [x] GO enrichment
- [x] KEGG enrichment
- [x] GSEA
- [x] Quarto / R Markdown reporting
- [x] Streamlit analytical dashboard
- [x] Downloadable analytical outputs
- [x] Local execution profile
- [x] Dataset and reference documentation
- [x] Git version control

### Engineering Roadmap

The following capabilities are planned and are **not presented as completed functionality**:

- [ ] Integrate GO / KEGG / GSEA directly into the Nextflow DAG
- [ ] SQLite analytical data layer
- [ ] Expanded automated testing
- [ ] GitHub Actions CI/CD
- [ ] Fully validated Docker execution
- [ ] SLURM/HPC execution profile
- [ ] Cloud-ready infrastructure configuration
- [ ] Structured logging and observability
- [ ] Expanded provenance tracking
- [ ] Kubernetes deployment demonstration
- [ ] Formal benchmark validation
- [ ] Production release packaging

---

## Documentation

Detailed project documentation is maintained separately from the README.

| Document | Purpose |
|---|---|
| [`dataset.md`](docs/dataset.md) | Dataset provenance and experimental design |
| [`methodology.md`](docs/methodology.md) | Analytical methodology |
| [`reference.md`](docs/reference.md) | Reference-resource specification |
| [`reproducibility.md`](docs/reproducibility.md) | Reproducibility strategy |
| [`benchmarking.md`](docs/benchmarking.md) | Benchmarking strategy |
| [`frozen_decisions.md`](docs/frozen_decisions.md) | Locked technical decisions |
| [`project_plan.md`](docs/project_plan.md) | Engineering roadmap |

---

## Scope & Limitations

RNAFlowX is currently an **engineering and bioinformatics portfolio platform**, not a clinical diagnostic workflow.

The four-sample dataset was deliberately selected to demonstrate complete end-to-end workflow execution on limited local hardware.

Because of the small biological sample size, biological findings should be interpreted as workflow-demonstration results rather than definitive experimental conclusions.

RNAFlowX is not intended for clinical decision-making.

---

## Future Architecture

The long-term engineering direction is:

```text
Input Data
    │
    ▼
Nextflow DSL2
    │
    ├── QC
    ├── Preprocessing
    ├── Quantification
    ├── Differential Expression
    ├── Functional Analysis
    └── Reporting
            │
            ▼
      Structured Results
            │
       ┌────┴────┐
       ▼         ▼
    Reports    Data Layer
                  │
                  ▼
             Dashboard
```

Future infrastructure work will focus on execution portability, automated validation, provenance, observability, and deployment rather than changing the core scientific workflow without justification.

---

## Author

**Sriram B**

B.Tech Biotechnology

**Bioinformatics • Workflow Engineering • Nextflow • Reproducible Computational Biology**

---

## License

This project is distributed under the terms defined in [`LICENSE`](LICENSE).

---

## Citation

Citation metadata for RNAFlowX is provided through [`CITATION.cff`](CITATION.cff).

