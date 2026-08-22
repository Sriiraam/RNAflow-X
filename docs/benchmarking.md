# RNAFlowX Benchmarking

## 1. Objective

RNAFlowX was benchmarked using a fresh Docker-based execution to measure pipeline runtime, task-level resource utilization, and execution reliability.

## 2. Benchmark Dataset

- Dataset: GSE342612
- Samples: 4 paired-end RNA-seq samples
- Comparison: Vehicle control vs 50 µM PFOS
- Input size: approximately 308 MB compressed FASTQ
- Organism: Homo sapiens
- Cell model: HMC3 human microglial cells

## 3. Benchmark Method

The benchmark was executed using:

```bash
NXF_WORK=work_benchmark \
./scripts/benchmark.sh docker results_benchmark

4. Overall Results
Metric	Result
Total tasks	19
Successful tasks	19
Failed tasks	0
Cached tasks	0
Success rate	100%
Nextflow runtime	17m 43s
Wall-clock runtime	17m 49.87s
CPU hours	0.4
Highest observed task RSS	~2.5 GB

5. Process-Level Performance
Process	Runtime	Peak RSS
FastQC Raw	22.7–58.0 s	203–247 MB
FastP	20.6–43.5 s	~1.2 GB
FastQC Trimmed	19.9–59.8 s	210–242 MB
Salmon	2m 10s–6m 49s	2.4–2.5 GB
tximport	41.0 s	558.5 MB
DESeq2	17.9 s	872.4 MB
MultiQC	10.2 s	199.2 MB

6. Bottleneck Analysis

Salmon quantification was the primary computational bottleneck.

Across the four samples, Salmon required approximately 2m 10s to 6m 49s per sample and approximately 2.4–2.5 GB peak resident memory.

The remaining workflow stages had substantially lower runtime and memory requirements.

7. Benchmark Artifacts

Benchmark evidence is stored under:

benchmark/
├── environment/
│   ├── container_info.txt
│   ├── software_versions.txt
│   └── system_info.txt
├── runs/
│   └── docker_local/
│       ├── execution_report.html
│       ├── execution_timeline.html
│       ├── execution_trace.txt
│       ├── system_metrics.txt
│       └── workflow_dag.html
└── summaries/
    ├── benchmark_summary.csv
    └── process_metrics.csv

8. Reproducibility

The benchmark records:

Hardware and operating-system information
CPU and memory information
Nextflow version
Java version
Docker version
RNAFlowX container information
Task-level execution metrics
Workflow timeline
Workflow DAG
System-level timing

9. Limitations

This benchmark represents one fresh execution using a deliberately small four-sample dataset on local hardware.

The results demonstrate reproducibility and local execution performance and should not be interpreted as large-scale HPC or cloud performance measurements.

Performance may vary with sequencing depth, sample count, hardware, storage performance, reference data, and execution environment.

SLURM/HPC and Azure Batch profiles are configuration-ready, but they have not been performance benchmarked.

10. Benchmark Status

Docker local benchmark: COMPLETED

Fresh execution: PASS
Tasks: 19/19 PASS
Failed tasks: 0
Cached tasks: 0
Benchmark artifacts: GENERATED
Environment provenance: GENERATED
Machine-readable summaries: GENERATED