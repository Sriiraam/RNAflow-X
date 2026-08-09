# RNAflow-X Benchmarking Plan

## 1. Objective

Benchmarking is used to quantify the computational performance of RNAflow-X rather than simply reporting that the pipeline completed successfully.

## 2. Metrics

The workflow will capture:

- Total pipeline execution time
- Per-process execution time
- CPU usage where available
- Peak memory usage
- Disk usage
- Number of input samples
- Input data size
- Output data size
- Number of completed processes
- Cache/resume behavior

## 3. Nextflow Reports

The following Nextflow execution reports will be generated:

- Timeline report
- Execution report
- Trace report
- DAG visualization

## 4. Benchmark Scenarios

At minimum, the project will evaluate:

### Benchmark 1 — Clean execution

Complete pipeline execution from an empty work directory.

### Benchmark 2 — Resume execution

Restart the workflow using Nextflow `-resume` after an intentional interruption or after completion.

The purpose is to demonstrate process caching and reproducibility.

### Benchmark 3 — Resource utilization

Record memory, CPU and storage behavior during execution.

## 5. Benchmark Context

Benchmark results must always be interpreted together with:

- Hardware specification
- Dataset size
- Number of samples
- Software versions
- Container versions
- Pipeline commit/version

This prevents benchmark numbers from being presented without computational context.

## 6. Benchmark Reproducibility

Each benchmark record should contain:

```text
Pipeline version
Git commit
Nextflow version
Container version
Operating environment
CPU
RAM
Storage
Dataset
Reference
Parameters
Execution date