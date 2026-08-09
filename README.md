## Frozen Project Configuration

RNAflow-X uses a deliberately constrained four-sample human bulk RNA-seq benchmark dataset to enable complete local execution without paid cloud infrastructure.

The dataset, reference, annotation, quantification strategy, execution constraints, reproducibility requirements, observability requirements, and benchmarking strategy are documented in:

- [Frozen Technical Decisions](docs/frozen_decisions.md)
- [Dataset](docs/dataset.md)
- [Methodology](docs/methodology.md)
- [Benchmarking](docs/benchmarking.md)
- [Reproducibility](docs/reproducibility.md)

The input sample manifest is:

```text
assets/samplesheet.csv