#!/usr/bin/env bash
set -euo pipefail

PROFILE="${1:-docker}"
OUTDIR="${2:-results_benchmark}"
WORKDIR="${3:-work_benchmark}"

BENCHMARK_DIR="benchmark"
RUN_DIR="${BENCHMARK_DIR}/runs/${PROFILE}_local"
ENV_DIR="${BENCHMARK_DIR}/environment"
SUMMARY_DIR="${BENCHMARK_DIR}/summaries"

mkdir -p "$RUN_DIR" "$ENV_DIR" "$SUMMARY_DIR"

echo "========================================"
echo " RNAFlowX Benchmark"
echo "========================================"
echo "Profile : $PROFILE"
echo "Output  : $OUTDIR"
echo "Work    : $WORKDIR"
echo "Started : $(date -u)"
echo

# ------------------------------------------------------------
# Environment provenance
# ------------------------------------------------------------

{
    echo "RNAFlowX Benchmark Environment"
    echo "Generated: $(date -u)"
    echo
    uname -a
    echo
    echo "CPU:"
    lscpu
    echo
    echo "Memory:"
    free -h
} > "$ENV_DIR/system_info.txt"

{
    echo "Nextflow:"
    nextflow -version
    echo
    echo "Docker:"
    docker --version
    echo
    echo "Java:"
    java -version 2>&1
} > "$ENV_DIR/software_versions.txt" 2>&1

{
    echo "Image: rnaflowx:1.0.0"
    docker image inspect rnaflowx:1.0.0 \
        --format 'ID={{.Id}} Size={{.Size}} Created={{.Created}}'
} > "$ENV_DIR/container_info.txt"

# ------------------------------------------------------------
# Fresh benchmark protection
# ------------------------------------------------------------

if [ -d "$WORKDIR" ]; then
    echo "ERROR: $WORKDIR already exists."
    echo "Benchmark aborted to prevent cache reuse."
    exit 1
fi

if [ -d "$OUTDIR" ]; then
    echo "ERROR: $OUTDIR already exists."
    echo "Benchmark aborted to protect previous results."
    exit 1
fi

# ------------------------------------------------------------
# Benchmark execution
# ------------------------------------------------------------

export NXF_WORK="$WORKDIR"

/usr/bin/time -v \
    -o "$RUN_DIR/system_metrics.txt" \
    nextflow run main.nf \
    -profile "$PROFILE" \
    --outdir "$OUTDIR"

# ------------------------------------------------------------
# Collect Nextflow benchmark artifacts
# ------------------------------------------------------------

cp "$OUTDIR/pipeline_info/execution_trace.txt" \
   "$RUN_DIR/execution_trace.txt"

cp "$OUTDIR/pipeline_info/execution_report.html" \
   "$RUN_DIR/execution_report.html"

cp "$OUTDIR/pipeline_info/execution_timeline.html" \
   "$RUN_DIR/execution_timeline.html"

cp "$OUTDIR/pipeline_info/workflow_dag.html" \
   "$RUN_DIR/workflow_dag.html"

# ------------------------------------------------------------
# Generate benchmark summaries
# ------------------------------------------------------------

python scripts/summarize_benchmark.py \
    "$RUN_DIR/execution_trace.txt" \
    "$SUMMARY_DIR"

echo
echo "========================================"
echo " RNAFlowX Benchmark Complete"
echo "========================================"
echo "Finished: $(date -u)"
echo
echo "Raw metrics : $RUN_DIR"
echo "Environment : $ENV_DIR"
echo "Summaries   : $SUMMARY_DIR"
