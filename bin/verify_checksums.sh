#!/usr/bin/env bash
set -euo pipefail

echo "Verifying FASTQ checksums..."
sha256sum -c checksums/fastq.sha256

echo
echo "Verifying tx2gene checksum..."
sha256sum -c checksums/tx2gene.sha256

echo
echo "Verifying Salmon index checksums..."
sha256sum -c checksums/salmon_index.sha256

echo
echo "All RNAFlowX checksums passed."
