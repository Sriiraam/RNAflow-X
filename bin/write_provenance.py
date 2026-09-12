#!/usr/bin/env python3

import argparse
import json
import subprocess
from pathlib import Path


def run(cmd):
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except Exception:
        return "unavailable"


parser = argparse.ArgumentParser()
parser.add_argument("--output", required=True)
parser.add_argument("--pipeline-version", required=True)
parser.add_argument("--params-json", required=True)
parser.add_argument("--run-name", required=True)
parser.add_argument("--session-id", required=True)
parser.add_argument("--start", required=True)
parser.add_argument("--complete", required=True)
parser.add_argument("--duration", required=True)
parser.add_argument("--success", required=True)
parser.add_argument("--profile", default="")
parser.add_argument("--command-line", default="")
parser.add_argument("--nextflow-version", required=True)
parser.add_argument("--container-image", required=True)
parser.add_argument("--samplesheet", required=True)
parser.add_argument("--salmon-index", required=True)
parser.add_argument("--tx2gene", required=True)
parser.add_argument("--metadata", required=True)
args = parser.parse_args()

try:
    resolved_params = json.loads(args.params_json)
except json.JSONDecodeError as exc:
    raise SystemExit(
        f"Invalid --params-json supplied to provenance writer: {exc}"
    )

container_id = run([
    "docker",
    "image",
    "inspect",
    args.container_image,
    "--format={{.Id}}"
])

git_commit = run(["git", "rev-parse", "HEAD"])
git_branch = run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
git_status = run(["git", "status", "--porcelain"])

data = {
    "pipeline": {
        "name": "RNAFlowX",
        "version": args.pipeline_version
    },
    "git": {
        "commit": git_commit,
        "branch": git_branch,
        "dirty": git_status not in ("", "unavailable")
    },
    "workflow": {
        "run_name": args.run_name,
        "session_id": args.session_id,
        "start": args.start,
        "complete": args.complete,
        "duration": args.duration,
        "success": args.success.lower() == "true",
        "profile": args.profile,
        "command_line": args.command_line,
        "nextflow_version": args.nextflow_version
    },
    "inputs": {
        "samplesheet": args.samplesheet,
        "salmon_index": args.salmon_index,
        "tx2gene": args.tx2gene,
        "metadata": args.metadata
    },
    "parameters": resolved_params,
    "container": {
        "image": args.container_image,
        "image_id": container_id
    },
    "checksums": {
        "fastq_manifest": "checksums/fastq.sha256",
        "tx2gene_manifest": "checksums/tx2gene.sha256",
        "salmon_index_manifest": "checksums/salmon_index.sha256"
    }
}

out = Path(args.output)
out.parent.mkdir(parents=True, exist_ok=True)

with out.open("w") as fh:
    json.dump(data, fh, indent=2)

print(f"Provenance written to {out}")
