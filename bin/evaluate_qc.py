#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def classify(value, pass_threshold, warn_threshold):
    if value >= pass_threshold:
        return "PASS"
    if value >= warn_threshold:
        return "WARN"
    return "FAIL"


parser = argparse.ArgumentParser()

parser.add_argument("--sample", required=True)
parser.add_argument("--fastp", required=True)
parser.add_argument("--salmon", required=True)

parser.add_argument("--q30-pass", type=float, required=True)
parser.add_argument("--q30-warn", type=float, required=True)

parser.add_argument("--retention-pass", type=float, required=True)
parser.add_argument("--retention-warn", type=float, required=True)

parser.add_argument("--mapping-pass", type=float, required=True)
parser.add_argument("--mapping-warn", type=float, required=True)

args = parser.parse_args()


with open(args.fastp) as fh:
    fastp = json.load(fh)

with open(args.salmon) as fh:
    salmon = json.load(fh)


summary = fastp.get("summary", {})

before = summary.get("before_filtering", {})
after = summary.get("after_filtering", {})


before_reads = before.get("total_reads")
after_reads = after.get("total_reads")

if before_reads is None or after_reads is None:
    raise ValueError("FASTP JSON missing total_reads")

if before_reads == 0:
    raise ValueError("FASTP before_filtering total_reads is zero")


retention_pct = (after_reads / before_reads) * 100.0


q30_rate = after.get("q30_rate")

if q30_rate is not None:
    q30_pct = q30_rate * 100.0
else:
    q30_bases = after.get("q30_bases")
    total_bases = after.get("total_bases")

    if q30_bases is None or total_bases in (None, 0):
        raise ValueError("FASTP JSON missing usable Q30 metrics")

    q30_pct = (q30_bases / total_bases) * 100.0


gc_content = after.get("gc_content")

if gc_content is None:
    gc_pct = None
else:
    gc_pct = gc_content * 100.0


mapping_pct = salmon.get("percent_mapped")

if mapping_pct is None:
    raise ValueError("Salmon meta_info.json missing percent_mapped")


q30_status = classify(
    q30_pct,
    args.q30_pass,
    args.q30_warn
)

retention_status = classify(
    retention_pct,
    args.retention_pass,
    args.retention_warn
)

mapping_status = classify(
    mapping_pct,
    args.mapping_pass,
    args.mapping_warn
)


statuses = [
    q30_status,
    retention_status,
    mapping_status
]

if "FAIL" in statuses:
    overall_status = "FAIL"
elif "WARN" in statuses:
    overall_status = "WARN"
else:
    overall_status = "PASS"


result = {
    "sample_id": args.sample,

    "metrics": {
        "q30_percent": round(q30_pct, 2),
        "read_retention_percent": round(retention_pct, 2),
        "gc_percent": None if gc_pct is None else round(gc_pct, 2),
        "salmon_mapping_percent": round(mapping_pct, 2)
    },

    "status": {
        "q30": q30_status,
        "read_retention": retention_status,
        "salmon_mapping": mapping_status,
        "overall": overall_status
    },

    "thresholds": {
        "q30": {
            "pass": args.q30_pass,
            "warn": args.q30_warn
        },
        "read_retention": {
            "pass": args.retention_pass,
            "warn": args.retention_warn
        },
        "salmon_mapping": {
            "pass": args.mapping_pass,
            "warn": args.mapping_warn
        }
    }
}


json_file = Path(f"{args.sample}_qc.json")

with json_file.open("w") as fh:
    json.dump(result, fh, indent=2)


tsv_file = Path(f"{args.sample}_qc.tsv")

with tsv_file.open("w", newline="") as fh:

    writer = csv.writer(fh, delimiter="\t")

    writer.writerow([
        "sample_id",
        "q30_percent",
        "q30_status",
        "read_retention_percent",
        "retention_status",
        "gc_percent",
        "salmon_mapping_percent",
        "mapping_status",
        "overall_status"
    ])

    writer.writerow([
        args.sample,
        f"{q30_pct:.2f}",
        q30_status,
        f"{retention_pct:.2f}",
        retention_status,
        "" if gc_pct is None else f"{gc_pct:.2f}",
        f"{mapping_pct:.2f}",
        mapping_status,
        overall_status
    ])


Path(f"{args.sample}_qc_status.txt").write_text(
    f"{args.sample}\t{overall_status}\n"
)

print(
    f"{args.sample}: "
    f"Q30={q30_pct:.2f}% ({q30_status}), "
    f"retention={retention_pct:.2f}% ({retention_status}), "
    f"mapping={mapping_pct:.2f}% ({mapping_status}), "
    f"overall={overall_status}"
)
