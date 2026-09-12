#!/usr/bin/env python3

import csv
import glob
import json


files = sorted(glob.glob("*_qc.tsv"))

if not files:
    raise SystemExit("No per-sample QC TSV files found")


rows = []

for file in files:
    with open(file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows.extend(reader)


fieldnames = [
    "sample_id",
    "q30_percent",
    "q30_status",
    "read_retention_percent",
    "retention_status",
    "gc_percent",
    "salmon_mapping_percent",
    "mapping_status",
    "overall_status"
]


with open("qc_summary.csv", "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(rows)


with open("qc_summary.json", "w") as fh:
    json.dump(rows, fh, indent=2)


statuses = [row["overall_status"] for row in rows]

if "FAIL" in statuses:
    overall = "FAIL"
elif "WARN" in statuses:
    overall = "WARN"
else:
    overall = "PASS"


with open("qc_status.txt", "w") as fh:
    fh.write(f"RNAFlowX QC STATUS: {overall}\n")
    fh.write(f"Samples evaluated: {len(rows)}\n")
    fh.write(f"PASS: {statuses.count('PASS')}\n")
    fh.write(f"WARN: {statuses.count('WARN')}\n")
    fh.write(f"FAIL: {statuses.count('FAIL')}\n")


print(
    f"RNAFlowX QC: {overall} | "
    f"PASS={statuses.count('PASS')} "
    f"WARN={statuses.count('WARN')} "
    f"FAIL={statuses.count('FAIL')}"
)
