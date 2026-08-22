#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: summarize_benchmark.py <trace.txt> <summary_dir>"
    )

trace_file = Path(sys.argv[1])
summary_dir = Path(sys.argv[2])
summary_dir.mkdir(parents=True, exist_ok=True)

if not trace_file.exists():
    raise SystemExit(f"Trace file not found: {trace_file}")

with trace_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

if not rows:
    raise SystemExit("Trace file contains no process records.")

# Keep useful process-level benchmark fields
wanted = [
    "name",
    "status",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
]

available = [c for c in wanted if c in rows[0]]

process_file = summary_dir / "process_metrics.csv"

with process_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=available)
    writer.writeheader()

    for row in rows:
        writer.writerow({c: row.get(c, "") for c in available})

# Overall task summary
total = len(rows)
completed = sum(
    1 for r in rows
    if r.get("status", "").upper() in {"COMPLETED", "CACHED"}
)
cached = sum(
    1 for r in rows
    if r.get("status", "").upper() == "CACHED"
)
failed = sum(
    1 for r in rows
    if r.get("status", "").upper() == "FAILED"
)

summary_file = summary_dir / "benchmark_summary.csv"

with summary_file.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["metric", "value"])
    writer.writerow(["total_tasks", total])
    writer.writerow(["completed_tasks", completed])
    writer.writerow(["cached_tasks", cached])
    writer.writerow(["failed_tasks", failed])

print(f"Created: {process_file}")
print(f"Created: {summary_file}")
