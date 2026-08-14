import csv
import gzip
import json
from pathlib import Path
import sys


prefix = Path(sys.argv[1])
comparison = Path(sys.argv[2])
expected_metric = "0.142857"
expected_ratio = "1.33333333333"
csv_suffixes = [
    "summary.csv",
    "extended.csv",
    "roc.all.csv.gz",
    "roc.Locations.SNP.csv.gz",
    "roc.Locations.SNP.PASS.csv.gz",
    "roc.Locations.INDEL.csv.gz",
    "roc.Locations.INDEL.PASS.csv.gz",
]
errors = []

for suffix in csv_suffixes:
    path = prefix.with_name(f"{prefix.name}.{suffix}")
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not any(row.get("METRIC.Recall") == expected_metric for row in rows):
        errors.append(f"{path.name}: missing METRIC.Recall={expected_metric}")
    if suffix in {"summary.csv", "extended.csv", "roc.all.csv.gz", "roc.Locations.SNP.csv.gz", "roc.Locations.SNP.PASS.csv.gz"}:
        if not any(expected_ratio in row.values() for row in rows):
            errors.append(f"{path.name}: missing ratio={expected_ratio}")

metrics_path = prefix.with_name(f"{prefix.name}.metrics.json.gz")
with gzip.open(metrics_path, "rt", encoding="utf-8") as handle:
    metrics = json.load(handle)["metrics"]

summary = next((table for table in metrics if table.get("id") == "summary.metrics"), None)
if summary is None:
    errors.append(f"{metrics_path.name}: missing table summary.metrics")
else:
    recall = next((item for item in summary["data"] if item.get("id") == "METRIC.Recall"), None)
    if recall is None or 0.14285699999999998 not in recall.get("values", []):
        errors.append(
            f"{metrics_path.name}: missing numeric METRIC.Recall=0.14285699999999998"
        )

location_ids = {
    "roc.Locations.SNP.PASS",
    "roc.Locations.INDEL",
    "roc.Locations.SNP",
    "roc.Locations.INDEL.PASS",
}
tables = {table["id"]: table for table in metrics if table.get("id") in location_ids}
for table_id in sorted(location_ids):
    table = tables.get(table_id)
    if table is None:
        errors.append(f"{metrics_path.name}: missing table {table_id}")
        continue
    column = next((item for item in table["data"] if item.get("id") == "Subset.Size"), None)
    if column is None:
        errors.append(f"{table_id}: missing Subset.Size")
        continue
    if column.get("type") != "string":
        errors.append(f"{table_id}: Subset.Size type={column.get('type')!r}, expected 'string'")
    if not all(type(value) is str for value in column.get("values", [])):
        errors.append(f"{table_id}: Subset.Size values are not all JSON strings")

if errors:
    raise SystemExit("\n".join(errors))

comparison.write_text(
    json.dumps({"case_id": "germline_output_parity", "lane": "happy", "ok": True}) + "\n",
    encoding="utf-8",
)
