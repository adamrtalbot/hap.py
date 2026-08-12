import csv
import json
from pathlib import Path
import sys


prefix = Path(sys.argv[1])
comparison = Path(sys.argv[2])
stats_path = prefix.with_name(f"{prefix.name}.stats.csv")
metrics_path = prefix.with_name(f"{prefix.name}.metrics.json")
bins = [
    "0.000000-0.200000",
    "0.200000-0.400000",
    "0.400000-0.600000",
    "0.600000-0.800000",
    "0.800000-1.000000",
]
expected_types = ["indels", "SNVs", "no-ALTs", "records", "MNPs", "others"] + [
    f"{prefix}.{interval}"
    for prefix in ["records", "SNVs", "indels"]
    for interval in bins
]
errors = []

with stats_path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
actual_types = [row["type"] for row in rows]
if actual_types != expected_types:
    errors.append(
        f"{stats_path.name}: type rows={actual_types!r}, expected {expected_types!r}"
    )

indels = next((row for row in rows if row["type"] == "indels"), None)
expected_indel_metrics = {
    "recall": "0.3333333333333333",
    "recall_upper": "0.8232639028687426",
    "recall2": "0.3333333333333333",
}
if indels is None:
    errors.append(f"{stats_path.name}: missing indels row")
else:
    for column, expected in expected_indel_metrics.items():
        actual = indels[column]
        if actual != expected:
            errors.append(
                f"{stats_path.name}: indels {column}={actual!r}, expected {expected!r}"
            )

with metrics_path.open(encoding="utf-8") as handle:
    data = json.load(handle)["metrics"][0]["data"]
for column in data:
    values = column.get("values", [])
    if len(values) != len(expected_types):
        errors.append(
            f"{metrics_path.name}: {column.get('id')} values={len(values)}, "
            f"expected {len(expected_types)}"
        )

if errors:
    raise SystemExit("\n".join(errors))

comparison.write_text(
    json.dumps({"case_id": "somatic_output_parity", "lane": "sompy", "ok": True})
    + "\n",
    encoding="utf-8",
)
