import gzip
import json
import pathlib
import sys


case_id, contract = sys.argv[1:]
differences = []
artifacts = sorted(path.name for path in pathlib.Path(".").glob("result*"))
required = {
    "result.extended.csv",
    "result.metrics.json.gz",
    "result.runinfo.json",
    "result.summary.csv",
    "result.vcf.gz",
    "result.vcf.gz.tbi",
}
missing = sorted(required - set(artifacts))
if missing:
    differences.append(
        {
            "kind": "artifact_set",
            "expected": sorted(required),
            "actual": artifacts,
            "reason": "required native vcfeval artifacts are missing",
        }
    )

stderr = pathlib.Path("stderr.log").read_text()
if contract == "deprecated_flags":
    warning = (
        "--engine-vcfeval-path and --engine-vcfeval-template "
        "are deprecated and ignored"
    )
    if warning not in stderr:
        differences.append(
            {
                "kind": "contract",
                "location": "/stderr",
                "expected": warning,
                "actual": stderr,
                "reason": "deprecated options did not emit the migration warning",
            }
        )
elif contract == "preserve_info" and pathlib.Path("result.vcf.gz").is_file():
    with gzip.open("result.vcf.gz", "rt", encoding="utf-8") as handle:
        records = [line for line in handle if line and not line.startswith("#")]
    decorated = any("RegionsExtent=" in line and "ctype=" in line for line in records)
    if not records or not decorated:
        differences.append(
            {
                "kind": "contract",
                "location": "/result.vcf.gz",
                "expected": "comparison metadata and RegionsExtent annotations",
                "actual": records[:3],
                "reason": "--preserve-info output was not decorated",
            }
        )

pathlib.Path("comparison.json").write_text(
    json.dumps(
        {
            "schema_version": 1,
            "lane": "happy",
            "case_id": case_id,
            "ok": not differences,
            "legacy_artifacts": [],
            "observed_artifacts": {"rust": artifacts},
            "differences": differences,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
