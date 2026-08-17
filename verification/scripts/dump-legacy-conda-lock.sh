#!/usr/bin/env bash
# Regenerate containers/happy-0.3.15.conda-lock.txt from the frozen Wave image
# used by HAPPY, PREPY, QFY, and VCFCHECK. The image ships no conda CLI, so the
# explicit lock is rebuilt from its conda-meta records.
set -euo pipefail

image=${LEGACY_IMAGE:-community.wave.seqera.io/library/happy-0.3.15:41c2102638513597}

docker run --rm --platform linux/amd64 --entrypoint /opt/conda/bin/python "$image" -c '
import glob, json
rows = []
for path in glob.glob("/opt/conda/conda-meta/*.json"):
    record = json.load(open(path))
    url = record.get("url") or ""
    md5 = record.get("md5") or ""
    if not url:
        continue
    rows.append((record.get("name", ""), url + ("#" + md5 if md5 else "")))
rows.sort()
print("@EXPLICIT")
for _, url in rows:
    print(url)
'
