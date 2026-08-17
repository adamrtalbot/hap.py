#!/usr/bin/env bash
# Regenerate an explicit conda lock from a pinned verification image. Both
# images ship no conda CLI, so the lock is rebuilt from their conda-meta
# records. The image is a required argument: a default would let the wrong
# image silently produce a governed artifact.
#
# Governed locks and the images they belong to are named in
# containers/*.conda-lock.txt headers and in verification/nextflow.config.
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $(basename "$0") <image>" >&2
    exit 2
fi

image=$1

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
