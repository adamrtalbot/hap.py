#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fixture_dir="$repo_dir/verification/assets/somatic_options"
hap_bin=${HAP_BIN:-"$repo_dir/target/release/hap"}
legacy_image=${LEGACY_IMAGE:-community.wave.seqera.io/library/hap.py_rtg-tools@sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6}
work_dir=$(mktemp -d "$repo_dir/verification/.somatic-options-oracle.XXXXXX")

cleanup() {
    if [[ ${KEEP_ORACLE_WORK:-0} == 1 ]]; then
        printf 'Oracle work retained at %s\n' "$work_dir"
        return
    fi
    find "$work_dir" -type f -delete
    find "$work_dir" -depth -type d -exec rmdir {} +
}
trap cleanup EXIT

command -v docker >/dev/null
command -v python3 >/dev/null
test -x "$hap_bin"

canonical_json() {
    python3 - "$1" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    payload = json.load(handle)
payload.pop("timestamp", None)
json.dump(payload, sys.stdout, sort_keys=True, separators=(",", ":"))
PY
}

compare_case() {
    local case_name=$1
    local legacy_dir="$work_dir/$case_name/legacy"
    local rust_dir="$work_dir/$case_name/rust"
    local legacy_names
    local rust_names

    legacy_names=$(cd "$legacy_dir" && find . -maxdepth 1 -type f -name 'result.*' -print | sort)
    rust_names=$(cd "$rust_dir" && find . -maxdepth 1 -type f -name 'result.*' -print | sort)
    test -n "$legacy_names"
    cmp <(printf '%s\n' "$legacy_names") <(printf '%s\n' "$rust_names")

    if [[ $case_name == explain_without_features ]]; then
        test -f "$legacy_dir/result.ambiclasses.csv"
        test -f "$legacy_dir/result.ambireasons.csv"
        test ! -e "$legacy_dir/result.features.csv"
    fi

    while IFS= read -r relative; do
        test -n "$relative" || continue
        if [[ $relative == ./result.metrics.json ]]; then
            cmp <(canonical_json "$legacy_dir/$relative") \
                <(canonical_json "$rust_dir/$relative")
        else
            cmp "$legacy_dir/$relative" "$rust_dir/$relative"
        fi
        printf 'PASS %s %s\n' "$case_name" "${relative#./}"
    done <<<"$legacy_names"
}

run_case() {
    local case_name=$1
    shift
    local legacy_dir="$work_dir/$case_name/legacy"
    local rust_dir="$work_dir/$case_name/rust"
    mkdir -p "$legacy_dir" "$rust_dir"

    docker run --rm --platform linux/amd64 \
        -v "$repo_dir:$repo_dir" \
        -w "$legacy_dir" \
        "$legacy_image" \
        bash -c '
            # The pinned image combines som.py with a pandas release that
            # removed this display-only option. Removing the dead formatting
            # call unlocks the legacy explanation serializer without changing
            # comparison, metrics, or CSV behavior.
            sed -i '\''/pandas.set_option("display.height"/d'\'' /opt/conda/bin/som.py
            exec som.py "$@"
        ' somatic-oracle \
        "$fixture_dir/truth.vcf" "$fixture_dir/query.vcf" -o result "$@"

    (
        cd "$rust_dir"
        "$hap_bin" somatic "$fixture_dir/truth.vcf" "$fixture_dir/query.vcf" \
            -o result "$@"
    )

    compare_case "$case_name"
}

missing_reference="$fixture_dir/missing.fa"

# Explanation tables do not depend on feature extraction. An explicit FP size
# also proves that an absent reference remains harmless on this path.
run_case explain_without_features \
    --reference "$missing_reference" \
    --fp-region-size 10 \
    --ambiguous "$fixture_dir/ambiguous.bed" \
    --explain_ambiguous --quiet

# A nonempty FP BED supplies the denominator, so som.py does not open the
# explicitly missing reference. The same BED is valid with its extra columns.
run_case fp_bed_without_reference \
    --reference "$missing_reference" \
    --false-positives "$fixture_dir/ambiguous.bed" \
    --quiet

# The multi-character -FN reporting switch must not suppress -N normalization.
run_case normalize_all_filtered_fn \
    --reference "$fixture_dir/reference.fa" \
    -P -N -FN --feature-table generic --fp-region-size 10 --quiet
