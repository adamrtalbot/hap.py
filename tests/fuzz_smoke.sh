#!/usr/bin/env bash
set -euo pipefail

runs=${1:-1000}
max_total_time=${FUZZ_MAX_TOTAL_TIME:-30}
smoke_root=$(mktemp -d)
trap 'rm -rf "$smoke_root"' EXIT

for target in vcf bcf fasta bed location; do
  mkdir "$smoke_root/$target"
  cp -R "fuzz/corpus/$target/." "$smoke_root/$target"
  if [[ "$target" == bcf ]]; then
    cp "verification/assets/fixtures/pre-matrix/input.bcf" \
      "$smoke_root/$target/seed-valid-native.bcf"
  fi
  cargo +nightly fuzz run "$target" "$smoke_root/$target" -- \
    -runs="$runs" -max_total_time="$max_total_time"
done
