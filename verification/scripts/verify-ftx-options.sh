#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fixture_dir="$repo_dir/verification/assets/fixtures/ftx-bam"
hap_bin=${HAP_BIN:-"$repo_dir/target/debug/hap"}
legacy_image=${LEGACY_IMAGE:-community.wave.seqera.io/library/hap.py_rtg-tools@sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6}
work_dir=$(mktemp -d "$repo_dir/verification/.ftx-options-oracle.XXXXXX")

cleanup() {
    find "$work_dir" -type f -delete
    rmdir "$work_dir"
}
trap cleanup EXIT

command -v docker >/dev/null
command -v bcftools >/dev/null
test -x "$hap_bin"

bcftools view -Ob -o "$work_dir/input.bcf" "$fixture_dir/input.vcf"

run_case() {
    local case_name=$1
    local input=$2
    shift 2

    docker run --rm --platform linux/amd64 \
        -v "$repo_dir:$repo_dir" \
        "$legacy_image" \
        ftx.py "$input" -o "$work_dir/$case_name.legacy" \
        --feature-label oracle "$@"

    "$hap_bin" ftx "$input" -o "$work_dir/$case_name.rust" \
        --feature-label oracle "$@"

    cmp "$work_dir/$case_name.legacy.csv" "$work_dir/$case_name.rust.csv"
    printf 'PASS %s\n' "$case_name"
}

run_case location "$fixture_dir/options.vcf" \
    --feature-table generic --reference "$fixture_dir/missing.fa" -l 1:7
run_case restrict_regions "$fixture_dir/options.vcf" \
    --feature-table generic --reference "$fixture_dir/missing.fa" -R "$fixture_dir/select.bed"
run_case target_regions "$fixture_dir/options.vcf" \
    --feature-table generic --reference "$fixture_dir/missing.fa" -T "$fixture_dir/select.bed"
run_case fix_chr "$fixture_dir/options.vcf" \
    --feature-table generic --reference "$fixture_dir/missing.fa" --fix-chr
run_case normalize "$fixture_dir/normalize.vcf" \
    --feature-table generic --reference "$fixture_dir/ref.fa" --normalize
run_case bcf_input "$work_dir/input.bcf" \
    --feature-table hcc.strelka.snv --reference "$fixture_dir/missing.fa"
run_case multi_bam "$fixture_dir/input.vcf" \
    --feature-table hcc.strelka.snv --reference "$fixture_dir/missing.fa" \
    --bam "$fixture_dir/reads.bam" --bam "$fixture_dir/reads2.bam"
run_case mutect_snv "$fixture_dir/mutect.vcf" \
    --feature-table hcc.mutect.snv --reference "$fixture_dir/missing.fa"
run_case varscan2_snv "$fixture_dir/varscan2.vcf" \
    --feature-table hcc.varscan2.snv --reference "$fixture_dir/missing.fa"
run_case pisces_snv "$fixture_dir/pisces.vcf" \
    --feature-table hcc.pisces.snv --reference "$fixture_dir/missing.fa"
