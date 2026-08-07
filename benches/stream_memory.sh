#!/usr/bin/env bash
set -euo pipefail

# Repeated peak-RSS benchmark for every primary variant workflow. This proves
# that these concrete command/options keep RSS growth below the configured
# allowance when record count grows 24x; it does not prove a universal memory
# bound for every input shape (dense comparison clusters and superloci are
# governed separately by documented active-window limits).
root=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/hap-rs-memory.XXXXXX")
reports=${HAP_MEMORY_REPORT_DIR:-"$root/target/stream-memory-reports"}
runs=${HAP_MEMORY_RUNS:-3}
mkdir -p "$reports"
trap 'rm -rf "$work"' EXIT

cargo build --locked --release --bin hap --manifest-path "$root/Cargo.toml"
hap_bin="$root/target/release/hap"

make_reference() {
  destination=$1
  bases=$2
  awk -v bases="$bases" 'BEGIN { print ">chr1"; for (i=0;i<bases;i++) printf "A"; print "" }' > "$destination"
}

make_vcf() {
  destination=$1
  records=$2
  filter_every=${3:-0}
  awk -v records="$records" -v filter_every="$filter_every" 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=chr1,length=250000000>"
    print "##FILTER=<ID=LowQual,Description=\"Synthetic filtered call\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    for (i=1;i<=records;i++) {
      filter=(filter_every && i%filter_every==0) ? "LowQual" : "PASS"
      print "chr1\t" i "\t.\tA\tC\t50\t" filter "\t.\tGT\t0/1"
    }
  }' > "$destination"
}

peak_rss() {
  report=$1
  shift
  if [[ $(uname -s) == Darwin ]]; then
    /usr/bin/time -l "$@" >/dev/null 2> "$report"
    awk '/maximum resident set size/ {print $1}' "$report"
  else
    /usr/bin/time -v "$@" >/dev/null 2> "$report"
    awk -F ': ' '/Maximum resident set size/ {print $2 * 1024}' "$report"
  fi
}

median() {
  sort -n | awk '{ values[NR]=$1 } END { print values[int((NR+1)/2)] }'
}

run_repeated() {
  workflow=$1
  size=$2
  shift 2
  rss_file="$work/$workflow.$size.rss"
  : > "$rss_file"
  for run in $(seq 1 "$runs"); do
    report="$reports/$workflow.$size.run-$run.time.txt"
    peak_rss "$report" "$@" >> "$rss_file"
  done
  median < "$rss_file"
}

reference="$work/reference.fa"
make_reference "$reference" 2500000

benchmark_size() {
  size=$1
  records=$2
  input="$work/$size.vcf"
  query="$work/$size.query.vcf"
  make_vcf "$input" "$records" 17
  make_vcf "$query" "$records" 0

  printf 'validate\t%s\n' "$(run_repeated validate "$size" "$hap_bin" validate "$input" -r "$reference" \
    --errors-bed "$work/$size.errors.bed" -o "$work/$size.validate.json")"
  printf 'preprocess\t%s\n' "$(run_repeated preprocess "$size" "$hap_bin" pre "$input" "$work/$size.pre.vcf.gz" \
    -r "$reference")"
  printf 'compare\t%s\n' "$(run_repeated compare "$size" "$hap_bin" germline "$input" "$query" -r "$reference" \
    -o "$work/$size.compare" --usefiltered-truth --no-roc --no-write-counts)"
  printf 'quantify_no_roc\t%s\n' "$(run_repeated quantify_no_roc "$size" "$hap_bin" quantify \
    "$work/$size.compare.vcf.gz" -o "$work/$size.qfy-no-roc" -r "$reference" \
    --type ga4gh --no-roc)"
  printf 'quantify_roc\t%s\n' "$(run_repeated quantify_roc "$size" "$hap_bin" quantify \
    "$work/$size.compare.vcf.gz" -o "$work/$size.qfy-roc" -r "$reference" \
    --type ga4gh --roc QQ)"
  printf 'somatic_features\t%s\n' "$(run_repeated somatic_features "$size" "$hap_bin" somatic "$input" "$query" \
    -r "$reference" -o "$work/$size.somatic" --feature-table generic --happy-stats)"
}

benchmark_size chromosome 100000 > "$work/chromosome.medians"
benchmark_size whole_genome 2400000 > "$work/whole-genome.medians"
paste "$work/chromosome.medians" "$work/whole-genome.medians" > "$reports/median-rss.tsv"

limit=$((64 * 1024 * 1024))
while read -r workflow chromosome repeated_workflow whole_genome; do
  [[ $workflow == "$repeated_workflow" ]]
  growth=$((whole_genome - chromosome))
  printf '%s\t%s\t%s\t%s\n' "$workflow" "$chromosome" "$whole_genome" "$growth"
  if (( growth > limit )); then
    printf 'peak RSS growth %s exceeds %s bytes\n' "$growth" "$limit" >&2
    exit 1
  fi
done < "$reports/median-rss.tsv"

printf 'reports=%s runs=%s growth_limit_bytes=%s\n' "$reports" "$runs" "$limit"
