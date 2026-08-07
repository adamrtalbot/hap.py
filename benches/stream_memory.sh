#!/usr/bin/env bash
set -euo pipefail

# Peak-RSS benchmark for the streaming validation path. It creates two inputs
# with identical record widths but a 24x record-count difference. Streaming is
# demonstrated when max RSS remains approximately flat while elapsed work and
# input bytes grow. macOS and GNU time use different RSS labels, so retain the
# complete reports as benchmark artifacts.
root=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/hap-rs-memory.XXXXXX")
trap 'rm -rf "$work"' EXIT

cargo build --locked --release --bin hap --manifest-path "$root/Cargo.toml"

make_vcf() {
  destination=$1
  records=$2
  awk -v records="$records" 'BEGIN {
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=chr1,length=250000000>"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    for (i = 1; i <= records; i++)
      print "chr1\t" i "\t.\tA\tC\t50\tPASS\t.\tGT\t0/1"
  }' > "$destination"
}

run_case() {
  name=$1
  records=$2
  input="$work/$name.vcf"
  output="$work/$name.json"
  make_vcf "$input" "$records"
  if [[ $(uname -s) == Darwin ]]; then
    if ! /usr/bin/time -l "$root/target/release/hap" validate "$input" -o "$output" \
      2> "$work/$name.time.txt"; then
      cat "$work/$name.time.txt" >&2
      return 1
    fi
  else
    if ! /usr/bin/time -v "$root/target/release/hap" validate "$input" -o "$output" \
      2> "$work/$name.time.txt"; then
      cat "$work/$name.time.txt" >&2
      return 1
    fi
  fi
  printf '%s input_bytes=%s\n' "$name" "$(wc -c < "$input")"
  cat "$work/$name.time.txt"
}

run_case chromosome 100000
run_case whole_genome 2400000

if [[ $(uname -s) == Darwin ]]; then
  chromosome_rss=$(awk '/maximum resident set size/ {print $1}' "$work/chromosome.time.txt")
  whole_genome_rss=$(awk '/maximum resident set size/ {print $1}' "$work/whole_genome.time.txt")
else
  chromosome_rss=$(awk -F ': ' '/Maximum resident set size/ {print $2 * 1024}' "$work/chromosome.time.txt")
  whole_genome_rss=$(awk -F ': ' '/Maximum resident set size/ {print $2 * 1024}' "$work/whole_genome.time.txt")
fi
growth=$((whole_genome_rss - chromosome_rss))
limit=$((64 * 1024 * 1024))
printf 'rss_growth_bytes=%s limit_bytes=%s\n' "$growth" "$limit"
if (( growth > limit )); then
  printf 'peak RSS grew beyond the streaming memory limit\n' >&2
  exit 1
fi
