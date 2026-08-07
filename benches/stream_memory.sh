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
  contigs=$2
  bases=$3
  awk -v contigs="$contigs" -v bases="$bases" 'BEGIN {
    for (chrom=1;chrom<=contigs;chrom++) {
      print ">chr" chrom
      for (i=0;i<bases;i++) printf "A"
      print ""
    }
  }' > "$destination"
}

make_vcf() {
  destination=$1
  records=$2
  filter_every=${3:-0}
  contigs=${4:-1}
  awk -v records="$records" -v filter_every="$filter_every" -v contigs="$contigs" 'BEGIN {
    print "##fileformat=VCFv4.2"
    for (chrom=1;chrom<=contigs;chrom++)
      print "##contig=<ID=chr" chrom ",length=10000000>"
    print "##FILTER=<ID=LowQual,Description=\"Synthetic filtered call\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"High-cardinality ROC score\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    emitted=0
    per_contig=int((records+contigs-1)/contigs)
    for (chrom=1;chrom<=contigs && emitted<records;chrom++) {
      for (pos=1;pos<=per_contig && emitted<records;pos++) {
        emitted++
        filter=(filter_every && emitted%filter_every==0) ? "LowQual" : "PASS"
        score=emitted % 100000
        print "chr" chrom "\t" pos * 100 "\t.\tA\tC\t" score "\t" filter "\t.\tGT:QQ\t0/1:" score
      }
    }
  }' > "$destination"
}

peak_rss() {
  report=$1
  shift
  if [[ $(uname -s) == Darwin ]]; then
    if ! /usr/bin/time -l "$@" >/dev/null 2> "$report"; then
      printf 'workflow failed; time report: %s\n' "$report" >&2
      cat "$report" >&2
      return 1
    fi
    rss=$(awk '/maximum resident set size/ {print $1}' "$report")
  else
    if ! /usr/bin/time -v "$@" >/dev/null 2> "$report"; then
      printf 'workflow failed; time report: %s\n' "$report" >&2
      cat "$report" >&2
      return 1
    fi
    rss=$(awk -F ': ' '/Maximum resident set size/ {printf "%.0f", $2 * 1024}' "$report")
  fi
  if [[ ! $rss =~ ^[0-9]+$ ]] || (( rss <= 0 )); then
    printf 'invalid peak RSS %q in %s\n' "$rss" "$report" >&2
    return 1
  fi
  printf '%s\n' "$rss"
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
    rss=$(peak_rss "$report" "$@") || return 1
    printf '%s\n' "$rss" >> "$rss_file"
  done
  result=$(median < "$rss_file")
  if [[ ! $result =~ ^[0-9]+$ ]] || (( result <= 0 )); then
    printf 'invalid median RSS for %s/%s: %q\n' "$workflow" "$size" "$result" >&2
    return 1
  fi
  printf '%s\n' "$result"
}

reference="$work/reference.fa"
make_reference "$reference" 24 10000000

benchmark_size() {
  size=$1
  records=$2
  contigs=$3
  input="$work/$size.vcf"
  query="$work/$size.query.vcf"
  make_vcf "$input" "$records" 17 "$contigs"
  make_vcf "$query" "$records" 0 "$contigs"

  measurement=$(run_repeated validate "$size" "$hap_bin" validate "$input" -r "$reference" \
    --errors-bed "$work/$size.errors.bed" -o "$work/$size.validate.json") || return 1
  printf 'validate\t%s\n' "$measurement"
  measurement=$(run_repeated preprocess "$size" "$hap_bin" pre "$input" "$work/$size.pre.vcf.gz" \
    -r "$reference") || return 1
  printf 'preprocess\t%s\n' "$measurement"
  measurement=$(run_repeated compare "$size" "$hap_bin" germline "$input" "$query" -r "$reference" \
    -o "$work/$size.compare" --usefiltered-truth --no-roc --no-write-counts) || return 1
  printf 'compare\t%s\n' "$measurement"
  measurement=$(run_repeated quantify_no_roc "$size" "$hap_bin" quantify \
    "$work/$size.compare.vcf.gz" -o "$work/$size.qfy-no-roc" -r "$reference" \
    --type ga4gh --no-roc) || return 1
  printf 'quantify_no_roc\t%s\n' "$measurement"
  measurement=$(run_repeated quantify_roc "$size" "$hap_bin" quantify \
    "$work/$size.compare.vcf.gz" -o "$work/$size.qfy-roc" -r "$reference" \
    --type ga4gh --roc QQ) || return 1
  printf 'quantify_roc\t%s\n' "$measurement"
  measurement=$(run_repeated somatic_features "$size" "$hap_bin" somatic "$input" "$query" \
    -r "$reference" -o "$work/$size.somatic" --feature-table generic --happy-stats) || return 1
  printf 'somatic_features\t%s\n' "$measurement"
}

benchmark_size chromosome 100000 1 > "$work/chromosome.medians"
benchmark_size whole_genome 2400000 24 > "$work/whole-genome.medians"
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
