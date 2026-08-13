# HG001 float rendering

This example case preserves the germline report-rendering discrepancy found by
comparing the public GIAB HG001 v4.2.1 truth set with the Illumina Platinum
Genomes NA12878 query:

- truth: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
- query: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/analysis/Illumina_PlatinumGenomes_NA12877_NA12878_09162015/hg38/2.0.1/NA12878/NA12878.vcf.gz`

The local 100-base reference and 28 VCF records reproduce the affected numeric
shapes: a `1/7` recall and an adjacent-double ratio. The ordinary HAPPY row runs
the pinned legacy implementation and hap-rs with `--no-roc`, then compares
`summary.csv`, `extended.csv`, the aggregate `roc.all.csv.gz`, and
`metrics.json.gz` through the unchanged paired-output comparator. Per-threshold
ROC tables stay disabled so their unrelated column-type behavior is not added
to this case.

The diagnosed cause is a dependent serialization stack: report cells require
Python 2's 12-significant-digit rendering after the pandas parser round trip,
metrics JSON must retain the parsed numeric value with Python's full float
representation, and notation decisions must use the legacy thresholds without
rounding adjacent `f64` values across them. Commits `eb5b1bd`, `5eddea6`, and
`fc0202e`, in that order, implement those corrections.
