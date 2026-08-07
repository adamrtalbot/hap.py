---
title: Troubleshooting
description: Common input, region, output, and comparison problems.
---

## Reference or contig errors

Check these inputs:

- the FASTA has a neighboring `.fai` index;
- each VCF, BED, and FASTA uses the same assembly;
- all contig names use or omit the `chr` prefix;
- reference alleles match the FASTA at each VCF position.

Use [`hap validate`](../../tools/validate/) to identify reference mismatches and
out-of-range records. Use the relevant `--fixchr` option for a consistent
prefix difference.

## Missing output directory

Create the parent directory before running commands that take a report prefix:

```bash
mkdir -p results
hap germline truth.vcf.gz query.vcf.gz \
  -r reference.fa -o results/sample
```

## Unexpected precision denominator

Pass the confident-region BED with `--false-positives`. `hap` labels unmatched
query calls outside those intervals as `UNK`.

## Restrict BED overlap

The xcmp engine rejects overlapping `-R, --restrict-regions` intervals. Merge
the BED first or use `-T, --target-regions` when overlaps are intentional.

## Different results between engines

Engines use different matching semantics. Record-level, haplotype-level,
vcfeval-style, and distance-based engines can classify a complex representation
in different ways. Record the engine in benchmarking methods and compare runs
that use the same engine.

## Verification failure

Open the lane's `comparison.json` and start with its first reported artifact and
location. A shared missing output also fails the gate. See
[Verification](../../project/verification/) for focused runs and allowed
metadata exclusions.
