# vcfcheck adversarial fixtures

These compact fixtures were authored for this parity audit and are released
under the repository license. They isolate stable behaviors of the pinned
Illumina hap.py v0.3.15 `vcfcheck` oracle:

- `missing-contig.vcf` exercises htslib VCF-to-BCF translation checks when a
  record contig has no header declaration.
- `sample-count-mismatch.vcf` has one declared sample and two record sample
  columns.
- `lowercase-x.vcf` distinguishes literal legacy `X`/`chrX` sex inference
  from case-insensitive contig matching.

All three files are synthetic; they contain no imported genomic data.
