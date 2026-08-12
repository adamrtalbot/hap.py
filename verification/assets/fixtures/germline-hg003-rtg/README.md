# HG003 RTG germline case

This verification row uses the public GIAB HG003 v4.2.1 truth set and the RTG
Ashkenazim Trio `family.merged.avr_0.1.vcf.gz` query, subsampled to `NA24149`.
The samplesheet retains the public source URLs. The files used to reproduce the
failure had these SHA-256 checksums:

- Query VCF: `fc58680f811d1cc3ac522756e2c26d436fcae2c24b0057d06c7eef128e0b8cfe`
- Query TBI: `ea8338ac9c7ee91f4e22df1c61585163b070746339b113acbddf40fb171b6352`
- Truth VCF: `a4f8fefa826f8d2c9457eacf6f4b2f6bfa29daef429d570ce04254c5e7b1121e`
- Truth TBI: `989305947e5af2e81ca8ad1087262a12cc87499b36163c49a7878dc9e18fcaa5`
- Confidence BED: `652afd3046705af3200f9c87c255fef11bb212dd76c75a19999c9b2df8a3180c`
- Uncompressed reference FASTA: `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`
- Uncompressed reference FAI: `502a1b8fb73ccd53285c28a0f12df90c818b4fe3de1e862ef47c593ef1a0a4b4`

The query and truth inputs pass through the same split, exact-deduplication,
sorting, normalization, and alternate-genotype filtering used by the public
VCF campaign before the legacy and Rust HAPPY processes run.
