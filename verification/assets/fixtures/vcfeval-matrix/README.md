# vcfeval matrix reference bundle

The verification fixture stores an RTG SDF in `ref.sdf.tar.gz`. RTG Tools 3.12.1
created the SDF from `ref.fa` with `rtg format -o ref.sdf ref.fa`. The archive
uses the `ref.sdf` directory as its root.

- RTG source tag: `3.12.1`
- RTG source commit: `32d4c2d2d340cb0288f68ddd008d29564bcfef13`
- Bundle SHA-256: `aacd186572c8d3726eb7ca2b731bb8d9e5df1f995ebd87b262de6045e4734b89`

Nextflow stages this bundle into `HAPPY_LEGACY` and gives `ref.fa` to
`HAPPY_RUST`. The Rust process has no bundle input.
