# vcfeval matrix reference bundle

`ref.sdf.tar.gz` is a verification-only RTG SDF generated from `ref.fa` with
RTG Tools 3.12.1 (`rtg format -o ref.sdf ref.fa`) and then archived with the
`ref.sdf` directory as the archive root.

- RTG source tag: `3.12.1`
- RTG source commit: `32d4c2d2d340cb0288f68ddd008d29564bcfef13`
- Bundle SHA-256: `aacd186572c8d3726eb7ca2b731bb8d9e5df1f995ebd87b262de6045e4734b89`

Nextflow stages this bundle only into `HAPPY_LEGACY`. `HAPPY_RUST` receives
`ref.fa` and never stages, opens, or extracts the bundle.
