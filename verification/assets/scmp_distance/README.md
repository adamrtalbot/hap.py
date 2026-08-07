# vcfeval reference bundle

`ref.sdf.tar.gz` is a verification-only RTG SDF generated from `ref.fa` with
RTG Tools 3.12.1 (`rtg format -o ref.sdf ref.fa`) and then archived with the
`ref.sdf` directory as the archive root.

- RTG source tag: `3.12.1`
- RTG source commit: `32d4c2d2d340cb0288f68ddd008d29564bcfef13`
- Bundle SHA-256: `e2b82260576fffaf4b7e576d14418c41329bf34a37ccaccfffdd2dee81b14935`

Nextflow stages this bundle only into `HAPPY_LEGACY`. `HAPPY_RUST` receives
`ref.fa` and never stages, opens, or extracts the bundle.
