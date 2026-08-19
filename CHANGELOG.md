# Changelog

## Unreleased

### Compatibility

- Restore legacy string typing for integer-looking `Subset.Size` values in
  germline location metrics JSON. This corrects hap-rs output that previously
  used JSON integers and does not change comparison normalization.
- Preserve legacy classified-row ordering and the narrow shared-insertion
  `hapfail` classification used by existing HAPPY artifacts.
- Apply the legacy preprocessing rule that removes FTX records when a sample
  calls the final `<NON_REF>` allele; uncalled and non-final shapes are retained.
- State the drop-in claim as a versioned pair: hap-rs 1.0.0 against one pinned
  legacy container digest and its conda lock, covering exit status and artifacts
  for every invocation the pinned parsers accept, observed on linux/amd64. No 0.x
  release carries the claim.
- Keep `--engine-vcfeval-path` and `--engine-vcfeval-template` supported. The
  scheduled removal in 1.0.0 is withdrawn; both options remain accepted and
  ignored, because the native engine reads the FASTA supplied with `--reference`
  and needs no SDF bundle.
- Return a non-zero status when `pre` or `quantify` receives an unknown option,
  rather than deferring the change to 1.0.0. Malformed invocations sit outside the
  covered invocation surface, so this narrows no compatibility claim.
- Report `germline` usage errors through the argument parser as well, replacing
  exit 1 and the help page hap-rs wrote to standard output. All six subcommands
  and their aliases now share one usage-error path: the parser's status, its
  message on standard error, and nothing on standard output.
- Take the reference from `--reference` alone. The `HG19` and `HGREF`
  environment variables and the `/opt/hap.py-data/hg19.fa` install path are no
  longer consulted, and `somatic` now requires `--reference` at parse time.
  This narrows the covered invocation surface: an environment-supplied
  reference was already outside it.
