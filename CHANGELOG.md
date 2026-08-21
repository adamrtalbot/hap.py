# Changelog

## Unreleased

### Compatibility

- Render every metric CSV float cell as the shortest round-trippable decimal,
  matching the pandas 0.24.2 `to_csv` build nf-core/variantbenchmarking runs:
  germline `summary.csv`, `extended.csv`, `*.roc.*.csv.gz`, somatic
  `*.stats.csv`, and the `ftx` feature tables. Values are unchanged; this
  replaces the earlier Python-2 twelve-significant-digit text.
- Type a germline location `Subset.Size` metrics-JSON column as `int64` when
  every cell is a bare integer, matching pandas 0.24.2 dtype inference; a column
  carrying a `%.6f` BED sum such as `141.000000` stays string-typed.
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
  longer consulted. This narrows the covered invocation surface: an
  environment-supplied reference was already outside it. `germline` and `pre`
  now report a missing reference instead of resolving one; `somatic` and `ftx`
  demand one only where they open it, so a somatic allele comparison with an
  FP BED still runs with no reference at all.
