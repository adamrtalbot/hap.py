# Close the map and measure what is left

The wayfinder map, [Wayfinder: decide whether hap-rs can replace the pinned
Legacy implementation](https://github.com/adamrtalbot/hap.py/issues/15), closes
at seven decisions, recorded in ADRs 0002 through 0008. Its destination asked
for "a ratified, versioned release-validation programme and verdict rule", and
a programme always has one more governance question, so the map generated
policy tickets faster than it retired fog. The direction it existed to find is
found: the claim, the reference authority, the invocation surface, the caller
boundary, the observed set, and the equivalence contract are all settled.

Fourteen open decision tickets remain. Eight of them are already answered by
artifacts in this repository and close by citation. Three questions remain
genuinely open, and each is answered by running the gate rather than by
interview. The map hands off to
[spec: define sufficient validation for hap-rs replacement](https://github.com/adamrtalbot/hap.py/issues/37),
which is rewritten to stand on its own rather than to point back here.

## Closed by citation

| Question | Answer, and where it already lives |
|---|---|
| Supported execution environments | `verification/README.md`. Parity is observed on linux/amd64, against two pinned digests and their conda locks, with authoritative runs in CI on `ubuntu-24.04`. |
| Campaign evidence inventory | ADR 0007 defines campaign evidence and excludes it from comparison. The harness writes `comparison.json`, `verification.json`, the streamed artifact records, `.command.*`, and captured output. |
| Evidence retention and audit | None is retained, because none is stored. `results/` is gitignored, there are no nf-test snapshots, and legacy and hap-rs run together in one execution and are diffed live, as ADR 0004 records. |
| Campaign sequencing and rerun rules | CI is the protocol. A rerun is a fresh paired execution, because nothing is carried between runs. |
| Release stop conditions | The parity gate. A committed parity case with an unapproved difference fails, and a failure stops the run under `errorStrategy = 'finish'`. |
| Residual-risk reporting and verdict | The exemption register is the residual-risk register, and it holds two permanent entries. The verdict is the gate's result over the committed corpus. |
| Independent scientific validity | Not required. The drop-in claim is agreement with the pinned Legacy implementation, not independent correctness, and CONTEXT.md already holds the two apart. |
| Concurrency and failure modes | Live product issues, not policy: the global publication lock and scratch placement are tracked as ordinary defects. |

## Answered by measurement

Three questions change what the harness contains rather than what a document
says, so each is resolved by observing the existing 154-comparison matrix and
tabling the gaps, not by grilling:

- **Coverage matrix, including option interactions.** What the committed rows
  already exercise across the six commands, input formats, engines, regions,
  and filtering modes, and which cells are deliberately absent.
- **Representative corpus, including public-data admission.** Which named
  inputs stand, and the provenance each one carries.
- **Deterministic replay and resource thresholds.** Derived from repeated gate
  runs, expressed as what may vary between runs.

## The rule this ADR sets

Never grill a question whose answer is measurable. Measure it, record the
measurement in `verification/README.md`, then decide once. Four of the seven
ADRs carry a correction added after a later session measured a premise the
earlier one had assumed, and that repair is where the hours went: ADR 0002
twice, on a digest that cannot run all six tools and on the pandas build the
rebuild took; ADR 0003, on `vcfcheck` accepting `--version`; ADR 0005, on
legacy consulting `TMPDIR`; and ADR 0006, on exit status where legacy
succeeds. Each one is a decision taken before its measurement.

An ADR states a decision and the rule it creates. Measurements belong in
`verification/README.md`, because an ADR that carries them becomes the input to
the next session and every session's prompt grows.

## Consequences

Eight decision tickets close with a citation apiece. Three merged tickets
remain, each resolved by a gate run. No further wayfinder session runs on this
map; work continues on the main flow from the rewritten spec.
