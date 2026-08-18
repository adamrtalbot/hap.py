# Keep the claim on the command line

The validated surface is the command line. No caller is inside the drop-in claim,
and no caller-specific behaviour is claimed. nf-core/variantbenchmarking 1.5.0 is
a confirmation run executed after the port is complete, to demonstrate hap-rs in a
real pipeline and surface edge cases the harness missed. It gates nothing and
settles no value.

## The caller boundary is empty

`0002-drop-in-claim-against-one-unpatched-legacy-image.md` ranks admitted caller
behaviour as a scope filter that can never settle a value. That rule stands with no
instances: there is no admitted-caller list to keep, and no second entry to
sign off. Every caller, whether a person at a shell, a Nextflow process, or any
other wrapper, gets the covered invocation surface from
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` and nothing more.

Nextflow issues command lines. It adds no interface of its own.

## The confirmation run

Pinned to the upstream tag `1.5.0`. The pin will move as the project progresses;
moving it is an ordinary edit, because nothing is gated on it.

Measured at that tag, the pipeline issues exactly one `hap.py` invocation and one
`som.py` invocation, and reaches them only through `params`:
`small_benchmark/main.nf` runs happy when `method` contains `happy` and
`analysis == germline`, and sompy when `method` contains `sompy` and
`analysis == somatic`.

| profile | genome | analysis | variant_type | method | invokes |
|---|---|---|---|---|---|
| `test_happy` | GRCh38 | germline | small | happy | happy |
| `test_ga4gh` | GRCh38 | germline | small | happy | happy, with `--engine=vcfeval --leftshift` |
| `test_full` | GRCh37 | germline | small | happy,rtgtools | happy |
| `test_full_somatic` | GRCh38 | somatic | indel | sompy,rtgtools | sompy |
| `test` | GRCh37 | germline | structural | svanalyzer | neither |

`test_full` and `test_full_somatic` are the real-context pair. No parameter sweep
belongs here; the coverage model owns which invocations earn validation cases.

What the pipeline depends on, measured at the tag: nine mandatory and two optional
output globs for happy, two mandatory and one optional for sompy; `bcftools view
-s TRUTH` and `-s QUERY` on the annotated VCF; four `bcftools filter -i
'FORMAT/BD="TP"|"FN"|"FP"'` expressions; a local module that splits
`features.csv` by TP and FP. It never scrapes a version, because it cannot: the
module hardcodes `val('0.3.15')`.

Two quirks at the pin. `happy/happy` maps its regions input to `-f`, which is
`--false-positives`, so a run supplying both regions and false positives emits
`-f` twice and argparse takes the last; `happy/sompy` uses `-R` correctly. The
`-f` spelling is what every parity observation to date exercised, so the existing
evidence covers the pinned invocation. Upstream fixed it to `-R` after 1.5.0,
which is the first reason the pin will move.

`test_ga4gh` overrides the container to an RTG-bearing image because the
biocontainer lacks `rtg`. hap-rs spawns no external process, so that override
becomes dead config.

## Absorbing the command-name change

Callers edit their call sites. Nothing ships to keep `hap.py` working.

Measured: no integration has absorbed the change yet. `testing/rust-parity/`
in the caller checkout shims it twice, once by baking `/usr/local/bin/hap.py`
and `som.py` wrappers into the image, which all sixteen local
`hap-rs-variantbenchmarking` images carry, and once by prepending a `native-bin`
directory holding the same two wrappers. Both are test scaffolding. A run through
a wrapper exercises the wrapper, so the confirmation run swaps in a local module
that calls `hap germline` and `hap somatic` directly. Running the legacy and
hap-rs modules side by side in one execution is preferred, because Nextflow's
concurrency then does the pairing.

That module swap is the caller-side migration diff. It answers whether a real
caller can absorb the change, in two lines you can read.

## The reference is an argument

The reference must be supplied on the command line. `HGREF` and `HG19` are not
supported, and the covered invocation surface excludes an environment-supplied
reference as an input form. It is not an exemption and the register stays at two
entries.

Measurement makes it a quirk to drop. All five Python tools print `set the environment variable 'HGREF' or
'HG19'` on every invocation including `--help`, yet `/opt/hap.py-data/` is absent
from the reference environment, so nothing sits behind the advertised default.
Legacy checks only that the named path exists. The pipeline always passes
`--reference`, so removing the fallback costs the one real caller nothing.

hap-rs currently implements the fallback at
`rust/src/application/preprocess/options.rs:173`,
`rust/src/application/compare/output.rs:634`, and
`rust/src/application/ftx/mod.rs:39`. It is removed.

GRCh37, hg19, and b37 stay supported by supplying the FASTA and its companions as
arguments. No reference data ships in any deployment bundle.

## Scratch

`TMPDIR` is not consulted. Where legacy offers a scratch option, hap-rs copies it:
measured, `hap.py` and `som.py` carry `--scratch-prefix` and `--keep-scratch`, while
`pre.py`, `qfy.py`, and `ftx.py` carry none. Nextflow passes `--scratch-prefix .` to
land scratch in the task working directory, where the task's disk accounting and
cleanup already reach.

Three commands therefore have no command-line way to place scratch, so "those three
create no scratch" is a requirement on hap-rs rather than a description of it. hap-rs
does not meet it today. It derives scratch from `TMPDIR` at
`rust/src/output.rs:638`, `rust/src/adapters/vcf.rs:923`, `rust/src/engines/roc.rs:47`
and `rust/src/application/roc_publication.rs:73`, and the last two sit on paths that
`germline` and `quantify` reach.

## The programmatic surface

Exit status is the only programmatic contract, which
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` already established
by putting message text and stream choice outside the claim. Nothing parses
standard output.

`--version` prints the version and exits 0 on every subcommand as ordinary tool
behaviour. No comparison is owed, and none is possible: measured, `hap.py
--version` and `vcfcheck --version` exit 0 printing an empty version string,
`som.py` and `ftx.py` do not advertise the option, and `pre.py --version` and
`qfy.py --version` exit 2 because the flag parses and the required positionals are
then missing. That last case needs no exemption, since a missing required argument
is a malformed invocation and already outside the surface. This also disposes of
the `validate` version wording in
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md`, which the same
measurement contradicts.

## No container

This repository ships no image and holds no recipe; installation is
`cargo install hap-rs`. Distribution goes through cargo and bioconda, from which
a Seqera Container follows. Which platforms hap-rs supports, and how they are
pinned, belongs to
[Define supported execution environments and the reproducibility boundary](https://github.com/adamrtalbot/hap.py/issues/29).

The confirmation run therefore compares a host hap-rs binary against
containerized legacy, which is what the existing harness config already does with
`container = null`. One observation is recorded for whoever builds an image later,
with no commitment attached: the reference image measures Entrypoint null with Cmd
`["/bin/bash"]`, and every measured caller runs a shell script rather than the bare
executable.

## Output discovery

Callers discover outputs by glob, and a glob matching nothing fails the task.
There is no closed-set rule: one VCF is produced, and if more than one ever
appears that is an issue to file and fix, not a contract to write. Avoiding
collisions between a tool's outputs and a pipeline's own patterns is the pipeline
developer's responsibility.

## Considered options

Admitting the pipeline as a gated caller was rejected. It would make release
depend on remote compute and public-data availability, and
`0001-separate-discovery-workflows-from-regressions.md` already separates a
full-scale workflow from the regressions that gate.

Admitting the two nf-core modules rather than the pipeline was rejected: the modules
are shared, so the boundary would reach nf-core pipelines nobody has measured.

Shipping legacy-named wrappers, in the image or anywhere else, was rejected;
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` bars them and callers
edit their call sites instead.

Keeping `HGREF` and `HG19` because legacy advertises them on every run was
rejected. The advertised default path is absent from the reference environment and
the one real caller always passes `--reference`.

Promoting `--version` or `validate --output-json` to a caller contract was
rejected as a dependency legacy never offered, and so one that no comparison can
validate.

## Consequences

The exemption register stays at two entries. Excluding an environment-supplied
reference narrows the covered surface rather than adding to the register, and
`hap validate --help` returning 0 stands as a deliberate divergence from a legacy
bug that nothing depends on.

Three hap-rs changes follow and belong outside this map: remove the `HGREF` and
`HG19` fallback, stop deriving scratch from `TMPDIR` in favour of
`--scratch-prefix`, and retire the wrapper scripts in the caller checkout's
parity harness.

[Define the exact hap-rs replacement compatibility claim](https://github.com/adamrtalbot/hap.py/issues/21)
inherits an empty caller boundary, so the claim it states is a command-line claim.
[Inventory product artifacts and observable side effects](https://github.com/adamrtalbot/hap.py/issues/22)
inventories artifacts with no completeness rule owed from here, and picks up
scratch placement as a filesystem side effect.
[Set release stop conditions and exception governance](https://github.com/adamrtalbot/hap.py/issues/35)
decides whether the confirmation run is named as a stop condition; nothing here
requires it.

Two residual risks stay open. A second artifact matching a caller's glob would
change that caller's output arity, and per-artifact comparison would not catch it;
file an issue if it happens. The confirmation run also lands after the port is
complete, so it finds its edge cases later than the harness would.
