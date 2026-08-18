# Bound the invocation surface to the pinned parsers

The supported invocation surface is `hap <subcommand>` for the six commands, and
it covers exactly the options, positionals, and input forms that the pinned
reference image's parsers accept. For a covered invocation the drop-in claim
reaches the exit status and the produced artifacts. Standard output and standard
error content, message wording, and which stream carries a message are outside
it. Malformed invocations are outside it entirely.

## Entry points

Legacy is six executables. hap-rs is one binary with six subcommands, no
`argv[0]` dispatch, and no legacy-named wrapper. Callers replace `hap.py` with
`hap germline`, `som.py` with `hap somatic`, `pre.py` with `hap pre`, `qfy.py`
with `hap quantify`, `ftx.py` with `hap ftx`, and `vcfcheck` with
`hap validate`. The aliases `compare`, `preprocess`, `prepy`, `ftxpy`, `qfy`,
and `vcfcheck` stay. That command-name change is a documented interface change
carrying a version, not an exemption, so it does not enter the register.

## Covered option set

Measured against the pinned image, the long-option counts are hap.py 62,
pre.py 26, qfy.py 27, som.py 41, ftx.py 12, vcfcheck 11. hap-rs accepts every one
of them. Short options match as well: every legacy short parses, with som.py's
`-FN` handled by the alias normalizer rather than advertised in `--help`. For five
of the six commands the sets are equal. Only `validate` differs: it lacks
`--version`, and it adds `--reference`, `--errors-bed`, `--regions`, `--targets`,
and `--output-json`, keeping `--output-file` as an alias and carrying `-r`, `-e`,
`-R`, and `-T` alongside.

`--force-interactive` looks added on three commands because legacy registers it
only when SGE is present, which the image lacks. Legacy tolerates it anyway:
`hap.py`, `pre.py`, and `qfy.py` each accept it and carry on to their input
checks, raising no unknown-argument error. Not a divergence.

Additions are declared extensions rather than claimed behavior. They must not
change the result of any legacy invocation, and no legacy comparison is owed for
them.

## Parser behavior that `--help` does not show

argparse expands unique long-option prefixes and consumes the next token as the
value, so `pre.py --refer in out` fails with too few arguments rather than
reading a reference. hap-rs emulates this for the five Python-derived commands
and deliberately not for `validate`, whose Boost parser has no such feature.
Boost accepts `on|off`, `yes|no`, `1|0`, and `true|false` for boolean options
while hap-rs accepts `true|false` only, so the other spellings are unsupported
even though legacy takes them, including the `--check-bcf-errors 1` form that
pre.py itself emits internally.

No legacy tool reads a VCF from standard input. `-` reaches htslib or bcftools
unchanged and fails in all four Python tools and in vcfcheck, so `-` is not a
supported input form and hap-rs owes nothing there.

## Malformed invocations

Unknown options, missing required arguments, extra positionals, and
non-canonical boolean spellings sit outside the surface. hap-rs documents what
it does in those cases, but none of it is a compatibility promise.

## How the surface is recorded

The six `--help` captures from the pinned image become a governed artifact under
`verification/`, alongside a written statement of the parser behavior described
above, which `--help` does not reveal. They regenerate from the image by script,
as `happy-0.3.15.conda-lock.txt` already does, so the recorded surface cannot
drift from the authority it describes.

## Considered options

Shipping six `argv[0]` shims would let existing scripts run untouched, which is
the most literal reading of drop-in. It was rejected because it moves six more
entry points into the validated surface and roughly doubles the invocation layer
of the coverage model, to buy an edit that call sites can make once.

Treating hap-rs additions as fully in scope was rejected because five `validate`
options have no legacy counterpart, so a comparison matrix has nothing to
compare them against. Restricting the surface to the intersection was rejected
in the other direction: it would make shipped options officially unsupported.

Putting stream discipline in the claim was tempting, because legacy hap.py
writes its whole 11 KB help to standard output at exit 1, and som.py and ftx.py
write BCFTOOLS error prose to standard output. It went out together with
byte-exact output comparison: reproducing argparse and Boost message text and
Python tracebacks, including line numbers inside `/opt/conda/bin/hap.py`, buys
diagnostics text that no result depends on.

Enumerating all 179 options by hand in the reference docs was rejected as a
second source of truth that can diverge from the image without anyone noticing.

Matching `vcfcheck --help`, which exits 1, was rejected. A non-zero exit for a
successful help request reads as a bug to anyone who has not read this decision.

## Consequences

The exemption register changes content while staying at two entries. Unknown
`pre` and `quantify` options returning zero is retired from it, because
malformed invocations are no longer claimed; the behavior and its 1.0.0
transition are now ordinary hap-rs documentation. `hap validate --help`
returning 0 where `vcfcheck` returns 1 joins the register in its place. The
register's length stays checkable at release, and this amends the register
recorded in `0002-drop-in-claim-against-one-unpatched-legacy-image.md`.

`hap validate --version` has to exist and return 0. Because standard output
content is outside the claim, this decision says nothing about what it prints,
and the version number itself is a provenance field either way.

The rows in the compatibility policy covering `germline` invalid-argument exits
and missing `pre` or `quantify` arguments describe hap-rs behavior from here on,
not agreement with legacy.

Whether any real caller can absorb the command-name change, or depends on
`vcfcheck --help` exiting 1, belongs to the downstream caller and integration
boundary. Which covered invocations earn validation cases belongs to the
coverage model, and the recorded surface is the denominator it needs. Report
families, index companions, and other observable side effects belong to the
artifact inventory.

One residual risk is unmeasured. The recorded surface is captured from `--help`,
so an option a legacy parser accepts without advertising would be missed.
Conditional registration is the known instance of that shape, and
`--force-interactive` is the one case found.

## Correction: `vcfcheck` does accept `--version`

The covered option set above says `validate` "lacks `--version`". Measured in the
pinned image, `vcfcheck --version` exits 0 and prints `vcfcheck version ` with an
empty version number, so the option is there and the sentence is wrong. Two of the
Python tools also behave differently from what a reading of `--help` suggests:

| command | `--version` |
|---|---|
| `hap.py` | exit 0, prints `Hap.py ` with an empty version |
| `vcfcheck` | exit 0, prints `vcfcheck version ` with an empty version |
| `pre.py`, `qfy.py` | advertised in `--help`, but exit 2, because the flag parses and the required positionals are then missing |
| `som.py`, `ftx.py` | not advertised, exit 2 |

No legacy tool yields a usable version string, which is why
nf-core/variantbenchmarking hardcodes `val('0.3.15')` in both its happy modules
with the comment that the tool provides no version on the CLI.

The consequence stated above, that `hap validate --version` has to exist and
return 0, is unchanged, and nothing in the register moves.
[Define the downstream caller and integration boundary](https://github.com/adamrtalbot/hap.py/issues/19)
then settled that `--version` prints the version and exits 0 on every subcommand as
ordinary tool behaviour, with no comparison owed and none possible, so the
divergence this correction exposes is not a compatibility question at all. The
`pre.py`/`qfy.py` exit 2 needs no exemption either: a missing required argument is
a malformed invocation, which this decision already places outside the surface.
