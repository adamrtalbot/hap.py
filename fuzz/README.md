# Parser fuzzing

The five targets exercise the native VCF, BCF, FASTA, BED, and location
parsers. Named corpus entries are tracked regression inputs and bounded CI
smoke starts from them. The BCF smoke corpus also includes the repository's
valid native BCF fixture, so the decoder is exercised beyond its magic header.

Run all targets with the same bounded gate used by CI:

```bash
tests/fuzz_smoke.sh 1000
```

Each target is bounded by both the requested run count and a 30-second wall
clock limit (`FUZZ_MAX_TOTAL_TIME` overrides the latter). CI uploads
`fuzz/artifacts` whenever a target fails.

For longer local campaigns, select a target and omit `-runs`:

```bash
cargo fuzz run vcf
```

When libFuzzer writes a reproducer under `fuzz/artifacts/<target>`:

1. Re-run it with `cargo fuzz run <target> <artifact>` and fix the defect.
2. Minimize it with `cargo fuzz tmin <target> <artifact>`, then copy the
   result into `fuzz/corpus/<target>` under a descriptive `seed-*.ext` name.
3. Run `tests/fuzz_smoke.sh 1` to prove every named seed reaches its native
   decoder without recreating the failure, and commit the seed with the fix.

The hash-named working corpus and `fuzz/artifacts` stay ignored. The current
named regressions retain truncated VCF fields, FASTA sequence before a header,
reversed BED coordinates, overflowing locations, truncated BCF header/record
cases, and a mutated BGZF block discovered by the bounded smoke.
Only the byte-exact valid BGZF fixture is decompressed by the harness; arbitrary
BCF mutations exercise the native uncompressed decoder without exposing the
dependency's unchecked malformed-block path. CI uploads the ignored
`fuzz/artifacts` tree on failure so a new reproducer is not lost.
