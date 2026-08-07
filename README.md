<p align="center">
  <img src="docs/src/assets/hap-rs-mark.svg" alt="hap-rs" width="180">
</p>

<h1 align="center">hap-rs</h1>

<p align="center">
  Haplotype-aware variant benchmarking in one native executable.
</p>

<p align="center">
  <a href="https://adamrtalbot.github.io/hap.py/">Documentation</a> ·
  <a href="https://adamrtalbot.github.io/hap.py/getting-started/quick-start/">Quick start</a> ·
  <a href="https://adamrtalbot.github.io/hap.py/project/verification/">Verification</a> ·
  <a href="CONTRIBUTING.md">Contributing</a>
</p>

---

`hap-rs` compares variant callsets, normalizes VCF data, extracts somatic
features, calculates benchmark metrics, and validates inputs. It packages the
established hap.py workflows as the single `hap` executable, implemented in
Rust and verified against pinned legacy releases.

## Benefits

- **One executable:** run all six commands without a Python, C, C++, Java, or
  RTG runtime.
- **Rust comparison engines:** use `vcfeval` with a FASTA reference instead of
  an SDF bundle.
- **Compatible interfaces:** keep the command options and report formats that
  existing workflows use.
- **Common genomics formats:** read VCF, BGZF-compressed VCF, and BCF with
  Tabix or CSI indexes.
- **Measured parity:** 153 comparisons cover all six commands.

## Commands

| Command | Use it to |
|---|---|
| [`hap germline`](https://adamrtalbot.github.io/hap.py/tools/germline/) | Compare germline truth and query callsets using haplotypes or another engine. |
| [`hap somatic`](https://adamrtalbot.github.io/hap.py/tools/somatic/) | Benchmark somatic calls by allele identity and genomic context. |
| [`hap pre`](https://adamrtalbot.github.io/hap.py/tools/pre/) | Normalize, filter, and transform a VCF or BCF. |
| [`hap ftx`](https://adamrtalbot.github.io/hap.py/tools/ftx/) | Extract caller-specific somatic features into a table. |
| [`hap quantify`](https://adamrtalbot.github.io/hap.py/tools/quantify/) | Produce stratified metrics from an annotated comparison VCF. |
| [`hap validate`](https://adamrtalbot.github.io/hap.py/tools/validate/) | Check VCF structure, alleles, samples, and reference consistency. |

## Quick start

> [!NOTE]
> `hap-rs` has no crates.io release yet. The first release will use this
> command.

```bash
cargo install hap-rs
```

Compare a query callset with a truth set:

```bash
hap germline truth.vcf.gz query.vcf.gz \
  --reference reference.fa \
  --false-positives confident.bed \
  --report-prefix results/sample
```

`hap` writes `results/sample.summary.csv` and
`results/sample.extended.csv`. See the
[quick start](https://adamrtalbot.github.io/hap.py/getting-started/quick-start/)
for input requirements and output details.

## Documentation

The [documentation website](https://adamrtalbot.github.io/hap.py/) covers
installation, each command, file formats, verification, and development.

## License

The Simplified BSD License covers `hap-rs`. See
[`LICENSE.txt`](LICENSE.txt) and [`THIRD_PARTY_LICENSES/`](THIRD_PARTY_LICENSES/)
for details. The machine-checkable
[`manifest.json`](THIRD_PARTY_LICENSES/manifest.json) maps imported source and
fixture files to immutable upstream revisions and packaged notices.
