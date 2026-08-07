---
title: Installation
description: Install hap-rs through Cargo.
---

## Requirements

- A supported operating system with a Rust toolchain
- Rust 1.89 or newer

The `hap` executable handles comparison without a Python, Java, or RTG runtime.

## Install with Cargo

```bash
cargo install hap-rs
```

:::caution[Cargo status]
`hap-rs` has no crates.io release yet. The first release will use this command.
Contributors can build the source with the
[contributing guide](../../project/contributing/).
:::

Check the installation:

```bash
hap --version
hap --help
```

## Upgrade

After the first crates.io release, install an update with:

```bash
cargo install hap-rs --force
```

`hap` writes reports to the path passed on the command line. Commands that
inspect reference alleles need a FASTA and its `.fai` index.
