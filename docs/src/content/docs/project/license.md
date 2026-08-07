---
title: License & Attribution
description: hap-rs licensing, lineage, and third-party notices.
---

## License

The Simplified BSD License covers `hap-rs`:

- Copyright © 2010–2015 Illumina, Inc.
- [`LICENSE.txt`](https://github.com/adamrtalbot/hap.py/blob/master/LICENSE.txt)
  lists the conditions for source and binary redistribution.
- The authors provide the software without warranty.

## Project lineage

`hap-rs` reimplements Illumina's hap.py toolkit in Rust. It keeps the benchmark
semantics and interfaces covered by the compatibility tests.

## Third-party components

Each Rust dependency keeps its own license. `Cargo.lock` records the version in
use.

RTG Tools behavior informed the Rust vcfeval implementation. The verification
lane runs RTG Tools 3.12.1 under its BSD 2-Clause license. See
[`THIRD_PARTY_LICENSES/rtg-tools.txt`](https://github.com/adamrtalbot/hap.py/blob/master/THIRD_PARTY_LICENSES/rtg-tools.txt).

License notices for imported fixtures sit beside the fixture or under
`THIRD_PARTY_LICENSES/`.
