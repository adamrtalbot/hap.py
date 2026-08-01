# Legacy oracle provenance

Verified on 2026-07-31 with the exact synthetic inputs below. Two independent
runs produced identical normalized outputs.

- Image tag: `community.wave.seqera.io/library/hap.py_rtg-tools:3d21f155636b42e0`
- Image digest: `sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6`
- Legacy package: `hap.py 0.3.15 py27hcb73b3d_0`
- Package archive SHA-256: `de3c93f3f8d6fc7c1e1d443562b127630bb888451a6598680204cf1707ad4c1e`
- `/opt/conda/bin/hap.py` SHA-256: `c2be55ae76d6f8e8efc91a701ddfda56fefa23fa1c511afd552ebe696cb38ff8`
- Repository revision at verification: `9a3ac919b8ee67bf734910756f85c01477f0604e`
- Upstream base revision: `84011695b2ff2406c16a335106db6831fb67fdfe`

This case is a local synthetic fixture and has no upstream source path; its exact
source is pinned by the input hashes below. The Conda package does not embed a
hap.py Git revision (`Haplo/version.py` has an empty `__version__`), so the
immutable image, package archive, and executable hashes above are the legacy
implementation pin.

Exact oracle invocation (with the repository mounted at the same absolute path):

```sh
hap.py /Users/adam.talbot/hap.py/tests/fixtures/synth-homopolymer-insertion/truth.vcf /Users/adam.talbot/hap.py/tests/fixtures/synth-homopolymer-insertion/query.vcf -r /Users/adam.talbot/hap.py/tests/fixtures/synth-homopolymer-insertion/ref.fa -o oracle --force-interactive -V -X
```

## SHA-256

```text
1a7c851527c4c32ac1eaef05e4e573e143326a9e1e5b6386f5469df55eb9cdca  truth.vcf
79757dd6a43515f9a1adbb0605057ab6239807a499ba41b29576f2f7eb7b40a6  query.vcf
819e6b81b9f0f6b9b23221c47912f3ddd771df596f8f475593b2be823cd1d21b  ref.fa
215c9f5370d7e6d9f2090f1ebbc67ef15c89de3030f2535a5e201d9b8c0ac3c0  ref.fa.fai
f676cfd8ba5f899cceca68a01e0e13dc9a839b350df20b55740eaecbf1aff269  expected/summary.csv
316504177c5dad1a408f2a8d20fc8cdd1280d1a5b6d6c679c4afad485910d73b  expected/extended.csv
cf2c6f90008ed197bf9891168b7d3018fb27cfd5a5d2e8eda01857a9a5357e5f  expected/vcf.normalized.vcf
```

`vcf.normalized.vcf` contains only non-header VCF records, matching the fixture
harness normalization.
