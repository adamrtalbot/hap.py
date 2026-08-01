# Legacy oracle provenance

Verified on 2026-07-31 with the exact fixture inputs below. Two independent
runs produced identical normalized outputs.

- Image tag: `community.wave.seqera.io/library/hap.py_rtg-tools:3d21f155636b42e0`
- Image digest: `sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6`
- Legacy package: `hap.py 0.3.15 py27hcb73b3d_0`
- Package archive SHA-256: `de3c93f3f8d6fc7c1e1d443562b127630bb888451a6598680204cf1707ad4c1e`
- `/opt/conda/bin/hap.py` SHA-256: `c2be55ae76d6f8e8efc91a701ddfda56fefa23fa1c511afd552ebe696cb38ff8`
- Upstream fixture revision: `84011695b2ff2406c16a335106db6831fb67fdfe`

The four input hashes match `example/haploid/{truth.vcf,query.vcf,test.fa,test.fa.fai}`
at the upstream fixture revision. The Conda package does not embed a hap.py Git
revision (`Haplo/version.py` has an empty `__version__`), so the immutable image,
package archive, and executable hashes above are the legacy implementation pin.

Exact oracle invocation (with the repository mounted at the same absolute path):

```sh
hap.py /Users/adam.talbot/hap.py/tests/fixtures/example-haploid/truth.vcf /Users/adam.talbot/hap.py/tests/fixtures/example-haploid/query.vcf -r /Users/adam.talbot/hap.py/tests/fixtures/example-haploid/ref.fa -o oracle --force-interactive -V -X
```

## SHA-256

```text
7652f6af2fd487d6ff8937be2b2c9229810406c6a84d35d0f7384e37f0feb3e1  truth.vcf
7c1db7cc8f53d7a04828953f963792b354b6a23c9d322accf6ee5551a145a2c1  query.vcf
d35561b3b8372b82ac86f9574a4176e3d5a45c30543b4c9c10284cb3915e8e92  ref.fa
be140264ae009a2f77dfe8bf1e49357f9cfae7deebf4ed11f34d14a9b65c4054  ref.fa.fai
a8566e40e57cbeb08886cb583a13e0ae4a919593fad8b730ff98ee60db8c1f95  expected/summary.csv
85115885ae9e5ce715c21df36cd6327fbf8477c6b63667342971aaa51c52ae6f  expected/extended.csv
573f936d623b482ef947ef4d212891656973187a458487c6c6c94cc162cbd9ac  expected/vcf.normalized.vcf
```

`vcf.normalized.vcf` contains only non-header VCF records, matching the fixture
harness normalization.
