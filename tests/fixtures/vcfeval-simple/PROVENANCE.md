# vcfeval-simple oracle provenance

The expected summary and VCF record were generated with the pinned legacy
image:

`community.wave.seqera.io/library/hap.py_rtg-tools@sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6`

The source command was:

```text
hap.py truth.vcf query.vcf -r ref.fa --engine vcfeval \
  --engine-vcfeval-path rtg --force-interactive --threads 1 \
  --no-json --no-roc -o oracle
```

`rtg-output.vcf` is the stable GA4GH handoff shape supplied by the pinned RTG
vcfeval run. The integration test stubs only the external RTG executable with
this saved handoff; Rust preprocessing, dispatch, GA4GH quantification, and
report generation remain live.
