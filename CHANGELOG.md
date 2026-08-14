# Changelog

## Unreleased

### Compatibility

- Restore legacy string typing for integer-looking `Subset.Size` values in
  germline location metrics JSON. This corrects hap-rs output that previously
  used JSON integers and does not change comparison normalization.
- Preserve legacy classified-row ordering and the narrow shared-insertion
  `hapfail` classification used by existing HAPPY artifacts.
- Apply the legacy preprocessing rule that removes FTX records when a sample
  calls the final `<NON_REF>` allele; uncalled and non-final shapes are retained.
