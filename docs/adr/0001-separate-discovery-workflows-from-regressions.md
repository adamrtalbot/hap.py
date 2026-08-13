# Separate full-data discovery from authoritative regressions

Full-scale public datasets are run from one shared public-data samplesheet as
Nextflow discovery workflows locally or on remote compute, not as part of the
ordinary nf-test gate. Each truth/query comparison remains one intact row;
HG001 discrepancies do not split the HG001 row. Their command and expected
comparison evidence are documented.

Every diagnosed defect is also reproduced by an extremely small example case
in the ordinary verification samplesheets and exercised by nf-test. Rows model
example scenarios rather than issue records, so one genuinely minimal scenario
may cover coupled defects. The public workflow confirms the combined fix on
realistic data, while nf-test provides affordable permanent regression proof.
