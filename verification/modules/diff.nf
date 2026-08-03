// verification/modules/diff.nf
//
// Pairs legacy and rust output sets for a single case and dispatches to
// the file-type-aware comparator in the staged Rust verifier. One task
// per (case, sample_id) pair runs in parallel.
//
// The input tuple is:
//   meta         - map with at least `id` and `case`
//   case_name    - string used in status output
//   prefix       - filename prefix passed to the comparator (e.g. 'result')
//   expected     - independently observed exact artifact names for this row
//   legacy_files - list of path objects from the *_LEGACY process
//   rust_files   - list of path objects from the *_RUST process
//
// The comparator writes `diff.log` (human payload) and `status.json`
// (machine-readable summary, one per case-sample).

process DIFF_OUTPUTS {
    tag { "${case_name}:${meta.id}" }
    publishDir { "${params.outdir}/${case_name}/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), val(case_name), val(prefix), val(expected_artifacts), path(legacy_files, stageAs: 'legacy/*'), path(rust_files, stageAs: 'rust/*')
    path verifier_bin, stageAs: 'verify-fixtures'

    output:
    path 'diff.log', emit: log
    tuple val(meta), val(case_name), path('status.json'), emit: status

    script:
    def expected_args = expected_artifacts
        .collect { artifact -> "--expected-artifact '${artifact}'" }
        .join(' ')
    """
    ./verify-fixtures compare-outputs \\
        --legacy-dir legacy \\
        --rust-dir rust \\
        --prefix ${prefix} \\
        ${expected_args} \\
        --case ${case_name} \\
        --sample ${meta.id} \\
        --report diff.log \\
        --status status.json
    """
}
