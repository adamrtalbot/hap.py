// verification/modules/report.nf
//
// Collects every per-case status.json produced by DIFF_OUTPUTS and
// renders a single markdown summary at ${outdir}/report.md. Parity
// gating (failing the run on any FAIL row) is enforced by the nf-test
// suite under verification/tests/, so this process always succeeds.

process REPORT {
    tag 'aggregate'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path status_files, stageAs: 'status???.json'
    path verifier_bin, stageAs: 'verify-fixtures'

    output:
    path 'report.md', emit: markdown
    path 'report.csv', emit: csv

    script:
    """
    ./verify-fixtures aggregate-report \\
        --image '${params.legacy_image}' \\
        --hap-bin '${params.hap_bin}' \\
        --markdown report.md \\
        --csv report.csv \\
        status???.json
    """
}
