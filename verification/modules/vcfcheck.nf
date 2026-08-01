// verification/modules/vcfcheck.nf
//
// Runs the legacy vcfcheck executable and the Rust validate compatibility
// command against the same upstream-pinned VCF.

process VCFCHECK_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/vcfcheck/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    vcfcheck ${args} ${input_vcf} \
        --output-file result.json
    """
}

process VCFCHECK_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/vcfcheck/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), val(args)
    path hap_bin

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    ./${hap_bin} validate ${args} ${input_vcf} \
        --output-json result.json
    """
}
