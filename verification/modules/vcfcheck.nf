// verification/modules/vcfcheck.nf
//
// Runs the legacy vcfcheck executable and the Rust validate compatibility
// command against the same upstream-pinned VCF.

process VCFCHECK_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/vcfcheck/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_indexes), val(input_mode), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    input_arg = input_mode == 'option' ? "--input-file ${input_vcf}" : input_vcf
    """
    vcfcheck ${args} ${input_arg} \
        --output-file result.json
    """
}

process VCFCHECK_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/vcfcheck/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_indexes), val(input_mode), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    input_arg = input_mode == 'option' ? "--input-file ${input_vcf}" : input_vcf
    """
    hap validate ${args} ${input_arg} \
        --output-json result.json
    """
}
