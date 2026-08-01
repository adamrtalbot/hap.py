// verification/modules/ftxpy.nf
//
// Mirrors nf-core/modules/happy/ftxpy/main.nf: feature-table extraction
// for a single VCF against a reference.

process FTXPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/ftxpy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(reference), path(reference_fai), path(bam), path(bam_index), val(has_bam), val(feature_table), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    def bam_args = has_bam ? "--bam ${bam}" : ''
    """
    ftx.py ${args} ${bam_args} ${input_vcf} \\
        -o result \\
        --reference ${reference} \\
        --feature-table ${feature_table}
    """
}

process FTXPY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/ftxpy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(reference), path(reference_fai), path(bam), path(bam_index), val(has_bam), val(feature_table), val(args)
    path hap_bin

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    def bam_args = has_bam ? "--bam ${bam}" : ''
    """
    ./${hap_bin} ftx ${args} ${bam_args} ${input_vcf} \\
        -o result \\
        --reference ${reference} \\
        --feature-table ${feature_table}
    """
}
