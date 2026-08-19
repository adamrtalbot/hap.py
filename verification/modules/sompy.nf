// verification/modules/sompy.nf
//
// Mirrors nf-core/modules/happy/sompy/main.nf.

process SOMPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/sompy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(reference), path(reference_fai), path(fp_bed), val(feature_table), val(args), path(bams), path(bam_indexes)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    def bam_args = bams.collect { bam -> "--bam ${bam}" }.join(' ')
    def reference_arg = reference ? "--reference ${reference}" : ''
    """
    som.py ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
        -o result \\
        ${reference_arg} \\
        --false-positives ${fp_bed} \\
        --feature-table ${feature_table}
    """
}

process SOMPY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/sompy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(reference), path(reference_fai), path(fp_bed), val(feature_table), val(args), path(bams), path(bam_indexes)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    def bam_args = bams.collect { bam -> "--bam ${bam}" }.join(' ')
    def reference_arg = reference ? "--reference ${reference}" : ''
    """
    hap somatic ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
        -o result \\
        ${reference_arg} \\
        --false-positives ${fp_bed} \\
        --feature-table ${feature_table}
    """
}
