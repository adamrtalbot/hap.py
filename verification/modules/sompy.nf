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
    """
    mkdir legacy_python_compat
    printf '%s\\n' \\
        'import pandas' \\
        '_set_option = pandas.set_option' \\
        'def set_option(*args, **kwargs):' \\
        '    if args and args[0] == "display.height":' \\
        '        return None' \\
        '    return _set_option(*args, **kwargs)' \\
        'pandas.set_option = set_option' \\
        > legacy_python_compat/sitecustomize.py

    PYTHONPATH=legacy_python_compat som.py ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
        -o result \\
        --reference ${reference} \\
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
    """
    hap somatic ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
        -o result \\
        --reference ${reference} \\
        --false-positives ${fp_bed} \\
        --feature-table ${feature_table}
    """
}
