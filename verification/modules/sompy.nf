// verification/modules/sompy.nf
//
// Mirrors nf-core/modules/happy/sompy/main.nf.

process SOMPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/sompy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf), path(query_vcf), path(reference), path(reference_fai), path(fp_bed), val(feature_table), val(args), path(bam), path(bam_index), val(has_bam)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    bam_args = has_bam ? "--bam ${bam}" : ''
    """
    som.py ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
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
    tuple val(meta), path(truth_vcf), path(query_vcf), path(reference), path(reference_fai), path(fp_bed), val(feature_table), val(args), path(bam), path(bam_index), val(has_bam)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    bam_args = has_bam ? "--bam ${bam}" : ''
    """
    hap somatic ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
        -o result \\
        --reference ${reference} \\
        --false-positives ${fp_bed} \\
        --feature-table ${feature_table}
    """
}
