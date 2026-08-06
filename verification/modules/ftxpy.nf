// verification/modules/ftxpy.nf
//
// Mirrors nf-core/modules/happy/ftxpy/main.nf: feature-table extraction
// for a single VCF against a reference.

process FTXPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/ftxpy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(reference), path(reference_fai), path(bams), path(bam_indexes), path(region_files), val(feature_table), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    def bam_args = bams.collect { bam -> "--bam ${bam}" }.join(' ')
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
    tuple val(meta), path(input_vcf), path(reference), path(reference_fai), path(bams), path(bam_indexes), path(region_files), val(feature_table), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    def bam_args = bams.collect { bam -> "--bam ${bam}" }.join(' ')
    """
    hap ftx ${args} ${bam_args} ${input_vcf} \\
        -o result \\
        --reference ${reference} \\
        --feature-table ${feature_table}
    """
}
