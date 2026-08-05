// verification/modules/qfy.nf
//
// Runs qfy.py and the Rust quantify compatibility command against the same
// xcmp-annotated VCF. QFY_ANNOTATE creates that shared input once from a
// bounded, upstream-pinned chr21 fixture using the reference comparator.

process QFY_ANNOTATE {
    tag { "${meta.id}" }
    container params.legacy_image

    input:
    tuple val(meta), path(truth_vcf), path(truth_tbi), path(query_vcf), path(query_tbi), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args)

    output:
    tuple val(meta), path('seed.vcf.gz'), path('seed.vcf.gz.tbi'), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args), emit: annotated

    script:
    """
    hap.py ${truth_vcf} ${query_vcf} \
        --reference ${reference} \
        --threads 1 \
        --engine xcmp \
        -l chr21:15000000-20000000 \
        --false-positives ${fp_bed} \
        -V \
        --no-json \
        --no-roc \
        -o seed
    """
}

process QFY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/qfy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_tbi), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    """
    qfy.py ${args} ${input_vcf} \
        --reference ${reference} \
        --false-positives ${fp_bed} \
        --threads ${task.cpus ?: 1} \
        --report-prefix result

    """
}

process QFY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/qfy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_tbi), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    """
    hap quantify ${args} ${input_vcf} \
        --reference ${reference} \
        --false-positives ${fp_bed} \
        --report-prefix result
    """
}
