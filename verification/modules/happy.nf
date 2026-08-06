// verification/modules/happy.nf
//
// Mirrors nf-core/modules/happy/happy/main.nf for the germline comparison.
// Two processes with identical input shape:
//   * HAPPY_LEGACY uses `hap.py` inside the pinned Wave container.
//   * HAPPY_RUST   uses `hap germline` from the locally-built binary.

process HAPPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/happy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(truth_indexes, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(query_indexes, stageAs: 'query/*'), path(reference), path(reference_fai), path(fp_bed), path(fp_indexes), path(stratification_files, stageAs: 'stratification/*'), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    """
    hap.py ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result

    """
}

process HAPPY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/happy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(truth_indexes, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(query_indexes, stageAs: 'query/*'), path(reference), path(reference_fai), path(fp_bed), path(fp_indexes), path(stratification_files, stageAs: 'stratification/*'), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    """
    hap germline ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result
    """
}
