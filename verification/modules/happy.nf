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
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(truth_indexes, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(query_indexes, stageAs: 'query/*'), path(reference), path(reference_indexes), path(fp_bed), path(fp_indexes), path(stratification_files, stageAs: 'stratification/*'), path(reference_sdf, stageAs: 'vcfeval-sdf/*'), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    def sdfSetup = reference_sdf
        ? "mkdir vcfeval-template.sdf && tar --no-same-owner --no-same-permissions --exclude='._*' -xzf vcfeval-sdf/*.tar.gz -C vcfeval-template.sdf --strip-components=1 && export RTG_JAVA_OPTS=-Xint"
        : ':'
    def sdfArgs = reference_sdf ? '--engine-vcfeval-template vcfeval-template.sdf' : ''
    """
    ${sdfSetup}
    hap.py ${args} ${sdfArgs} ${truth_vcf} ${query_vcf} \\
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
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(truth_indexes, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(query_indexes, stageAs: 'query/*'), path(reference), path(reference_indexes), path(fp_bed), path(fp_indexes), path(stratification_files, stageAs: 'stratification/*'), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.log', hidden: true, emit: runlogs

    script:
    def rustReference = reference.name.endsWith('.gz') ? 'rust-reference.fa' : reference
    def referenceSetup = reference.name.endsWith('.gz')
        ? "gzip -dc ${reference} > ${rustReference}"
        : ':'
    """
    ${referenceSetup}
    hap germline ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${rustReference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result
    """
}
