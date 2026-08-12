// verification/modules/happy.nf
//
// Mirrors nf-core/modules/happy/happy/main.nf for the germline comparison.
// Two processes with identical input shape:
//   * HAPPY_LEGACY uses `hap.py` inside the pinned Wave container.
//   * HAPPY_RUST   uses `hap germline` from the locally-built binary.

process HAPPY_PREPARE {
    tag { "${meta.id}" }

    container 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'

    input:
    tuple val(meta), path(truth_vcf, stageAs: 'truth/*'), path(truth_indexes, stageAs: 'truth/*'), path(query_vcf, stageAs: 'query/*'), path(query_indexes, stageAs: 'query/*'), path(reference), path(reference_indexes), path(fp_bed), path(fp_indexes), path(stratification_files, stageAs: 'stratification/*'), path(reference_sdf, stageAs: 'vcfeval-sdf/*'), val(args)

    output:
    tuple val(meta), path('truth.prepared.vcf.gz'), path('truth.prepared.vcf.gz.tbi'), path('query.prepared.vcf.gz'), path('query.prepared.vcf.gz.tbi'), path(reference), path(reference_indexes), path(fp_bed), path(fp_indexes), path(stratification_files), path(reference_sdf), val(args), emit: rows

    script:
    """
    set -o pipefail

    bcftools norm --fasta-ref ${reference} --multiallelics -any --check-ref w --output-type u ${truth_vcf} \
        | bcftools norm --fasta-ref ${reference} --rm-dup exact --check-ref w --output-type u \
        | bcftools sort --output-type u \
        | bcftools norm --fasta-ref ${reference} --check-ref w --output-type u \
        | bcftools view --include 'GT="alt"' --output-type z --write-index=tbi --output truth.prepared.vcf.gz

    bcftools sort --output-type u ${query_vcf} \
        | bcftools view --samples '${meta.query_sample}' --output-type u \
        | bcftools norm --fasta-ref ${reference} --multiallelics -any --check-ref w --output-type u \
        | bcftools norm --fasta-ref ${reference} --rm-dup exact --check-ref w --output-type u \
        | bcftools sort --output-type u \
        | bcftools norm --fasta-ref ${reference} --check-ref w --output-type u \
        | bcftools view --include 'GT="alt"' --output-type z --write-index=tbi --output query.prepared.vcf.gz
    """
}

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
    """
    hap germline ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result
    """
}

process HAPPY_RUST_CONTRACT {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/happy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'
    publishDir { "${params.outdir}/happy/${meta.id}" }, mode: 'copy', pattern: 'comparison.json'

    input:
    tuple val(meta), path(truth_vcf), path(query_vcf), path(reference), path(reference_fai), path(fp_bed), path(contract_validator), val(args), val(contract)

    output:
    tuple val(meta), val('happy'), path('comparison.json'), emit: comparison
    path 'result*', emit: artifacts

    script:
    """
    hap germline ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result \\
        2> stderr.log

    python3 ${contract_validator} '${meta.id}' '${contract}'
    """
}
