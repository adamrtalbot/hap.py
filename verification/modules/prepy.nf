// verification/modules/prepy.nf
//
// Mirrors nf-core/modules/happy/prepy/main.nf: normalisation and
// preprocessing of a single VCF against a reference, optionally
// restricted to a regions BED.

process PREPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/prepy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_tbi), path(reference), path(reference_fai), path(regions_bed), path(regions_bed_tbi), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    pre.py ${args} ${input_vcf} result.vcf.gz \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        -R ${regions_bed}
    """
}

process PREPY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/prepy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_tbi), path(reference), path(reference_fai), path(regions_bed), path(regions_bed_tbi), val(args)
    path hap_bin

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    ./${hap_bin} pre ${args} ${input_vcf} result.vcf.gz \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        -R ${regions_bed}
    """
}
