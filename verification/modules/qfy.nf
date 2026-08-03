// verification/modules/qfy.nf
//
// Runs qfy.py and the Rust quantify compatibility command against the same
// xcmp-annotated VCF. QFY_ANNOTATE creates that shared input once from a
// bounded, upstream-pinned chr21 fixture using the oracle comparator.

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
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    qfy.py ${args} ${input_vcf} \
        --reference ${reference} \
        --false-positives ${fp_bed} \
        --threads ${task.cpus ?: 1} \
        --report-prefix result

    # The legacy threaded reporter can exit successfully after dropping rows
    # and appending a malformed aggregate row. Fail the task so the isolated,
    # bounded retry policy regenerates a structurally valid oracle report.
    gzip -cd result.roc.all.csv.gz | awk -F',' '
        NR == 1 {
            expected_fields = NF
            if (expected_fields != 65 && expected_fields != 71) exit 1
            if (\$1 != "Type" || \$2 != "Subtype" || \$3 != "Subset" ||
                \$4 != "Filter" || \$5 != "Genotype" || \$6 != "QQ.Field" ||
                \$7 != "QQ") exit 1
            next
        }
        {
            if (NF != expected_fields) exit 1
            if (\$1 != "SNP" && \$1 != "INDEL") exit 1
            if (\$1 == "" || \$2 == "" || \$3 == "" || \$4 == "" ||
                \$5 == "" || \$6 == "" || \$7 == "") exit 1
            rows++
        }
        END { if (rows == 0) exit 1 }
    '
    """
}

process QFY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/qfy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(input_vcf), path(input_tbi), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args)
    path hap_bin

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    ./${hap_bin} quantify ${args} ${input_vcf} \
        --reference ${reference} \
        --false-positives ${fp_bed} \
        --report-prefix result
    """
}
