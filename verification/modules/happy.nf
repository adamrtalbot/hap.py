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
    tuple val(meta), path(truth_vcf), path(truth_tbi), path(query_vcf), path(query_tbi), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    hap.py ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result

    # hap.py can occasionally exit successfully after appending a malformed
    # aggregate ROC row. Turn that upstream data race into a task failure so
    # the retry policy in nextflow.config regenerates the oracle artifacts.
    gzip -cd result.roc.all.csv.gz | awk -F',' '
        NR == 1 {
            if (\$1 != "Type" || \$2 != "Subtype" || \$3 != "Subset" ||
                \$4 != "Filter" || \$5 != "Genotype" || \$6 != "QQ.Field" ||
                \$7 != "QQ") exit 1
            next
        }
        {
            if (\$1 == "" || \$2 == "" || \$3 == "" || \$4 == "" ||
                \$5 == "" || \$6 == "" || \$7 == "") exit 1
            rows++
        }
        END { if (rows == 0) exit 1 }
    '
    """
}

process HAPPY_RUST {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/happy/${meta.id}/rust" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf), path(truth_tbi), path(query_vcf), path(query_tbi), path(reference), path(reference_fai), path(fp_bed), path(fp_bed_tbi), val(args)
    path hap_bin

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    """
    ./${hap_bin} germline ${args} ${truth_vcf} ${query_vcf} \\
        --reference ${reference} \\
        --threads ${task.cpus ?: 1} \\
        --false-positives ${fp_bed} \\
        -o result
    """
}
