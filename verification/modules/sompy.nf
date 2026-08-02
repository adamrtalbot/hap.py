// verification/modules/sompy.nf
//
// Mirrors nf-core/modules/happy/sompy/main.nf. The `-P` and `--count-unk`
// flags match the shape exercised by `rust/src/verification.rs` so the
// two harnesses stay comparable.

process SOMPY_LEGACY {
    tag { "${meta.id}" }
    publishDir { "${params.outdir}/sompy/${meta.id}/legacy" }, mode: 'copy', pattern: 'result*'

    input:
    tuple val(meta), path(truth_vcf), path(query_vcf), path(reference), path(reference_fai), path(fp_bed), val(feature_table), val(args), path(bam), path(bam_index), val(has_bam)

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    bam_args = has_bam ? "--bam ${bam}" : ''
    """
    # The pinned image's pandas removed this display-only option. Explanation
    # output reaches the stale call, so patch a private task copy while leaving
    # the immutable oracle installation untouched.
    cp /opt/conda/bin/som.py legacy-som.py
    sed -i '/pandas.set_option("display.height"/d' legacy-som.py

    # Execute the patched bytes while preserving the installed oracle's
    # logical script path for metadata, command lines, imports, and tracebacks.
    PYTHONPATH=/opt/conda/lib/python27 python -c 'import sys; sys.argv[0]="/opt/conda/bin/som.py"; exec compile(open("legacy-som.py", "rb").read(), sys.argv[0], "exec") in {"__name__":"__main__", "__file__":sys.argv[0]}' ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
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
    path hap_bin

    output:
    tuple val(meta), path('result*'), emit: outputs
    path '.command.{log,sh}', emit: runlogs, optional: true

    script:
    bam_args = has_bam ? "--bam ${bam}" : ''
    """
    # if bed is compressed, uncompress here

    ./${hap_bin} somatic ${args} ${bam_args} ${truth_vcf} ${query_vcf} \\
        -o result \\
        --reference ${reference} \\
        --false-positives ${fp_bed} \\
        --feature-table ${feature_table}
    """
}
