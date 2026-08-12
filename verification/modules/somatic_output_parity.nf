process SOMATIC_OUTPUT_PARITY {
    publishDir { params.outdir }, mode: 'copy'

    input:
    path truth_vcf
    path query_vcf
    path reference
    path reference_fai
    path confident_bed
    path validator

    output:
    path 'comparison.json'
    path 'result*'

    script:
    """
    hap somatic ${truth_vcf} ${query_vcf} \
        --reference ${reference} \
        --false-positives ${confident_bed} \
        --feature-table hcc.strelka.indel \
        --bin-afs \
        -P \
        -o result

    UV_CACHE_DIR=.uv-cache uv run --no-project --offline python ${validator} result comparison.json
    """
}
