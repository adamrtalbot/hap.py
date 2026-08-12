#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMATIC_OUTPUT_PARITY } from './modules/somatic_output_parity'

def fixture(String name) {
    file("${projectDir}/assets/fixtures/somatic-output-parity/${name}")
}

workflow {
    SOMATIC_OUTPUT_PARITY(
        fixture('truth.vcf'),
        fixture('query.vcf'),
        file("${projectDir}/assets/fixtures/ftx-bam/ref.fa"),
        file("${projectDir}/assets/fixtures/ftx-bam/ref.fa.fai"),
        fixture('confident.bed'),
        file("${projectDir}/scripts/validate_somatic_output_parity.py"),
    )
}
