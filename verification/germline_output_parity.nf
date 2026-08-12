#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GERMLINE_OUTPUT_PARITY } from './modules/germline_output_parity'

def fixture(String name) {
    file("${projectDir}/assets/fixtures/germline-output-parity/${name}")
}

workflow {
    GERMLINE_OUTPUT_PARITY(
        fixture('truth.vcf'),
        fixture('query.vcf'),
        fixture('ref.fa'),
        fixture('ref.fa.fai'),
        fixture('confident.bed'),
        file("${projectDir}/scripts/validate_germline_output_parity.py"),
    )
}
