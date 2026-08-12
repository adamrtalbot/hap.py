#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMPY_LEGACY } from './modules/sompy'
include { SOMPY_RUST } from './modules/sompy'
include { DIFF_OUTPUTS } from './modules/diff'

def fixture(String name) {
    file("${projectDir}/assets/fixtures/somatic-af-type-parity/${name}")
}

workflow {
    inputs = channel.of(
        tuple(
            [id: 'somatic_af_type_parity', case_name: 'sompy'],
            fixture('truth.vcf'),
            fixture('query.vcf'),
            file("${projectDir}/assets/fixtures/ftx-bam/ref.fa"),
            file("${projectDir}/assets/fixtures/ftx-bam/ref.fa.fai"),
            fixture('confident.bed'),
            'hcc.strelka.indel',
            '-P --bin-afs',
            [],
            [],
        )
    )

    SOMPY_LEGACY(inputs)
    SOMPY_RUST(inputs)

    paired = SOMPY_LEGACY.out.outputs
        .join(SOMPY_RUST.out.outputs, by: 0)
        .map { meta, legacy_files, rust_files ->
            tuple(meta, 'sompy', 'result', legacy_files, rust_files)
        }
    DIFF_OUTPUTS(paired)
}
