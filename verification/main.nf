#!/usr/bin/env nextflow
//
// verification/main.nf
//
// Concordance pipeline: runs every supported subcommand through both the
// legacy hap.py oracle (inside the pinned Wave container) and the local
// `hap` binary (no container), then diffs outputs per case and writes an
// aggregated markdown report.
//
// Input fixtures are fetched directly with Nextflow's native `file()`
// helper from GitHub permalinks pinned to the upstream Illumina/hap.py
// commit declared in `params.fixture_base`. No local staging step
// is required: each URL resolves to a cached file on first access.
//
// Each row of each samplesheet is one independent Nextflow task per
// side, so adding a row fans the matrix out in parallel.

nextflow.enable.dsl = 2

include { HAPPY_LEGACY ; HAPPY_RUST } from './modules/happy'
include { SOMPY_LEGACY ; SOMPY_RUST } from './modules/sompy'
include { PREPY_LEGACY ; PREPY_RUST } from './modules/prepy'
include { FTXPY_LEGACY ; FTXPY_RUST } from './modules/ftxpy'
include { QFY_ANNOTATE ; QFY_LEGACY ; QFY_RUST } from './modules/qfy'
include { VCFCHECK_LEGACY ; VCFCHECK_RUST } from './modules/vcfcheck'
include { DIFF_OUTPUTS as DIFF_HAPPY    } from './modules/diff'
include { DIFF_OUTPUTS as DIFF_SOMPY    } from './modules/diff'
include { DIFF_OUTPUTS as DIFF_PREPY    } from './modules/diff'
include { DIFF_OUTPUTS as DIFF_FTXPY    } from './modules/diff'
include { DIFF_OUTPUTS as DIFF_QFY      } from './modules/diff'
include { DIFF_OUTPUTS as DIFF_VCFCHECK } from './modules/diff'
include { REPORT                        } from './modules/report'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

def selected_cases() {
    params.cases.toString().tokenize(',').collect { entry -> entry.trim() }
}

// Resolve a repo-relative path against params.fixture_base (URL or
// local path). `file()` will download-and-cache URLs, or symlink
// local paths, into each task's work directory.
def fixture(rel) {
    def path = rel.toString()
    if (path.startsWith('local:')) {
        return file("${projectDir}/${path.substring('local:'.length())}")
    }
    file("${params.fixture_base}/${path}")
}

// Nextflow path inputs are intentionally uniform across compressed/indexed
// and plain-text fixtures.  The process commands never reference the staged
// index variables directly: htslib-compatible tools discover a real sibling
// index by name, while plain-text inputs receive this harmless marker solely
// to keep the tuple shape stable.
def fixtureIndex(rel, marker) {
    def path = rel.toString()
    if (path.endsWith('.gz') || path.endsWith('.bcf')) {
        def suffix = path.endsWith('.bcf') ? '.csi' : '.tbi'
        return fixture("${path}${suffix}")
    }
    file("${projectDir}/assets/no-index-${marker}")
}

def samples(samplesheet, transform) {
    channel.fromPath(samplesheet).splitCsv(header: true).map(transform)
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    hap = file(params.hap_bin, checkIfExists: true)
    verifier = file(params.verify_bin, checkIfExists: true)
    cases = selected_cases()

    statuses = channel.empty()

    // -----------------------------------------------------------------------
    // happy: hap.py germline comparison
    // -----------------------------------------------------------------------
    if (cases.contains('happy')) {
        happy_in = samples(params.happy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'happy']
            tuple(
                meta,
                fixture(row.truth_vcf),
                fixtureIndex(row.truth_vcf, 'truth'),
                fixture(row.query_vcf),
                fixtureIndex(row.query_vcf, 'query'),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.fp_bed),
                fixtureIndex(row.fp_bed, 'fp'),
                (row.args ?: '').toString(),
            )
        }

        HAPPY_LEGACY(happy_in)
        HAPPY_RUST(happy_in, hap)

        happy_pair = HAPPY_LEGACY.out.outputs
            .join(HAPPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'happy', 'result', legacy_files, rust_files)
            }
        DIFF_HAPPY(happy_pair, verifier)
        statuses = statuses.mix(DIFF_HAPPY.out.status)
    }

    // -----------------------------------------------------------------------
    // sompy: som.py somatic comparison
    // -----------------------------------------------------------------------
    if (cases.contains('sompy')) {
        sompy_in = samples(params.sompy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'sompy']
            tuple(
                meta,
                fixture(row.truth_vcf),
                fixture(row.query_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.fp_bed),
                (row.feature_table ?: 'generic').toString(),
                (row.args ?: '').toString(),
                row.bam ? fixture(row.bam) : file("${projectDir}/assets/no-index-bam"),
                row.bam ? fixture("${row.bam}.bai") : file("${projectDir}/assets/no-index-bam-index"),
                row.bam ? true : false,
            )
        }

        SOMPY_LEGACY(sompy_in)
        SOMPY_RUST(sompy_in, hap)

        sompy_pair = SOMPY_LEGACY.out.outputs
            .join(SOMPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'sompy', 'result', legacy_files, rust_files)
            }
        DIFF_SOMPY(sompy_pair, verifier)
        statuses = statuses.mix(DIFF_SOMPY.out.status)
    }

    // -----------------------------------------------------------------------
    // prepy: pre.py preprocessing / normalisation
    // -----------------------------------------------------------------------
    if (cases.contains('prepy')) {
        prepy_in = samples(params.prepy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'prepy']
            tuple(
                meta,
                fixture(row.input_vcf),
                fixtureIndex(row.input_vcf, 'input'),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.regions_bed),
                fixtureIndex(row.regions_bed, 'regions'),
                (row.args ?: '').toString(),
            )
        }

        PREPY_LEGACY(prepy_in)
        PREPY_RUST(prepy_in, hap)

        prepy_pair = PREPY_LEGACY.out.outputs
            .join(PREPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'prepy', 'result', legacy_files, rust_files)
            }
        DIFF_PREPY(prepy_pair, verifier)
        statuses = statuses.mix(DIFF_PREPY.out.status)
    }

    // -----------------------------------------------------------------------
    // ftxpy: ftx.py feature extraction
    // -----------------------------------------------------------------------
    if (cases.contains('ftxpy')) {
        ftxpy_in = samples(params.ftxpy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'ftxpy']
            tuple(
                meta,
                fixture(row.input_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                row.bam ? fixture(row.bam) : file("${projectDir}/assets/no-index-bam"),
                row.bam ? fixture("${row.bam}.bai") : file("${projectDir}/assets/no-index-bam-index"),
                row.bam ? true : false,
                (row.feature_table ?: 'generic').toString(),
                (row.args ?: '').toString(),
            )
        }

        FTXPY_LEGACY(ftxpy_in)
        FTXPY_RUST(ftxpy_in, hap)

        ftxpy_pair = FTXPY_LEGACY.out.outputs
            .join(FTXPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'ftxpy', 'result', legacy_files, rust_files)
            }
        DIFF_FTXPY(ftxpy_pair, verifier)
        statuses = statuses.mix(DIFF_FTXPY.out.status)
    }

    // -----------------------------------------------------------------------
    // qfy: quantify an already xcmp-annotated comparison VCF
    // -----------------------------------------------------------------------
    if (cases.contains('qfy')) {
        qfy_in = samples(params.qfy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'qfy']
            tuple(
                meta,
                fixture(row.truth_vcf),
                fixtureIndex(row.truth_vcf, 'truth'),
                fixture(row.query_vcf),
                fixtureIndex(row.query_vcf, 'query'),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.fp_bed),
                fixtureIndex(row.fp_bed, 'fp'),
                (row.args ?: '').toString(),
            )
        }

        QFY_ANNOTATE(qfy_in)
        QFY_LEGACY(QFY_ANNOTATE.out.annotated)
        QFY_RUST(QFY_ANNOTATE.out.annotated, hap)

        qfy_pair = QFY_LEGACY.out.outputs
            .join(QFY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'qfy', 'result', legacy_files, rust_files)
            }
        DIFF_QFY(qfy_pair, verifier)
        statuses = statuses.mix(DIFF_QFY.out.status)
    }

    // -----------------------------------------------------------------------
    // vcfcheck: VCF validation and summary counts
    // -----------------------------------------------------------------------
    if (cases.contains('vcfcheck')) {
        vcfcheck_in = samples(params.vcfcheck_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'vcfcheck']
            tuple(
                meta,
                fixture(row.input_vcf),
                (row.args ?: '').toString(),
            )
        }

        VCFCHECK_LEGACY(vcfcheck_in)
        VCFCHECK_RUST(vcfcheck_in, hap)

        vcfcheck_pair = VCFCHECK_LEGACY.out.outputs
            .join(VCFCHECK_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'vcfcheck', 'result', legacy_files, rust_files)
            }
        DIFF_VCFCHECK(vcfcheck_pair, verifier)
        statuses = statuses.mix(DIFF_VCFCHECK.out.status)
    }

    // -----------------------------------------------------------------------
    // Aggregate all per-case statuses into a single markdown report.
    //
    // `statuses` carries tuples of shape [meta, case_name, status.json].
    // REPORT only needs the JSON payloads; the meta/case_name stay in
    // the per-case publishDir layout.
    // -----------------------------------------------------------------------
    status_jsons = statuses
        .map { _meta, _case, json -> json }
        .collect()

    REPORT(status_jsons, verifier)
}
