#!/usr/bin/env nextflow
//
// verification/main.nf
//
// Concordance pipeline: runs every supported subcommand through both the
// legacy hap.py reference (inside the pinned Wave container) and the local
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

include { HAPPY_LEGACY ; HAPPY_RUST ; HAPPY_RUST_CONTRACT } from './modules/happy'
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

// Stage real sibling indexes when the fixture format uses one. Plain-text
// inputs use an empty collection so no placeholder file is required.
def fixtureIndexes(rel) {
    def path = rel.toString()
    if (path.endsWith('.gz') || path.endsWith('.bcf')) {
        def suffix = path.endsWith('.bcf') ? '.csi' : '.tbi'
        return [fixture("${path}${suffix}")]
    }
    []
}

def samples(samplesheet, transform) {
    channel.fromPath(samplesheet).splitCsv(header: true).map(transform)
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    cases = selected_cases()

    statuses = channel.empty()

    // -----------------------------------------------------------------------
    // happy: hap.py germline comparison
    // -----------------------------------------------------------------------
    if (cases.contains('happy')) {
        happy_rows = samples(params.happy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'happy']
            def additionalStratificationBeds = (row.stratification_beds ?: '')
                .toString()
                .tokenize(';')
            def stratificationFiles = ([row.stratification_tsv, row.stratification_bed] + additionalStratificationBeds)
                .findAll { path -> path }
                .collect { path -> fixture(path) }
            tuple(
                meta,
                fixture(row.truth_vcf),
                fixtureIndexes(row.truth_vcf),
                fixture(row.query_vcf),
                fixtureIndexes(row.query_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.fp_bed),
                fixtureIndexes(row.fp_bed),
                stratificationFiles,
                row.reference_sdf ? [fixture(row.reference_sdf)] : [],
                (row.args ?: '').toString(),
            )
        }

        happy_legacy_in = happy_rows.map { meta, truth, truthIndexes, query, queryIndexes, reference, referenceFai, fpBed, fpIndexes, stratificationFiles, referenceSdf, args ->
            tuple(
                meta,
                truth,
                truthIndexes,
                query,
                queryIndexes,
                reference,
                referenceFai,
                fpBed,
                fpIndexes,
                stratificationFiles,
                referenceSdf,
                args,
            )
        }
        happy_rust_in = happy_rows.map { meta, truth, truthIndexes, query, queryIndexes, reference, referenceFai, fpBed, fpIndexes, stratificationFiles, _referenceSdf, args ->
            tuple(
                meta,
                truth,
                truthIndexes,
                query,
                queryIndexes,
                reference,
                referenceFai,
                fpBed,
                fpIndexes,
                stratificationFiles,
                args,
            )
        }

        HAPPY_LEGACY(happy_legacy_in)
        HAPPY_RUST(happy_rust_in)

        happy_pair = HAPPY_LEGACY.out.outputs
            .join(HAPPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'happy', 'result', legacy_files, rust_files)
            }
        DIFF_HAPPY(happy_pair)
        statuses = statuses.mix(DIFF_HAPPY.out.comparison)

        happy_contracts = channel.of(
            tuple(
                [id: 'matrix_vcfeval_deprecated_flags', case_name: 'happy'],
                fixture('local:assets/fixtures/vcfeval-matrix/truth.vcf'),
                fixture('local:assets/fixtures/vcfeval-matrix/query.vcf'),
                fixture('local:assets/fixtures/vcfeval-matrix/ref.fa'),
                fixture('local:assets/fixtures/vcfeval-matrix/ref.fa.fai'),
                fixture('local:assets/fixtures/vcfeval-matrix/confident.bed'),
                fixture('local:scripts/validate_vcfeval_contract.py'),
                '--engine vcfeval --no-leftshift -D --no-adjust-conf-regions --engine-vcfeval-path /definitely/missing/rtg --engine-vcfeval-template /definitely/missing/template.sdf',
                'deprecated_flags',
            ),
            tuple(
                [id: 'matrix_vcfeval_preserve_info', case_name: 'happy'],
                fixture('local:assets/fixtures/vcfeval-matrix/truth.vcf'),
                fixture('local:assets/fixtures/vcfeval-matrix/query.vcf'),
                fixture('local:assets/fixtures/vcfeval-matrix/ref.fa'),
                fixture('local:assets/fixtures/vcfeval-matrix/ref.fa.fai'),
                fixture('local:assets/fixtures/vcfeval-matrix/confident.bed'),
                fixture('local:scripts/validate_vcfeval_contract.py'),
                '--engine vcfeval --no-leftshift -D --no-adjust-conf-regions --preserve-info',
                'preserve_info',
            ),
        )
        HAPPY_RUST_CONTRACT(happy_contracts)
        statuses = statuses.mix(HAPPY_RUST_CONTRACT.out.comparison)
    }

    // -----------------------------------------------------------------------
    // sompy: som.py somatic comparison
    // -----------------------------------------------------------------------
    if (cases.contains('sompy')) {
        sompy_in = samples(params.sompy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'sompy']
            def bamPaths = [row.bam, row.bam2].findAll { path -> path }
            tuple(
                meta,
                fixture(row.truth_vcf),
                fixture(row.query_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.fp_bed),
                (row.feature_table ?: 'generic').toString(),
                (row.args ?: '').toString(),
                bamPaths.collect { path -> fixture(path) },
                bamPaths.collect { path -> fixture("${path}.bai") },
            )
        }

        SOMPY_LEGACY(sompy_in)
        SOMPY_RUST(sompy_in)

        sompy_pair = SOMPY_LEGACY.out.outputs
            .join(SOMPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'sompy', 'result', legacy_files, rust_files)
            }
        DIFF_SOMPY(sompy_pair)
        statuses = statuses.mix(DIFF_SOMPY.out.comparison)
    }

    // -----------------------------------------------------------------------
    // prepy: pre.py preprocessing / normalisation
    // -----------------------------------------------------------------------
    if (cases.contains('prepy')) {
        prepy_in = samples(params.prepy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'prepy', cpus: (row.cpus ?: 1).toString().toInteger()]
            tuple(
                meta,
                fixture(row.input_vcf),
                fixtureIndexes(row.input_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.regions_bed),
                fixtureIndexes(row.regions_bed),
                (row.args ?: '').toString(),
            )
        }

        PREPY_LEGACY(prepy_in)
        PREPY_RUST(prepy_in)

        prepy_pair = PREPY_LEGACY.out.outputs
            .join(PREPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'prepy', 'result', legacy_files, rust_files)
            }
        DIFF_PREPY(prepy_pair)
        statuses = statuses.mix(DIFF_PREPY.out.comparison)
    }

    // -----------------------------------------------------------------------
    // ftxpy: ftx.py feature extraction
    // -----------------------------------------------------------------------
    if (cases.contains('ftxpy')) {
        ftxpy_in = samples(params.ftxpy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'ftxpy']
            def bamPaths = [row.bam, row.bam2].findAll { path -> path }
            tuple(
                meta,
                fixture(row.input_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                bamPaths.collect { path -> fixture(path) },
                bamPaths.collect { path -> fixture("${path}.bai") },
                row.regions_bed ? [fixture(row.regions_bed)] : [],
                (row.feature_table ?: 'generic').toString(),
                (row.args ?: '').toString(),
            )
        }

        FTXPY_LEGACY(ftxpy_in)
        FTXPY_RUST(ftxpy_in)

        ftxpy_pair = FTXPY_LEGACY.out.outputs
            .join(FTXPY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'ftxpy', 'result', legacy_files, rust_files)
            }
        DIFF_FTXPY(ftxpy_pair)
        statuses = statuses.mix(DIFF_FTXPY.out.comparison)
    }

    // -----------------------------------------------------------------------
    // qfy: quantify an already xcmp-annotated comparison VCF
    // -----------------------------------------------------------------------
    if (cases.contains('qfy')) {
        qfy_in = samples(params.qfy_samplesheet) { row ->
            def meta = [id: row.sample_id, case_name: 'qfy']
            def stratificationFiles = [row.stratification_tsv, row.stratification_bed]
                .findAll { path -> path }
                .collect { path -> fixture(path) }
            tuple(
                meta,
                fixture(row.truth_vcf),
                fixtureIndexes(row.truth_vcf),
                fixture(row.query_vcf),
                fixtureIndexes(row.query_vcf),
                fixture(row.reference),
                fixture("${row.reference}.fai"),
                fixture(row.fp_bed),
                fixtureIndexes(row.fp_bed),
                stratificationFiles,
                (row.args ?: '').toString(),
            )
        }

        QFY_ANNOTATE(qfy_in)
        QFY_LEGACY(QFY_ANNOTATE.out.annotated)
        QFY_RUST(QFY_ANNOTATE.out.annotated)

        qfy_pair = QFY_LEGACY.out.outputs
            .join(QFY_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'qfy', 'result', legacy_files, rust_files)
            }
        DIFF_QFY(qfy_pair)
        statuses = statuses.mix(DIFF_QFY.out.comparison)
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
                fixtureIndexes(row.input_vcf),
                (row.input_mode ?: 'positional').toString(),
                (row.args ?: '').toString(),
            )
        }

        VCFCHECK_LEGACY(vcfcheck_in)
        VCFCHECK_RUST(vcfcheck_in)

        vcfcheck_pair = VCFCHECK_LEGACY.out.outputs
            .join(VCFCHECK_RUST.out.outputs, by: 0)
            .map { meta, legacy_files, rust_files ->
                tuple(meta, 'vcfcheck', 'result', legacy_files, rust_files)
            }
        DIFF_VCFCHECK(vcfcheck_pair)
        statuses = statuses.mix(DIFF_VCFCHECK.out.comparison)
    }

    // -----------------------------------------------------------------------
    // Aggregate all per-case statuses into a single markdown report.
    //
    // `statuses` carries tuples of shape [meta, case_name, status.json].
    // REPORT only needs the JSON payloads; the meta/case_name stay in
    // the per-case publishDir layout.
    // -----------------------------------------------------------------------
    comparison_jsons = statuses
        .map { _meta, _case, json -> json }
        .collect()

    REPORT(comparison_jsons)
}
