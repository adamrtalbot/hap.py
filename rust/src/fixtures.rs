use crate::cli::SomaticGtMode;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug)]
pub struct CompareFixtureCase {
    pub id: &'static str,
    pub pass_only: bool,
    pub fp_bed: Option<&'static str>,
    pub restrict_bed: Option<&'static str>,
}

#[derive(Clone, Copy, Debug)]
pub struct PreprocessFixtureCase {
    pub id: &'static str,
    pub input: &'static str,
    pub reference: &'static str,
    pub output: &'static str,
    pub pass_only: bool,
    pub regions_bed: Option<&'static str>,
    pub locations: Option<&'static str>,
    pub fixchr: Option<bool>,
    pub set_gt: Option<SomaticGtMode>,
    pub somatic: bool,
    pub filter_nonref: bool,
    pub convert_gvcf_to_vcf: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct QuantifyFixtureCase {
    pub id: &'static str,
    pub input_vcf: &'static str,
    pub reference: &'static str,
}

#[derive(Clone, Copy, Debug)]
pub struct SomaticFixtureCase {
    pub id: &'static str,
    pub truth: &'static str,
    pub query: &'static str,
    pub reference: &'static str,
    pub fp_bed: Option<&'static str>,
    pub count_unk: bool,
    pub include_nonpass: bool,
}

#[derive(Clone, Copy, Debug)]
pub enum ValidateFixtureKind {
    SummaryJson,
    ErrorsBed,
}

#[derive(Clone, Copy, Debug)]
pub struct ValidateFixtureCase {
    pub id: &'static str,
    pub input: &'static str,
    pub reference: Option<&'static str>,
    pub kind: ValidateFixtureKind,
}

pub const COMPARE_FIXTURE_CASES: &[CompareFixtureCase] = &[
    CompareFixtureCase {
        id: "example-haploid",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "numeric-chr-prefix",
        pass_only: false,
        fp_bed: Some("fp.bed"),
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-snp-match",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-snp-mismatch",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-mnv-vs-split",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-homopolymer-insertion",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-empty-query",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-pass-only-filtered-query",
        pass_only: true,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-auto-fixchr",
        pass_only: false,
        fp_bed: None,
        restrict_bed: None,
    },
    CompareFixtureCase {
        id: "synth-region-restricted",
        pass_only: false,
        fp_bed: None,
        restrict_bed: Some("restrict.bed"),
    },
];

pub const PREPROCESS_FIXTURE_CASES: &[PreprocessFixtureCase] = &[
    PreprocessFixtureCase {
        id: "preprocess-pass-fixchr",
        input: "input.vcf",
        reference: "ref.fa",
        output: "processed.vcf.gz",
        pass_only: true,
        regions_bed: Some("restrict.bed"),
        locations: None,
        fixchr: Some(true),
        set_gt: None,
        somatic: false,
        filter_nonref: false,
        convert_gvcf_to_vcf: false,
    },
    PreprocessFixtureCase {
        id: "preprocess-gvcf-somatic",
        input: "input.vcf",
        reference: "ref.fa",
        output: "processed.vcf.gz",
        pass_only: false,
        regions_bed: None,
        locations: None,
        fixchr: Some(true),
        set_gt: Some(SomaticGtMode::Het),
        somatic: false,
        filter_nonref: true,
        convert_gvcf_to_vcf: true,
    },
];

pub const QUANTIFY_FIXTURE_CASES: &[QuantifyFixtureCase] = &[QuantifyFixtureCase {
    id: "quantify-simple",
    input_vcf: "annotated.vcf",
    reference: "ref.fa",
}];

pub const SOMATIC_FIXTURE_CASES: &[SomaticFixtureCase] = &[
    SomaticFixtureCase {
        id: "somatic-simple",
        truth: "truth.vcf",
        query: "query.vcf",
        reference: "ref.fa",
        fp_bed: None,
        count_unk: false,
        include_nonpass: true,
    },
    SomaticFixtureCase {
        id: "somatic-fp-unk",
        truth: "truth.vcf",
        query: "query.vcf",
        reference: "ref.fa",
        fp_bed: Some("fp.bed"),
        count_unk: true,
        include_nonpass: true,
    },
];

pub const VALIDATE_FIXTURE_CASES: &[ValidateFixtureCase] = &[
    ValidateFixtureCase {
        id: "validate-summary",
        input: "input.vcf",
        reference: None,
        kind: ValidateFixtureKind::SummaryJson,
    },
    ValidateFixtureCase {
        id: "validate-errors",
        input: "input.vcf",
        reference: Some("ref.fa"),
        kind: ValidateFixtureKind::ErrorsBed,
    },
];

pub fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn fixtures_root() -> PathBuf {
    repo_root().join("tests/fixtures")
}

impl CompareFixtureCase {
    pub fn fixture_dir(&self) -> PathBuf {
        fixtures_root().join(self.id)
    }

    pub fn expected_dir(&self) -> PathBuf {
        self.fixture_dir().join("expected")
    }

    pub fn reference_path(&self) -> PathBuf {
        self.fixture_dir().join("ref.fa")
    }

    pub fn truth_path(&self) -> PathBuf {
        self.fixture_dir().join("truth.vcf")
    }

    pub fn query_path(&self) -> PathBuf {
        self.fixture_dir().join("query.vcf")
    }

    pub fn fp_bed_path(&self) -> Option<PathBuf> {
        self.fp_bed.map(|path| self.fixture_dir().join(path))
    }

    pub fn restrict_bed_path(&self) -> Option<PathBuf> {
        self.restrict_bed.map(|path| self.fixture_dir().join(path))
    }
}

impl PreprocessFixtureCase {
    pub fn fixture_dir(&self) -> PathBuf {
        fixtures_root().join(self.id)
    }

    pub fn expected_dir(&self) -> PathBuf {
        self.fixture_dir().join("expected")
    }

    pub fn input_path(&self) -> PathBuf {
        self.fixture_dir().join(self.input)
    }

    pub fn reference_path(&self) -> PathBuf {
        self.fixture_dir().join(self.reference)
    }

    pub fn output_path(&self) -> PathBuf {
        self.fixture_dir().join(self.output)
    }

    pub fn regions_bed_path(&self) -> Option<PathBuf> {
        self.regions_bed.map(|path| self.fixture_dir().join(path))
    }
}

impl QuantifyFixtureCase {
    pub fn fixture_dir(&self) -> PathBuf {
        fixtures_root().join(self.id)
    }

    pub fn expected_dir(&self) -> PathBuf {
        self.fixture_dir().join("expected")
    }

    pub fn input_vcf_path(&self) -> PathBuf {
        self.fixture_dir().join(self.input_vcf)
    }

    pub fn reference_path(&self) -> PathBuf {
        self.fixture_dir().join(self.reference)
    }
}

impl SomaticFixtureCase {
    pub fn fixture_dir(&self) -> PathBuf {
        fixtures_root().join(self.id)
    }

    pub fn expected_dir(&self) -> PathBuf {
        self.fixture_dir().join("expected")
    }

    pub fn truth_path(&self) -> PathBuf {
        self.fixture_dir().join(self.truth)
    }

    pub fn query_path(&self) -> PathBuf {
        self.fixture_dir().join(self.query)
    }

    pub fn reference_path(&self) -> PathBuf {
        self.fixture_dir().join(self.reference)
    }

    pub fn fp_bed_path(&self) -> Option<PathBuf> {
        self.fp_bed.map(|path| self.fixture_dir().join(path))
    }
}

impl ValidateFixtureCase {
    pub fn fixture_dir(&self) -> PathBuf {
        fixtures_root().join(self.id)
    }

    pub fn expected_dir(&self) -> PathBuf {
        self.fixture_dir().join("expected")
    }

    pub fn input_path(&self) -> PathBuf {
        self.fixture_dir().join(self.input)
    }

    pub fn reference_path(&self) -> Option<PathBuf> {
        self.reference.map(|path| self.fixture_dir().join(path))
    }
}
