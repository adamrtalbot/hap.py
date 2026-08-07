/// Counts accumulated for one side and classification of a benchmark.
#[derive(Clone, Debug, Default)]
pub(crate) struct CountsBucket {
    pub(crate) total: usize,
    pub(crate) ti: usize,
    pub(crate) tv: usize,
    pub(crate) het: usize,
    pub(crate) homalt: usize,
}

/// Truth and query counts for one variant grouping.
#[derive(Clone, Debug, Default)]
pub(crate) struct TypeCounts {
    pub(crate) truth_total: CountsBucket,
    pub(crate) truth_tp: CountsBucket,
    pub(crate) truth_fn: CountsBucket,
    pub(crate) query_total: CountsBucket,
    pub(crate) query_tp: CountsBucket,
    pub(crate) query_fp: CountsBucket,
    pub(crate) query_unk: CountsBucket,
}

/// A comparison record plus the domain facts needed by report engines.
#[derive(Clone, Debug)]
pub(crate) struct AnnotatedRow {
    pub(crate) sort_key: (String, usize, usize, usize),
    pub(crate) record: super::RawVcfRecord,
    /// Whether the originating query variant was PASS-filtered.
    pub(crate) query_pass: bool,
    /// False-positive subclass (`gt` or `al`) when this is an FP row.
    pub(crate) fp_class: Option<&'static str>,
    /// Legacy xcmp block-level comparison classification.
    pub(crate) xcmp_ctype: Option<&'static str>,
    pub(crate) xcmp_hap_match: bool,
}
