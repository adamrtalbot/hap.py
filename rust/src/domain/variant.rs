//! Neutral variant record shared by codecs and comparison engines.

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PrimitiveIdentity {
    /// One-based inclusive reference coordinates before final VCF padding.
    /// Insertions use the empty interval immediately before `start`, so
    /// `end == start - 1`.
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) alt: String,
}

#[derive(Clone, Debug)]
pub(crate) struct RawVcfRecord {
    pub(crate) chrom: String,
    pub(crate) pos: usize,
    pub(crate) id: String,
    pub(crate) ref_allele: String,
    pub(crate) alt_allele: String,
    pub(crate) qual: String,
    pub(crate) filter: String,
    pub(crate) info: String,
    pub(crate) format: Option<String>,
    pub(crate) samples: Vec<String>,
    /// Internal preprocessing provenance. Mixed edits retain their source
    /// anchor and sort separately from ordinary insertion primitives during
    /// cross-record location aggregation. This field is never serialized.
    pub(crate) mixed_edit_primitive: bool,
    /// Internal edit before final VCF padding and later left shifts. Legacy
    /// removes exact duplicate edits before padding, while distinct edits may
    /// later serialize to the same REF/ALT spelling. This field exists only
    /// until location aggregation and is never published.
    pub(crate) primitive_identity: Option<PrimitiveIdentity>,
}
