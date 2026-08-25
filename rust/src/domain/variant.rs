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

/// Legacy edit-classification bits for one REF/ALT pair: bit 1 substitution,
/// bit 2 insertion, bit 4 deletion. Symbolic alternates map to insertion or
/// deletion by their `<DEL...>` spelling.
pub(crate) fn allele_edit_bits(reference: &str, alternate: &str) -> u8 {
    if alternate.starts_with('<') {
        return if alternate.starts_with("<DEL") { 4 } else { 2 };
    }
    let ref_bytes = reference.as_bytes();
    let alt_bytes = alternate.as_bytes();
    let prefix = ref_bytes
        .iter()
        .zip(alt_bytes)
        .take_while(|(left, right)| left == right)
        .count();
    let suffix_limit = (ref_bytes.len() - prefix).min(alt_bytes.len() - prefix);
    let suffix = (0..suffix_limit)
        .take_while(|offset| {
            ref_bytes[ref_bytes.len() - 1 - offset] == alt_bytes[alt_bytes.len() - 1 - offset]
        })
        .count();
    let ref_remaining = ref_bytes.len() - prefix - suffix;
    let alt_remaining = alt_bytes.len() - prefix - suffix;
    match (ref_remaining, alt_remaining) {
        (0, 0) => 0,
        (0, _) => 2,
        (_, 0) => 4,
        (left, right) if left == right => 1,
        (left, right) if left < right => 1 | 2,
        _ => 1 | 4,
    }
}

/// Legacy short name for a `allele_edit_bits` value combined with the optional
/// reference-overlap flag (bit 0x08).
pub(crate) fn legacy_type_bits(bits: u8) -> &'static str {
    const NAMES: [&str; 16] = [
        "nc", "s", "i", "si", "d", "sd", "id", "sid", "r", "rs", "ri", "rsi", "rd", "rsd", "rid",
        "rsid",
    ];
    NAMES[usize::from(bits & 0x0f)]
}
