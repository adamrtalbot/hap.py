//! Neutral variant record shared by codecs and comparison engines.

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
}

impl RawVcfRecord {
    pub(crate) fn sample_values_contain(&self, needle: &str) -> bool {
        self.samples.iter().any(|sample| sample.contains(needle))
    }

    pub(crate) fn replace_sample_values(&mut self, from: &str, to: &str) {
        for sample in &mut self.samples {
            *sample = sample.replace(from, to);
        }
    }
}
