//! Allele-frequency binning and legacy metric rounding.

use super::MAX_AF_BINS;
use anyhow::{Context, Result, bail};

pub(super) fn rounded_metric(numerator: usize, denominator: usize) -> Option<f64> {
    (denominator != 0).then(|| round_four(numerator as f64 / denominator as f64))
}

pub(super) fn round_four(value: f64) -> f64 {
    (value * 10_000.0).round() / 10_000.0
}

pub(super) fn parse_af_bins(raw: &str) -> Result<Vec<(f64, f64)>> {
    let fields = raw.split(',').collect::<Vec<_>>();
    if fields.is_empty() || fields.iter().any(|field| field.trim().is_empty()) {
        bail!("AF bin size list must not be empty");
    }
    let bins = fields
        .into_iter()
        .map(|part| {
            part.parse::<f64>()
                .with_context(|| format!("failed to parse AF bin size '{part}'"))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut out = Vec::new();
    let mut start: f64 = 0.0;
    let mut idx = 0usize;
    for _ in 0..MAX_AF_BINS {
        if !start.is_finite() || start >= 1.0 || bins.is_empty() {
            return Ok(out);
        }
        let mut end = start + bins[idx];
        if end >= 1.0 {
            end = 1.000_000_01;
        }
        if start >= end {
            return Ok(out);
        }
        out.push((start, end));
        start = end;
        idx = (idx + 1) % bins.len();
    }
    if !start.is_finite() || start >= 1.0 {
        Ok(out)
    } else {
        bail!("AF bin sizes produce more than {MAX_AF_BINS} bins")
    }
}

fn format_af_bound(value: f64) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else {
        format!("{value:.6}")
    }
}

pub(super) fn format_af_interval(start: f64, end: f64) -> String {
    format!("{}-{}", format_af_bound(start), format_af_bound(end))
}

#[cfg(test)]
mod tests {
    use super::{format_af_interval, parse_af_bins};

    #[test]
    fn bins_include_one_using_the_legacy_epsilon() {
        assert_eq!(parse_af_bins("inf").unwrap(), vec![(0.0, 1.000_000_01)]);
        assert_eq!(format_af_interval(0.0, 1.000_000_01), "0.000000-1.000000");
    }
}
