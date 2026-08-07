//! Neutral data carried between ROC contribution and accumulation stages.

use crate::domain::CountsBucket;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) struct RowKey {
    pub(super) ty: String,
    pub(super) subtype: String,
    pub(super) subset: String,
    pub(super) filter: String,
    pub(super) genotype: String,
    pub(super) qq_field: String,
}

impl RowKey {
    #[cfg(test)]
    pub(super) fn new(ty: &str, subtype: &str, subset: &str, filter: &str) -> Self {
        Self::new_with_qq_field(ty, subtype, subset, filter, "QUAL")
    }

    pub(super) fn new_with_qq_field(
        ty: &str,
        subtype: &str,
        subset: &str,
        filter: &str,
        qq_field: &str,
    ) -> Self {
        Self {
            ty: ty.to_string(),
            subtype: subtype.to_string(),
            subset: subset.to_string(),
            filter: filter.to_string(),
            genotype: "*".to_string(),
            qq_field: qq_field.to_string(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct Cumul {
    pub(super) truth_tp: CountsBucket,
    pub(super) truth_fn: CountsBucket,
    pub(super) query_tp: CountsBucket,
    pub(super) query_fp: CountsBucket,
    pub(super) query_unk: CountsBucket,
    pub(super) fp_gt: usize,
    pub(super) fp_al: usize,
}

impl Cumul {
    pub(super) fn add(&mut self, other: &Self) {
        add_bucket(&mut self.truth_tp, &other.truth_tp);
        add_bucket(&mut self.truth_fn, &other.truth_fn);
        add_bucket(&mut self.query_tp, &other.query_tp);
        add_bucket(&mut self.query_fp, &other.query_fp);
        add_bucket(&mut self.query_unk, &other.query_unk);
        self.fp_gt += other.fp_gt;
        self.fp_al += other.fp_al;
    }

    pub(super) fn truth_total(&self) -> CountsBucket {
        sum_buckets(&self.truth_tp, &self.truth_fn)
    }

    pub(super) fn query_total(&self) -> CountsBucket {
        sum_buckets(
            &sum_buckets(&self.query_tp, &self.query_fp),
            &self.query_unk,
        )
    }
}

fn add_bucket(dst: &mut CountsBucket, src: &CountsBucket) {
    dst.total += src.total;
    dst.ti += src.ti;
    dst.tv += src.tv;
    dst.het += src.het;
    dst.homalt += src.homalt;
}

fn sum_buckets(left: &CountsBucket, right: &CountsBucket) -> CountsBucket {
    CountsBucket {
        total: left.total + right.total,
        ti: left.ti + right.ti,
        tv: left.tv + right.tv,
        het: left.het + right.het,
        homalt: left.homalt + right.homalt,
    }
}

pub(super) struct Sample<'a> {
    pub(super) format_keys: &'a [&'a str],
    pub(super) parts: &'a [&'a str],
    pub(super) gt: Option<&'a str>,
    pub(super) bd: Option<&'a str>,
    pub(super) bi: Option<&'a str>,
    pub(super) bvt: Option<&'a str>,
    pub(super) blt: Option<&'a str>,
    pub(super) qq: Option<&'a str>,
}

impl<'a> Sample<'a> {
    pub(super) fn new(format_keys: &'a [&'a str], parts: &'a [&'a str]) -> Self {
        let lookup = |name: &str| -> Option<&'a str> {
            format_keys
                .iter()
                .position(|key| *key == name)
                .and_then(|index| parts.get(index).copied())
        };
        Self {
            format_keys,
            parts,
            gt: lookup("GT"),
            bd: lookup("BD"),
            bi: lookup("BI"),
            bvt: lookup("BVT"),
            blt: lookup("BLT"),
            qq: lookup("QQ"),
        }
    }

    pub(super) fn variant_type(&self) -> Option<&'a str> {
        self.bvt.filter(|value| matches!(*value, "SNP" | "INDEL"))
    }

    pub(super) fn roc_qq(&self) -> Option<f64> {
        self.qq.and_then(|raw| raw.parse::<f64>().ok())
    }

    pub(super) fn roc_value(&self, field: &str, record_qual: &str, info: &str) -> Option<f64> {
        if field == "QUAL" || field == "QQ" {
            return self.roc_qq();
        }
        if field == "." {
            return None;
        }
        if let Some(value) = info.split(';').find_map(|entry| {
            entry
                .split_once('=')
                .filter(|(key, _)| *key == field)
                .map(|(_, value)| value)
        }) {
            return parse_roc_number(value);
        }
        self.format_keys
            .iter()
            .position(|key| *key == field)
            .and_then(|index| self.parts.get(index).copied())
            .and_then(parse_roc_number)
            .or_else(|| {
                (field == "QUAL")
                    .then(|| parse_roc_number(record_qual))
                    .flatten()
            })
    }
}

fn parse_roc_number(raw: &str) -> Option<f64> {
    raw.split(',')
        .next()
        .filter(|value| !matches!(*value, "" | "."))
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite())
}
