//! Bounded disk-backed storage for comparison rows and source metadata.

use super::output::DecorationIndex;
use super::{AnnotatedRow, MAX_CLUSTER_VARIANTS};
use crate::adapters::vcf::{self, ValidatedVcfReader, ValidatedVcfRecord, VariantKey};
use crate::domain::{FpClass, RawVcfRecord, SortKey, XcmpCtype};
use anyhow::{Context, Result, bail};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const COMPARISON_ROW_CHUNK: usize = 65_536;
const COMPARISON_MERGE_FAN_IN: usize = 32;
type ComparisonSortKey = (
    String,
    usize,
    u8,
    usize,
    usize,
    String,
    String,
    u8,
    String,
    usize,
);

fn legacy_same_key_rank(samples: &[String]) -> u8 {
    if samples
        .first()
        .is_some_and(|sample| !sample.starts_with("./.") && !sample.contains("NOCALL:nocall"))
    {
        return 0;
    }
    let query_gt = samples
        .get(1)
        .and_then(|sample| sample.split(':').next())
        .unwrap_or(".");
    let alleles = query_gt.split(['/', '|']).collect::<Vec<_>>();
    if alleles.len() == 2 && alleles[0] != "0" && alleles[0] != "." && alleles[0] == alleles[1] {
        1
    } else {
        2
    }
}

pub(super) struct ComparisonRowSpool {
    buffer: Vec<(ComparisonSortKey, AnnotatedRow)>,
    chunks: Vec<tempfile::TempPath>,
    serial: usize,
}

impl ComparisonRowSpool {
    pub(super) fn new() -> Self {
        Self {
            buffer: Vec::new(),
            chunks: Vec::new(),
            serial: 0,
        }
    }

    pub(super) fn push(
        &mut self,
        row: AnnotatedRow,
        filtered_truth_match: bool,
        sort_line: String,
    ) -> Result<()> {
        let raw = row.record.raw();
        // Pinned classified VCFs first group filtered truth matches, then keep
        // the side/type ranks assigned by row construction before comparing
        // allele spelling. The HG001 graph-order rule therefore remains the
        // leading rank while ordinary truth/query row precedence is retained.
        let key = (
            row.sort_key.chrom.clone(),
            row.sort_key.pos,
            u8::from(!filtered_truth_match),
            row.sort_key.side_rank,
            row.sort_key.type_rank,
            raw.ref_allele.clone(),
            raw.alt_allele.clone(),
            legacy_same_key_rank(&raw.samples),
            sort_line,
            self.serial,
        );
        self.serial = self
            .serial
            .checked_add(1)
            .context("comparison row serial overflow")?;
        self.buffer.push((key, row));
        if self.buffer.len() >= COMPARISON_ROW_CHUNK {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer.sort_by(|left, right| left.0.cmp(&right.0));
        let mut chunk =
            tempfile::NamedTempFile::new().context("failed to create comparison sort chunk")?;
        {
            let mut writer = BufWriter::new(chunk.as_file_mut());
            for (key, row) in self.buffer.drain(..) {
                write_comparison_spool_row(&mut writer, &key, &row)?;
            }
            writer.flush()?;
        }
        self.chunks.push(chunk.into_temp_path());
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<ComparisonRowFile> {
        self.flush()?;
        let chunks = collapse_comparison_chunks(self.chunks)?;
        let mut output =
            tempfile::NamedTempFile::new().context("failed to create ordered comparison spool")?;
        {
            let mut writer = BufWriter::new(output.as_file_mut());
            let mut merge = ComparisonRowMerge::open(chunks)?;
            while let Some(entry) = merge.next_keyed() {
                let (key, row) = entry?;
                write_comparison_spool_row(&mut writer, &key, &row)?;
            }
            writer.flush()?;
        }
        Ok(ComparisonRowFile {
            path: output.into_temp_path(),
        })
    }
}

pub(super) struct ComparisonRowFile {
    path: tempfile::TempPath,
}

impl ComparisonRowFile {
    pub(super) fn rows(&self) -> Result<ComparisonRowReader> {
        Ok(ComparisonRowReader {
            lines: BufReader::new(File::open(&self.path)?).lines(),
        })
    }
}

pub(super) struct ComparisonRowReader {
    lines: std::io::Lines<BufReader<File>>,
}

impl Iterator for ComparisonRowReader {
    type Item = Result<AnnotatedRow>;

    fn next(&mut self) -> Option<Self::Item> {
        self.lines.next().map(|line| {
            let line = line?;
            parse_comparison_spool_row(&line).map(|(_, row)| row)
        })
    }
}

struct ComparisonRowMerge {
    _chunks: Vec<tempfile::TempPath>,
    readers: Vec<std::io::Lines<BufReader<File>>>,
    current: Vec<Option<(ComparisonSortKey, AnnotatedRow)>>,
    heap: BinaryHeap<Reverse<(ComparisonSortKey, usize)>>,
}

impl ComparisonRowMerge {
    fn open(chunks: Vec<tempfile::TempPath>) -> Result<Self> {
        let mut readers = Vec::with_capacity(chunks.len());
        for chunk in &chunks {
            readers.push(BufReader::new(File::open(chunk)?).lines());
        }
        let mut merge = Self {
            current: (0..readers.len()).map(|_| None).collect(),
            readers,
            heap: BinaryHeap::new(),
            _chunks: chunks,
        };
        for index in 0..merge.readers.len() {
            merge.advance(index)?;
        }
        Ok(merge)
    }

    fn advance(&mut self, index: usize) -> Result<()> {
        let Some(line) = self.readers[index].next() else {
            return Ok(());
        };
        let entry = parse_comparison_spool_row(&line?)?;
        self.heap.push(Reverse((entry.0.clone(), index)));
        self.current[index] = Some(entry);
        Ok(())
    }

    fn next_keyed(&mut self) -> Option<Result<(ComparisonSortKey, AnnotatedRow)>> {
        let Reverse((_, index)) = self.heap.pop()?;
        let entry = self.current[index]
            .take()
            .expect("comparison merge heap entry has a current row");
        if let Err(error) = self.advance(index) {
            return Some(Err(error));
        }
        Some(Ok(entry))
    }
}

fn write_comparison_spool_row(
    writer: &mut dyn Write,
    key: &ComparisonSortKey,
    row: &AnnotatedRow,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        key.0,
        key.1,
        key.2,
        key.3,
        key.4,
        key.5,
        key.6,
        key.7,
        key.9,
        hex_encode(key.8.as_bytes()),
        u8::from(row.query_pass),
        row.fp_class.map_or(".", FpClass::as_str),
        row.xcmp_ctype.map_or(".", XcmpCtype::as_str),
        u8::from(row.xcmp_hap_match),
        row.record.raw().to_line()
    )?;
    Ok(())
}

fn parse_comparison_spool_row(line: &str) -> Result<(ComparisonSortKey, AnnotatedRow)> {
    let mut fields = line.splitn(15, '\t');
    let chrom = fields
        .next()
        .context("comparison spool lacks chromosome")?
        .to_string();
    let pos = fields
        .next()
        .context("comparison spool lacks position")?
        .parse()?;
    let filtered = fields
        .next()
        .context("comparison spool lacks filtered rank")?
        .parse()?;
    let row_side_rank = fields
        .next()
        .context("comparison spool lacks row side rank")?
        .parse()?;
    let row_type_rank = fields
        .next()
        .context("comparison spool lacks row type rank")?
        .parse()?;
    let reference = fields
        .next()
        .context("comparison spool lacks reference sort key")?
        .to_string();
    let alternate = fields
        .next()
        .context("comparison spool lacks alternate sort key")?
        .to_string();
    let same_key_rank = fields
        .next()
        .context("comparison spool lacks same-key rank")?
        .parse()?;
    let serial = fields
        .next()
        .context("comparison spool lacks serial")?
        .parse()?;
    let sort_line = String::from_utf8(hex_decode(
        fields
            .next()
            .context("comparison spool lacks original sort row")?,
    )?)
    .context("comparison spool sort row is not UTF-8")?;
    let query_pass = fields.next() == Some("1");
    let fp_class = match fields.next() {
        Some("gt") => Some(FpClass::Gt),
        Some("al") => Some(FpClass::Al),
        _ => None,
    };
    let xcmp_ctype = match fields.next() {
        Some(".") | None => None,
        Some("simple:match") => Some(XcmpCtype::SimpleMatch),
        Some("simple:mismatch") => Some(XcmpCtype::SimpleMismatch),
        Some("hap:match") => Some(XcmpCtype::HapMatch),
        Some("hap:mismatch") => Some(XcmpCtype::HapMismatch),
        Some("hapfail:mismatch") => Some(XcmpCtype::HapfailMismatch),
        Some(value) => bail!("comparison spool contains unknown XCMP context {value}"),
    };
    let xcmp_hap_match = fields.next() == Some("1");
    let row_line = fields.next().context("comparison spool lacks VCF row")?;
    let record = RawVcfRecord::from_line(row_line, std::path::Path::new("comparison-row.spool"))?;
    let key = (
        chrom.clone(),
        pos,
        filtered,
        row_side_rank,
        row_type_rank,
        reference,
        alternate,
        same_key_rank,
        sort_line,
        serial,
    );
    Ok((
        key,
        AnnotatedRow {
            sort_key: SortKey::new(chrom, pos, row_side_rank, row_type_rank),
            record: record.into(),
            query_pass,
            fp_class,
            xcmp_ctype,
            xcmp_hap_match,
        },
    ))
}

fn hex_encode(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(HEX[(byte >> 4) as usize] as char);
        output.push(HEX[(byte & 0x0f) as usize] as char);
    }
    output
}

fn hex_decode(value: &str) -> Result<Vec<u8>> {
    if !value.len().is_multiple_of(2) {
        bail!("comparison spool sort row has odd hex length");
    }
    value
        .as_bytes()
        .chunks_exact(2)
        .map(|pair| {
            let digit = |byte: u8| match byte {
                b'0'..=b'9' => Ok(byte - b'0'),
                b'a'..=b'f' => Ok(byte - b'a' + 10),
                _ => bail!("comparison spool sort row contains invalid hex"),
            };
            Ok((digit(pair[0])? << 4) | digit(pair[1])?)
        })
        .collect()
}

fn collapse_comparison_chunks(
    mut chunks: Vec<tempfile::TempPath>,
) -> Result<Vec<tempfile::TempPath>> {
    while chunks.len() > COMPARISON_MERGE_FAN_IN {
        let mut merged = Vec::with_capacity(chunks.len().div_ceil(COMPARISON_MERGE_FAN_IN));
        let mut remaining = chunks.into_iter();
        loop {
            let batch = remaining
                .by_ref()
                .take(COMPARISON_MERGE_FAN_IN)
                .collect::<Vec<_>>();
            if batch.is_empty() {
                break;
            }
            let mut output = tempfile::NamedTempFile::new()
                .context("failed to create comparison merge chunk")?;
            {
                let mut writer = BufWriter::new(output.as_file_mut());
                let mut merge = ComparisonRowMerge::open(batch)?;
                while let Some(entry) = merge.next_keyed() {
                    let (key, row) = entry?;
                    write_comparison_spool_row(&mut writer, &key, &row)?;
                }
                writer.flush()?;
            }
            merged.push(output.into_temp_path());
        }
        chunks = merged;
    }
    Ok(chunks)
}

pub(super) struct ComparisonContigSpool {
    path: tempfile::TempPath,
}

pub(super) fn spool_comparison_contigs(
    path: &Path,
) -> Result<BTreeMap<String, ComparisonContigSpool>> {
    let mut spools = BTreeMap::new();
    let mut active: Option<(String, tempfile::NamedTempFile)> = None;
    for record in vcf::open_validated_vcf(path)? {
        let record = record?;
        if active
            .as_ref()
            .is_none_or(|(chrom, _)| chrom != &record.raw().chrom)
        {
            if let Some((chrom, mut file)) = active.take() {
                file.as_file_mut().flush()?;
                if spools
                    .insert(
                        chrom.clone(),
                        ComparisonContigSpool {
                            path: file.into_temp_path(),
                        },
                    )
                    .is_some()
                {
                    bail!("comparison records for chromosome {chrom} are not contiguous");
                }
            }
            active = Some((
                record.raw().chrom.clone(),
                tempfile::NamedTempFile::new()
                    .context("failed to create comparison metadata spool")?,
            ));
        }
        writeln!(
            active
                .as_mut()
                .expect("comparison metadata spool was just created")
                .1
                .as_file_mut(),
            "{}",
            record.raw().to_line()
        )?;
    }
    if let Some((chrom, mut file)) = active {
        file.as_file_mut().flush()?;
        spools.insert(
            chrom,
            ComparisonContigSpool {
                path: file.into_temp_path(),
            },
        );
    }
    Ok(spools)
}

const COMPARISON_METADATA_LOOKBEHIND: usize = 1_024;
const MAX_COMPARISON_METADATA_RETAINED_RECORDS: usize = MAX_CLUSTER_VARIANTS * 2;

pub(super) struct ComparisonMetadataCursor {
    reader: ValidatedVcfReader,
    retained: VecDeque<ValidatedVcfRecord>,
    pending: Option<ValidatedVcfRecord>,
    furthest_start: usize,
}

pub(super) struct ActiveComparisonMetadata {
    pub(super) chrom: String,
    pub(super) truth: Option<ComparisonMetadataCursor>,
    pub(super) query: Option<ComparisonMetadataCursor>,
}

/// Collection settings for one metadata sweep: whether to keep INFO, which QQ
/// field drives ROC decorations, and whether to record filtered truth keys.
pub(super) struct CollectOptions<'a> {
    pub(super) preserve_info: bool,
    pub(super) roc_field: &'a str,
    pub(super) collect_filtered: bool,
}

impl ComparisonMetadataCursor {
    pub(super) fn open(spool: &ComparisonContigSpool) -> Result<Self> {
        Ok(Self {
            reader: vcf::open_validated_vcf(&spool.path)?,
            retained: VecDeque::new(),
            pending: None,
            furthest_start: 0,
        })
    }

    pub(super) fn collect(
        &mut self,
        chrom: &str,
        start: usize,
        end: usize,
        filtered_truth_keys: &mut BTreeSet<VariantKey>,
        decorations: &mut DecorationIndex,
        options: CollectOptions<'_>,
    ) -> Result<()> {
        let CollectOptions {
            preserve_info,
            roc_field,
            collect_filtered,
        } = options;
        if start.saturating_add(COMPARISON_METADATA_LOOKBEHIND) < self.furthest_start {
            bail!(
                "comparison metadata range for {chrom}:{start}-{end} moved more than {} bases behind the streaming cursor",
                COMPARISON_METADATA_LOOKBEHIND
            );
        }
        self.furthest_start = self.furthest_start.max(start);
        let retain_from = self
            .furthest_start
            .saturating_sub(COMPARISON_METADATA_LOOKBEHIND);
        while self
            .retained
            .front()
            .is_some_and(|record| record.raw().pos < retain_from)
        {
            self.retained.pop_front();
        }

        loop {
            let record = match self.pending.take() {
                Some(record) => record,
                None => match self.reader.next() {
                    Some(record) => record?,
                    None => break,
                },
            };
            if record.raw().pos > end {
                self.pending = Some(record);
                break;
            }
            if self.retained.len() >= MAX_COMPARISON_METADATA_RETAINED_RECORDS {
                bail!(
                    "comparison metadata window for {chrom}:{start}-{end} exceeds the {} retained-record resource limit",
                    MAX_COMPARISON_METADATA_RETAINED_RECORDS
                );
            }
            self.retained.push_back(record);
        }

        for record in self
            .retained
            .iter()
            .map(ValidatedVcfRecord::raw)
            .filter(|record| record.pos >= start && record.pos <= end)
        {
            if collect_filtered && !record.is_pass() {
                filtered_truth_keys.insert(VariantKey {
                    chrom: chrom.to_string(),
                    pos: record.pos,
                    ref_allele: record.ref_allele.clone(),
                    alt_allele: record.alt_allele.clone(),
                });
            }
            decorations.observe(record, preserve_info, roc_field);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{ComparisonRowSpool, legacy_same_key_rank};
    use crate::application::compare::AnnotatedRow;
    use crate::domain::{RawVcfRecord, SortKey};

    fn row(
        reference: &str,
        alternate: &str,
        side_rank: usize,
        samples: Vec<String>,
    ) -> AnnotatedRow {
        let raw = RawVcfRecord {
            chrom: "chr1".to_string(),
            pos: 25,
            id: ".".to_string(),
            ref_allele: reference.to_string(),
            alt_allele: alternate.to_string(),
            qual: "60".to_string(),
            filter: ".".to_string(),
            info: ".".to_string(),
            format: Some("GT:BD:BK:BI:BVT:BLT:QQ".to_string()),
            samples,
            mixed_edit_primitive: false,
            primitive_identity: None,
        };
        AnnotatedRow {
            sort_key: SortKey::new("chr1".to_string(), 25, 1, side_rank),
            record: raw.into(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }
    }

    #[test]
    fn duplicate_keys_order_truth_then_query_homalt_then_query_het() {
        let truth_only = vec![
            "0/1:UNK:lm:i6_15:INDEL:het:.".to_string(),
            "./.:.:.:.:NOCALL:nocall:0".to_string(),
        ];
        let query_homalt = vec![
            "./.:.:.:.:NOCALL:nocall:.".to_string(),
            "1/1:UNK:lm:i6_15:INDEL:homalt:0".to_string(),
        ];
        let query_het = vec![
            "./.:.:.:.:NOCALL:nocall:.".to_string(),
            "1/0:UNK:.:ti:SNP:het:0".to_string(),
        ];

        assert_eq!(legacy_same_key_rank(&truth_only), 0);
        assert_eq!(legacy_same_key_rank(&query_homalt), 1);
        assert_eq!(legacy_same_key_rank(&query_het), 2);
    }

    #[test]
    fn production_spool_preserves_truth_before_query_row_rank() {
        let query = row(
            "A",
            "G",
            1,
            vec![
                "./.:.:.:.:NOCALL:nocall:.".to_string(),
                "0/1:TP:gm:ti:SNP:het:55".to_string(),
            ],
        );
        let truth = row(
            "ACG",
            "GTA",
            0,
            vec![
                "0/1:TP:gm:ti:SNP:het:45".to_string(),
                "./.:.:.:.:NOCALL:nocall:0".to_string(),
            ],
        );
        let mut spool = ComparisonRowSpool::new();
        spool
            .push(query, false, "query".to_string())
            .expect("query row spools");
        spool
            .push(truth, false, "truth".to_string())
            .expect("truth row spools");

        let rows = spool
            .finish()
            .expect("spool finishes")
            .rows()
            .expect("spool opens")
            .collect::<Result<Vec<_>, _>>()
            .expect("rows decode");

        assert_eq!(rows[0].sort_key.type_rank, 0);
        assert_eq!(rows[0].record.raw().ref_allele, "ACG");
        assert_eq!(rows[1].sort_key.type_rank, 1);
        assert_eq!(rows[1].record.raw().ref_allele, "A");
    }
}
