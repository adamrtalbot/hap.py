//! Bounded-memory preprocessing spools and external record ordering.

use super::{LEFT_SHIFT_WINDOW, SymbolicDeletionMaterialization};
use crate::adapters::vcf::{ValidatedVcfReader, ValidatedVcfRecord};
use crate::domain::{PrimitiveIdentity, QueryProvenance, RawVcfRecord};
use anyhow::{Context, Result};
use std::cmp::Reverse;
use std::collections::VecDeque;
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

const PREPROCESS_SORT_CHUNK_RECORDS: usize = 65_536;
const PREPROCESS_SORT_MERGE_FAN_IN: usize = 32;
const PREPROCESS_SORT_CHUNK_VERSION: &str = "v2";

fn write_preprocess_sort_record(
    writer: &mut impl Write,
    rank: usize,
    pos: usize,
    serial: usize,
    record: &RawVcfRecord,
) -> Result<()> {
    let (identity_present, identity_start, identity_end, identity_alt) =
        match record.primitive_identity.as_ref() {
            Some(identity) => (1, identity.start, identity.end, identity.alt.as_str()),
            None => (0, 0, 0, "."),
        };
    writeln!(
        writer,
        "{PREPROCESS_SORT_CHUNK_VERSION}\t{rank}\t{pos}\t{serial}\t{}\t{identity_present}\t{identity_start}\t{identity_end}\t{identity_alt}\t{}",
        u8::from(record.mixed_edit_primitive),
        record.to_line()
    )?;
    Ok(())
}

pub(super) struct PreparedRecordSpool {
    path: tempfile::TempPath,
    offsets: Vec<u64>,
    contigs: Vec<String>,
}

impl PreparedRecordSpool {
    pub(super) fn reader(&self) -> Result<PreparedRecordReader<'_>> {
        Ok(PreparedRecordReader {
            reader: BufReader::new(File::open(&self.path)?),
            offsets: &self.offsets,
            next_index: None,
        })
    }

    pub(super) fn len(&self) -> Result<u64> {
        Ok(std::fs::metadata(&self.path)?.len())
    }

    pub(super) fn record_count(&self) -> usize {
        self.offsets.len()
    }

    pub(super) fn contigs(&self) -> &[String] {
        &self.contigs
    }
}

pub(super) struct PreparedRecordSpoolWriter {
    writer: BufWriter<tempfile::NamedTempFile>,
    offsets: Vec<u64>,
    encoded: Vec<u8>,
    bytes_written: u64,
    contig_set: HashSet<String>,
    contigs: Vec<String>,
}

impl PreparedRecordSpoolWriter {
    pub(super) fn new() -> Result<Self> {
        Ok(Self {
            writer: BufWriter::new(
                tempfile::NamedTempFile::new().context("failed to create prepared-record spool")?,
            ),
            offsets: Vec::new(),
            encoded: Vec::new(),
            bytes_written: 0,
            contig_set: HashSet::new(),
            contigs: Vec::new(),
        })
    }

    pub(super) fn push(
        &mut self,
        record: &RawVcfRecord,
        symbolic_deletion: Option<SymbolicDeletionMaterialization>,
    ) -> Result<()> {
        if self.contig_set.insert(record.chrom.clone()) {
            self.contigs.push(record.chrom.clone());
        }
        let marker = match symbolic_deletion {
            None => 0,
            Some(SymbolicDeletionMaterialization::LeadingAnchor) => 1,
            Some(SymbolicDeletionMaterialization::ContigStart) => 2,
        };
        self.encoded.clear();
        self.encoded.push(marker);
        encode_raw_record(&mut self.encoded, record)?;
        self.offsets.push(self.bytes_written);
        self.writer.write_all(&self.encoded)?;
        self.bytes_written = self
            .bytes_written
            .checked_add(u64::try_from(self.encoded.len())?)
            .context("prepared-record spool size overflow")?;
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<PreparedRecordSpool> {
        self.writer.flush()?;
        let file = self
            .writer
            .into_inner()
            .map_err(|error| error.into_error())?;
        Ok(PreparedRecordSpool {
            path: file.into_temp_path(),
            offsets: self.offsets,
            contigs: self.contigs,
        })
    }
}

pub(super) struct PreparedRecordReader<'a> {
    reader: BufReader<File>,
    offsets: &'a [u64],
    next_index: Option<usize>,
}

impl PreparedRecordReader<'_> {
    pub(super) fn read_at(
        &mut self,
        index: usize,
    ) -> Result<(RawVcfRecord, Option<SymbolicDeletionMaterialization>)> {
        let offset = self
            .offsets
            .get(index)
            .copied()
            .context("prepared-record spool index is out of bounds")?;
        if self.next_index != Some(index) {
            self.reader.seek(SeekFrom::Start(offset))?;
        }
        let marker = read_u8(&mut self.reader)?;
        let symbolic_deletion = match marker {
            0 => None,
            1 => Some(SymbolicDeletionMaterialization::LeadingAnchor),
            2 => Some(SymbolicDeletionMaterialization::ContigStart),
            _ => anyhow::bail!("invalid prepared-record spool marker {marker}"),
        };
        let record = decode_raw_record(&mut self.reader)?;
        self.next_index = index.checked_add(1);
        Ok((record, symbolic_deletion))
    }
}

fn encode_raw_record(output: &mut Vec<u8>, record: &RawVcfRecord) -> Result<()> {
    encode_string(output, &record.chrom)?;
    output.extend_from_slice(&u64::try_from(record.pos)?.to_le_bytes());
    encode_string(output, &record.id)?;
    encode_string(output, &record.ref_allele)?;
    encode_string(output, &record.alt_allele)?;
    encode_string(output, &record.qual)?;
    encode_string(output, &record.filter)?;
    encode_string(output, &record.info)?;
    match &record.format {
        Some(format) => {
            output.push(1);
            encode_string(output, format)?;
        }
        None => output.push(0),
    }
    output.extend_from_slice(&u32::try_from(record.samples.len())?.to_le_bytes());
    for sample in &record.samples {
        encode_string(output, sample)?;
    }
    Ok(())
}

fn encode_string(output: &mut Vec<u8>, value: &str) -> Result<()> {
    output.extend_from_slice(&u32::try_from(value.len())?.to_le_bytes());
    output.extend_from_slice(value.as_bytes());
    Ok(())
}

fn decode_raw_record(reader: &mut impl Read) -> Result<RawVcfRecord> {
    let chrom = decode_string(reader)?;
    let pos = usize::try_from(read_u64(reader)?)?;
    let id = decode_string(reader)?;
    let ref_allele = decode_string(reader)?;
    let alt_allele = decode_string(reader)?;
    let qual = decode_string(reader)?;
    let filter = decode_string(reader)?;
    let info = decode_string(reader)?;
    let format = match read_u8(reader)? {
        0 => None,
        1 => Some(decode_string(reader)?),
        marker => anyhow::bail!("invalid prepared-record FORMAT marker {marker}"),
    };
    let sample_count = usize::try_from(read_u32(reader)?)?;
    let mut samples = Vec::with_capacity(sample_count);
    for _ in 0..sample_count {
        samples.push(decode_string(reader)?);
    }
    Ok(RawVcfRecord {
        chrom,
        pos,
        id,
        ref_allele,
        alt_allele,
        qual,
        filter,
        info,
        format,
        samples,
        mixed_edit_primitive: false,
        primitive_identity: None,
    })
}

fn decode_string(reader: &mut impl Read) -> Result<String> {
    let len = usize::try_from(read_u32(reader)?)?;
    let mut value = vec![0; len];
    reader.read_exact(&mut value)?;
    String::from_utf8(value).context("prepared-record spool field is not UTF-8")
}

fn read_u8(reader: &mut impl Read) -> Result<u8> {
    let mut bytes = [0; 1];
    reader.read_exact(&mut bytes)?;
    Ok(bytes[0])
}

fn read_u32(reader: &mut impl Read) -> Result<u32> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(reader: &mut impl Read) -> Result<u64> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

#[derive(Default)]
struct StreamPositionState {
    maximum_position: usize,
    ordinals: BTreeMap<usize, usize>,
}

pub(super) struct PreprocessSpool {
    sorted: bool,
    unsorted: tempfile::NamedTempFile,
    chunks: Vec<tempfile::TempPath>,
    buffer: Vec<(usize, usize, usize, RawVcfRecord)>,
    stream_count: usize,
    stream_position_states: Vec<HashMap<String, StreamPositionState>>,
    contig_ranks: HashMap<String, usize>,
    emitted_contig_set: HashSet<String>,
    pub(super) emitted_contigs: Vec<String>,
    next_rank: usize,
    pub(super) serial: usize,
}

impl PreprocessSpool {
    pub(super) fn new(
        sorted: bool,
        stream_count: usize,
        declared_contigs: &[String],
    ) -> Result<Self> {
        let contig_ranks = declared_contigs
            .iter()
            .enumerate()
            .map(|(rank, contig)| (contig.clone(), rank))
            .collect::<HashMap<_, _>>();
        let next_rank = contig_ranks.len();
        Ok(Self {
            sorted,
            unsorted: tempfile::NamedTempFile::new()
                .context("failed to create preprocess spool")?,
            chunks: Vec::new(),
            buffer: Vec::new(),
            stream_count: stream_count.max(1),
            stream_position_states: (0..stream_count.max(1)).map(|_| HashMap::new()).collect(),
            contig_ranks,
            emitted_contig_set: HashSet::new(),
            emitted_contigs: Vec::new(),
            next_rank,
            serial: 0,
        })
    }

    pub(super) fn push(&mut self, record: ValidatedVcfRecord, stream_index: usize) -> Result<()> {
        let record = record.into_raw();
        if self.emitted_contig_set.insert(record.chrom.clone()) {
            self.emitted_contigs.push(record.chrom.clone());
        }
        self.seed_contig_rank(&record.chrom);
        if !self.sorted {
            writeln!(self.unsorted.as_file_mut(), "{}", record.to_line())?;
            self.serial += 1;
            return Ok(());
        }
        let rank = if let Some(rank) = self.contig_ranks.get(record.chrom.as_str()) {
            *rank
        } else {
            let rank = self.next_rank;
            self.next_rank += 1;
            self.contig_ranks.insert(record.chrom.clone(), rank);
            rank
        };
        let stream_states = self
            .stream_position_states
            .get_mut(stream_index)
            .context("preprocess stream index is out of bounds")?;
        if !stream_states.contains_key(record.chrom.as_str()) {
            stream_states.insert(record.chrom.clone(), StreamPositionState::default());
        }
        let position_state = stream_states
            .get_mut(record.chrom.as_str())
            .context("preprocess stream contig state was just created")?;
        position_state.maximum_position = position_state.maximum_position.max(record.pos);
        let minimum_retained_position = position_state
            .maximum_position
            .saturating_sub(LEFT_SHIFT_WINDOW);
        while position_state
            .ordinals
            .first_key_value()
            .is_some_and(|(position, _)| *position < minimum_retained_position)
        {
            position_state.ordinals.pop_first();
        }
        let ordinal = position_state.ordinals.entry(record.pos).or_default();
        let tie_break = ordinal
            .checked_mul(self.stream_count)
            .and_then(|value| value.checked_add(stream_index))
            .context("preprocess stream ordering overflow")?;
        *ordinal += 1;
        self.buffer.push((rank, record.pos, tie_break, record));
        self.serial += 1;
        if self.buffer.len() >= PREPROCESS_SORT_CHUNK_RECORDS {
            self.flush_chunk()?;
        }
        Ok(())
    }

    pub(super) fn seed_contig_rank(&mut self, contig: &str) {
        if !self.contig_ranks.contains_key(contig) {
            let rank = self.next_rank;
            self.next_rank += 1;
            self.contig_ranks.insert(contig.to_string(), rank);
        }
    }

    pub(super) fn sort_emitted_contigs(&mut self) {
        self.emitted_contigs
            .sort_by_key(|contig| self.contig_ranks.get(contig).copied().unwrap_or(usize::MAX));
    }

    #[cfg(test)]
    pub(super) fn retained_position_count(&self) -> usize {
        self.stream_position_states
            .iter()
            .flat_map(HashMap::values)
            .map(|state| state.ordinals.len())
            .sum()
    }

    fn flush_chunk(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer
            .sort_by_key(|(rank, pos, serial, _)| (*rank, *pos, *serial));
        let mut chunk =
            tempfile::NamedTempFile::new().context("failed to create preprocess sort chunk")?;
        for (rank, pos, serial, record) in self.buffer.drain(..) {
            write_preprocess_sort_record(chunk.as_file_mut(), rank, pos, serial, &record)?;
        }
        chunk.as_file_mut().flush()?;
        self.chunks.push(chunk.into_temp_path());
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<PreprocessRecords> {
        if !self.sorted {
            self.unsorted.as_file_mut().flush()?;
            let path = self.unsorted.into_temp_path();
            let reader = crate::adapters::vcf::open_validated_vcf(&path)?;
            return Ok(PreprocessRecords::Unsorted {
                _path: path,
                reader,
            });
        }
        self.flush_chunk()?;
        Ok(PreprocessRecords::Sorted(ExternalRecordMerge::new(
            self.chunks,
            self.stream_count,
        )?))
    }
}

pub(super) enum PreprocessRecords {
    Unsorted {
        _path: tempfile::TempPath,
        reader: ValidatedVcfReader,
    },
    Sorted(ExternalRecordMerge),
}

pub(super) struct LocationAggregatedRecords {
    inner: PreprocessRecords,
    enabled: bool,
    pending: Option<Result<ValidatedVcfRecord>>,
    ready: VecDeque<Result<ValidatedVcfRecord>>,
}

impl LocationAggregatedRecords {
    pub(super) fn new(inner: PreprocessRecords, enabled: bool) -> Self {
        Self {
            inner,
            enabled,
            pending: None,
            ready: VecDeque::new(),
        }
    }
}

impl Iterator for LocationAggregatedRecords {
    type Item = Result<ValidatedVcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(record) = self.ready.pop_front() {
            return Some(record);
        }
        if !self.enabled {
            return self.pending.take().or_else(|| self.inner.next());
        }

        let first = match self.pending.take().or_else(|| self.inner.next())? {
            Ok(record) => record,
            Err(error) => return Some(Err(error)),
        };
        let chrom = first.raw().chrom.clone();
        let pos = first.raw().pos;
        let mut group = vec![(first.provenance().source_index(), first.into_raw())];
        loop {
            match self.inner.next() {
                Some(Ok(record)) if record.raw().chrom == chrom && record.raw().pos == pos => {
                    group.push((record.provenance().source_index(), record.into_raw()));
                }
                Some(record) => {
                    self.pending = Some(record);
                    break;
                }
                None => break,
            }
        }

        let mut streams = BTreeMap::<usize, Vec<RawVcfRecord>>::new();
        for (stream, record) in group {
            streams.entry(stream.unwrap_or(0)).or_default().push(record);
        }
        let streams = streams
            .into_values()
            .map(crate::engines::variant_pipeline::aggregate_location_records)
            .collect::<Vec<_>>();
        // Legacy merges independent block/location streams round-wise at an
        // equal normalized position: first record from every stream, then
        // the second from every stream. Keeping stream provenance through
        // the indexed spool prevents duplicate location copies from being
        // mistaken for two alleles of one diploid call.
        let maximum_stream_records = streams.iter().map(Vec::len).max().unwrap_or(0);
        for record_index in 0..maximum_stream_records {
            for stream in &streams {
                if let Some(record) = stream.get(record_index) {
                    self.ready.push_back(ValidatedVcfRecord::try_from_raw(
                        record.clone(),
                        QueryProvenance::Unavailable,
                    ));
                }
            }
        }
        self.ready.pop_front()
    }
}

impl Iterator for PreprocessRecords {
    type Item = Result<ValidatedVcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Unsorted { reader, .. } => reader.next(),
            Self::Sorted(reader) => reader.next(),
        }
    }
}

type SortKey = (usize, usize, usize, usize);
type KeyedRecord = ((usize, usize, usize), ValidatedVcfRecord);

pub(super) struct ExternalRecordMerge {
    _chunks: Vec<tempfile::TempPath>,
    readers: Vec<BufReader<File>>,
    current: Vec<Option<KeyedRecord>>,
    heap: BinaryHeap<Reverse<SortKey>>,
    stream_count: usize,
}

impl ExternalRecordMerge {
    fn new(chunks: Vec<tempfile::TempPath>, stream_count: usize) -> Result<Self> {
        Self::open(
            collapse_preprocess_chunks(chunks, stream_count)?,
            stream_count,
        )
    }

    fn open(chunks: Vec<tempfile::TempPath>, stream_count: usize) -> Result<Self> {
        let readers = chunks
            .iter()
            .map(File::open)
            .map(|file| file.map(BufReader::new))
            .collect::<std::io::Result<Vec<_>>>()?;
        let current = vec![None; readers.len()];
        let mut merge = Self {
            _chunks: chunks,
            readers,
            current,
            heap: BinaryHeap::new(),
            stream_count: stream_count.max(1),
        };
        for index in 0..merge.readers.len() {
            merge.read_next(index)?;
        }
        Ok(merge)
    }

    fn read_next(&mut self, index: usize) -> Result<()> {
        let mut line = String::new();
        if self.readers[index].read_line(&mut line)? == 0 {
            self.current[index] = None;
            return Ok(());
        }
        let line = line.trim_end_matches(['\r', '\n']);
        let mut fields = line.splitn(10, '\t');
        let version = fields
            .next()
            .context("preprocess sort chunk lacks version")?;
        if version != PREPROCESS_SORT_CHUNK_VERSION {
            anyhow::bail!("unsupported preprocess sort chunk version {version}");
        }
        let rank = fields
            .next()
            .context("preprocess sort chunk lacks rank")?
            .parse()?;
        let pos = fields
            .next()
            .context("preprocess sort chunk lacks position")?
            .parse()?;
        let serial = fields
            .next()
            .context("preprocess sort chunk lacks serial")?
            .parse()?;
        let mixed_edit_primitive = match fields
            .next()
            .context("preprocess sort chunk lacks mixed-edit marker")?
        {
            "0" => false,
            "1" => true,
            marker => anyhow::bail!("invalid preprocess mixed-edit marker {marker}"),
        };
        let identity_present = fields
            .next()
            .context("preprocess sort chunk lacks primitive-identity marker")?;
        let identity_start = fields
            .next()
            .context("preprocess sort chunk lacks primitive-identity start")?;
        let identity_end = fields
            .next()
            .context("preprocess sort chunk lacks primitive-identity end")?;
        let identity_alt = fields
            .next()
            .context("preprocess sort chunk lacks primitive-identity ALT")?;
        let primitive_identity = match identity_present {
            "0" => None,
            "1" => Some(PrimitiveIdentity {
                start: identity_start.parse()?,
                end: identity_end.parse()?,
                alt: identity_alt.to_string(),
            }),
            marker => anyhow::bail!("invalid preprocess primitive-identity marker {marker}"),
        };
        let mut raw = RawVcfRecord::from_line(
            fields
                .next()
                .context("preprocess sort chunk lacks record")?,
            Path::new("preprocess-sort-chunk"),
        )?;
        raw.mixed_edit_primitive = mixed_edit_primitive;
        raw.primitive_identity = primitive_identity;
        let provenance = QueryProvenance::source(serial % self.stream_count, self.stream_count)?;
        let record = ValidatedVcfRecord::try_from_raw(raw, provenance)?;
        self.current[index] = Some(((rank, pos, serial), record));
        self.heap.push(Reverse((rank, pos, serial, index)));
        Ok(())
    }

    fn next_keyed(&mut self) -> Option<Result<KeyedRecord>> {
        let Reverse((_, _, _, index)) = self.heap.pop()?;
        let Some(entry) = self.current[index].take() else {
            return Some(Err(anyhow::anyhow!(
                "preprocess merge heap referenced an exhausted chunk"
            )));
        };
        if let Err(error) = self.read_next(index) {
            return Some(Err(error));
        }
        Some(Ok(entry))
    }
}

impl Iterator for ExternalRecordMerge {
    type Item = Result<ValidatedVcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_keyed()
            .map(|entry| entry.map(|(_, record)| record))
    }
}

fn collapse_preprocess_chunks(
    mut chunks: Vec<tempfile::TempPath>,
    stream_count: usize,
) -> Result<Vec<tempfile::TempPath>> {
    while chunks.len() > PREPROCESS_SORT_MERGE_FAN_IN {
        let mut merged = Vec::with_capacity(chunks.len().div_ceil(PREPROCESS_SORT_MERGE_FAN_IN));
        let mut remaining = chunks.into_iter();
        loop {
            let batch = remaining
                .by_ref()
                .take(PREPROCESS_SORT_MERGE_FAN_IN)
                .collect::<Vec<_>>();
            if batch.is_empty() {
                break;
            }
            let mut output = tempfile::NamedTempFile::new()
                .context("failed to create preprocess merge chunk")?;
            {
                let mut writer = BufWriter::new(output.as_file_mut());
                let mut merge = ExternalRecordMerge::open(batch, stream_count)?;
                while let Some(entry) = merge.next_keyed() {
                    let ((rank, pos, serial), record) = entry?;
                    write_preprocess_sort_record(&mut writer, rank, pos, serial, record.raw())?;
                }
                writer.flush()?;
            }
            merged.push(output.into_temp_path());
        }
        chunks = merged;
    }
    Ok(chunks)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(pos: usize, identity: Option<PrimitiveIdentity>) -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr1".to_string(),
            pos,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            info: ".".to_string(),
            format: Some("GT".to_string()),
            samples: vec!["0/1".to_string()],
            mixed_edit_primitive: false,
            primitive_identity: identity,
        }
    }

    #[test]
    fn legacy_only_primitive_identity_survives_multi_pass_chunk_rewrite() -> Result<()> {
        let mut chunks = Vec::new();
        for serial in 0..=PREPROCESS_SORT_MERGE_FAN_IN {
            let identity = (serial % 2 == 0).then(|| PrimitiveIdentity {
                start: serial + 1,
                end: serial,
                alt: format!("I{serial}"),
            });
            let mut chunk = tempfile::NamedTempFile::new()?;
            write_preprocess_sort_record(
                chunk.as_file_mut(),
                0,
                serial + 1,
                serial,
                &record(serial + 1, identity),
            )?;
            chunk.as_file_mut().flush()?;
            chunks.push(chunk.into_temp_path());
        }

        let collapsed = collapse_preprocess_chunks(chunks, 1)?;
        assert_eq!(collapsed.len(), 2, "33 chunks require one rewrite pass");
        let observed = ExternalRecordMerge::open(collapsed, 1)?
            .map(|entry| entry.map(|record| record.into_raw().primitive_identity))
            .collect::<Result<Vec<_>>>()?;

        assert_eq!(observed.len(), PREPROCESS_SORT_MERGE_FAN_IN + 1);
        for (serial, identity) in observed.into_iter().enumerate() {
            if serial % 2 == 0 {
                assert_eq!(
                    identity,
                    Some(PrimitiveIdentity {
                        start: serial + 1,
                        end: serial,
                        alt: format!("I{serial}"),
                    })
                );
            } else {
                assert!(identity.is_none());
            }
        }
        Ok(())
    }
}
