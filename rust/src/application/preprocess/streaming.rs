//! Bounded-memory preprocessing spools and external record ordering.

use super::LEFT_SHIFT_WINDOW;
use crate::adapters::vcf::{ValidatedVcfReader, ValidatedVcfRecord};
use crate::domain::{QueryProvenance, RawVcfRecord};
use anyhow::{Context, Result};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const PREPROCESS_SORT_CHUNK_RECORDS: usize = 65_536;
const PREPROCESS_SORT_MERGE_FAN_IN: usize = 32;

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
    pub(super) fn new(sorted: bool, stream_count: usize) -> Result<Self> {
        Ok(Self {
            sorted,
            unsorted: tempfile::NamedTempFile::new()
                .context("failed to create preprocess spool")?,
            chunks: Vec::new(),
            buffer: Vec::new(),
            stream_count: stream_count.max(1),
            stream_position_states: (0..stream_count.max(1)).map(|_| HashMap::new()).collect(),
            contig_ranks: HashMap::new(),
            emitted_contig_set: HashSet::new(),
            emitted_contigs: Vec::new(),
            next_rank: 0,
            serial: 0,
        })
    }

    pub(super) fn push(&mut self, record: ValidatedVcfRecord, stream_index: usize) -> Result<()> {
        let record = record.into_raw();
        if !self.emitted_contig_set.contains(record.chrom.as_str()) {
            self.emitted_contig_set.insert(record.chrom.clone());
            self.emitted_contigs.push(record.chrom.clone());
        }
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
            writeln!(
                chunk.as_file_mut(),
                "{rank}\t{pos}\t{serial}\t{}",
                record.to_line()
            )?;
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
}

impl ExternalRecordMerge {
    fn new(chunks: Vec<tempfile::TempPath>) -> Result<Self> {
        Self::open(collapse_preprocess_chunks(chunks)?)
    }

    fn open(chunks: Vec<tempfile::TempPath>) -> Result<Self> {
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
        let mut fields = line.splitn(4, '\t');
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
        let raw = RawVcfRecord::from_line(
            fields
                .next()
                .context("preprocess sort chunk lacks record")?,
            Path::new("preprocess-sort-chunk"),
        )?;
        let record = ValidatedVcfRecord::try_from_raw(raw, QueryProvenance::Unavailable)?;
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
                let mut merge = ExternalRecordMerge::open(batch)?;
                while let Some(entry) = merge.next_keyed() {
                    let ((rank, pos, serial), record) = entry?;
                    writeln!(
                        writer,
                        "{rank}\t{pos}\t{serial}\t{}",
                        record.raw().to_line()
                    )?;
                }
                writer.flush()?;
            }
            merged.push(output.into_temp_path());
        }
        chunks = merged;
    }
    Ok(chunks)
}
