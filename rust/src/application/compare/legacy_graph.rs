use super::genotype::parse_gt_alleles;
use crate::adapters::vcf::Variant;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum Color {
    Grey,
    Black,
    Red,
    Blue,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Kind {
    Source,
    Alternative,
    HomRef,
    Sink,
}

#[derive(Clone, Debug)]
struct Node {
    start: isize,
    end: isize,
    color: Color,
    kind: Kind,
    edit: Option<Edit>,
}

#[derive(Clone, Debug)]
struct Edit {
    start: isize,
    end: isize,
    alt: String,
}

#[derive(Clone, Debug)]
struct PathState {
    node: usize,
    color: Color,
    mask: u64,
    homs: usize,
    edits: Vec<Edit>,
    sequences_seen: BTreeSet<String>,
}

fn trimmed_edit(variant: &Variant, allele: usize) -> Option<Edit> {
    if allele == 0 {
        return None;
    }
    let alt = variant.key.alt_allele.split(',').nth(allele - 1)?;
    let reference = variant.key.ref_allele.as_bytes();
    let alternate = alt.as_bytes();
    let mut prefix = 0;
    while prefix < reference.len()
        && prefix < alternate.len()
        && reference[prefix].eq_ignore_ascii_case(&alternate[prefix])
    {
        prefix += 1;
    }
    let mut suffix = 0;
    while suffix < reference.len().saturating_sub(prefix)
        && suffix < alternate.len().saturating_sub(prefix)
        && reference[reference.len() - 1 - suffix]
            .eq_ignore_ascii_case(&alternate[alternate.len() - 1 - suffix])
    {
        suffix += 1;
    }
    let start = variant.key.pos as isize + prefix as isize;
    let end = variant.key.pos as isize + reference.len() as isize - suffix as isize - 1;
    Some(Edit {
        start,
        end,
        alt: String::from_utf8(alternate[prefix..alternate.len() - suffix].to_vec()).ok()?,
    })
}

/// Smallest trimmed edit start over the GT-selected alleles — the coordinate
/// `GraphReference::makeGraph` records after `trimLeft`/`trimRight`. Records
/// with no selected alternate allele keep their VCF position.
pub(super) fn selected_edit_start(variant: &Variant) -> isize {
    parse_gt_alleles(&variant.gt)
        .into_iter()
        .filter(|allele| *allele > 0)
        .filter_map(|allele| trimmed_edit(variant, allele))
        .map(|edit| edit.start)
        .min()
        .unwrap_or(variant.key.pos as isize)
}

fn variant_nodes(variant: &Variant) -> Vec<Node> {
    let alleles = parse_gt_alleles(&variant.gt);
    if alleles.len() != 2
        || alleles
            .iter()
            .any(|allele| *allele > variant.key.alt_allele.split(',').count())
    {
        return Vec::new();
    }
    if alleles[0] == alleles[1] {
        return trimmed_edit(variant, alleles[0])
            .map(|edit| {
                vec![Node {
                    start: edit.start,
                    end: edit.end,
                    color: Color::Black,
                    kind: Kind::Alternative,
                    edit: Some(edit),
                }]
            })
            .unwrap_or_default();
    }

    let colors = if variant.gt.contains('|') {
        [Color::Red, Color::Blue]
    } else {
        [Color::Grey, Color::Grey]
    };
    let mut nodes = Vec::with_capacity(2);
    for (index, allele) in alleles.iter().copied().enumerate() {
        if let Some(edit) = trimmed_edit(variant, allele) {
            nodes.push(Node {
                start: edit.start,
                end: edit.end,
                color: colors[index],
                kind: Kind::Alternative,
                edit: Some(edit),
            });
        } else {
            let other = alleles[index ^ 1];
            let Some(other_edit) = trimmed_edit(variant, other) else {
                continue;
            };
            let alt_start = other_edit.start;
            let alt_end = other_edit.end;
            let start = if alt_start > alt_end {
                alt_start.max(alt_end)
            } else {
                alt_start
            };
            nodes.push(Node {
                start,
                end: start - 1,
                color: colors[index],
                kind: Kind::HomRef,
                edit: None,
            });
        }
    }
    nodes
}

fn compatible(left: &Node, right: &Node) -> bool {
    left.kind == Kind::Source
        || (left.end < right.start
            && (left.color <= Color::Black
                || right.color <= Color::Black
                || left.color == right.color))
}

fn build_graph(variants: &[Variant]) -> (Vec<Node>, Vec<BTreeSet<usize>>) {
    let start = variants
        .first()
        .map(|variant| variant.key.pos as isize)
        .unwrap_or(1);
    let mut nodes = vec![Node {
        start: start - 1,
        end: start - 1,
        color: Color::Grey,
        kind: Kind::Source,
        edit: None,
    }];
    let mut adjacency = vec![BTreeSet::new()];
    let mut previous = vec![0usize];

    for variant in variants {
        let current_nodes = variant_nodes(variant);
        if current_nodes.is_empty() {
            continue;
        }
        let current = current_nodes
            .into_iter()
            .map(|node| {
                let index = nodes.len();
                nodes.push(node);
                adjacency.push(BTreeSet::new());
                index
            })
            .collect::<Vec<_>>();
        let mut next_previous = current.clone();
        let mut has_predecessor = vec![false; current.len()];
        for predecessor in previous.iter().copied() {
            let mut superseded = false;
            for (slot, current_index) in current.iter().copied().enumerate() {
                if compatible(&nodes[predecessor], &nodes[current_index]) {
                    adjacency[predecessor].insert(current_index);
                    has_predecessor[slot] = true;
                    superseded = true;
                }
            }
            if !superseded {
                next_previous.push(predecessor);
            }
        }
        previous = next_previous;

        for (slot, current_index) in current.iter().copied().enumerate() {
            if has_predecessor[slot] {
                continue;
            }
            if let Some(predecessor) = (0..current_index)
                .rev()
                .find(|candidate| compatible(&nodes[*candidate], &nodes[current_index]))
            {
                adjacency[predecessor].insert(current_index);
            } else {
                adjacency[0].insert(current_index);
            }
        }
    }

    let max_end = previous
        .iter()
        .map(|index| nodes[*index].end)
        .max()
        .unwrap_or(start);
    let sink = nodes.len();
    nodes.push(Node {
        start: max_end,
        end: max_end,
        color: Color::Grey,
        kind: Kind::Sink,
        edit: None,
    });
    adjacency.push(BTreeSet::new());
    for predecessor in previous {
        adjacency[predecessor].insert(sink);
    }
    (nodes, adjacency)
}

fn path_sequence(
    reference: &str,
    region_start: usize,
    region_end: usize,
    edits: &[Edit],
) -> Option<String> {
    let mut output = String::new();
    let mut cursor = region_start as isize;
    let region_end = region_end as isize;
    for edit in edits {
        if edit.start < cursor || edit.start > region_end + 1 || edit.end > region_end {
            return None;
        }
        if cursor < edit.start {
            output.push_str(
                std::str::from_utf8(
                    &reference.as_bytes()[(cursor - 1) as usize..(edit.start - 1) as usize],
                )
                .ok()?,
            );
        }
        output.push_str(&edit.alt);
        cursor = edit.end + 1;
    }
    if cursor <= region_end {
        output.push_str(
            std::str::from_utf8(&reference.as_bytes()[(cursor - 1) as usize..region_end as usize])
                .ok()?,
        );
    }
    Some(output.to_ascii_uppercase())
}

pub(super) fn signatures(
    variants: &[Variant],
    reference: &str,
    region_start: usize,
    region_end: usize,
    max_paths: usize,
) -> Option<BTreeSet<String>> {
    let unphased_hets = variants
        .iter()
        .filter(|variant| {
            let alleles = parse_gt_alleles(&variant.gt);
            !variant.gt.contains('|') && alleles.len() == 2 && alleles[0] != alleles[1]
        })
        .count();
    if unphased_hets >= usize::BITS as usize || (1usize << unphased_hets) > max_paths {
        return None;
    }

    let (nodes, adjacency) = build_graph(variants);
    let sink = nodes.len() - 1;
    let mut node_masks = vec![0u64; nodes.len()];
    let mut next_mask = 1u64;
    let mut homs = 0usize;
    for (index, node) in nodes.iter().enumerate() {
        if node.kind != Kind::Alternative {
            continue;
        }
        if node.color == Color::Black {
            homs += 1;
        } else {
            if next_mask == 0 {
                return None;
            }
            node_masks[index] = next_mask;
            next_mask <<= 1;
        }
    }
    let all_hets = next_mask.wrapping_sub(1);
    let reference_sequence =
        std::str::from_utf8(&reference.as_bytes()[region_start - 1..region_end])
            .ok()?
            .to_ascii_uppercase();
    let mut queue = VecDeque::from([PathState {
        node: 0,
        color: Color::Grey,
        mask: 0,
        homs: 0,
        edits: Vec::new(),
        sequences_seen: BTreeSet::from([reference_sequence.clone()]),
    }]);
    let mut paths = Vec::new();
    while let Some(state) = queue.pop_front() {
        for next in adjacency[state.node].iter().copied() {
            let node = &nodes[next];
            if node.color > Color::Black && state.color > Color::Black && node.color != state.color
            {
                continue;
            }
            let mut next_state = PathState {
                node: next,
                color: state.color.max(node.color),
                mask: state.mask | node_masks[next],
                homs: state.homs
                    + usize::from(node.kind == Kind::Alternative && node.color == Color::Black),
                edits: state.edits.clone(),
                sequences_seen: state.sequences_seen.clone(),
            };
            if let Some(edit) = &node.edit {
                next_state.edits.push(edit.clone());
                let Some(sequence) =
                    path_sequence(reference, region_start, region_end, &next_state.edits)
                else {
                    return None;
                };
                if next_state.mask != 0 && next_state.sequences_seen.contains(&sequence) {
                    continue;
                }
                if next_state.mask != 0 {
                    next_state.sequences_seen.insert(sequence);
                }
            }
            if next == sink {
                if next_state.homs == homs {
                    let Some(sequence) =
                        path_sequence(reference, region_start, region_end, &next_state.edits)
                    else {
                        return None;
                    };
                    paths.push((next_state.mask, sequence));
                    if paths.len() >= max_paths {
                        break;
                    }
                }
            } else {
                queue.push_back(next_state);
            }
        }
        if paths.len() >= max_paths {
            break;
        }
    }
    if paths.is_empty() {
        if all_hets == 0 && homs == 0 {
            return Some(BTreeSet::from([format!(
                "homref:{reference_sequence}|{reference_sequence}"
            )]));
        }
        return None;
    }
    if all_hets == 0 {
        return Some(
            paths
                .into_iter()
                .map(|(_, sequence)| format!("hom:{sequence}|{sequence}"))
                .collect(),
        );
    }
    let mut paths_by_mask = BTreeMap::<u64, BTreeSet<String>>::new();
    for (mask, sequence) in paths {
        paths_by_mask.entry(mask).or_default().insert(sequence);
    }
    let mut signatures = BTreeSet::new();
    for (left_mask, left_sequences) in &paths_by_mask {
        let opposite = (!left_mask) & all_hets;
        if let Some(right_sequences) = paths_by_mask.get(&opposite) {
            for left_sequence in left_sequences {
                for right_sequence in right_sequences {
                    let mut pair = [left_sequence, right_sequence];
                    pair.sort();
                    let kind = if pair[0] == &reference_sequence || pair[1] == &reference_sequence {
                        "het"
                    } else {
                        "hetalt"
                    };
                    signatures.insert(format!("{kind}:{}|{}", pair[0], pair[1]));
                }
            }
        }
    }
    (!signatures.is_empty()).then_some(signatures)
}
