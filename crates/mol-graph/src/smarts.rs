use std::collections::HashMap;

use crate::error::MolError;
use crate::mol::{BondOrder, MolAtom, MolGraph};

/// Compiled clean-room SMARTS query graph supporting atom/bond primitives, boolean OR `,`,
/// negation `!`, ring constraints `R`/`rN`, degree `D`, total H `H`, and charge.
#[derive(Debug, Clone)]
pub struct SmartsQuery {
    nodes: Vec<SmartsNode>,
    edges: Vec<SmartsEdge>,
    adj: Vec<Vec<(usize, usize)>>,
}

#[derive(Debug, Clone)]
struct SmartsNode {
    alternatives: Vec<Vec<AtomPred>>, // comma-separated OR clauses of AND primitives
}

#[derive(Debug, Clone)]
enum AtomPred {
    Negated(Box<AtomPred>),
    Any,
    AnyAromatic,
    AnyAliphatic,
    Element { sym: String, aromatic: Option<bool> },
    AtomicNum(u8),
    InRing(Option<bool>),
    RingSize(usize),
    TotalH(u8),
    Degree(usize),
    Charge(i8),
}

#[derive(Debug, Clone)]
struct SmartsEdge {
    cond: BondCond,
}

#[derive(Debug, Clone, Copy)]
enum BondCond {
    Any,
    SingleOrAromatic,
    Order(BondOrder),
    InRing,
}

impl SmartsQuery {
    /// Parse a SMARTS pattern string cleanly.
    pub fn parse(pattern: &str) -> Result<Self, MolError> {
        let s = pattern.trim();
        if s.is_empty() {
            return Err(MolError::InvalidSmarts {
                pos: 0,
                msg: "empty SMARTS pattern".to_string(),
            });
        }
        let bytes = s.as_bytes();
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut adj: Vec<Vec<(usize, usize)>> = Vec::new();
        let mut stack: Vec<usize> = Vec::new();
        let mut current: Option<usize> = None;
        let mut pending_bond: Option<BondCond> = None;
        let mut rings: HashMap<u16, (usize, Option<BondCond>, usize)> = HashMap::new();
        let mut i = 0usize;

        while i < bytes.len() {
            let ch = bytes[i];
            match ch {
                b'(' => {
                    let cur = current.ok_or(MolError::InvalidSmarts {
                        pos: i,
                        msg: "branch without preceding atom".to_string(),
                    })?;
                    stack.push(cur);
                    i += 1;
                }
                b')' => {
                    current = Some(stack.pop().ok_or(MolError::InvalidSmarts {
                        pos: i,
                        msg: "unmatched ')'".to_string(),
                    })?);
                    i += 1;
                }
                b'-' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Single));
                    i += 1;
                }
                b'=' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Double));
                    i += 1;
                }
                b'#' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Triple));
                    i += 1;
                }
                b':' => {
                    pending_bond = Some(BondCond::Order(BondOrder::Aromatic));
                    i += 1;
                }
                b'~' => {
                    pending_bond = Some(BondCond::Any);
                    i += 1;
                }
                b'@' => {
                    pending_bond = Some(BondCond::InRing);
                    i += 1;
                }
                b'0'..=b'9' => {
                    let cur = current.ok_or(MolError::InvalidSmarts {
                        pos: i,
                        msg: "ring digit without atom".to_string(),
                    })?;
                    let rnum = u16::from(ch - b'0');
                    if let Some((other, bo1, _)) = rings.remove(&rnum) {
                        let cond = bo1.or(pending_bond).unwrap_or(BondCond::SingleOrAromatic);
                        let eid = edges.len();
                        edges.push(SmartsEdge { cond });
                        adj[other].push((cur, eid));
                        adj[cur].push((other, eid));
                    } else {
                        rings.insert(rnum, (cur, pending_bond, i));
                    }
                    pending_bond = None;
                    i += 1;
                }
                b'[' => {
                    let close = bytes[i + 1..]
                        .iter()
                        .position(|&b| b == b']')
                        .map(|o| i + 1 + o)
                        .ok_or(MolError::InvalidSmarts {
                            pos: i,
                            msg: "unclosed '['".to_string(),
                        })?;
                    let inner = std::str::from_utf8(&bytes[i + 1..close]).map_err(|_| {
                        MolError::InvalidSmarts {
                            pos: i,
                            msg: "invalid UTF-8".to_string(),
                        }
                    })?;
                    let node = parse_smarts_bracket(inner, i + 1)?;
                    let nid = nodes.len();
                    nodes.push(node);
                    adj.push(Vec::new());
                    if let Some(prev) = current {
                        let cond = pending_bond.unwrap_or(BondCond::SingleOrAromatic);
                        let eid = edges.len();
                        edges.push(SmartsEdge { cond });
                        adj[prev].push((nid, eid));
                        adj[nid].push((prev, eid));
                    }
                    current = Some(nid);
                    pending_bond = None;
                    i = close + 1;
                }
                b'*' | b'A' | b'a' | b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I'
                | b'c' | b'n' | b'o' | b's' | b'p' => {
                    let (pred, step) = if ch == b'C' && bytes.get(i + 1) == Some(&b'l') {
                        (
                            AtomPred::Element {
                                sym: "Cl".to_string(),
                                aromatic: Some(false),
                            },
                            2,
                        )
                    } else if ch == b'B' && bytes.get(i + 1) == Some(&b'r') {
                        (
                            AtomPred::Element {
                                sym: "Br".to_string(),
                                aromatic: Some(false),
                            },
                            2,
                        )
                    } else {
                        let p = match ch {
                            b'*' => AtomPred::Any,
                            b'A' => AtomPred::AnyAliphatic,
                            b'a' => AtomPred::AnyAromatic,
                            b'c' => AtomPred::Element {
                                sym: "C".to_string(),
                                aromatic: Some(true),
                            },
                            b'n' => AtomPred::Element {
                                sym: "N".to_string(),
                                aromatic: Some(true),
                            },
                            b'o' => AtomPred::Element {
                                sym: "O".to_string(),
                                aromatic: Some(true),
                            },
                            b's' => AtomPred::Element {
                                sym: "S".to_string(),
                                aromatic: Some(true),
                            },
                            b'p' => AtomPred::Element {
                                sym: "P".to_string(),
                                aromatic: Some(true),
                            },
                            _ => AtomPred::Element {
                                sym: (ch as char).to_string(),
                                aromatic: Some(false),
                            },
                        };
                        (p, 1)
                    };
                    let nid = nodes.len();
                    nodes.push(SmartsNode {
                        alternatives: vec![vec![pred]],
                    });
                    adj.push(Vec::new());
                    if let Some(prev) = current {
                        let cond = pending_bond.unwrap_or(BondCond::SingleOrAromatic);
                        let eid = edges.len();
                        edges.push(SmartsEdge { cond });
                        adj[prev].push((nid, eid));
                        adj[nid].push((prev, eid));
                    }
                    current = Some(nid);
                    pending_bond = None;
                    i += step;
                }
                _ => {
                    return Err(MolError::InvalidSmarts {
                        pos: i,
                        msg: format!("unsupported SMARTS token '{}'", ch as char),
                    });
                }
            }
        }

        if let Some((rnum, (_, _, pos))) = rings.into_iter().next() {
            return Err(MolError::InvalidSmarts {
                pos,
                msg: format!("unclosed SMARTS ring {rnum}"),
            });
        }
        Ok(Self { nodes, edges, adj })
    }

    /// True if `mol` contains at least one match of this SMARTS pattern.
    pub fn matches(&self, mol: &MolGraph) -> bool {
        !self.find_matches(mol, true).is_empty()
    }

    /// Return atom mappings `[target_atom_idx_for_q0, ...]` matching this pattern.
    pub fn find_all_matches(&self, mol: &MolGraph) -> Vec<Vec<usize>> {
        self.find_matches(mol, false)
    }

    fn find_matches(&self, mol: &MolGraph, stop_at_first: bool) -> Vec<Vec<usize>> {
        if self.nodes.is_empty() || self.nodes.len() > mol.atoms.len() {
            return Vec::new();
        }
        let mut results = Vec::new();
        let mut mapping = vec![usize::MAX; self.nodes.len()];
        let mut used = vec![false; mol.atoms.len()];

        for start in 0..mol.atoms.len() {
            if !node_matches(&self.nodes[0], &mol.atoms[start], mol) {
                continue;
            }
            mapping[0] = start;
            used[start] = true;
            self.backtrack(1, mol, &mut mapping, &mut used, &mut results, stop_at_first);
            used[start] = false;
            mapping[0] = usize::MAX;
            if stop_at_first && !results.is_empty() {
                break;
            }
        }
        results
    }

    fn backtrack(
        &self,
        q_idx: usize,
        mol: &MolGraph,
        mapping: &mut [usize],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
        stop_at_first: bool,
    ) {
        if q_idx == self.nodes.len() {
            results.push(mapping.to_vec());
            return;
        }
        // Restrict candidate target atoms using already-mapped neighbors of q_idx
        let prev_req: Vec<(usize, BondCond)> = self.adj[q_idx]
            .iter()
            .filter_map(|&(qn, eid)| {
                if qn < q_idx {
                    Some((mapping[qn], self.edges[eid].cond))
                } else {
                    None
                }
            })
            .collect();

        let candidates: Vec<usize> = if let Some(&(first_t, _)) = prev_req.first() {
            mol.adj[first_t].iter().map(|&(nbr, _)| nbr).collect()
        } else {
            (0..mol.atoms.len()).collect()
        };

        for t_idx in candidates {
            if used[t_idx] || !node_matches(&self.nodes[q_idx], &mol.atoms[t_idx], mol) {
                continue;
            }
            let all_edges_ok = prev_req.iter().all(|&(t_nbr, cond)| {
                mol.adj[t_idx].iter().any(|&(cand_nbr, bid)| {
                    cand_nbr == t_nbr
                        && bond_matches(cond, mol.bonds[bid].order, mol.bonds[bid].in_ring)
                })
            });
            if !all_edges_ok {
                continue;
            }
            mapping[q_idx] = t_idx;
            used[t_idx] = true;
            self.backtrack(q_idx + 1, mol, mapping, used, results, stop_at_first);
            used[t_idx] = false;
            mapping[q_idx] = usize::MAX;
            if stop_at_first && !results.is_empty() {
                return;
            }
        }
    }
}

fn parse_smarts_bracket(inner: &str, base_pos: usize) -> Result<SmartsNode, MolError> {
    let mut alternatives = Vec::new();
    for part in inner.split(',') {
        let mut preds = Vec::new();
        let bytes = part.as_bytes();
        let mut i = 0usize;
        while i < bytes.len() {
            let neg = bytes[i] == b'!';
            if neg {
                i += 1;
            }
            if bytes.get(i) == Some(&b';') || bytes.get(i) == Some(&b'&') {
                i += 1;
                continue;
            }
            let (pred, next_i) = parse_single_pred(bytes, i, base_pos + i)?;
            preds.push(if neg {
                AtomPred::Negated(Box::new(pred))
            } else {
                pred
            });
            i = next_i;
            while bytes.get(i) == Some(&b';') || bytes.get(i) == Some(&b'&') {
                i += 1;
            }
        }
        alternatives.push(preds);
    }
    Ok(SmartsNode { alternatives })
}

fn parse_single_pred(
    bytes: &[u8],
    pos: usize,
    err_pos: usize,
) -> Result<(AtomPred, usize), MolError> {
    if pos >= bytes.len() {
        return Err(MolError::InvalidSmarts {
            pos: err_pos,
            msg: "dangling negation or operator".to_string(),
        });
    }
    let ch = bytes[pos];
    match ch {
        b'*' => Ok((AtomPred::Any, pos + 1)),
        b'A' => Ok((AtomPred::AnyAliphatic, pos + 1)),
        b'a' => Ok((AtomPred::AnyAromatic, pos + 1)),
        b'#' => {
            let (z, next_i) = read_u32(bytes, pos + 1).ok_or(MolError::InvalidSmarts {
                pos: err_pos,
                msg: "expected atomic number after '#'".to_string(),
            })?;
            Ok((AtomPred::AtomicNum(z as u8), next_i))
        }
        b'R' => {
            if let Some((_n, next_i)) = read_u32(bytes, pos + 1) {
                Ok((AtomPred::InRing(Some(true)), next_i))
            } else {
                Ok((AtomPred::InRing(Some(true)), pos + 1))
            }
        }
        b'r' => {
            if let Some((n, next_i)) = read_u32(bytes, pos + 1) {
                Ok((AtomPred::RingSize(n as usize), next_i))
            } else {
                Ok((AtomPred::InRing(Some(true)), pos + 1))
            }
        }
        b'D' => {
            let (d, next_i) = read_u32(bytes, pos + 1).unwrap_or((1, pos + 1));
            Ok((AtomPred::Degree(d as usize), next_i))
        }
        b'H' => {
            let (h, next_i) = read_u32(bytes, pos + 1).unwrap_or((1, pos + 1));
            Ok((AtomPred::TotalH(h as u8), next_i))
        }
        b'+' => {
            let (c, next_i) = read_u32(bytes, pos + 1).unwrap_or((1, pos + 1));
            Ok((AtomPred::Charge(c as i8), next_i))
        }
        b'-' => {
            let (c, next_i) = read_u32(bytes, pos + 1).unwrap_or((1, pos + 1));
            Ok((AtomPred::Charge(-(c as i8)), next_i))
        }
        _ => {
            if ch.is_ascii_alphabetic() {
                let arom = Some(ch.is_ascii_lowercase());
                let sym = if let Some(&sec) = bytes.get(pos + 1)
                    && sec.is_ascii_lowercase()
                    && crate::element::element_by_symbol(&format!(
                        "{}{}",
                        ch.to_ascii_uppercase() as char,
                        sec as char
                    ))
                    .is_some()
                {
                    (
                        format!("{}{}", ch.to_ascii_uppercase() as char, sec as char),
                        pos + 2,
                    )
                } else {
                    ((ch.to_ascii_uppercase() as char).to_string(), pos + 1)
                };
                Ok((
                    AtomPred::Element {
                        sym: sym.0,
                        aromatic: arom,
                    },
                    sym.1,
                ))
            } else {
                Err(MolError::InvalidSmarts {
                    pos: err_pos,
                    msg: format!("unknown SMARTS atom primitive '{}'", ch as char),
                })
            }
        }
    }
}

fn read_u32(bytes: &[u8], start: usize) -> Option<(u32, usize)> {
    let mut i = start;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == start {
        return None;
    }
    let val = std::str::from_utf8(&bytes[start..i])
        .ok()?
        .parse::<u32>()
        .ok()?;
    Some((val, i))
}

fn node_matches(node: &SmartsNode, atom: &MolAtom, mol: &MolGraph) -> bool {
    node.alternatives
        .iter()
        .any(|conj| conj.iter().all(|p| pred_matches(p, atom, mol)))
}

fn pred_matches(pred: &AtomPred, atom: &MolAtom, mol: &MolGraph) -> bool {
    match pred {
        AtomPred::Negated(inner) => !pred_matches(inner, atom, mol),
        AtomPred::Any => true,
        AtomPred::AnyAromatic => atom.aromatic,
        AtomPred::AnyAliphatic => !atom.aromatic,
        AtomPred::Element { sym, aromatic } => {
            atom.symbol == *sym && aromatic.is_none_or(|ar| atom.aromatic == ar)
        }
        AtomPred::AtomicNum(z) => atom.atomic_number == *z,
        AtomPred::InRing(flag) => flag.map_or(atom.in_ring, |req| atom.in_ring == req),
        AtomPred::RingSize(sz) => atom.ring_sizes.contains(sz),
        AtomPred::TotalH(h) => atom.total_h() == *h,
        AtomPred::Degree(d) => mol.adj[atom.id].len() == *d,
        AtomPred::Charge(c) => atom.formal_charge == *c,
    }
}

fn bond_matches(cond: BondCond, order: BondOrder, in_ring: bool) -> bool {
    match cond {
        BondCond::Any => true,
        BondCond::SingleOrAromatic => {
            matches!(
                order,
                BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative
            )
        }
        BondCond::Order(req) => order == req,
        BondCond::InRing => in_ring,
    }
}
