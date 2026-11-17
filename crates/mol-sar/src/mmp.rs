use std::collections::{BTreeMap, BTreeSet, VecDeque};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use mol_graph::{BondOrder, BondStereo, MolAtom, MolGraph};

use crate::table::SarRow;

/// Detected Matched Molecular Pair (`R1 -> R2` around a shared anchored core `[*]...`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MatchedMolecularPair {
    pub mol_a_name: String,
    pub mol_b_name: String,
    pub core_smiles: String,
    pub r_group_a: String,
    pub r_group_b: String,
    pub delta_pic50: Option<f64>,
    pub delta_clogp: Option<f64>,
    pub delta_tpsa: Option<f64>,
    pub delta_mw: f64,
}

/// Detected Activity Cliff pair ranked by Structure-Activity Landscape Index (`SALI`),
/// excluding identical/duplicate structures.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ActivityCliff {
    pub row_a_idx: usize,
    pub row_b_idx: usize,
    pub mol_a_name: String,
    pub mol_b_name: String,
    pub tanimoto: f64,
    pub delta_pic50: f64,
    pub sali: f64,
}

/// Discordant duplicate compound pair detected via canonical `structural_key`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DuplicatePair {
    pub structural_key: String,
    pub mol_a_name: String,
    pub mol_b_name: String,
    pub delta_pic50: Option<f64>,
}

/// Detect duplicate structures (including Kekulé vs aromatic representations).
pub fn find_duplicate_compounds(rows: &[SarRow]) -> Vec<DuplicatePair> {
    let mut by_key: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (idx, r) in rows.iter().enumerate() {
        by_key.entry(&r.structural_key).or_default().push(idx);
    }
    let mut out = Vec::new();
    for (key, idxs) in by_key {
        if idxs.len() < 2 {
            continue;
        }
        for i in 0..idxs.len() {
            for j in (i + 1)..idxs.len() {
                let ra = &rows[idxs[i]];
                let rb = &rows[idxs[j]];
                let dp = match (ra.pic50, rb.pic50) {
                    (Some(pa), Some(pb)) => Some((pa - pb).abs()),
                    _ => None,
                };
                out.push(DuplicatePair {
                    structural_key: key.to_string(),
                    mol_a_name: ra.name.clone(),
                    mol_b_name: rb.name.clone(),
                    delta_pic50: dp,
                });
            }
        }
    }
    out
}

/// Detect Activity Cliffs across unique molecules (excluding duplicates) deterministically.
pub fn find_activity_cliffs(rows: &[SarRow], min_tanimoto: f64) -> Vec<ActivityCliff> {
    let mut cliffs: Vec<ActivityCliff> = (0..rows.len())
        .into_par_iter()
        .flat_map(|i| {
            let mut local = Vec::new();
            let Some(p_i) = rows[i].pic50 else {
                return local;
            };
            for j in (i + 1)..rows.len() {
                if rows[i].structural_key == rows[j].structural_key {
                    // Never rank duplicate compounds as an Activity Cliff
                    continue;
                }
                let Some(p_j) = rows[j].pic50 else {
                    continue;
                };
                let sim = rows[i].fp.tanimoto(&rows[j].fp);
                if sim >= min_tanimoto {
                    let dp = (p_i - p_j).abs();
                    if dp > 0.3 {
                        let denom = (1.0 - sim).max(0.01);
                        local.push(ActivityCliff {
                            row_a_idx: i,
                            row_b_idx: j,
                            mol_a_name: rows[i].name.clone(),
                            mol_b_name: rows[j].name.clone(),
                            tanimoto: sim,
                            delta_pic50: dp,
                            sali: dp / denom,
                        });
                    }
                }
            }
            local
        })
        .collect();

    cliffs.sort_by(|a, b| {
        b.sali
            .total_cmp(&a.sali)
            .then_with(|| a.mol_a_name.cmp(&b.mol_a_name))
            .then_with(|| a.mol_b_name.cmp(&b.mol_b_name))
    });
    cliffs.truncate(100);
    cliffs
}

/// Reproducible single-cut Matched Molecular Pairs (`BTreeMap` + attachment dummy `[*]`),
/// rejecting duplicate molecules and position-isomer false positives.
pub fn find_matched_molecular_pairs(rows: &[SarRow]) -> Vec<MatchedMolecularPair> {
    let mut index: BTreeMap<String, Vec<(usize, String)>> = BTreeMap::new();

    for (r_idx, row) in rows.iter().enumerate() {
        for (core, rgrp) in single_bond_cuts_anchored(&row.mol) {
            index.entry(core).or_default().push((r_idx, rgrp));
        }
    }

    let mut seen_pairs: BTreeSet<(usize, usize, String)> = BTreeSet::new();
    let mut out = Vec::new();
    for (core, hits) in index {
        if hits.len() < 2 {
            continue;
        }
        for i in 0..hits.len() {
            for j in (i + 1)..hits.len() {
                let (a_idx, ref ra) = hits[i];
                let (b_idx, ref rb) = hits[j];
                if a_idx == b_idx
                    || ra == rb
                    || rows[a_idx].structural_key == rows[b_idx].structural_key
                {
                    continue;
                }
                let pair_key = if a_idx < b_idx {
                    (a_idx, b_idx, core.clone())
                } else {
                    (b_idx, a_idx, core.clone())
                };
                if !seen_pairs.insert(pair_key) {
                    continue;
                }
                let row_a = &rows[a_idx];
                let row_b = &rows[b_idx];
                let delta_pic50 = match (row_a.pic50, row_b.pic50) {
                    (Some(pa), Some(pb)) => Some(pb - pa),
                    _ => None,
                };
                let (delta_clogp, delta_tpsa) = match (&row_a.admet, &row_b.admet) {
                    (Some(aa), Some(ab)) => (
                        Some(ab.consensus_logp - aa.consensus_logp),
                        Some(ab.tpsa - aa.tpsa),
                    ),
                    _ => (None, None),
                };
                out.push(MatchedMolecularPair {
                    mol_a_name: row_a.name.clone(),
                    mol_b_name: row_b.name.clone(),
                    core_smiles: core.clone(),
                    r_group_a: ra.clone(),
                    r_group_b: rb.clone(),
                    delta_pic50,
                    delta_clogp,
                    delta_tpsa,
                    delta_mw: row_b.mw - row_a.mw,
                });
            }
        }
    }

    out.sort_by(|a, b| {
        let abs_a = a.delta_pic50.map_or(0.0, f64::abs);
        let abs_b = b.delta_pic50.map_or(0.0, f64::abs);
        abs_b
            .total_cmp(&abs_a)
            .then_with(|| a.core_smiles.cmp(&b.core_smiles))
            .then_with(|| a.r_group_a.cmp(&b.r_group_a))
            .then_with(|| a.r_group_b.cmp(&b.r_group_b))
            .then_with(|| a.mol_a_name.cmp(&b.mol_a_name))
            .then_with(|| a.mol_b_name.cmp(&b.mol_b_name))
    });
    out.truncate(150);
    out
}

fn single_bond_cuts_anchored(mol: &MolGraph) -> Vec<(String, String)> {
    let total_heavy = mol.heavy_atom_count();
    if total_heavy < 6 {
        return Vec::new();
    }
    let mut res = Vec::new();
    for bid in 0..mol.bonds.len() {
        let b = &mol.bonds[bid];
        if b.in_ring || b.order != BondOrder::Single {
            continue;
        }
        let comp_u = component_without_bond(mol, b.u, bid);
        let comp_v: BTreeSet<usize> = (0..mol.atoms.len())
            .filter(|id| !comp_u.contains(id))
            .collect();
        if comp_u.is_empty() || comp_v.is_empty() {
            continue;
        }
        let (core_nodes, core_cut, r_nodes, r_cut) = if comp_u.len() >= comp_v.len() {
            (&comp_u, b.u, &comp_v, b.v)
        } else {
            (&comp_v, b.v, &comp_u, b.u)
        };
        if core_nodes.len() < 5 || r_nodes.len() > (total_heavy / 2).max(1) {
            continue;
        }
        if let (Some(c_smi), Some(r_smi)) = (
            induced_canonical_with_anchor(mol, core_nodes, core_cut),
            induced_canonical_with_anchor(mol, r_nodes, r_cut),
        ) {
            res.push((c_smi, r_smi));
        }
    }
    res
}

fn component_without_bond(mol: &MolGraph, start: usize, ignored_bid: usize) -> BTreeSet<usize> {
    let mut comp = BTreeSet::new();
    let mut q = VecDeque::new();
    comp.insert(start);
    q.push_back(start);
    while let Some(u) = q.pop_front() {
        for &(v, bid) in &mol.adj[u] {
            if bid != ignored_bid && comp.insert(v) {
                q.push_back(v);
            }
        }
    }
    comp
}

fn induced_canonical_with_anchor(
    mol: &MolGraph,
    nodes: &BTreeSet<usize>,
    cut_atom: usize,
) -> Option<String> {
    let mut sub = MolGraph::new();
    let mut remap = vec![None; mol.atoms.len()];
    for &old in nodes {
        let mut a = mol.atoms[old].clone();
        a.coords_2d = None;
        a.chiral = mol_graph::ChiralTag::None;
        a.chiral_h_slot = None;
        let nid = sub.add_atom(a);
        remap[old] = Some(nid);
    }
    for b in &mol.bonds {
        if let (Some(nu), Some(nv)) = (remap[b.u], remap[b.v]) {
            let _ = sub.add_bond(nu, nv, b.order, b.stereo);
        }
    }
    let anchor_id = remap[cut_atom]?;
    let dummy = MolAtom::new(0, "*").ok()?;
    let dummy_id = sub.add_atom(dummy);
    let _ = sub.add_bond(anchor_id, dummy_id, BondOrder::Single, BondStereo::None);
    sub.perceive().ok()?;
    Some(sub.canonical_smiles())
}
