use std::collections::{BTreeSet, HashMap, VecDeque};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use mol_graph::MolGraph;

use crate::table::SarRow;

/// Detected Matched Molecular Pair (`R1 -> R2` around a shared core).
#[derive(Debug, Clone, Serialize, Deserialize)]
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

/// Detected Activity Cliff pair ranked by Structure-Activity Landscape Index (`SALI`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ActivityCliff {
    pub row_a_idx: usize,
    pub row_b_idx: usize,
    pub mol_a_name: String,
    pub mol_b_name: String,
    pub tanimoto: f64,
    pub delta_pic50: f64,
    pub sali: f64,
}

/// Detect Activity Cliffs across rows with numeric `pIC50` values.
pub fn find_activity_cliffs(rows: &[SarRow], min_tanimoto: f64) -> Vec<ActivityCliff> {
    let mut cliffs: Vec<ActivityCliff> = (0..rows.len())
        .into_par_iter()
        .flat_map(|i| {
            let mut local = Vec::new();
            let Some(p_i) = rows[i].pic50 else {
                return local;
            };
            for j in (i + 1)..rows.len() {
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

    cliffs.sort_by(|a, b| b.sali.total_cmp(&a.sali));
    cliffs.truncate(100);
    cliffs
}

/// Clean-room single-cut Matched Molecular Pairs (`core` >= 50% heavy atoms + terminal `R-group`).
pub fn find_matched_molecular_pairs(rows: &[SarRow]) -> Vec<MatchedMolecularPair> {
    // Map canonical core SMILES -> list of (row_idx, r_group_smiles)
    let mut index: HashMap<String, Vec<(usize, String)>> = HashMap::new();

    for (r_idx, row) in rows.iter().enumerate() {
        let cuts = single_bond_cuts(&row.mol);
        for (core, rgrp) in cuts {
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
                if a_idx == b_idx || ra == rb {
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
    out.truncate(150);
    out
}

fn single_bond_cuts(mol: &MolGraph) -> Vec<(String, String)> {
    let total_heavy = mol.heavy_atom_count();
    if total_heavy < 6 {
        return Vec::new();
    }
    let mut res = Vec::new();
    for bid in 0..mol.bonds.len() {
        let b = &mol.bonds[bid];
        if b.in_ring {
            continue;
        }
        let comp_u = component_without_bond(mol, b.u, bid);
        let comp_v: BTreeSet<usize> = (0..mol.atoms.len())
            .filter(|id| !comp_u.contains(id))
            .collect();
        if comp_u.is_empty() || comp_v.is_empty() {
            continue;
        }
        // Core must be the larger scaffold, R-group at most 40% of molecule
        let (core_nodes, r_nodes) = if comp_u.len() >= comp_v.len() {
            (&comp_u, &comp_v)
        } else {
            (&comp_v, &comp_u)
        };
        if core_nodes.len() < 5 || r_nodes.len() > (total_heavy / 2).max(1) {
            continue;
        }
        if let (Some(c_smi), Some(r_smi)) = (
            induced_canonical(mol, core_nodes),
            induced_canonical(mol, r_nodes),
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

fn induced_canonical(mol: &MolGraph, nodes: &BTreeSet<usize>) -> Option<String> {
    let mut sub = MolGraph::new();
    let mut remap = vec![None; mol.atoms.len()];
    for &old in nodes {
        let mut a = mol.atoms[old].clone();
        a.coords_2d = None;
        let nid = sub.add_atom(a);
        remap[old] = Some(nid);
    }
    for b in &mol.bonds {
        if let (Some(nu), Some(nv)) = (remap[b.u], remap[b.v]) {
            let _ = sub.add_bond(nu, nv, b.order, b.stereo);
        }
    }
    sub.perceive().ok()?;
    Some(sub.canonical_smiles())
}
