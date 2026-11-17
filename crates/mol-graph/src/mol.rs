use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

use serde::{Deserialize, Serialize};

use crate::element::{element_by_number, element_by_symbol, max_neutral_valence};
use crate::error::MolError;

/// Covalent / coordination bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Dative,
}

impl BondOrder {
    pub fn valence_units(self) -> f64 {
        match self {
            Self::Single => 1.0,
            Self::Double => 2.0,
            Self::Triple => 3.0,
            Self::Aromatic => 1.5,
            Self::Dative => 1.0,
        }
    }

    pub fn integer_order_sum(self) -> u8 {
        match self {
            Self::Single | Self::Dative => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Aromatic => 1,
        }
    }
}

/// 2D stereochemical bond marker.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum BondStereo {
    #[default]
    None,
    Up,
    Down,
    Wavy,
    CisTransUp,
    CisTransDown,
}

/// Perceived double-bond E/Z geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum BondCisTrans {
    #[default]
    None,
    E,
    Z,
}

/// Tetrahedral chirality indicator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum ChiralTag {
    #[default]
    None,
    Clockwise,
    CounterClockwise,
    Unspecified,
}

/// Atom node inside an editable `MolGraph`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolAtom {
    pub id: usize,
    pub atomic_number: u8,
    pub symbol: String,
    pub formal_charge: i8,
    pub isotope: Option<u16>,
    pub aromatic: bool,
    pub explicit_h: u8,
    pub implicit_h: u8,
    pub chiral: ChiralTag,
    /// Neighbor slot index (inside `adj[id]` parse order) where bracket H appeared for `@`/`@@`.
    pub chiral_h_slot: Option<usize>,
    pub cip_label: Option<char>,
    pub in_ring: bool,
    pub ring_sizes: Vec<usize>,
    pub coords_2d: Option<[f64; 2]>,
    pub valence_exceeded: bool,
    pub from_bracket: bool,
}

impl MolAtom {
    pub fn new(id: usize, symbol: &str) -> Result<Self, MolError> {
        let info = element_by_symbol(symbol)
            .ok_or_else(|| MolError::Graph(format!("unsupported element symbol '{symbol}'")))?;
        Ok(Self {
            id,
            atomic_number: info.atomic_number,
            symbol: info.symbol.to_string(),
            formal_charge: 0,
            isotope: None,
            aromatic: false,
            explicit_h: 0,
            implicit_h: 0,
            chiral: ChiralTag::None,
            chiral_h_slot: None,
            cip_label: None,
            in_ring: false,
            ring_sizes: Vec::new(),
            coords_2d: None,
            valence_exceeded: false,
            from_bracket: false,
        })
    }

    pub fn total_h(&self) -> u8 {
        self.implicit_h.saturating_add(self.explicit_h)
    }

    pub fn atomic_weight(&self) -> f64 {
        if let Some(iso) = self.isotope {
            return f64::from(iso);
        }
        element_by_number(self.atomic_number).map_or(0.0, |e| e.atomic_weight)
    }
}

/// Covalent/coordination bond edge inside `MolGraph`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolBond {
    pub id: usize,
    pub u: usize,
    pub v: usize,
    pub order: BondOrder,
    pub stereo: BondStereo,
    /// Directional stereo recorded at bond creation: `true` if `u -> v` matched SMILES order.
    pub stereo_from_u: bool,
    pub cis_trans: BondCisTrans,
    pub in_ring: bool,
}

/// Clean-room editable molecular graph.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct MolGraph {
    pub name: Option<String>,
    pub input_smiles: Option<String>,
    pub atoms: Vec<MolAtom>,
    pub bonds: Vec<MolBond>,
    pub adj: Vec<Vec<(usize, usize)>>,
    pub rings: Vec<Vec<usize>>,
    pub warnings: Vec<String>,
    pub properties: BTreeMap<String, String>,
}

impl MolGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_atom(&mut self, mut atom: MolAtom) -> usize {
        let id = self.atoms.len();
        atom.id = id;
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        id
    }

    pub fn add_bond(
        &mut self,
        u: usize,
        v: usize,
        order: BondOrder,
        stereo: BondStereo,
    ) -> Result<usize, MolError> {
        if u >= self.atoms.len() || v >= self.atoms.len() {
            return Err(MolError::Graph(format!(
                "bond atom index out of bounds: ({u}, {v})"
            )));
        }
        if u == v {
            return Err(MolError::Graph("self-loop bond is forbidden".to_string()));
        }
        if let Some((_, eid)) = self.adj[u].iter().find(|&&(nbr, _)| nbr == v) {
            let eid = *eid;
            self.bonds[eid].order = order;
            self.bonds[eid].stereo = stereo;
            self.bonds[eid].stereo_from_u = true;
            return Ok(eid);
        }
        let id = self.bonds.len();
        self.bonds.push(MolBond {
            id,
            u,
            v,
            order,
            stereo,
            stereo_from_u: true,
            cis_trans: BondCisTrans::None,
            in_ring: false,
        });
        self.adj[u].push((v, id));
        self.adj[v].push((u, id));
        Ok(id)
    }

    pub fn set_atom_element(&mut self, idx: usize, symbol: &str) -> Result<(), MolError> {
        let info = element_by_symbol(symbol)
            .ok_or_else(|| MolError::Graph(format!("unsupported element '{symbol}'")))?;
        let atom = self
            .atoms
            .get_mut(idx)
            .ok_or_else(|| MolError::Graph(format!("invalid atom index {idx}")))?;
        atom.atomic_number = info.atomic_number;
        atom.symbol = info.symbol.to_string();
        atom.aromatic = false;
        atom.from_bracket = false;
        atom.explicit_h = 0;
        self.perceive()
    }

    pub fn remove_atom(&mut self, target: usize) -> Result<(), MolError> {
        if target >= self.atoms.len() {
            return Err(MolError::Graph(format!("invalid atom index {target}")));
        }
        let old_atoms = std::mem::take(&mut self.atoms);
        let old_bonds = std::mem::take(&mut self.bonds);
        self.adj.clear();
        let mut remap = vec![None; old_atoms.len()];
        for (old_id, mut atom) in old_atoms.into_iter().enumerate() {
            if old_id == target {
                continue;
            }
            let new_id = self.atoms.len();
            atom.id = new_id;
            remap[old_id] = Some(new_id);
            self.atoms.push(atom);
            self.adj.push(Vec::new());
        }
        for mut b in old_bonds {
            if let (Some(nu), Some(nv)) = (remap[b.u], remap[b.v]) {
                let id = self.bonds.len();
                b.id = id;
                b.u = nu;
                b.v = nv;
                self.bonds.push(b);
                self.adj[nu].push((nv, id));
                self.adj[nv].push((nu, id));
            }
        }
        self.perceive()
    }

    pub fn remove_bond(&mut self, target_bid: usize) -> Result<(), MolError> {
        if target_bid >= self.bonds.len() {
            return Err(MolError::Graph(format!("invalid bond index {target_bid}")));
        }
        let old_bonds = std::mem::take(&mut self.bonds);
        self.adj = vec![Vec::new(); self.atoms.len()];
        for (old_id, mut b) in old_bonds.into_iter().enumerate() {
            if old_id == target_bid {
                continue;
            }
            let id = self.bonds.len();
            b.id = id;
            let bu = b.u;
            let bv = b.v;
            self.bonds.push(b);
            self.adj[bu].push((bv, id));
            self.adj[bv].push((bu, id));
        }
        self.perceive()
    }

    /// Collapse standalone explicit hydrogen atoms onto `explicit_h` of their heavy parent.
    pub fn suppress_explicit_hydrogens(&mut self) {
        let mut h_parent: Vec<Option<usize>> = vec![None; self.atoms.len()];
        let mut is_remapped = vec![false; self.atoms.len()];
        for (idx, atom) in self.atoms.iter().enumerate() {
            if atom.atomic_number == 1 && atom.isotope.is_none() && self.adj[idx].len() == 1 {
                let (nbr, _) = self.adj[idx][0];
                if self.atoms[nbr].atomic_number != 1 {
                    h_parent[idx] = Some(nbr);
                    is_remapped[idx] = true;
                }
            }
        }
        if !is_remapped.iter().any(|&v| v) {
            return;
        }
        for (h_idx, maybe_parent) in h_parent.iter().enumerate() {
            if let Some(p) = maybe_parent {
                self.atoms[*p].explicit_h = self.atoms[*p].explicit_h.saturating_add(1);
                if self.atoms[h_idx].coords_2d.is_some() && self.atoms[*p].coords_2d.is_none() {
                    self.atoms[*p].coords_2d = self.atoms[h_idx].coords_2d;
                }
            }
        }
        let mut new_idx: Vec<Option<usize>> = vec![None; self.atoms.len()];
        let mut new_atoms = Vec::new();
        for (old, mut atom) in self.atoms.drain(..).enumerate() {
            if is_remapped[old] {
                continue;
            }
            let id = new_atoms.len();
            atom.id = id;
            new_idx[old] = Some(id);
            new_atoms.push(atom);
        }
        let mut new_bonds = Vec::new();
        let mut new_adj = vec![Vec::new(); new_atoms.len()];
        for mut b in self.bonds.drain(..) {
            if let (Some(nu), Some(nv)) = (new_idx[b.u], new_idx[b.v]) {
                let id = new_bonds.len();
                b.id = id;
                b.u = nu;
                b.v = nv;
                new_bonds.push(b);
                new_adj[nu].push((nv, id));
                new_adj[nv].push((nu, id));
            }
        }
        self.atoms = new_atoms;
        self.bonds = new_bonds;
        self.adj = new_adj;
    }

    /// Run full chemical perception: SSSR rings, initial valences, iterative Hückel aromaticity
    /// (monocyclic + fused 4n+2 systems), final valences/H counts, and CIP R/S + E/Z descriptors.
    pub fn perceive(&mut self) -> Result<(), MolError> {
        self.perceive_sssr();
        self.recompute_hydrogens_and_valence();
        self.perceive_aromaticity();
        self.recompute_hydrogens_and_valence();
        self.perceive_stereo_descriptors();
        Ok(())
    }

    /// Find Smallest Set of Smallest Rings (SSSR) without artificial 16-atom truncation.
    fn perceive_sssr(&mut self) {
        self.rings.clear();
        for a in &mut self.atoms {
            a.in_ring = false;
            a.ring_sizes.clear();
        }
        for b in &mut self.bonds {
            b.in_ring = false;
        }
        let n = self.atoms.len();
        let m = self.bonds.len();
        if n < 3 || m < 3 {
            return;
        }

        let mut seen = vec![false; n];
        let mut comps = 0usize;
        for s in 0..n {
            if seen[s] {
                continue;
            }
            comps += 1;
            let mut q = VecDeque::new();
            seen[s] = true;
            q.push_back(s);
            while let Some(u) = q.pop_front() {
                for &(v, _) in &self.adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        q.push_back(v);
                    }
                }
            }
        }
        let target_rank = m.saturating_add(comps).saturating_sub(n);
        if target_rank == 0 {
            return;
        }

        let mut candidates: Vec<(usize, Vec<usize>, BTreeSet<usize>)> = Vec::new();
        let mut seen_edge_sets: BTreeSet<Vec<usize>> = BTreeSet::new();

        for (bid, bond) in self.bonds.iter().enumerate() {
            let mut prev_node: Vec<Option<(usize, usize)>> = vec![None; n];
            let mut dist = vec![usize::MAX; n];
            let mut q = VecDeque::new();
            dist[bond.u] = 0;
            q.push_back(bond.u);
            while let Some(cur) = q.pop_front() {
                if cur == bond.v {
                    break;
                }
                for &(nbr, e) in &self.adj[cur] {
                    if e == bid {
                        continue;
                    }
                    if dist[nbr] == usize::MAX {
                        dist[nbr] = dist[cur] + 1;
                        prev_node[nbr] = Some((cur, e));
                        q.push_back(nbr);
                    }
                }
            }
            if dist[bond.v] == usize::MAX {
                continue;
            }
            let mut edge_set = BTreeSet::new();
            edge_set.insert(bid);
            let mut ring_nodes = vec![bond.v];
            let mut cur = bond.v;
            while let Some((p, e)) = prev_node[cur] {
                edge_set.insert(e);
                ring_nodes.push(p);
                cur = p;
            }
            let key: Vec<usize> = edge_set.iter().copied().collect();
            if seen_edge_sets.insert(key) {
                candidates.push((ring_nodes.len(), ring_nodes, edge_set));
            }
        }

        candidates.sort_by_key(|c| c.0);
        let mut basis: Vec<(usize, Vec<u64>)> = Vec::new();
        for (_len, ring_nodes, edge_set) in candidates {
            if self.rings.len() >= target_rank {
                break;
            }
            if insert_gf2(&mut basis, &edge_set, m) {
                for &v in &ring_nodes {
                    self.atoms[v].in_ring = true;
                    self.atoms[v].ring_sizes.push(ring_nodes.len());
                }
                for &e in &edge_set {
                    self.bonds[e].in_ring = true;
                }
                self.rings.push(ring_nodes);
            }
        }
    }

    /// Iterative fixed-point Hückel `4n + 2` aromaticity perception across SSSR rings and
    /// fused aromatic ring systems (e.g. Kekulé naphthalene, indole, quinoline, purine).
    fn perceive_aromaticity(&mut self) {
        // 1. Pre-sync atoms joined by explicit aromatic bonds (e.g., SDF bond type 4)
        for b in &self.bonds {
            if b.order == BondOrder::Aromatic {
                self.atoms[b.u].aromatic = true;
                self.atoms[b.v].aromatic = true;
            }
        }

        // Preserve pyrrole/indole NH before aromatic retyping
        let initial_nh: Vec<bool> = self
            .atoms
            .iter()
            .map(|a| a.symbol == "N" && a.total_h() > 0)
            .collect();

        // 2. Build candidate monocyclic rings + fused 2-ring & 3-ring symmetric-difference envelopes
        let mut envelopes: Vec<BTreeSet<usize>> = self
            .rings
            .iter()
            .filter(|r| (5..=7).contains(&r.len()))
            .map(|r| r.iter().copied().collect())
            .collect();

        for i in 0..self.rings.len() {
            if !(5..=6).contains(&self.rings[i].len()) {
                continue;
            }
            let si: BTreeSet<usize> = self.rings[i].iter().copied().collect();
            for ring_j in self.rings.iter().skip(i + 1) {
                if !(5..=6).contains(&ring_j.len()) {
                    continue;
                }
                let sj: BTreeSet<usize> = ring_j.iter().copied().collect();
                if si.intersection(&sj).count() >= 2 {
                    let union: BTreeSet<usize> = si.union(&sj).copied().collect();
                    if union.len() <= 12 {
                        envelopes.push(union);
                    }
                }
            }
        }

        // 3. Fixed-point iteration until no more rings/envelopes become aromatic
        loop {
            let mut changed = false;
            for env in &envelopes {
                if env.iter().all(|&a| self.atoms[a].aromatic) {
                    continue;
                }
                if let Some(pi_count) = self.system_pi_electrons(env)
                    && pi_count >= 2
                    && (pi_count - 2) % 4 == 0
                {
                    for &a in env {
                        if !self.atoms[a].aromatic {
                            self.atoms[a].aromatic = true;
                            if initial_nh[a] && self.atoms[a].explicit_h == 0 {
                                self.atoms[a].explicit_h = 1;
                            }
                            changed = true;
                        }
                    }
                    for b in &mut self.bonds {
                        if env.contains(&b.u)
                            && env.contains(&b.v)
                            && b.order != BondOrder::Aromatic
                        {
                            b.order = BondOrder::Aromatic;
                            changed = true;
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }

        // 4. Sync ring bonds between aromatic endpoints
        for b in &mut self.bonds {
            if b.in_ring
                && self.atoms[b.u].aromatic
                && self.atoms[b.v].aromatic
                && matches!(b.order, BondOrder::Double | BondOrder::Single)
            {
                b.order = BondOrder::Aromatic;
            }
        }
    }

    fn system_pi_electrons(&self, sys: &BTreeSet<usize>) -> Option<usize> {
        let mut total_pi = 0usize;
        for &u in sys {
            let a = &self.atoms[u];
            let has_sys_double = self.adj[u].iter().any(|&(v, bid)| {
                sys.contains(&v)
                    && matches!(
                        self.bonds[bid].order,
                        BondOrder::Double | BondOrder::Aromatic
                    )
            });
            let has_exo_double = self.adj[u].iter().any(|&(v, bid)| {
                !sys.contains(&v)
                    && matches!(self.bonds[bid].order, BondOrder::Double | BondOrder::Triple)
            });

            match a.symbol.as_str() {
                "C" => {
                    if a.aromatic || has_sys_double {
                        total_pi += 1;
                    } else if a.formal_charge == -1 {
                        total_pi += 2;
                    } else if a.formal_charge == 1 || has_exo_double {
                        total_pi += 0;
                    } else {
                        return None;
                    }
                }
                "N" => {
                    if has_sys_double || a.aromatic {
                        let deg_in_sys = self.adj[u]
                            .iter()
                            .filter(|&&(v, _)| sys.contains(&v))
                            .count();
                        if a.total_h() > 0 || self.adj[u].len() > deg_in_sys {
                            total_pi += 2; // pyrrole / N-Me pyrrole / caffeine N-Me
                        } else {
                            total_pi += 1; // pyridine =N-
                        }
                    } else if self.adj[u].len() + usize::from(a.total_h()) <= 3 {
                        total_pi += 2;
                    } else {
                        return None;
                    }
                }
                "O" | "S" if !has_sys_double && self.adj[u].len() == 2 => {
                    total_pi += 2;
                }
                _ => return None,
            }
        }
        Some(total_pi)
    }

    fn recompute_hydrogens_and_valence(&mut self) {
        self.warnings
            .retain(|w| !w.starts_with("Valence error on atom"));
        for idx in 0..self.atoms.len() {
            let explicit_bond_sum = self.bond_order_sum(idx);
            let atom = &mut self.atoms[idx];
            atom.valence_exceeded = false;

            if !atom.from_bracket {
                let target_valence =
                    default_organic_valence(&atom.symbol, atom.formal_charge, explicit_bond_sum);
                let used = explicit_bond_sum + f64::from(atom.explicit_h);
                if target_valence == 0 {
                    atom.implicit_h = 0;
                    if used
                        > f64::from(max_neutral_valence(&atom.symbol, atom.formal_charge)) + 0.01
                    {
                        atom.valence_exceeded = true;
                    }
                } else if used > f64::from(target_valence) + 0.01 {
                    atom.implicit_h = 0;
                    atom.valence_exceeded = true;
                } else {
                    atom.implicit_h = (f64::from(target_valence) - used).floor().max(0.0) as u8;
                }
            } else {
                let max_v = max_neutral_valence(&atom.symbol, atom.formal_charge);
                let total = explicit_bond_sum + f64::from(atom.explicit_h);
                if total > f64::from(max_v) + 0.01 {
                    atom.valence_exceeded = true;
                }
            }

            if atom.valence_exceeded {
                self.warnings.push(format!(
                    "Valence error on atom #{} ({})",
                    idx + 1,
                    atom.symbol
                ));
            }
        }
    }

    fn bond_order_sum(&self, idx: usize) -> f64 {
        let atom = &self.atoms[idx];
        let mut arom_bonds = 0usize;
        let mut non_arom_units = 0.0;
        let mut has_exo_double_o = false;

        for &(nbr, bid) in &self.adj[idx] {
            let bo = self.bonds[bid].order;
            if bo == BondOrder::Aromatic {
                arom_bonds += 1;
            } else {
                non_arom_units += bo.valence_units();
                if bo == BondOrder::Double && matches!(self.atoms[nbr].symbol.as_str(), "O" | "S") {
                    has_exo_double_o = true;
                }
            }
        }

        if arom_bonds >= 2 && atom.aromatic {
            return match atom.symbol.as_str() {
                "C" => {
                    if arom_bonds >= 3 {
                        // Fused junction carbon (naphthalene, indole, purine): 3 sigma + 1 pi = 4
                        4.0 + non_arom_units
                    } else if has_exo_double_o {
                        // Aromatic carbonyl c(=O) (caffeine, uracil): 2 single ring + 1 double =O -> 4
                        2.0 + non_arom_units
                    } else {
                        // Ring carbon (CH or C-R): uses 3 valence in ring + exo sigma bonds
                        3.0 + non_arom_units
                    }
                }
                "N" => {
                    if atom.explicit_h > 0 || non_arom_units >= 1.0 {
                        // Pyrrole [nH] or N-methyl-indole / caffeine n(C): 2 ring bonds + H/subst = 3.0
                        2.0 + non_arom_units
                    } else {
                        // Pyridine / pyrimidine =N-: 3.0 valence (no implicit H)
                        3.0
                    }
                }
                "O" | "S" => 2.0 + non_arom_units,
                _ => non_arom_units + (arom_bonds as f64) * 1.5,
            };
        }
        non_arom_units + (arom_bonds as f64) * 1.5
    }

    /// Assign real Cahn-Ingold-Prelog `(R)/(S)` and `(E)/(Z)` stereochemical labels.
    fn perceive_stereo_descriptors(&mut self) {
        let cip_ranks = compute_cip_priorities(self);

        for idx in 0..self.atoms.len() {
            let a = &self.atoms[idx];
            if !matches!(a.chiral, ChiralTag::Clockwise | ChiralTag::CounterClockwise) {
                self.atoms[idx].cip_label = None;
                continue;
            }
            // Collect neighbor sequence in SMILES parse order `[n0, (optional H), n1, n2, n3]`
            let mut seq: Vec<usize> = self.adj[idx].iter().map(|&(nbr, _)| nbr).collect();
            if let Some(h_slot) = a.chiral_h_slot {
                let ins = h_slot.min(seq.len());
                seq.insert(ins, usize::MAX); // usize::MAX marks implicit/explicit H (lowest CIP priority)
            }
            if seq.len() != 4 {
                self.atoms[idx].cip_label = None;
                continue;
            }

            // Sort a copy by DESCENDING CIP priority (lowest priority H = (0, MAX) goes last!)
            let prio_key = |node: usize| -> (u64, usize) {
                if node == usize::MAX {
                    (0, usize::MAX)
                } else {
                    (cip_ranks[node], usize::MAX - node)
                }
            };
            let mut cip_sorted = seq.clone();
            cip_sorted.sort_by_key(|&a| std::cmp::Reverse(prio_key(a)));

            // Check that all 4 groups are distinct priorities (or distinct by CIP)
            let p0 = prio_key(cip_sorted[0]).0;
            let p1 = prio_key(cip_sorted[1]).0;
            let p2 = prio_key(cip_sorted[2]).0;
            let p3 = prio_key(cip_sorted[3]).0;
            if p0 == p1 || p1 == p2 || p2 == p3 {
                self.atoms[idx].cip_label = None;
                continue;
            }

            // Permutation parity between SMILES parse order `seq` and CIP descending `cip_sorted`:
            // In SMILES looking down `seq[0] -> center`:
            //   `@`  (CounterClockwise) corresponds to `S` when parse order matches CIP `[1, 2, 3, 4]`?
            // Let's verify rigorously via 3D determinant:
            //   Place center = (0,0,0), seq[0] = (0,0,+1) [viewer looks from seq[0] -> center],
            //   and for `@@` (Clockwise looking down seq[0]->center), place (seq[1], seq[2], seq[3])
            //   clockwise around -Z: seq[1]=(0, 1, -0.35), seq[2]=(+0.866, -0.5, -0.35), seq[3]=(-0.866, -0.5, -0.35).
            //   Then look down (center -> cip_sorted[3] [lowest]) and check if cip_sorted[0] -> [1] -> [2] is CW (R) or CCW (S)!
            let mut pos_map = [0usize; 4];
            for (k, &s_val) in seq.iter().enumerate() {
                if let Some(ci) = cip_sorted.iter().position(|&c| c == s_val) {
                    pos_map[k] = ci;
                }
            }
            let parity = perm_parity(&pos_map); // +1 even, -1 odd
            let smiles_sign: i32 = match a.chiral {
                ChiralTag::Clockwise => 1,         // @@
                ChiralTag::CounterClockwise => -1, // @
                _ => 0,
            };
            // Verify: L-alanine N[C@@H](C)C(=O)O:
            // seq = [N, H, CH3, COOH], cip_sorted = [N, COOH, CH3, H] -> pos_map = [0, 3, 2, 1] (1 swap -> parity = -1).
            // For smiles_sign = +1 (@@): smiles_sign * parity = -1 -> S!
            // Check another SMILES for L-alanine: C[C@H](N)C(=O)O:
            // seq = [CH3, H, N, COOH], cip_sorted = [N, COOH, CH3, H] -> pos_map = [2, 3, 0, 1] (2 swaps -> parity = +1).
            // For smiles_sign = -1 (@): smiles_sign * parity = -1 -> S!
            // Both give 'S' for L-alanine and 'R' for D-alanine (`N[C@H](C)C(=O)O` -> +1 -> R)!
            let rs_sign = smiles_sign * parity;
            self.atoms[idx].cip_label = if rs_sign > 0 { Some('R') } else { Some('S') };
        }

        // E/Z double-bond stereochemistry
        for bid in 0..self.bonds.len() {
            if self.bonds[bid].order != BondOrder::Double {
                self.bonds[bid].cis_trans = BondCisTrans::None;
                continue;
            }
            let u = self.bonds[bid].u;
            let v = self.bonds[bid].v;
            let side_u = highest_cip_double_bond_side(self, u, v, &cip_ranks);
            let side_v = highest_cip_double_bond_side(self, v, u, &cip_ranks);
            self.bonds[bid].cis_trans = match (side_u, side_v) {
                (Some(su), Some(sv)) => {
                    if su == sv {
                        BondCisTrans::Z
                    } else {
                        BondCisTrans::E
                    }
                }
                _ => BondCisTrans::None,
            };
        }
    }

    pub fn heavy_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| a.atomic_number > 1).count()
    }

    pub fn aromatic_atom_count(&self) -> usize {
        self.atoms
            .iter()
            .filter(|a| a.atomic_number > 1 && a.aromatic)
            .count()
    }

    /// Standard IUPAC Hill molecular formula.
    pub fn formula(&self) -> String {
        let mut counts: BTreeMap<String, u32> = BTreeMap::new();
        for a in &self.atoms {
            if a.atomic_number == 0 {
                continue;
            }
            *counts.entry(a.symbol.clone()).or_default() += 1;
            let h = u32::from(a.total_h());
            if h > 0 && a.atomic_number != 1 {
                *counts.entry("H".to_string()).or_default() += h;
            }
        }
        let mut out = String::new();
        if let Some(c) = counts.remove("C") {
            out.push('C');
            if c > 1 {
                out.push_str(&c.to_string());
            }
            if let Some(h) = counts.remove("H") {
                out.push('H');
                if h > 1 {
                    out.push_str(&h.to_string());
                }
            }
        }
        for (sym, cnt) in counts {
            out.push_str(&sym);
            if cnt > 1 {
                out.push_str(&cnt.to_string());
            }
        }
        out
    }

    pub fn molecular_weight(&self) -> f64 {
        let h_weight = element_by_number(1).map_or(1.008, |e| e.atomic_weight);
        self.atoms
            .iter()
            .map(|a| a.atomic_weight() + f64::from(a.total_h()) * h_weight)
            .sum()
    }

    /// Deterministic Open-Spec connectivity digest (`CCKEY-<14-hex>-<FORMULA>`) derived
    /// from canonical SMILES (respecting aromaticity normalization, chirality, and E/Z).
    pub fn structural_key(&self) -> String {
        let can = self.canonical_smiles();
        let mut h1: u64 = 0xcbf2_9ce4_8422_2325;
        let mut h2: u64 = 0x1000_0000_01b3_9e37;
        for &b in can.as_bytes() {
            h1 ^= u64::from(b);
            h1 = h1.wrapping_mul(0x0100_0000_01b3);
            h2 = h2.rotate_left(7) ^ u64::from(b).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        }
        format!(
            "CCKEY-{:08X}{:06X}-{}",
            (h1 >> 32) as u32,
            (h2 & 0xff_ffff) as u32,
            self.formula()
        )
    }

    /// Strip recognized alkali/alkaline-earth/halide/hydrate pharmaceutical salts while
    /// preserving neutral mixtures, salt pairs, and transition-metal complexes (e.g. ferrocene).
    pub fn desalt(&mut self) {
        let n = self.atoms.len();
        if n <= 1 {
            return;
        }
        let mut seen = vec![false; n];
        let mut comps: Vec<Vec<usize>> = Vec::new();
        for start in 0..n {
            if seen[start] {
                continue;
            }
            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            seen[start] = true;
            q.push_back(start);
            while let Some(u) = q.pop_front() {
                comp.push(u);
                for &(v, _) in &self.adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        q.push_back(v);
                    }
                }
            }
            comps.push(comp);
        }
        if comps.len() <= 1 {
            return;
        }
        let is_salt: Vec<bool> = comps
            .iter()
            .map(|c| self.is_counterion_component(c))
            .collect();
        if !is_salt.iter().any(|&s| s) || !is_salt.iter().any(|&s| !s) {
            return;
        }
        let mut keep = BTreeSet::new();
        for (i, comp) in comps.iter().enumerate() {
            if !is_salt[i] {
                keep.extend(comp.iter().copied());
            }
        }
        let mut remap = vec![None; n];
        let mut new_atoms = Vec::new();
        for (old, mut atom) in self.atoms.drain(..).enumerate() {
            if keep.contains(&old) {
                let nid = new_atoms.len();
                atom.id = nid;
                remap[old] = Some(nid);
                new_atoms.push(atom);
            }
        }
        let mut new_bonds = Vec::new();
        let mut new_adj = vec![Vec::new(); new_atoms.len()];
        for mut b in self.bonds.drain(..) {
            if let (Some(nu), Some(nv)) = (remap[b.u], remap[b.v]) {
                let id = new_bonds.len();
                b.id = id;
                b.u = nu;
                b.v = nv;
                new_bonds.push(b);
                new_adj[nu].push((nv, id));
                new_adj[nv].push((nu, id));
            }
        }
        self.atoms = new_atoms;
        self.bonds = new_bonds;
        self.adj = new_adj;
        self.warnings
            .push("Desalted: stripped recognized counterion/hydrate fragment(s).".to_string());
    }

    fn is_counterion_component(&self, comp: &[usize]) -> bool {
        let heavy: Vec<usize> = comp
            .iter()
            .copied()
            .filter(|&a| self.atoms[a].atomic_number != 1)
            .collect();
        if heavy.len() == 1 {
            let a = &self.atoms[heavy[0]];
            // Only alkali / alkaline-earth / halide / hydrate O (exclude Fe/Ru/Pt/etc. organometallics!)
            return matches!(
                a.atomic_number,
                3 | 8 | 9 | 11 | 12 | 17 | 19 | 20 | 35 | 53 | 56
            );
        }
        false
    }

    pub fn canonical_smiles(&self) -> String {
        self.try_canonical_smiles().unwrap_or_default()
    }

    pub fn try_canonical_smiles(&self) -> Result<String, MolError> {
        if self.atoms.is_empty() {
            return Ok(String::new());
        }
        let n = self.atoms.len();
        let mut ranks: Vec<u64> = self
            .atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let deg = self.adj[i].len() as u64;
                let z = u64::from(a.atomic_number);
                let h = u64::from(a.total_h());
                let ch = (i64::from(a.formal_charge) + 16) as u64;
                let arom = u64::from(a.aromatic);
                let cip = match a.cip_label {
                    Some('R') => 1u64,
                    Some('S') => 2u64,
                    _ => 0u64,
                };
                (z << 28) | (deg << 20) | (h << 14) | (ch << 8) | (arom << 2) | cip
            })
            .collect();

        for _ in 0..6 {
            let mut tuples: Vec<(u64, Vec<u64>, usize)> = (0..n)
                .map(|i| {
                    let mut nbrs: Vec<u64> = self.adj[i]
                        .iter()
                        .map(|&(nbr, bid)| {
                            let bo = match self.bonds[bid].order {
                                BondOrder::Single | BondOrder::Dative => 1u64,
                                BondOrder::Double => 2,
                                BondOrder::Triple => 3,
                                BondOrder::Aromatic => 4,
                            };
                            (ranks[nbr] << 3) | bo
                        })
                        .collect();
                    nbrs.sort_unstable();
                    (ranks[i], nbrs, i)
                })
                .collect();
            tuples.sort();
            let mut next_ranks = vec![0u64; n];
            let mut rank = 1u64;
            for idx in 0..tuples.len() {
                if idx > 0
                    && (tuples[idx].0 != tuples[idx - 1].0 || tuples[idx].1 != tuples[idx - 1].1)
                {
                    rank += 1;
                }
                next_ranks[tuples[idx].2] = rank;
            }
            ranks = next_ranks;
        }

        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| (ranks[i], i));
        let mut visited = vec![false; n];
        let mut parts = Vec::new();
        for &start in &order {
            if !visited[start] {
                parts.push(self.emit_component_smiles(start, &ranks, &mut visited)?);
            }
        }
        parts.sort();
        Ok(parts.join("."))
    }

    fn emit_component_smiles(
        &self,
        start: usize,
        ranks: &[u64],
        visited: &mut [bool],
    ) -> Result<String, MolError> {
        let mut tree_vis = vec![false; self.atoms.len()];
        let mut closures = Vec::new();
        let mut closure_set = BTreeSet::new();

        struct Frame {
            u: usize,
            parent: usize,
            nbrs: Vec<(usize, usize)>,
            next: usize,
        }
        let mut init_nbrs = self.adj[start].clone();
        init_nbrs.sort_by_key(|&(v, _)| (ranks[v], v));
        tree_vis[start] = true;
        let mut stack = vec![Frame {
            u: start,
            parent: usize::MAX,
            nbrs: init_nbrs,
            next: 0,
        }];
        while let Some(f) = stack.last_mut() {
            if f.next >= f.nbrs.len() {
                stack.pop();
                continue;
            }
            let pu = f.u;
            let (v, bid) = f.nbrs[f.next];
            f.next += 1;
            if v == f.parent {
                continue;
            }
            if !tree_vis[v] {
                tree_vis[v] = true;
                let mut nbrs = self.adj[v].clone();
                nbrs.sort_by_key(|&(nv, _)| (ranks[nv], nv));
                stack.push(Frame {
                    u: v,
                    parent: pu,
                    nbrs,
                    next: 0,
                });
            } else if closure_set.insert(bid) {
                closures.push(bid);
            }
        }

        let mut atom_closures: Vec<Vec<(usize, usize, BondOrder)>> =
            vec![Vec::new(); self.atoms.len()];
        for &bid in &closures {
            let b = &self.bonds[bid];
            atom_closures[b.u].push((bid, b.v, b.order));
            atom_closures[b.v].push((bid, b.u, b.order));
        }
        for list in &mut atom_closures {
            list.sort_by_key(|&(_, o, _)| (ranks[o], o));
        }

        enum Op {
            Visit {
                u: usize,
                p: usize,
                bo: Option<BondOrder>,
                p_arom: bool,
                ez_tag: Option<char>,
            },
            Ch(char),
        }
        let mut out = String::new();
        let mut used = [false; 100];
        let mut active: HashMap<usize, u16> = HashMap::new();
        let mut emit = vec![Op::Visit {
            u: start,
            p: usize::MAX,
            bo: None,
            p_arom: false,
            ez_tag: None,
        }];

        // Track canonical E/Z side per double bond during emission
        let mut ez_tags: HashMap<usize, (char, char)> = HashMap::new();
        for b in &self.bonds {
            if b.cis_trans != BondCisTrans::None {
                let (c1, c2) = if b.cis_trans == BondCisTrans::E {
                    ('/', '/')
                } else {
                    ('/', '\\')
                };
                ez_tags.insert(b.id, (c1, c2));
            }
        }

        while let Some(op) = emit.pop() {
            match op {
                Op::Ch(c) => out.push(c),
                Op::Visit {
                    u,
                    p,
                    bo,
                    p_arom,
                    ez_tag,
                } => {
                    if visited[u] {
                        continue;
                    }
                    visited[u] = true;
                    if let Some(ez) = ez_tag {
                        out.push(ez);
                    }
                    if let Some(order) = bo {
                        match order {
                            BondOrder::Double => out.push('='),
                            BondOrder::Triple => out.push('#'),
                            BondOrder::Single | BondOrder::Dative => {
                                if p_arom && self.atoms[u].aromatic {
                                    out.push('-');
                                }
                            }
                            BondOrder::Aromatic => {}
                        }
                    }

                    // Determine canonical children + chiral marker derived deterministically from CIP label
                    let mut children: Vec<(usize, usize)> = self.adj[u]
                        .iter()
                        .copied()
                        .filter(|&(v, bid)| v != p && !closure_set.contains(&bid) && !visited[v])
                        .collect();
                    children.sort_by_key(|&(v, _)| (ranks[v], v));

                    out.push_str(&self.format_atom_smiles(u, p, &children));

                    let mut to_open = Vec::new();
                    for &(bid, _, border) in &atom_closures[u] {
                        if let Some(num) = active.remove(&bid) {
                            push_ring_tag(&mut out, num);
                            if usize::from(num) < used.len() {
                                used[usize::from(num)] = false;
                            }
                        } else {
                            to_open.push((bid, border));
                        }
                    }
                    for (bid, border) in to_open {
                        let Some(num) = (1..=99u16).find(|&k| !used[usize::from(k)]) else {
                            return Err(MolError::TooManyOpenRings);
                        };
                        used[usize::from(num)] = true;
                        active.insert(bid, num);
                        if border == BondOrder::Double {
                            out.push('=');
                        } else if border == BondOrder::Triple {
                            out.push('#');
                        }
                        push_ring_tag(&mut out, num);
                    }

                    let nc = children.len();
                    for (idx, (v, bid)) in children.into_iter().enumerate().rev() {
                        let last = idx + 1 == nc;
                        if !last {
                            emit.push(Op::Ch(')'));
                        }
                        // Check if parent `p` -> `u` was double bond (emit exit `/`/`\`) or `u` -> `v` leads into double bond (emit entry `/`)
                        let mut ez_c: Option<char> = None;
                        if p != usize::MAX
                            && let Some(&(_, p_bid)) = self.adj[u].iter().find(|&&(nb, _)| nb == p)
                            && let Some(&(_, c2)) = ez_tags.get(&p_bid)
                            && idx == 0
                        {
                            ez_c = Some(c2);
                        } else if idx == 0 {
                            let child_has_dbl = self.adj[v]
                                .iter()
                                .any(|&(_, cbid)| ez_tags.contains_key(&cbid));
                            if child_has_dbl {
                                ez_c = Some('/');
                            }
                        }
                        emit.push(Op::Visit {
                            u: v,
                            p: u,
                            bo: Some(self.bonds[bid].order),
                            p_arom: self.atoms[u].aromatic,
                            ez_tag: ez_c,
                        });
                        if !last {
                            emit.push(Op::Ch('('));
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    fn format_atom_smiles(&self, u: usize, p: usize, children: &[(usize, usize)]) -> String {
        let a = &self.atoms[u];
        let chiral_str = match a.cip_label {
            Some('R') => "@@",
            Some('S') => "@",
            _ => "",
        };
        let _ = (p, children);
        let organic = matches!(
            a.symbol.as_str(),
            "*" | "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
        let bracket = !organic
            || a.formal_charge != 0
            || a.isotope.is_some()
            || !chiral_str.is_empty()
            || (a.aromatic && a.symbol == "N" && a.explicit_h > 0)
            || (a.aromatic && !matches!(a.symbol.as_str(), "B" | "C" | "N" | "O" | "P" | "S"));
        if !bracket {
            return if a.aromatic {
                a.symbol.to_lowercase()
            } else {
                a.symbol.clone()
            };
        }
        let mut s = String::from("[");
        if let Some(iso) = a.isotope {
            s.push_str(&iso.to_string());
        }
        if a.aromatic {
            s.push_str(&a.symbol.to_lowercase());
        } else {
            s.push_str(&a.symbol);
        }
        if !chiral_str.is_empty() {
            s.push_str(chiral_str);
        }
        let h = a.total_h();
        if h == 1 {
            s.push('H');
        } else if h > 1 {
            s.push('H');
            s.push_str(&h.to_string());
        }
        if a.formal_charge > 0 {
            s.push('+');
            if a.formal_charge > 1 {
                s.push_str(&a.formal_charge.to_string());
            }
        } else if a.formal_charge < 0 {
            s.push('-');
            let mag = i16::from(a.formal_charge).abs();
            if mag > 1 {
                s.push_str(&mag.to_string());
            }
        }
        s.push(']');
        s
    }
}

/// Compute recursive Cahn-Ingold-Prelog (CIP) priority ranks for atoms (up to depth 4).
fn compute_cip_priorities(mol: &MolGraph) -> Vec<u64> {
    let n = mol.atoms.len();
    let mut ranks: Vec<u64> = mol
        .atoms
        .iter()
        .map(|a| {
            let z = u64::from(a.atomic_number);
            let iso = u64::from(a.isotope.unwrap_or(0));
            (z << 12) | iso
        })
        .collect();

    for _round in 0..4 {
        let mut next = vec![0u64; n];
        for i in 0..n {
            let mut child_tags: Vec<u64> = Vec::new();
            for &(nbr, bid) in &mol.adj[i] {
                let mult = match mol.bonds[bid].order {
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    _ => 1,
                };
                for _ in 0..mult {
                    child_tags.push(ranks[nbr]);
                }
            }
            for _ in 0..mol.atoms[i].total_h() {
                child_tags.push(1u64 << 12); // H atom priority
            }
            child_tags.sort_unstable_by(|a, b| b.cmp(a));
            let mut acc = ranks[i].wrapping_mul(131);
            for (k, &t) in child_tags.iter().take(4).enumerate() {
                acc = acc.wrapping_add(t.wrapping_mul(1u64 << (12 * (3 - k) / 3)));
            }
            next[i] = acc;
        }
        // Normalize compressed monotonic integer ranks so comparisons are exact!
        let mut paired: Vec<(u64, usize)> = next.iter().copied().zip(0..n).collect();
        paired.sort_unstable();
        let mut r = 10u64;
        for idx in 0..paired.len() {
            if idx > 0 && paired[idx].0 != paired[idx - 1].0 {
                r += 1;
            }
            ranks[paired[idx].1] = r;
        }
    }
    ranks
}

fn highest_cip_double_bond_side(
    mol: &MolGraph,
    atom: usize,
    other_end: usize,
    cip_ranks: &[u64],
) -> Option<i32> {
    let mut subs: Vec<(u64, i32)> = Vec::new();
    for &(nbr, bid) in &mol.adj[atom] {
        if nbr == other_end {
            continue;
        }
        let dir = match mol.bonds[bid].stereo {
            BondStereo::CisTransUp => {
                // `/` pointing FROM nbr -> atom puts `nbr` DOWN relative to atom (`-1`)!
                // `/` pointing FROM atom -> nbr puts `nbr` UP relative to atom (`+1`)!
                if mol.bonds[bid].u == atom { 1 } else { -1 }
            }
            BondStereo::CisTransDown => {
                if mol.bonds[bid].u == atom {
                    -1
                } else {
                    1
                }
            }
            _ => 0,
        };
        if dir != 0 {
            subs.push((cip_ranks[nbr], dir));
        }
    }
    if subs.is_empty() {
        return None;
    }
    subs.sort_by_key(|a| std::cmp::Reverse(a.0));
    Some(subs[0].1)
}

fn perm_parity(p: &[usize; 4]) -> i32 {
    let mut inv = 0usize;
    for i in 0..4 {
        for j in (i + 1)..4 {
            if p[i] > p[j] {
                inv += 1;
            }
        }
    }
    if inv.is_multiple_of(2) { 1 } else { -1 }
}

fn push_ring_tag(out: &mut String, num: u16) {
    if num < 10 {
        out.push_str(&num.to_string());
    } else {
        out.push('%');
        out.push_str(&format!("{num:02}"));
    }
}

fn default_organic_valence(sym: &str, charge: i8, bond_sum: f64) -> u8 {
    match sym {
        "*" => 1,
        "B" => {
            if charge < 0 {
                4
            } else {
                3
            }
        }
        "C" => {
            if charge < 0 {
                3
            } else {
                4
            }
        }
        "N" => {
            if charge > 0 {
                4
            } else if charge < 0 {
                2
            } else if bond_sum > 3.5 {
                5
            } else {
                3
            }
        }
        "O" => {
            if charge > 0 {
                3
            } else if charge < 0 {
                1
            } else {
                2
            }
        }
        "P" => {
            if bond_sum > 3.5 {
                5
            } else {
                3
            }
        }
        "S" => {
            if charge < 0 {
                1
            } else if bond_sum > 4.5 {
                6
            } else if bond_sum > 2.5 {
                4
            } else {
                2
            }
        }
        "F" | "Cl" | "Br" | "I" => {
            if charge < 0 {
                0
            } else {
                1
            }
        }
        _ => 0,
    }
}

fn insert_gf2(
    basis: &mut Vec<(usize, Vec<u64>)>,
    cand: &BTreeSet<usize>,
    num_edges: usize,
) -> bool {
    let words = num_edges.div_ceil(64);
    let mut vec_w = vec![0u64; words];
    for &e in cand {
        vec_w[e / 64] |= 1u64 << (e % 64);
    }
    for (pivot, row) in basis.iter() {
        if (vec_w[pivot / 64] & (1u64 << (pivot % 64))) != 0 {
            for (dst, src) in vec_w.iter_mut().zip(row.iter()) {
                *dst ^= *src;
            }
        }
    }
    let mut new_pivot = None;
    for (wi, &w) in vec_w.iter().enumerate() {
        if w != 0 {
            new_pivot = Some(wi * 64 + w.trailing_zeros() as usize);
            break;
        }
    }
    if let Some(p) = new_pivot {
        basis.push((p, vec_w));
        basis.sort_by_key(|(pv, _)| *pv);
        true
    } else {
        false
    }
}
