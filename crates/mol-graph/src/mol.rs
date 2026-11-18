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

/// Optional PDB / biomolecular residue metadata for peptide & complex export.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PdbAtomMeta {
    pub atom_name: String,
    pub res_name: String,
    pub res_seq: u16,
    pub chain: char,
    pub is_hetatm: bool,
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
    /// Exact OpenSMILES / 3D-derived 4-ligand neighbor sequence (`usize::MAX` marks bracket/suppressed H).
    #[serde(default)]
    pub chiral_seq: Vec<usize>,
    pub cip_label: Option<char>,
    pub in_ring: bool,
    pub ring_sizes: Vec<usize>,
    pub coords_2d: Option<[f64; 2]>,
    #[serde(default)]
    pub coords_3d: Option<[f64; 3]>,
    #[serde(default)]
    pub pdb_meta: Option<PdbAtomMeta>,
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
            chiral_seq: Vec::new(),
            cip_label: None,
            in_ring: false,
            ring_sizes: Vec::new(),
            coords_2d: None,
            coords_3d: None,
            pdb_meta: None,
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
                let pa = &mut self.atoms[*p];
                for (slot_i, nid) in pa.chiral_seq.iter_mut().enumerate() {
                    if *nid == h_idx {
                        *nid = usize::MAX;
                        if pa.chiral_h_slot.is_none() {
                            pa.chiral_h_slot = Some(slot_i);
                        }
                    }
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
        for atom in &mut new_atoms {
            atom.chiral_seq = atom
                .chiral_seq
                .iter()
                .filter_map(|&nid| {
                    if nid == usize::MAX {
                        Some(usize::MAX)
                    } else {
                        new_idx.get(nid).copied().flatten()
                    }
                })
                .collect();
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

        // 1b. Infer unambiguous aromatic NH protons (e.g. pyrrole, indole, uracil from type-4 SDFs)
        self.infer_unambiguous_aromatic_nh();

        // Preserve pyrrole/indole/uracil NH before aromatic retyping
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

        // Build fused 2-ring envelopes strictly via shared bond incidence (O(R) instead of O(R^2))
        let mut edge_to_rings: HashMap<usize, Vec<usize>> = HashMap::new();
        for (ri, r) in self.rings.iter().enumerate() {
            if !(5..=6).contains(&r.len()) {
                continue;
            }
            for w in 0..r.len() {
                let u = r[w];
                let v = r[(w + 1) % r.len()];
                if let Some(&(_, bid)) = self.adj[u].iter().find(|&&(nb, _)| nb == v) {
                    edge_to_rings.entry(bid).or_default().push(ri);
                }
            }
        }
        let mut seen_pairs: BTreeSet<(usize, usize)> = BTreeSet::new();
        for rlist in edge_to_rings.values() {
            for a_i in 0..rlist.len() {
                for b_i in (a_i + 1)..rlist.len() {
                    let i = rlist[a_i].min(rlist[b_i]);
                    let j = rlist[a_i].max(rlist[b_i]);
                    if seen_pairs.insert((i, j)) {
                        let mut union: BTreeSet<usize> = self.rings[i].iter().copied().collect();
                        union.extend(self.rings[j].iter().copied());
                        if union.len() <= 12 {
                            envelopes.push(union);
                        }
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

    /// Resolve unambiguous aromatic `[nH]` vs pyridine `=n-` hydrogens from Hückel `4n + 2`
    /// electron counts when reading pre-aromatized rings (such as bond-type-4 SDF files without explicit H).
    fn infer_unambiguous_aromatic_nh(&mut self) {
        let rings = self.rings.clone();
        for ring in rings {
            if !(5..=7).contains(&ring.len()) {
                continue;
            }
            let ring_set: BTreeSet<usize> = ring.iter().copied().collect();
            let all_arom = ring.iter().all(|&u| self.atoms[u].aromatic)
                && self.bonds.iter().all(|b| {
                    !(ring_set.contains(&b.u) && ring_set.contains(&b.v))
                        || b.order == BondOrder::Aromatic
                });
            if !all_arom {
                continue;
            }

            let mut fixed_pi = 0usize;
            let mut unknown_n: Vec<usize> = Vec::new();
            let mut valid = true;
            for &u in &ring {
                let a = &self.atoms[u];
                let has_exo_double = self.adj[u].iter().any(|&(v, bid)| {
                    !ring_set.contains(&v)
                        && matches!(self.bonds[bid].order, BondOrder::Double | BondOrder::Triple)
                });
                match a.symbol.as_str() {
                    "C" => {
                        if !has_exo_double {
                            fixed_pi += 1;
                        }
                    }
                    "O" | "S" => {
                        if !has_exo_double && self.adj[u].len() == 2 {
                            fixed_pi += 2;
                        } else {
                            valid = false;
                            break;
                        }
                    }
                    "N" => {
                        if a.formal_charge > 0 {
                            fixed_pi += 1;
                        } else if a.explicit_h > 0 || self.adj[u].len() > 2 || a.formal_charge < 0 {
                            fixed_pi += 2;
                        } else {
                            unknown_n.push(u);
                        }
                    }
                    _ => {
                        valid = false;
                        break;
                    }
                }
            }
            if !valid || unknown_n.is_empty() || unknown_n.len() > 4 {
                continue;
            }

            let k = unknown_n.len();
            let mut huckel_masks: Vec<usize> = Vec::new();
            for mask in 0..(1usize << k) {
                let mut total = fixed_pi;
                for i in 0..k {
                    total += if (mask & (1 << i)) != 0 { 2 } else { 1 };
                }
                if total >= 2 && (total - 2).is_multiple_of(4) {
                    huckel_masks.push(mask);
                }
            }
            if huckel_masks.len() == 1 {
                let mask = huckel_masks[0];
                for (i, &u) in unknown_n.iter().enumerate() {
                    if (mask & (1 << i)) != 0 {
                        self.atoms[u].explicit_h = 1;
                    }
                }
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
                    if has_exo_double {
                        total_pi += 0;
                    } else if a.aromatic || has_sys_double {
                        total_pi += 1;
                    } else if a.formal_charge == -1 {
                        total_pi += 2;
                    } else if a.formal_charge == 1 {
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

    /// Effective 4-neighbor tetrahedral sequence `[n0, n1, n2, n3]` matching `ChiralTag`.
    pub(crate) fn effective_chiral_seq(&self, idx: usize) -> Option<Vec<usize>> {
        let a = self.atoms.get(idx)?;
        if !a.chiral_seq.is_empty() && a.chiral_seq.len() == 4 {
            // Verify every neighbor in adj[idx] is present once
            let ok = self.adj[idx]
                .iter()
                .all(|&(nbr, _)| a.chiral_seq.contains(&nbr));
            if ok {
                return Some(a.chiral_seq.clone());
            }
        }
        let mut seq: Vec<usize> = self.adj[idx].iter().map(|&(nbr, _)| nbr).collect();
        if let Some(h_slot) = a.chiral_h_slot {
            let ins = h_slot.min(seq.len());
            seq.insert(ins, usize::MAX);
        } else if a.total_h() == 1 && seq.len() == 3 {
            seq.insert(1.min(seq.len()), usize::MAX);
        }
        (seq.len() == 4).then_some(seq)
    }

    /// Assign real multi-sphere Cahn-Ingold-Prelog `(R)/(S)` and `(E)/(Z)` stereochemical labels.
    fn perceive_stereo_descriptors(&mut self) {
        for idx in 0..self.atoms.len() {
            let chiral = self.atoms[idx].chiral;
            if !matches!(chiral, ChiralTag::Clockwise | ChiralTag::CounterClockwise) {
                self.atoms[idx].cip_label = None;
                continue;
            }
            let Some(seq) = self.effective_chiral_seq(idx) else {
                self.atoms[idx].cip_label = None;
                continue;
            };

            let mut cip_sorted = seq.clone();
            cip_sorted.sort_by(|&a, &b| cip_compare_branches(self, idx, b, a));

            let has_tie = cip_sorted
                .windows(2)
                .any(|w| cip_compare_branches(self, idx, w[0], w[1]) == std::cmp::Ordering::Equal);
            if has_tie {
                self.atoms[idx].cip_label = None;
                continue;
            }

            let mut pos_map = [0usize; 4];
            for (k, &s_val) in seq.iter().enumerate() {
                if let Some(ci) = cip_sorted.iter().position(|&c| c == s_val) {
                    pos_map[k] = ci;
                }
            }
            let parity = perm_parity(&pos_map);
            let smiles_sign: i32 = match chiral {
                ChiralTag::Clockwise => 1,         // @@
                ChiralTag::CounterClockwise => -1, // @
                _ => 0,
            };
            let rs_sign = smiles_sign * parity;
            self.atoms[idx].cip_label = if rs_sign > 0 { Some('R') } else { Some('S') };
        }

        // E/Z double-bond stereochemistry using full CIP ligand ranking (SMILES / or 2D CTAB geometry)
        for bid in 0..self.bonds.len() {
            if self.bonds[bid].order != BondOrder::Double {
                self.bonds[bid].cis_trans = BondCisTrans::None;
                continue;
            }
            let u = self.bonds[bid].u;
            let v = self.bonds[bid].v;
            let side_u = highest_cip_double_bond_nbr(self, u, v).map(|(_, s)| s);
            let side_v = highest_cip_double_bond_nbr(self, v, u).map(|(_, s)| s);
            self.bonds[bid].cis_trans = match (side_u, side_v) {
                (Some(su), Some(sv)) => {
                    if su == sv {
                        BondCisTrans::Z
                    } else {
                        BondCisTrans::E
                    }
                }
                _ if !self.bonds[bid].in_ring && self.bonds[bid].stereo != BondStereo::Wavy => {
                    geometric_ez_from_2d(self, u, v)
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
        for atom in &mut new_atoms {
            atom.chiral_seq = atom
                .chiral_seq
                .iter()
                .filter_map(|&nid| {
                    if nid == usize::MAX {
                        Some(usize::MAX)
                    } else {
                        remap.get(nid).copied().flatten()
                    }
                })
                .collect();
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
                let has_chiral = u64::from(a.cip_label.is_some());
                (z << 28) | (deg << 20) | (h << 14) | (ch << 8) | (arom << 2) | has_chiral
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
                            let ez = match self.bonds[bid].cis_trans {
                                BondCisTrans::E => 1u64,
                                BondCisTrans::Z => 2,
                                BondCisTrans::None => 0,
                            };
                            (ranks[nbr] << 5) | (bo << 2) | ez
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

        // Track entry substituent & incoming vertical side per E/Z double bond during emission:
        // map: bond_id -> (entry_substituent_atom_id, entry_side_relative_to_double_atom)
        let mut ez_entry: HashMap<usize, (usize, i32)> = HashMap::new();

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

                    // Determine canonical children and exact SMILES text neighbor order `emit_seq`
                    let mut children: Vec<(usize, usize)> = self.adj[u]
                        .iter()
                        .copied()
                        .filter(|&(v, bid)| v != p && !closure_set.contains(&bid) && !visited[v])
                        .collect();
                    children.sort_by_key(|&(v, _)| (ranks[v], v));

                    let mut emit_seq: Vec<usize> = Vec::with_capacity(4);
                    if p != usize::MAX {
                        emit_seq.push(p);
                    }
                    if self.atoms[u].total_h() == 1 {
                        emit_seq.push(usize::MAX);
                    }
                    for &(bid, other, _) in &atom_closures[u] {
                        if active.contains_key(&bid) {
                            emit_seq.push(other);
                        }
                    }
                    for &(bid, other, _) in &atom_closures[u] {
                        if !active.contains_key(&bid) {
                            emit_seq.push(other);
                        }
                    }
                    for &(v, _) in &children {
                        emit_seq.push(v);
                    }

                    out.push_str(&self.format_atom_smiles(u, &emit_seq));

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
                    let root_ez_pair = if p == usize::MAX {
                        let dbl = children.iter().position(|&(_, cbid)| {
                            matches!(
                                self.bonds[cbid].cis_trans,
                                BondCisTrans::E | BondCisTrans::Z
                            )
                        });
                        dbl.and_then(|di| {
                            children
                                .iter()
                                .enumerate()
                                .position(|(ci, _)| ci != di)
                                .map(|si| (si, children[di].1))
                        })
                    } else {
                        None
                    };
                    let p_ez_info = (p != usize::MAX)
                        .then(|| {
                            self.adj[u]
                                .iter()
                                .find(|&&(nb, _)| nb == p)
                                .map(|&(_, p_bid)| p_bid)
                        })
                        .flatten()
                        .filter(|&p_bid| {
                            matches!(
                                self.bonds[p_bid].cis_trans,
                                BondCisTrans::E | BondCisTrans::Z
                            )
                        })
                        .and_then(|p_bid| {
                            ez_entry
                                .get(&p_bid)
                                .copied()
                                .map(|(entry_sub, entry_side)| (p_bid, entry_sub, entry_side))
                        });

                    let next_dbls: Vec<Option<usize>> = children
                        .iter()
                        .map(|&(v, _)| {
                            self.adj[v]
                                .iter()
                                .find(|&&(nb, cbid)| {
                                    nb != u
                                        && matches!(
                                            self.bonds[cbid].cis_trans,
                                            BondCisTrans::E | BondCisTrans::Z
                                        )
                                        && !ez_entry.contains_key(&cbid)
                                })
                                .map(|&(_, cbid)| cbid)
                        })
                        .collect();

                    // Prefer exiting `p_bid` along a child that simultaneously enters the next
                    // conjugated E/Z bond so trisubstituted branched alkenes never drop slashes.
                    let preferred_exit_idx = if p_ez_info.is_some() {
                        next_dbls
                            .iter()
                            .position(Option::is_some)
                            .or((nc > 0).then_some(0))
                    } else {
                        None
                    };

                    for (idx, (v, bid)) in children.into_iter().enumerate().rev() {
                        let last = idx + 1 == nc;
                        if !last {
                            emit.push(Op::Ch(')'));
                        }
                        let mut ez_c: Option<char> = None;
                        let next_dbl = next_dbls[idx];
                        if let Some((p_bid, entry_sub, entry_side)) = p_ez_info
                            && (preferred_exit_idx == Some(idx) || next_dbl.is_some())
                        {
                            let want_z = self.bonds[p_bid].cis_trans == BondCisTrans::Z;
                            let hi_p = highest_cip_end_nbr(self, p, u);
                            let hi_u = highest_cip_end_nbr(self, u, p);
                            let entry_is_hi = hi_p == Some(entry_sub);
                            let exit_is_hi = hi_u == Some(v);
                            let hi_entry_side = if entry_is_hi { entry_side } else { -entry_side };
                            let hi_exit_side = if want_z {
                                hi_entry_side
                            } else {
                                -hi_entry_side
                            };
                            let exit_side = if exit_is_hi {
                                hi_exit_side
                            } else {
                                -hi_exit_side
                            };
                            ez_c = Some(if exit_side > 0 { '/' } else { '\\' });
                        } else if let Some((si, d_bid)) = root_ez_pair
                            && idx == si
                        {
                            ez_c = Some('/');
                            ez_entry.insert(d_bid, (v, 1));
                        }

                        if let Some(cbid) = next_dbl
                            && !ez_entry.contains_key(&cbid)
                        {
                            let ch = ez_c.unwrap_or('/');
                            ez_c = Some(ch);
                            // Bond `u -> v` labeled `/` puts `u` DOWN (-1) relative to `v`; `\` puts `u` UP (+1)
                            let entry_side_into_v = if ch == '/' { -1 } else { 1 };
                            ez_entry.insert(cbid, (u, entry_side_into_v));
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

    fn format_atom_smiles(&self, u: usize, emit_seq: &[usize]) -> String {
        let a = &self.atoms[u];
        let chiral_str = match (a.chiral, self.effective_chiral_seq(u)) {
            (ChiralTag::Clockwise | ChiralTag::CounterClockwise, Some(orig_seq))
                if emit_seq.len() == 4 =>
            {
                let orig_sign: i32 = if a.chiral == ChiralTag::Clockwise {
                    1
                } else {
                    -1
                };
                let mut perm = [0usize; 4];
                let mut valid = true;
                for (k, &ov) in orig_seq.iter().enumerate() {
                    if let Some(ei) = emit_seq.iter().position(|&ev| ev == ov) {
                        perm[k] = ei;
                    } else {
                        valid = false;
                        break;
                    }
                }
                if valid {
                    if orig_sign * perm_parity(&perm) > 0 {
                        "@@"
                    } else {
                        "@"
                    }
                } else {
                    ""
                }
            }
            _ => "",
        };
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

/// Compare two ligand branches rooted at `root_a` and `root_b` attached to `anchor` using
/// directed breadth-first IUPAC Cahn-Ingold-Prelog expansion (multiple bonds and cycle closures
/// spawn leaf phantom duplicate atoms, scaled in halves so aromatic bonds = 3 vs single = 2).
fn cip_compare_branches(
    mol: &MolGraph,
    anchor: usize,
    root_a: usize,
    root_b: usize,
) -> std::cmp::Ordering {
    let sig_a = cip_bfs_spheres(mol, anchor, root_a, 12);
    let sig_b = cip_bfs_spheres(mol, anchor, root_b, 12);
    for d in 0..sig_a.len().max(sig_b.len()) {
        let sa = sig_a.get(d).map(Vec::as_slice).unwrap_or(&[]);
        let sb = sig_b.get(d).map(Vec::as_slice).unwrap_or(&[]);
        if let ord @ (std::cmp::Ordering::Less | std::cmp::Ordering::Greater) = sa.cmp(sb) {
            return ord;
        }
    }
    std::cmp::Ordering::Equal
}

fn cip_atom_key(mol: &MolGraph, node: usize) -> (u16, u16) {
    if node == usize::MAX {
        (1, 0)
    } else {
        let a = &mol.atoms[node];
        (u16::from(a.atomic_number), a.isotope.unwrap_or(0))
    }
}

fn bond_weight_halves(order: BondOrder) -> usize {
    match order {
        BondOrder::Double => 4,
        BondOrder::Triple => 6,
        BondOrder::Aromatic => 3,
        _ => 2,
    }
}

fn cip_bfs_spheres(
    mol: &MolGraph,
    anchor: usize,
    root: usize,
    max_depth: usize,
) -> Vec<Vec<(u16, u16)>> {
    let mut spheres = Vec::new();
    let r_key = cip_atom_key(mol, root);
    let root_mult = if root == usize::MAX {
        2
    } else {
        mol.adj[anchor]
            .iter()
            .find(|&&(nb, _)| nb == root)
            .map_or(2, |&(_, bid)| bond_weight_halves(mol.bonds[bid].order))
    };
    let mut s0 = Vec::with_capacity(root_mult);
    for _ in 0..root_mult {
        s0.push(r_key);
    }
    s0.sort_unstable_by(|a, b| b.cmp(a));
    spheres.push(s0);

    if root == usize::MAX {
        return spheres;
    }

    // Linear-time visited-deduplicated breadth-first spheres (`O(N + M)`, immune to cyclic DoS)
    let n = mol.atoms.len();
    let mut visited = vec![false; n];
    visited[anchor] = true;
    visited[root] = true;
    let mut frontier: Vec<(usize, usize)> = vec![(root, anchor)];
    for _depth in 1..max_depth {
        let mut next_frontier: Vec<(usize, usize)> = Vec::new();
        let mut bag: Vec<(u16, u16)> = Vec::new();
        for (cur, parent) in frontier {
            for &(nbr, bid) in &mol.adj[cur] {
                if nbr == parent {
                    continue;
                }
                let w = bond_weight_halves(mol.bonds[bid].order);
                let key = cip_atom_key(mol, nbr);
                for _ in 0..w {
                    bag.push(key);
                }
                if !visited[nbr] {
                    visited[nbr] = true;
                    next_frontier.push((nbr, cur));
                }
            }
            for _ in 0..mol.atoms[cur].total_h() {
                bag.push((1, 0));
                bag.push((1, 0));
            }
        }
        if bag.is_empty() && next_frontier.is_empty() {
            break;
        }
        bag.sort_unstable_by(|a, b| b.cmp(a));
        spheres.push(bag);
        if next_frontier.is_empty() {
            break;
        }
        frontier = next_frontier;
    }
    spheres
}

pub fn highest_cip_end_nbr(mol: &MolGraph, atom: usize, other_end: usize) -> Option<usize> {
    let mut subs: Vec<usize> = mol.adj[atom]
        .iter()
        .map(|&(nb, _)| nb)
        .filter(|&nb| nb != other_end)
        .collect();
    if subs.is_empty() {
        return None;
    }
    if subs.len() == 1 {
        return Some(subs[0]);
    }
    if subs.len() == 2 {
        return match cip_compare_branches(mol, atom, subs[0], subs[1]) {
            std::cmp::Ordering::Greater => Some(subs[0]),
            std::cmp::Ordering::Less => Some(subs[1]),
            std::cmp::Ordering::Equal => None,
        };
    }
    subs.sort_by(|&a, &b| cip_compare_branches(mol, atom, b, a));
    (cip_compare_branches(mol, atom, subs[0], subs[1]) != std::cmp::Ordering::Equal)
        .then_some(subs[0])
}

fn geometric_ez_from_2d(mol: &MolGraph, u: usize, v: usize) -> BondCisTrans {
    let (Some(pu), Some(pv)) = (mol.atoms[u].coords_2d, mol.atoms[v].coords_2d) else {
        return BondCisTrans::None;
    };
    let dx = pv[0] - pu[0];
    let dy = pv[1] - pu[1];
    let len = dx.hypot(dy);
    if len < 1e-4 {
        return BondCisTrans::None;
    }
    let nx = -dy / len;
    let ny = dx / len;
    let Some(hi_u) = highest_cip_end_nbr(mol, u, v) else {
        return BondCisTrans::None;
    };
    let Some(hi_v) = highest_cip_end_nbr(mol, v, u) else {
        return BondCisTrans::None;
    };
    let (Some(phu), Some(phv)) = (mol.atoms[hi_u].coords_2d, mol.atoms[hi_v].coords_2d) else {
        return BondCisTrans::None;
    };
    let su = (phu[0] - pu[0]) * nx + (phu[1] - pu[1]) * ny;
    let sv = (phv[0] - pv[0]) * nx + (phv[1] - pv[1]) * ny;
    if su.abs() < 1e-3 || sv.abs() < 1e-3 {
        return BondCisTrans::None;
    }
    if (su > 0.0) == (sv > 0.0) {
        BondCisTrans::Z
    } else {
        BondCisTrans::E
    }
}

fn highest_cip_double_bond_nbr(
    mol: &MolGraph,
    atom: usize,
    other_end: usize,
) -> Option<(usize, i32)> {
    let mut subs: Vec<(usize, i32)> = Vec::new();
    for &(nbr, bid) in &mol.adj[atom] {
        if nbr == other_end {
            continue;
        }
        let dir = match mol.bonds[bid].stereo {
            BondStereo::CisTransUp => {
                if mol.bonds[bid].u == atom {
                    1
                } else {
                    -1
                }
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
        subs.push((nbr, dir));
    }
    if subs.is_empty() || subs.iter().all(|&(_, d)| d == 0) {
        return None;
    }
    if subs.len() == 2 {
        let (n0, d0) = subs[0];
        let (n1, d1) = subs[1];
        let eff0 = if d0 != 0 { d0 } else { -d1 };
        let eff1 = if d1 != 0 { d1 } else { -d0 };
        return match cip_compare_branches(mol, atom, n0, n1) {
            std::cmp::Ordering::Greater => Some((n0, eff0)),
            std::cmp::Ordering::Less => Some((n1, eff1)),
            std::cmp::Ordering::Equal => None,
        };
    }
    let (n0, d0) = subs[0];
    (d0 != 0).then_some((n0, d0))
}

pub(crate) fn perm_parity(p: &[usize; 4]) -> i32 {
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
