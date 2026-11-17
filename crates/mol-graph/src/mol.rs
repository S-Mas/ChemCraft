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
    /// Numeric bond-order contribution for valence accounting.
    pub fn valence_units(self) -> f64 {
        match self {
            Self::Single => 1.0,
            Self::Double => 2.0,
            Self::Triple => 3.0,
            Self::Aromatic => 1.5,
            Self::Dative => 1.0,
        }
    }

    /// Integer valence estimate used for implicit-H and hypervalency checks.
    pub fn integer_order_sum(self) -> u8 {
        match self {
            Self::Single | Self::Dative => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Aromatic => 1, // resolved via aromatic pi-credit
        }
    }
}

/// 2D stereochemical bond drawing marker (ChemDraw/VectorCraft wedge/hash/wavy).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum BondStereo {
    #[default]
    None,
    /// Solid wedged bond pointing toward viewer from `u` to `v`.
    Up,
    /// Hashed wedged bond pointing away from viewer from `u` to `v`.
    Down,
    /// Wavy / unspecified stereocenter bond.
    Wavy,
    /// SMILES `/` directional hint for E/Z double-bond geometry.
    CisTransUp,
    /// SMILES `\` directional hint for E/Z double-bond geometry.
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
    Clockwise,        // @@ / R candidate
    CounterClockwise, // @ / S candidate
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
    pub cip_label: Option<char>, // 'R' or 'S'
    pub in_ring: bool,
    pub ring_sizes: Vec<usize>,
    pub coords_2d: Option<[f64; 2]>,
    pub valence_exceeded: bool,
    pub from_bracket: bool,
}

impl MolAtom {
    /// Construct an atom with standard element metadata.
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
            cip_label: None,
            in_ring: false,
            ring_sizes: Vec::new(),
            coords_2d: None,
            valence_exceeded: false,
            from_bracket: false,
        })
    }

    /// Total hydrogens attached (`implicit_h + explicit_h`).
    pub fn total_h(&self) -> u8 {
        self.implicit_h.saturating_add(self.explicit_h)
    }

    /// Atomic weight (IUPAC standard unless an isotope is set).
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

    /// Add an atom node and return its assigned zero-based index.
    pub fn add_atom(&mut self, mut atom: MolAtom) -> usize {
        let id = self.atoms.len();
        atom.id = id;
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        id
    }

    /// Add a bond between atoms `u` and `v`. Updates existing bond if present.
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
            return Ok(eid);
        }
        let id = self.bonds.len();
        self.bonds.push(MolBond {
            id,
            u,
            v,
            order,
            stereo,
            cis_trans: BondCisTrans::None,
            in_ring: false,
        });
        self.adj[u].push((v, id));
        self.adj[v].push((u, id));
        Ok(id)
    }

    /// Interactive hotkey operation: change element at atom `idx` and re-perceive valences.
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

    /// Remove an atom and its incident bonds cleanly, reindexing remaining atoms and bonds.
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

    /// Remove a bond by index cleanly and re-perceive molecular state.
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

    /// Run full chemical perception: SSSR rings, Hückel aromaticity, implicit H,
    /// valence error diagnostics, and E/Z + R/S stereochemical hints.
    pub fn perceive(&mut self) -> Result<(), MolError> {
        self.perceive_sssr();
        self.perceive_aromaticity();
        self.recompute_hydrogens_and_valence();
        self.perceive_stereo_descriptors();
        Ok(())
    }

    /// Find Smallest Set of Smallest Rings (SSSR) using a clean-room Horton-style
    /// minimum-cycle basis over GF(2) edge space.
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

        // Connected components count -> expected cycle rank nu = m - n + c
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

        // Gather candidate shortest cycles through each edge and vertex (Horton basis).
        let mut candidates: Vec<(usize, Vec<usize>, BTreeSet<usize>)> = Vec::new();
        let mut seen_edge_sets: BTreeSet<Vec<usize>> = BTreeSet::new();

        for (bid, bond) in self.bonds.iter().enumerate() {
            // Shortest path between bond.u and bond.v ignoring edge `bid`
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
            if dist[bond.v] == usize::MAX || dist[bond.v] + 1 > 16 {
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

    /// Clean-room Hückel `4n + 2` aromaticity perception across SSSR rings + fused 6/5 systems.
    fn perceive_aromaticity(&mut self) {
        for b in &self.bonds {
            if b.order == BondOrder::Aromatic {
                self.atoms[b.u].aromatic = true;
                self.atoms[b.v].aromatic = true;
            }
        }
        for ring in &self.rings {
            if ring.len() < 5 || ring.len() > 7 {
                continue;
            }
            let already_arom = ring.iter().all(|&a| self.atoms[a].aromatic);
            if already_arom {
                continue;
            }
            if let Some(pi_count) = self.ring_pi_electrons(ring)
                && pi_count >= 2
                && (pi_count - 2) % 4 == 0
            {
                let ring_set: BTreeSet<usize> = ring.iter().copied().collect();
                for &a in ring {
                    self.atoms[a].aromatic = true;
                }
                for b in &mut self.bonds {
                    if ring_set.contains(&b.u) && ring_set.contains(&b.v) {
                        b.order = BondOrder::Aromatic;
                    }
                }
            }
        }

        // Sync aromatic bonds if endpoints are aromatic in the same ring
        let n_bonds = self.bonds.len();
        for i in 0..n_bonds {
            let b = &self.bonds[i];
            if b.in_ring
                && self.atoms[b.u].aromatic
                && self.atoms[b.v].aromatic
                && b.order == BondOrder::Double
            {
                self.bonds[i].order = BondOrder::Aromatic;
            }
        }
    }

    fn ring_pi_electrons(&self, ring: &[usize]) -> Option<usize> {
        let mut total_pi = 0usize;
        let ring_set: BTreeSet<usize> = ring.iter().copied().collect();
        for &u in ring {
            let a = &self.atoms[u];
            if a.aromatic {
                // Already marked aromatic from SMILES lowercase
                let has_double_in_ring = self.adj[u].iter().any(|&(v, bid)| {
                    ring_set.contains(&v) && self.bonds[bid].order == BondOrder::Double
                });
                if has_double_in_ring {
                    total_pi += 1;
                    continue;
                }
            }
            let has_ring_double = self.adj[u].iter().any(|&(v, bid)| {
                ring_set.contains(&v) && self.bonds[bid].order == BondOrder::Double
            });
            let has_exo_double = self.adj[u].iter().any(|&(v, bid)| {
                !ring_set.contains(&v)
                    && matches!(self.bonds[bid].order, BondOrder::Double | BondOrder::Triple)
            });

            match a.symbol.as_str() {
                "C" => {
                    if has_ring_double || (has_exo_double && a.formal_charge == 0) {
                        // Exocyclic C=O (e.g. uracil/caffeine) donates 0 pi if carbonyl oxygen withdraws
                        if has_ring_double {
                            total_pi += 1;
                        } else {
                            return None;
                        }
                    } else if a.formal_charge == 1 {
                        total_pi += 0;
                    } else if a.formal_charge == -1 {
                        total_pi += 2;
                    } else {
                        return None;
                    }
                }
                "N" => {
                    if has_ring_double {
                        total_pi += 1;
                    } else if self.adj[u].len() + usize::from(a.total_h()) <= 3 {
                        total_pi += 2; // pyrrole / indole NH
                    } else {
                        return None;
                    }
                }
                "O" | "S" if !has_ring_double && self.adj[u].len() == 2 => {
                    total_pi += 2; // furan / thiophene
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
                if used > f64::from(target_valence) + 0.01 {
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
        let mut raw = 0.0;
        let mut arom_count = 0usize;
        for &(_, bid) in &self.adj[idx] {
            let bo = self.bonds[bid].order;
            if bo == BondOrder::Aromatic {
                arom_count += 1;
            }
            raw += bo.valence_units();
        }
        if arom_count >= 2 && atom.aromatic {
            // Two aromatic bonds represent ~3.0 valence units cleanly unless pyrrole NH / furan O
            let has_exocyclic = self.adj[idx]
                .iter()
                .any(|&(_, bid)| self.bonds[bid].order != BondOrder::Aromatic);
            if matches!(atom.symbol.as_str(), "N") && (has_exocyclic || atom.explicit_h > 0) {
                raw = 2.0 + if has_exocyclic { 1.0 } else { 0.0 };
            } else if matches!(atom.symbol.as_str(), "O" | "S") && !has_exocyclic {
                raw = 2.0;
            } else {
                raw = raw.round();
                if atom.symbol == "C" && arom_count == 2 && !has_exocyclic {
                    raw = 3.0;
                } else if atom.symbol == "C" && arom_count >= 2 && has_exocyclic {
                    raw = 4.0;
                }
            }
        }
        raw
    }

    fn perceive_stereo_descriptors(&mut self) {
        // Assign lightweight E/Z and R/S labels when explicit stereo markers exist
        for a in &mut self.atoms {
            a.cip_label = match a.chiral {
                ChiralTag::Clockwise => Some('R'),
                ChiralTag::CounterClockwise => Some('S'),
                ChiralTag::None | ChiralTag::Unspecified => None,
            };
        }
        for i in 0..self.bonds.len() {
            if self.bonds[i].order != BondOrder::Double {
                self.bonds[i].cis_trans = BondCisTrans::None;
                continue;
            }
            let u = self.bonds[i].u;
            let v = self.bonds[i].v;
            let dir_u = self.adj[u]
                .iter()
                .find_map(|&(_, bid)| match self.bonds[bid].stereo {
                    BondStereo::CisTransUp => Some(1),
                    BondStereo::CisTransDown => Some(-1),
                    _ => None,
                });
            let dir_v = self.adj[v]
                .iter()
                .find_map(|&(_, bid)| match self.bonds[bid].stereo {
                    BondStereo::CisTransUp => Some(1),
                    BondStereo::CisTransDown => Some(-1),
                    _ => None,
                });
            self.bonds[i].cis_trans = match (dir_u, dir_v) {
                (Some(a), Some(b)) => {
                    if a == b {
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
        self.atoms.iter().filter(|a| a.atomic_number != 1).count()
    }

    pub fn aromatic_atom_count(&self) -> usize {
        self.atoms
            .iter()
            .filter(|a| a.atomic_number != 1 && a.aromatic)
            .count()
    }

    /// Standard IUPAC Hill molecular formula.
    pub fn formula(&self) -> String {
        let mut counts: BTreeMap<String, u32> = BTreeMap::new();
        for a in &self.atoms {
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

    /// IUPAC average molecular weight (g/mol).
    pub fn molecular_weight(&self) -> f64 {
        let h_weight = element_by_number(1).map_or(1.008, |e| e.atomic_weight);
        self.atoms
            .iter()
            .map(|a| a.atomic_weight() + f64::from(a.total_h()) * h_weight)
            .sum()
    }

    /// Deterministic Open-Spec connectivity digest (`CCKEY-<14-hex>-<FORMULA>`)
    /// derived strictly from canonical graph topology (no proprietary InChI binary).
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

    /// Strip recognized counterions/hydrates while preserving neutral mixtures and salt pairs.
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
            return matches!(
                a.atomic_number,
                3 | 8 | 9 | 11 | 12 | 13 | 17 | 19 | 20 | 25 | 26 | 29 | 30 | 35 | 53 | 56
            );
        }
        false
    }

    /// Deterministic canonical SMILES (infallible convenience wrapper).
    pub fn canonical_smiles(&self) -> String {
        self.try_canonical_smiles().unwrap_or_default()
    }

    /// Deterministic canonical SMILES, failing cleanly if >99 ring closures are open concurrently.
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
                (z << 24) | (deg << 16) | (h << 10) | (ch << 4) | arom
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
        }];

        while let Some(op) = emit.pop() {
            match op {
                Op::Ch(c) => out.push(c),
                Op::Visit { u, p, bo, p_arom } => {
                    if visited[u] {
                        continue;
                    }
                    visited[u] = true;
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
                    out.push_str(&self.format_atom_smiles(u));

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

                    let mut children: Vec<(usize, usize)> = self.adj[u]
                        .iter()
                        .copied()
                        .filter(|&(v, bid)| v != p && !closure_set.contains(&bid) && !visited[v])
                        .collect();
                    children.sort_by_key(|&(v, _)| (ranks[v], v));
                    let nc = children.len();
                    for (idx, (v, bid)) in children.into_iter().enumerate().rev() {
                        let last = idx + 1 == nc;
                        if !last {
                            emit.push(Op::Ch(')'));
                        }
                        emit.push(Op::Visit {
                            u: v,
                            p: u,
                            bo: Some(self.bonds[bid].order),
                            p_arom: self.atoms[u].aromatic,
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

    fn format_atom_smiles(&self, u: usize) -> String {
        let a = &self.atoms[u];
        let organic = matches!(
            a.symbol.as_str(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
        let bracket = !organic
            || a.formal_charge != 0
            || a.isotope.is_some()
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
            } else if bond_sum > 3.5 {
                5
            } else {
                3
            }
        }
        "O" => {
            if charge > 0 {
                3
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
            if bond_sum > 4.5 {
                6
            } else if bond_sum > 2.5 {
                4
            } else {
                2
            }
        }
        "F" | "Cl" | "Br" | "I" => 1,
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
