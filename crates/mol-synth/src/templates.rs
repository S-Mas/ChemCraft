#![forbid(unsafe_code)]

use std::collections::{HashSet, VecDeque};

use mol_graph::{BondOrder, BondStereo, MolAtom, MolGraph, parse_smiles};
use serde::{Deserialize, Serialize};

/// A single chemically validated retrosynthetic disconnection (Target <= Precursors).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Disconnection {
    pub template_id: String,
    pub reaction_name: String,
    pub reaction_class: String,
    /// Precursor canonical SMILES (1 for intramolecular/FGI, 2 for convergent coupling).
    pub precursor_smiles: Vec<String>,
    /// Precursor parsed graphs (not serialized).
    #[serde(skip)]
    pub precursors: Vec<MolGraph>,
    pub reagents: String,
    pub solvent: String,
    pub temperature: String,
    /// Typical medicinal-chemistry isolated yield fraction in `[0.30, 0.96]`.
    pub est_yield: f64,
    /// Step difficulty / preference score (lower is preferred in Retro* search).
    pub step_cost: f64,
    pub workup_notes: String,
    /// Chemoselectivity and protecting-group alerts for this specific substrate pair.
    pub warnings: Vec<String>,
}

/// Functional group cap attached to a cleaved radical atom during retrosynthetic disconnection.
#[derive(Debug, Clone, Copy)]
enum CapGroup {
    NoneH,
    AromNH,
    Oh,
    Br,
    Cl,
    I,
    BoronicAcid,
    CarbonylO,
    Cyano,
    Azide,
    Nitro,
}

fn attach_cap(mol: &mut MolGraph, idx: usize, cap: CapGroup) -> bool {
    match cap {
        CapGroup::NoneH => {
            if let Some(a) = mol.atoms.get_mut(idx) {
                if a.atomic_number == 7 && a.aromatic && a.ring_sizes.contains(&5) {
                    a.from_bracket = true;
                    a.explicit_h = 1;
                    a.implicit_h = 0;
                } else {
                    a.from_bracket = false;
                    a.explicit_h = 0;
                }
            }
            true
        }
        CapGroup::AromNH => {
            if let Some(a) = mol.atoms.get_mut(idx) {
                a.from_bracket = true;
                a.explicit_h = 1;
                a.implicit_h = 0;
            }
            true
        }
        CapGroup::Oh => add_single_atom(mol, idx, "O", BondOrder::Single),
        CapGroup::Br => add_single_atom(mol, idx, "Br", BondOrder::Single),
        CapGroup::Cl => add_single_atom(mol, idx, "Cl", BondOrder::Single),
        CapGroup::I => add_single_atom(mol, idx, "I", BondOrder::Single),
        CapGroup::CarbonylO => add_single_atom(mol, idx, "O", BondOrder::Double),
        CapGroup::BoronicAcid => {
            let Ok(b_atom) = MolAtom::new(mol.atoms.len(), "B") else {
                return false;
            };
            let b_id = mol.add_atom(b_atom);
            if mol
                .add_bond(idx, b_id, BondOrder::Single, BondStereo::None)
                .is_err()
            {
                return false;
            }
            add_single_atom(mol, b_id, "O", BondOrder::Single)
                && add_single_atom(mol, b_id, "O", BondOrder::Single)
        }
        CapGroup::Cyano => {
            let Ok(c_atom) = MolAtom::new(mol.atoms.len(), "C") else {
                return false;
            };
            let c_id = mol.add_atom(c_atom);
            if mol
                .add_bond(idx, c_id, BondOrder::Single, BondStereo::None)
                .is_err()
            {
                return false;
            }
            add_single_atom(mol, c_id, "N", BondOrder::Triple)
        }
        CapGroup::Azide => {
            let Ok(mut n1) = MolAtom::new(mol.atoms.len(), "N") else {
                return false;
            };
            n1.from_bracket = false;
            let n1_id = mol.add_atom(n1);
            let Ok(mut n2) = MolAtom::new(mol.atoms.len(), "N") else {
                return false;
            };
            n2.formal_charge = 1;
            n2.from_bracket = true;
            let n2_id = mol.add_atom(n2);
            let Ok(mut n3) = MolAtom::new(mol.atoms.len(), "N") else {
                return false;
            };
            n3.formal_charge = -1;
            n3.from_bracket = true;
            let n3_id = mol.add_atom(n3);
            mol.add_bond(idx, n1_id, BondOrder::Single, BondStereo::None)
                .is_ok()
                && mol
                    .add_bond(n1_id, n2_id, BondOrder::Double, BondStereo::None)
                    .is_ok()
                && mol
                    .add_bond(n2_id, n3_id, BondOrder::Double, BondStereo::None)
                    .is_ok()
        }
        CapGroup::Nitro => {
            let Ok(mut n) = MolAtom::new(mol.atoms.len(), "N") else {
                return false;
            };
            n.formal_charge = 1;
            n.from_bracket = true;
            let n_id = mol.add_atom(n);
            let Ok(o1) = MolAtom::new(mol.atoms.len(), "O") else {
                return false;
            };
            let o1_id = mol.add_atom(o1);
            let Ok(mut o2) = MolAtom::new(mol.atoms.len(), "O") else {
                return false;
            };
            o2.formal_charge = -1;
            o2.from_bracket = true;
            let o2_id = mol.add_atom(o2);
            mol.add_bond(idx, n_id, BondOrder::Single, BondStereo::None)
                .is_ok()
                && mol
                    .add_bond(n_id, o1_id, BondOrder::Double, BondStereo::None)
                    .is_ok()
                && mol
                    .add_bond(n_id, o2_id, BondOrder::Single, BondStereo::None)
                    .is_ok()
        }
    }
}

fn add_single_atom(mol: &mut MolGraph, idx: usize, sym: &str, order: BondOrder) -> bool {
    let Ok(atom) = MolAtom::new(mol.atoms.len(), sym) else {
        return false;
    };
    let nid = mol.add_atom(atom);
    mol.add_bond(idx, nid, order, BondStereo::None).is_ok()
}

/// Validates and canonicalizes a generated precursor graph by full perception + SMILES round-trip.
fn finalize_precursor(mut g: MolGraph) -> Option<(MolGraph, String)> {
    if g.perceive().is_err() {
        return None;
    }
    if g.atoms.is_empty() || g.atoms.iter().any(|a| a.valence_exceeded) {
        return None;
    }
    let smi = g.try_canonical_smiles().ok()?;
    if smi.is_empty() || smi.contains('.') {
        return None;
    }
    // Round-trip re-parse guarantees clean perception, stereochemistry & zero broken aromaticity
    let reparsed = parse_smiles(&smi).ok()?;
    if reparsed.atoms.iter().any(|a| a.valence_exceeded) {
        return None;
    }
    let canon = reparsed.canonical_smiles();
    Some((reparsed, canon))
}

/// Extracts a clean re-indexed subgraph for the atom set `keep`, returning `(subgraph, old_to_new)`.
fn extract_subgraph(
    mol: &MolGraph,
    keep: &HashSet<usize>,
    skip_bid: Option<usize>,
) -> (MolGraph, Vec<Option<usize>>) {
    let mut sub = MolGraph::new();
    let mut map = vec![None; mol.atoms.len()];
    for (old_i, a) in mol.atoms.iter().enumerate() {
        if keep.contains(&old_i) {
            let mut copy = a.clone();
            copy.coords_2d = None;
            copy.coords_3d = None;
            copy.chiral_seq.clear();
            copy.chiral_h_slot = None;
            // Clear bracket unless explicit charge/isotope or aromatic pyrrole NH
            if copy.formal_charge == 0
                && copy.isotope.is_none()
                && !(copy.aromatic && copy.atomic_number == 7 && copy.explicit_h > 0)
            {
                copy.from_bracket = false;
                copy.explicit_h = 0;
            }
            let new_i = sub.add_atom(copy);
            map[old_i] = Some(new_i);
        }
    }
    for (bid, b) in mol.bonds.iter().enumerate() {
        if Some(bid) == skip_bid {
            continue;
        }
        if let (Some(nu), Some(nv)) = (map[b.u], map[b.v]) {
            let _ = sub.add_bond(nu, nv, b.order, BondStereo::None);
        }
    }
    (sub, map)
}

/// Result of cleaving a single bond `(u, v)`.
enum BondCleavage {
    /// Intermolecular split into two disconnected fragments `(mol_u, new_u, mol_v, new_v)`.
    TwoFragments(MolGraph, usize, MolGraph, usize),
    /// Intramolecular ring-opening leaving one connected fragment `(mol_all, new_u, new_v)`.
    OneRingOpened(MolGraph, usize, usize),
}

fn cleave_bond(mol: &MolGraph, bid: usize, u: usize, v: usize) -> Option<BondCleavage> {
    let mut comp_u = HashSet::new();
    let mut q = VecDeque::new();
    comp_u.insert(u);
    q.push_back(u);
    while let Some(cur) = q.pop_front() {
        for &(nb, nb_bid) in &mol.adj[cur] {
            if nb_bid == bid {
                continue;
            }
            if comp_u.insert(nb) {
                q.push_back(nb);
            }
        }
    }

    if comp_u.contains(&v) {
        // Still connected -> intramolecular ring-opening (only allow for macrocycles/rings >= 7)
        let keep_all: HashSet<usize> = (0..mol.atoms.len()).collect();
        let (g, map) = extract_subgraph(mol, &keep_all, Some(bid));
        Some(BondCleavage::OneRingOpened(g, map[u]?, map[v]?))
    } else {
        let mut comp_v = HashSet::new();
        for i in 0..mol.atoms.len() {
            if !comp_u.contains(&i) {
                comp_v.insert(i);
            }
        }
        if comp_u.is_empty() || comp_v.is_empty() {
            return None;
        }
        let (gu, map_u) = extract_subgraph(mol, &comp_u, Some(bid));
        let (gv, map_v) = extract_subgraph(mol, &comp_v, Some(bid));
        Some(BondCleavage::TwoFragments(gu, map_u[u]?, gv, map_v[v]?))
    }
}

/// Helper predicates on `MolGraph`
fn has_double_bonded_o(mol: &MolGraph, c_idx: usize) -> Option<usize> {
    mol.adj.get(c_idx)?.iter().find_map(|&(nb, bid)| {
        if mol.atoms.get(nb).is_some_and(|a| a.atomic_number == 8)
            && mol
                .bonds
                .get(bid)
                .is_some_and(|b| b.order == BondOrder::Double)
        {
            Some(nb)
        } else {
            None
        }
    })
}

fn count_double_bonded_o(mol: &MolGraph, idx: usize) -> usize {
    mol.adj.get(idx).map_or(0, |nbrs| {
        nbrs.iter()
            .filter(|&&(nb, bid)| {
                mol.atoms.get(nb).is_some_and(|a| a.atomic_number == 8)
                    && mol
                        .bonds
                        .get(bid)
                        .is_some_and(|b| b.order == BondOrder::Double)
            })
            .count()
    })
}

fn is_electron_deficient_hetaryl_carbon(mol: &MolGraph, c_idx: usize) -> bool {
    let Some(a) = mol.atoms.get(c_idx) else {
        return false;
    };
    if a.atomic_number != 6 || !a.aromatic {
        return false;
    }
    // Check if adjacent to an aromatic ring nitrogen (2/4-position of pyridine, pyrimidine, purine, quinazoline)
    mol.rings.iter().any(|ring| {
        ring.len() == 6
            && ring.contains(&c_idx)
            && ring.iter().any(|&ri| {
                mol.atoms
                    .get(ri)
                    .is_some_and(|x| x.atomic_number == 7 && x.aromatic)
            })
    }) && mol.adj[c_idx].iter().any(|&(nb, _)| {
        mol.atoms
            .get(nb)
            .is_some_and(|x| x.atomic_number == 7 && x.aromatic)
    })
}

fn is_carbonyl_carbon(mol: &MolGraph, idx: usize) -> bool {
    mol.atoms
        .get(idx)
        .is_some_and(|a| a.atomic_number == 6 && !a.aromatic)
        && has_double_bonded_o(mol, idx).is_some()
}

fn is_sulfonyl_sulfur(mol: &MolGraph, idx: usize) -> bool {
    mol.atoms
        .get(idx)
        .is_some_and(|a| a.atomic_number == 16 && !a.aromatic)
        && count_double_bonded_o(mol, idx) >= 2
}

fn is_sp3_carbon(mol: &MolGraph, idx: usize) -> bool {
    let Some(a) = mol.atoms.get(idx) else {
        return false;
    };
    if a.atomic_number != 6 || a.aromatic {
        return false;
    }
    mol.adj[idx].iter().all(|&(_, bid)| {
        mol.bonds
            .get(bid)
            .is_some_and(|b| b.order == BondOrder::Single)
    })
}

fn is_alkene_carbon(mol: &MolGraph, idx: usize) -> bool {
    let Some(a) = mol.atoms.get(idx) else {
        return false;
    };
    if a.atomic_number != 6 || a.aromatic {
        return false;
    }
    mol.adj[idx].iter().any(|&(nb, bid)| {
        mol.atoms
            .get(nb)
            .is_some_and(|x| x.atomic_number == 6 && !x.aromatic)
            && mol
                .bonds
                .get(bid)
                .is_some_and(|b| b.order == BondOrder::Double)
    })
}

fn is_alkyne_carbon(mol: &MolGraph, idx: usize) -> bool {
    let Some(a) = mol.atoms.get(idx) else {
        return false;
    };
    if a.atomic_number != 6 || a.aromatic {
        return false;
    }
    mol.adj[idx].iter().any(|&(nb, bid)| {
        mol.atoms.get(nb).is_some_and(|x| x.atomic_number == 6)
            && mol
                .bonds
                .get(bid)
                .is_some_and(|b| b.order == BondOrder::Triple)
    })
}

fn is_macrocycle_bond(mol: &MolGraph, u: usize, v: usize) -> bool {
    // Check if every ring containing both u and v has size >= 8 (macrocycle, not a 5/6-ring!)
    let shared_rings: Vec<usize> = mol
        .rings
        .iter()
        .filter(|r| r.contains(&u) && r.contains(&v))
        .map(Vec::len)
        .collect();
    !shared_rings.is_empty() && shared_rings.iter().all(|&sz| sz >= 8)
}

/// Checks competing functional groups on precursors and emits actionable protecting-group advice.
fn audit_chemoselectivity(reaction_class: &str, precursors: &[MolGraph]) -> Vec<String> {
    let mut warnings = Vec::new();
    let mut free_nh2_count = 0usize;
    let mut free_cooh_count = 0usize;
    let mut free_phenol_count = 0usize;
    let mut nitro_or_azide = false;

    for g in precursors {
        for (i, a) in g.atoms.iter().enumerate() {
            if a.atomic_number == 7 && !a.aromatic && a.total_h() >= 1 {
                let adj_co = g.adj[i]
                    .iter()
                    .any(|&(nb, _)| is_carbonyl_carbon(g, nb) || is_sulfonyl_sulfur(g, nb));
                if !adj_co {
                    free_nh2_count += 1;
                }
            }
            if is_carbonyl_carbon(g, i) {
                let has_oh = g.adj[i].iter().any(|&(nb, bid)| {
                    g.atoms
                        .get(nb)
                        .is_some_and(|x| x.atomic_number == 8 && x.total_h() >= 1)
                        && g.bonds
                            .get(bid)
                            .is_some_and(|b| b.order == BondOrder::Single)
                });
                if has_oh {
                    free_cooh_count += 1;
                }
            }
            if a.atomic_number == 8
                && a.total_h() >= 1
                && g.adj[i]
                    .iter()
                    .any(|&(nb, _)| g.atoms.get(nb).is_some_and(|x| x.aromatic))
            {
                free_phenol_count += 1;
            }
            if a.atomic_number == 7 && a.formal_charge > 0 {
                nitro_or_azide = true;
            }
        }
    }

    match reaction_class {
        "Amide Coupling" => {
            if free_nh2_count >= 2 {
                warnings.push("Multiple nucleophilic amines detected: protect competing amine as Boc (Boc2O, Et3N) and deprotect with TFA/DCM".to_string());
            }
            if free_cooh_count >= 2 {
                warnings.push("Multiple carboxylic acids detected: protect non-reacting acid as methyl/t-Bu ester".to_string());
            }
        }
        "Williamson / Mitsunobu" | "N-Alkylation" => {
            if free_nh2_count >= 2 || (free_nh2_count >= 1 && free_phenol_count >= 1) {
                warnings.push("Competing N-H / O-H nucleophiles: use Boc/TBS orthogonal protection for selective alkylation".to_string());
            }
        }
        "Suzuki-Miyaura" | "Buchwald-Hartwig" | "Sonogashira" => {
            if free_cooh_count >= 1 {
                warnings.push("Free carboxylic acid present: ensure excess base (K2CO3/Cs2CO3) or use methyl ester during Pd coupling".to_string());
            }
        }
        "Grignard Addition" => {
            if free_nh2_count >= 1 || free_cooh_count >= 1 || free_phenol_count >= 1 {
                warnings.push("Protic groups (-OH/-NH/-COOH) quench organomagnesium reagents: protect as TBS/Boc first".to_string());
            }
        }
        "Nitro Reduction"
            if nitro_or_azide
                && precursors.iter().any(|g| {
                    g.atoms
                        .iter()
                        .any(|a| a.atomic_number == 35 || a.atomic_number == 53)
                }) =>
        {
            warnings.push("Aryl Br/I present: use Fe/NH4Cl or SnCl2 instead of H2/Pd-C to avoid hydrodehalogenation".to_string());
        }
        _ => {}
    }

    warnings
}

#[allow(clippy::too_many_arguments)]
fn push_disconnection(
    out: &mut Vec<Disconnection>,
    seen: &mut HashSet<String>,
    target_key: &str,
    template_id: &str,
    reaction_name: &str,
    reaction_class: &str,
    raw_precursors: Vec<MolGraph>,
    reagents: &str,
    solvent: &str,
    temperature: &str,
    est_yield: f64,
    workup_notes: &str,
) {
    let mut final_mols = Vec::with_capacity(raw_precursors.len());
    let mut final_smis = Vec::with_capacity(raw_precursors.len());
    for g in raw_precursors {
        let Some((pg, smi)) = finalize_precursor(g) else {
            return;
        };
        if pg.structural_key() == target_key {
            return; // Reject no-op self-loops
        }
        final_mols.push(pg);
        final_smis.push(smi);
    }

    let mut dedup_smis = final_smis.clone();
    dedup_smis.sort();
    let sig = format!("{template_id}::{}", dedup_smis.join("."));
    if !seen.insert(sig) {
        return;
    }

    let warnings = audit_chemoselectivity(reaction_class, &final_mols);
    let pg_pen = 0.22 * (warnings.len() as f64);
    let step_cost = (-est_yield.clamp(0.15, 0.96).ln()) + pg_pen;

    out.push(Disconnection {
        template_id: template_id.to_string(),
        reaction_name: reaction_name.to_string(),
        reaction_class: reaction_class.to_string(),
        precursor_smiles: final_smis,
        precursors: final_mols,
        reagents: reagents.to_string(),
        solvent: solvent.to_string(),
        temperature: temperature.to_string(),
        est_yield,
        step_cost,
        workup_notes: workup_notes.to_string(),
        warnings,
    });
}

/// Applies all 25+ curated medicinal-chemistry retrosynthetic disconnection templates directly on
/// the molecular graph `mol`, returning chemically validated, deduplicated precursor sets.
#[must_use]
pub fn enumerate_disconnections(mol: &MolGraph) -> Vec<Disconnection> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    let target_key = mol.structural_key();
    let total_ha = mol.heavy_atom_count();
    if total_ha < 4 {
        return out;
    }

    // --- PART A: Bond-cleavage disconnections (Intermolecular + Macrocyclization) ---
    for (bid, bond) in mol.bonds.iter().enumerate() {
        // Exocyclic biaryl bonds in canonical SMILES (e.g. `c1ccccc1c2ccccc2`) may carry
        // `BondOrder::Aromatic` with `in_ring == false`; treat them as single exocyclic C-C bonds.
        let is_single_or_exocyclic_biaryl =
            bond.order == BondOrder::Single || (bond.order == BondOrder::Aromatic && !bond.in_ring);
        if !is_single_or_exocyclic_biaryl {
            continue;
        }
        let u = bond.u;
        let v = bond.v;
        let is_macro = bond.in_ring && is_macrocycle_bond(mol, u, v);
        if bond.in_ring && !is_macro {
            continue; // Handled in Part B (ring-forming cycloadditions)
        }

        // Check both orientations `(a, b)` = `(u, v)` and `(v, u)`
        for (a, b) in [(u, v), (v, u)] {
            let za = mol.atoms[a].atomic_number;
            let zb = mol.atoms[b].atomic_number;

            // 1. Amide Coupling / Urea / Carbamate: C(=O) - N
            if is_carbonyl_carbon(mol, a)
                && zb == 7
                && !mol.atoms[b].aromatic
                && mol.atoms[b].formal_charge == 0
            {
                let other_n = mol.adj[a]
                    .iter()
                    .any(|&(nb, _)| nb != b && mol.atoms[nb].atomic_number == 7);
                let other_o_single = mol.adj[a].iter().any(|&(nb, bid2)| {
                    nb != b
                        && mol.atoms[nb].atomic_number == 8
                        && mol
                            .bonds
                            .get(bid2)
                            .is_some_and(|x| x.order == BondOrder::Single)
                });
                if !other_n && !other_o_single {
                    if let Some(cleave) = cleave_bond(mol, bid, a, b) {
                        match cleave {
                            BondCleavage::TwoFragments(mut g_ac, ac_i, mut g_am, am_i) => {
                                if attach_cap(&mut g_ac, ac_i, CapGroup::Oh)
                                    && attach_cap(&mut g_am, am_i, CapGroup::NoneH)
                                {
                                    push_disconnection(
                                        &mut out,
                                        &mut seen,
                                        &target_key,
                                        "AMIDE_HATU",
                                        "Amide Coupling (HATU / DIPEA)",
                                        "Amide Coupling",
                                        vec![g_ac, g_am],
                                        "HATU (1.2 eq), DIPEA (2.5 eq)",
                                        "DMF",
                                        "0 -> 25 °C, 4 h",
                                        0.88,
                                        "Quench with H2O, extract EtOAc, wash sat. NaHCO3 & brine, silica flash",
                                    );
                                }
                            }
                            BondCleavage::OneRingOpened(mut g, ac_i, am_i) => {
                                if attach_cap(&mut g, ac_i, CapGroup::Oh)
                                    && attach_cap(&mut g, am_i, CapGroup::NoneH)
                                {
                                    push_disconnection(
                                        &mut out,
                                        &mut seen,
                                        &target_key,
                                        "MACROLACTAM",
                                        "Macrolactamization (HATU, High Dilution)",
                                        "Amide Coupling",
                                        vec![g],
                                        "HATU (1.3 eq), DIPEA (3.0 eq)",
                                        "DMF (2 mM high dilution)",
                                        "25 °C, 16 h",
                                        0.68,
                                        "Slow syringe-pump addition, concentrate, reverse-phase or silica prep",
                                    );
                                }
                            }
                        }
                    }
                } else if other_n && !bond.in_ring {
                    // Urea: R1-NH-C(=O)-NR2R3 -> R1-NH2 + HNR2R3 via CDI
                    let Some((other_n_idx, other_bid)) = mol.adj[a]
                        .iter()
                        .copied()
                        .find(|&(nb, _)| nb != b && mol.atoms[nb].atomic_number == 7)
                    else {
                        continue;
                    };
                    if let Some(BondCleavage::TwoFragments(_, _, mut g_b, bi)) =
                        cleave_bond(mol, bid, a, b)
                        && let Some(BondCleavage::TwoFragments(_, _, mut g_n2, n2i)) =
                            cleave_bond(mol, other_bid, a, other_n_idx)
                        && attach_cap(&mut g_b, bi, CapGroup::NoneH)
                        && attach_cap(&mut g_n2, n2i, CapGroup::NoneH)
                    {
                        push_disconnection(
                            &mut out,
                            &mut seen,
                            &target_key,
                            "UREA_CDI",
                            "Urea Formation (CDI / Triphosgene)",
                            "Urea Formation",
                            vec![g_n2, g_b],
                            "CDI (1.15 eq), DIPEA (2.0 eq)",
                            "DCM / THF",
                            "0 -> 25 °C, 6 h",
                            0.84,
                            "Activate first amine with CDI (1 h), add second amine, wash 1M HCl/NaHCO3",
                        );
                    }
                } else if other_o_single && !bond.in_ring {
                    // Carbamate: R1-O-C(=O)-NR2R3 -> R1-OH + HNR2R3 via CDI
                    let Some((o_idx, o_bid)) = mol.adj[a].iter().copied().find(|&(nb, bid2)| {
                        nb != b
                            && mol.atoms[nb].atomic_number == 8
                            && mol
                                .bonds
                                .get(bid2)
                                .is_some_and(|x| x.order == BondOrder::Single)
                    }) else {
                        continue;
                    };
                    if let Some(BondCleavage::TwoFragments(_, _, mut g_am, bi)) =
                        cleave_bond(mol, bid, a, b)
                        && let Some(BondCleavage::TwoFragments(_, _, mut g_al, oi)) =
                            cleave_bond(mol, o_bid, a, o_idx)
                        && attach_cap(&mut g_am, bi, CapGroup::NoneH)
                        && attach_cap(&mut g_al, oi, CapGroup::NoneH)
                    {
                        push_disconnection(
                            &mut out,
                            &mut seen,
                            &target_key,
                            "CARBAMATE_CDI",
                            "Carbamate Formation (CDI)",
                            "Carbamate Formation",
                            vec![g_al, g_am],
                            "CDI (1.2 eq), DMAP (0.1 eq)",
                            "THF / MeCN",
                            "25 -> 50 °C, 8 h",
                            0.82,
                            "Form imidazolide from alcohol + CDI, add amine, aqueous workup",
                        );
                    }
                }
            }

            // 2. Esterification: C(=O) - O - R (non-acid, non-carbamate)
            if is_carbonyl_carbon(mol, a)
                && zb == 8
                && !mol.atoms[b].aromatic
                && mol.atoms[b].total_h() == 0
                && !mol.adj[a]
                    .iter()
                    .any(|&(nb, _)| mol.atoms[nb].atomic_number == 7)
                && let Some(cleave) = cleave_bond(mol, bid, a, b)
            {
                match cleave {
                    BondCleavage::TwoFragments(mut g_ac, ac_i, mut g_al, al_i) => {
                        if attach_cap(&mut g_ac, ac_i, CapGroup::Oh)
                            && attach_cap(&mut g_al, al_i, CapGroup::NoneH)
                        {
                            push_disconnection(
                                &mut out,
                                &mut seen,
                                &target_key,
                                "ESTER_EDC",
                                "Steglich Esterification (EDC / DMAP)",
                                "Esterification",
                                vec![g_ac, g_al],
                                "EDC·HCl (1.25 eq), DMAP (0.15 eq)",
                                "DCM",
                                "0 -> 25 °C, 6 h",
                                0.86,
                                "Wash 1M HCl, sat. NaHCO3, brine; dry Na2SO4 and purify",
                            );
                        }
                    }
                    BondCleavage::OneRingOpened(mut g, ac_i, al_i) => {
                        if attach_cap(&mut g, ac_i, CapGroup::Oh)
                            && attach_cap(&mut g, al_i, CapGroup::NoneH)
                        {
                            push_disconnection(
                                &mut out,
                                &mut seen,
                                &target_key,
                                "MACROLACTONE",
                                "Yamaguchi Macrolactonization",
                                "Esterification",
                                vec![g],
                                "2,4,6-Trichlorobenzoyl chloride, Et3N, DMAP",
                                "Toluene (high dilution)",
                                "25 -> 80 °C, 12 h",
                                0.65,
                                "Slow addition into DMAP/toluene, quench NaHCO3",
                            );
                        }
                    }
                }
            }

            // 3. Sulfonamide Coupling: S(=O)2 - N
            if is_sulfonyl_sulfur(mol, a)
                && zb == 7
                && !mol.atoms[b].aromatic
                && let Some(BondCleavage::TwoFragments(mut g_s, si, mut g_n, ni)) =
                    cleave_bond(mol, bid, a, b)
                && attach_cap(&mut g_s, si, CapGroup::Cl)
                && attach_cap(&mut g_n, ni, CapGroup::NoneH)
            {
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "SULFONAMIDE",
                    "Sulfonamide Formation (R-SO2Cl + Amine)",
                    "Sulfonamide Coupling",
                    vec![g_s, g_n],
                    "Et3N (2.0 eq) or Pyridine",
                    "DCM",
                    "0 -> 25 °C, 3 h",
                    0.90,
                    "Add sulfonyl chloride dropwise at 0 °C, wash 1M HCl & brine",
                );
            }

            // 4. SNAr Heteroaryl Amination & Buchwald-Hartwig C-N Coupling: Ar - N
            if za == 6
                && mol.atoms[a].aromatic
                && zb == 7
                && !mol.atoms[b].aromatic
                && mol.atoms[b].formal_charge == 0
                && !bond.in_ring
            {
                let n_is_amide_or_oxidized = mol.adj[b].iter().any(|&(nb, _)| {
                    is_carbonyl_carbon(mol, nb)
                        || is_sulfonyl_sulfur(mol, nb)
                        || mol.atoms[nb].atomic_number == 8
                });
                if !n_is_amide_or_oxidized
                    && let Some(BondCleavage::TwoFragments(mut g_ar, ari, mut g_am, ami)) =
                        cleave_bond(mol, bid, a, b)
                {
                    if is_electron_deficient_hetaryl_carbon(mol, a) {
                        if attach_cap(&mut g_ar, ari, CapGroup::Cl)
                            && attach_cap(&mut g_am, ami, CapGroup::NoneH)
                        {
                            push_disconnection(
                                &mut out,
                                &mut seen,
                                &target_key,
                                "SNAR_AMINE",
                                "SNAr Heteroaryl Amination (HetAr-Cl + Amine)",
                                "SNAr Amination",
                                vec![g_ar, g_am],
                                "DIPEA (2.5 eq) or K2CO3",
                                "DMSO or NMP",
                                "90–110 °C, 8 h",
                                0.87,
                                "Dilute with H2O, filter precipitate or extract EtOAc, flash column",
                            );
                        }
                    } else if attach_cap(&mut g_ar, ari, CapGroup::Br)
                        && attach_cap(&mut g_am, ami, CapGroup::NoneH)
                    {
                        push_disconnection(
                            &mut out,
                            &mut seen,
                            &target_key,
                            "BUCHWALD_HARTWIG",
                            "Buchwald-Hartwig C-N Coupling (Pd / XPhos)",
                            "Buchwald-Hartwig",
                            vec![g_ar, g_am],
                            "Pd2(dba)3 (2 mol%), XPhos (4 mol%), Cs2CO3 (2.0 eq)",
                            "1,4-Dioxane (degassed N2)",
                            "100 °C, 12 h",
                            0.80,
                            "Filter through Celite (EtOAc), concentrate, silica flash",
                        );
                    }
                }
            }

            // 5. SNAr Heteroaryl Etherification & Williamson Ether: Ar/R - O - C(sp3)
            if za == 8 && !mol.atoms[a].aromatic && zb == 6 && is_sp3_carbon(mol, b) {
                let other_is_carbonyl = mol.adj[a]
                    .iter()
                    .any(|&(nb, _)| nb != b && is_carbonyl_carbon(mol, nb));
                if !other_is_carbonyl && let Some(cleave) = cleave_bond(mol, bid, a, b) {
                    match cleave {
                        BondCleavage::TwoFragments(mut g_o, oi, mut g_c, ci) => {
                            if attach_cap(&mut g_o, oi, CapGroup::NoneH)
                                && attach_cap(&mut g_c, ci, CapGroup::Br)
                            {
                                let ar_o = mol.adj[a]
                                    .iter()
                                    .any(|&(nb, _)| nb != b && mol.atoms[nb].aromatic);
                                let (tid, rname, reag, yld) = if ar_o {
                                    (
                                        "WILLIAMSON_PHENOL",
                                        "Williamson Aryl Ether Synthesis (Phenol + R-Br)",
                                        "K2CO3 (2.0 eq) or Cs2CO3",
                                        0.87,
                                    )
                                } else {
                                    (
                                        "WILLIAMSON_ALKYL",
                                        "Williamson Ether Synthesis (Alcohol + R-Br)",
                                        "NaH (1.2 eq, 60% oil)",
                                        0.79,
                                    )
                                };
                                push_disconnection(
                                    &mut out,
                                    &mut seen,
                                    &target_key,
                                    tid,
                                    rname,
                                    "Williamson / Mitsunobu",
                                    vec![g_o, g_c],
                                    reag,
                                    "DMF (anhydrous)",
                                    "0 -> 60 °C, 6 h",
                                    yld,
                                    "Quench sat. NH4Cl, extract EtOAc, wash H2O/brine, dry Na2SO4",
                                );
                            }
                        }
                        BondCleavage::OneRingOpened(mut g, oi, ci) => {
                            if attach_cap(&mut g, oi, CapGroup::NoneH)
                                && attach_cap(&mut g, ci, CapGroup::Br)
                            {
                                push_disconnection(
                                    &mut out,
                                    &mut seen,
                                    &target_key,
                                    "MACRO_ETHER",
                                    "Intramolecular Williamson Macroetherification",
                                    "Williamson / Mitsunobu",
                                    vec![g],
                                    "Cs2CO3 (2.5 eq)",
                                    "DMF (3 mM high dilution)",
                                    "65 °C, 16 h",
                                    0.72,
                                    "High-dilution syringe addition to suppress dimerization, silica flash",
                                );
                            }
                        }
                    }
                }
            }

            // 5b. SNAr Heteroaryl-O-R: HetAr(N) - O
            if is_electron_deficient_hetaryl_carbon(mol, a)
                && zb == 8
                && !mol.atoms[b].aromatic
                && !bond.in_ring
                && let Some(BondCleavage::TwoFragments(mut g_ar, ari, mut g_o, oi)) =
                    cleave_bond(mol, bid, a, b)
                && attach_cap(&mut g_ar, ari, CapGroup::Cl)
                && attach_cap(&mut g_o, oi, CapGroup::NoneH)
            {
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "SNAR_ETHER",
                    "SNAr Heteroaryl Etherification (HetAr-Cl + ROH)",
                    "SNAr Etherification",
                    vec![g_ar, g_o],
                    "Cs2CO3 (2.0 eq) or NaH (1.2 eq)",
                    "DMF or THF",
                    "60–80 °C, 6 h",
                    0.85,
                    "Quench H2O, extract EtOAc, wash brine, silica flash",
                );
            }

            // 6. Heteroaryl / Indole / Pyrazole / Imidazole / Maleimide N-Alkylation: HetN - C(sp3)
            if za == 7
                && is_sp3_carbon(mol, b)
                && (mol.atoms[a].aromatic
                    || (mol.atoms[a].in_ring
                        && mol.adj[a]
                            .iter()
                            .any(|&(nb, _)| is_carbonyl_carbon(mol, nb))))
            {
                let cap_n = if mol.atoms[a].aromatic {
                    CapGroup::AromNH
                } else {
                    CapGroup::NoneH
                };
                if let Some(cleave) = cleave_bond(mol, bid, a, b) {
                    match cleave {
                        BondCleavage::TwoFragments(mut g_n, ni, mut g_c, ci) => {
                            if attach_cap(&mut g_n, ni, cap_n)
                                && attach_cap(&mut g_c, ci, CapGroup::Br)
                            {
                                push_disconnection(
                                    &mut out,
                                    &mut seen,
                                    &target_key,
                                    "HET_N_ALKYLATION",
                                    "Heterocycle / Indole N-Alkylation (Het-NH + R-Br)",
                                    "N-Alkylation",
                                    vec![g_n, g_c],
                                    "Cs2CO3 (2.0 eq) or NaH (1.2 eq)",
                                    "DMF (anhydrous)",
                                    "0 -> 60 °C, 8 h",
                                    0.84,
                                    "Quench ice/NH4Cl, extract EtOAc, wash brine, silica flash",
                                );
                            }
                        }
                        BondCleavage::OneRingOpened(mut g, ni, ci) => {
                            if attach_cap(&mut g, ni, cap_n) && attach_cap(&mut g, ci, CapGroup::Br)
                            {
                                push_disconnection(
                                    &mut out,
                                    &mut seen,
                                    &target_key,
                                    "MACRO_N_ALKYLATION",
                                    "Intramolecular N-Alkylation Macrocyclization (Cs2CO3)",
                                    "N-Alkylation",
                                    vec![g],
                                    "Cs2CO3 (2.5 eq), KI (0.2 eq)",
                                    "DMF (2–5 mM high dilution)",
                                    "60–75 °C, 16 h",
                                    0.74,
                                    "High-dilution macrocyclization (Ruboxistaurin/bisindolylmaleimide protocol), flash column",
                                );
                            }
                        }
                    }
                }
            }

            // 7. Reductive Amination: N(non-amide, non-arom) - CH2/CH(R)
            if za == 7
                && !mol.atoms[a].aromatic
                && mol.atoms[a].formal_charge == 0
                && is_sp3_carbon(mol, b)
                && mol.atoms[b].total_h() >= 1
                && !bond.in_ring
            {
                let is_amide_or_sulf = mol.adj[a].iter().any(|&(nb, _)| {
                    is_carbonyl_carbon(mol, nb)
                        || is_sulfonyl_sulfur(mol, nb)
                        || mol.atoms[nb].atomic_number == 8
                });
                if !is_amide_or_sulf
                    && let Some(BondCleavage::TwoFragments(mut g_am, ami, mut g_c, ci)) =
                        cleave_bond(mol, bid, a, b)
                    && attach_cap(&mut g_am, ami, CapGroup::NoneH)
                    && attach_cap(&mut g_c, ci, CapGroup::CarbonylO)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "REDUCTIVE_AMINATION",
                        "Reductive Amination (NaBH(OAc)3, Aldehyde/Ketone + Amine)",
                        "Reductive Amination",
                        vec![g_c, g_am],
                        "NaBH(OAc)3 (1.5 eq), AcOH (1.0 eq)",
                        "DCE or THF",
                        "25 °C, 12 h",
                        0.83,
                        "Quench sat. NaHCO3, extract DCM, dry Na2SO4, silica flash",
                    );
                }
            }

            // 8. Suzuki-Miyaura Biaryl Coupling: Ar1 - Ar2 (only check once per bond `a < b`!)
            if a < b
                && za == 6
                && zb == 6
                && mol.atoms[a].aromatic
                && mol.atoms[b].aromatic
                && !bond.in_ring
            {
                // Polarity 1: Ar1-Br + Ar2-B(OH)2
                if let Some(BondCleavage::TwoFragments(mut g1, i1, mut g2, i2)) =
                    cleave_bond(mol, bid, a, b)
                    && attach_cap(&mut g1, i1, CapGroup::Br)
                    && attach_cap(&mut g2, i2, CapGroup::BoronicAcid)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "SUZUKI_BIARYL",
                        "Suzuki-Miyaura Cross-Coupling (Ar-Br + Ar'-B(OH)2)",
                        "Suzuki-Miyaura",
                        vec![g1, g2],
                        "Pd(dppf)Cl2 (3 mol%), K2CO3 (2.5 eq)",
                        "1,4-Dioxane / H2O (4:1, N2)",
                        "90 °C, 8 h",
                        0.85,
                        "Filter Celite, partition EtOAc/H2O, wash brine, silica flash",
                    );
                }
                // Polarity 2: Ar1-B(OH)2 + Ar2-Br
                if let Some(BondCleavage::TwoFragments(mut g1, i1, mut g2, i2)) =
                    cleave_bond(mol, bid, a, b)
                    && attach_cap(&mut g1, i1, CapGroup::BoronicAcid)
                    && attach_cap(&mut g2, i2, CapGroup::Br)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "SUZUKI_BIARYL_REV",
                        "Suzuki-Miyaura Cross-Coupling (Ar-B(OH)2 + Ar'-Br)",
                        "Suzuki-Miyaura",
                        vec![g1, g2],
                        "Pd(dppf)Cl2 (3 mol%), K2CO3 (2.5 eq)",
                        "1,4-Dioxane / H2O (4:1, N2)",
                        "90 °C, 8 h",
                        0.85,
                        "Filter Celite, partition EtOAc/H2O, wash brine, silica flash",
                    );
                }
            }

            // 9. Sonogashira Coupling: Ar - C#C-R
            if za == 6
                && mol.atoms[a].aromatic
                && is_alkyne_carbon(mol, b)
                && !bond.in_ring
                && let Some(BondCleavage::TwoFragments(mut g_ar, ari, mut g_alk, alki)) =
                    cleave_bond(mol, bid, a, b)
                && attach_cap(&mut g_ar, ari, CapGroup::I)
                && attach_cap(&mut g_alk, alki, CapGroup::NoneH)
            {
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "SONOGASHIRA",
                    "Sonogashira Coupling (Ar-I + Terminal Alkyne)",
                    "Sonogashira",
                    vec![g_ar, g_alk],
                    "Pd(PPh3)2Cl2 (2 mol%), CuI (4 mol%), Et3N (3.0 eq)",
                    "THF / DMF (degassed)",
                    "25 -> 60 °C, 6 h",
                    0.83,
                    "Filter Celite, wash NH4Cl & brine, flash chromatography",
                );
            }

            // 10. Heck Alkenylation: Ar - CH=CH-R
            if za == 6
                && mol.atoms[a].aromatic
                && is_alkene_carbon(mol, b)
                && !bond.in_ring
                && let Some(BondCleavage::TwoFragments(mut g_ar, ari, mut g_alk, alki)) =
                    cleave_bond(mol, bid, a, b)
                && attach_cap(&mut g_ar, ari, CapGroup::Br)
                && attach_cap(&mut g_alk, alki, CapGroup::NoneH)
            {
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "HECK",
                    "Mizoroki-Heck Alkenylation (Ar-Br + Alkene)",
                    "Heck Coupling",
                    vec![g_ar, g_alk],
                    "Pd(OAc)2 (2 mol%), P(o-tol)3 (4 mol%), Et3N (2.5 eq)",
                    "DMF or MeCN",
                    "100 °C, 12 h",
                    0.78,
                    "Cool, filter Celite, extract EtOAc, flash column",
                );
            }

            // 11. Thioether Alkylation: S(divalent) - C(sp3)
            if za == 16
                && !mol.atoms[a].aromatic
                && count_double_bonded_o(mol, a) == 0
                && is_sp3_carbon(mol, b)
                && !bond.in_ring
                && let Some(BondCleavage::TwoFragments(mut g_s, si, mut g_c, ci)) =
                    cleave_bond(mol, bid, a, b)
                && attach_cap(&mut g_s, si, CapGroup::NoneH)
                && attach_cap(&mut g_c, ci, CapGroup::Br)
            {
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "THIOETHER_ALKYL",
                    "Thioether S-Alkylation (R-SH + R'-Br)",
                    "Thioether Alkylation",
                    vec![g_s, g_c],
                    "K2CO3 (2.0 eq)",
                    "DMF or Acetone",
                    "25 -> 50 °C, 4 h",
                    0.89,
                    "Filter salts, partition EtOAc/H2O, wash brine",
                );
            }

            // 12. Grignard / Organolithium Addition: Ar - CH(OH) - R
            if za == 6 && mol.atoms[a].aromatic && zb == 6 && is_sp3_carbon(mol, b) && !bond.in_ring
            {
                let oh_nbr = mol.adj[b].iter().copied().find(|&(nb, bid2)| {
                    mol.atoms[nb].atomic_number == 8
                        && mol.atoms[nb].total_h() >= 1
                        && mol
                            .bonds
                            .get(bid2)
                            .is_some_and(|x| x.order == BondOrder::Single)
                });
                if let Some((oh_idx, oh_bid)) = oh_nbr
                    && let Some(BondCleavage::TwoFragments(mut g_ar, ari, mut g_c, ci)) =
                        cleave_bond(mol, bid, a, b)
                {
                    // Convert C-OH in g_c into C=O by upgrading the C-O bond order
                    let _ = (oh_idx, oh_bid);
                    if let Some((_, c_o_bid)) = g_c.adj[ci].iter().copied().find(|&(nb, _)| {
                        g_c.atoms[nb].atomic_number == 8 && g_c.atoms[nb].total_h() >= 1
                    }) {
                        if let Some(bo) = g_c.bonds.get_mut(c_o_bid) {
                            bo.order = BondOrder::Double;
                        }
                        if attach_cap(&mut g_ar, ari, CapGroup::Br) {
                            push_disconnection(
                                &mut out,
                                &mut seen,
                                &target_key,
                                "GRIGNARD_ADD",
                                "Grignard Addition (Ar-MgBr from Ar-Br + Aldehyde/Ketone)",
                                "Grignard Addition",
                                vec![g_ar, g_c],
                                "Mg turnings (1.2 eq) or i-PrMgCl·LiCl (1.1 eq)",
                                "THF (anhydrous)",
                                "-78 -> 0 °C, 3 h",
                                0.77,
                                "Quench sat. NH4Cl at 0 °C, extract EtOAc, silica flash",
                            );
                        }
                    }
                }
            }

            // 13. Friedel-Crafts / Weinreb Ketone Acylation: Ar - C(=O) - R
            if za == 6 && mol.atoms[a].aromatic && is_carbonyl_carbon(mol, b) && !bond.in_ring {
                let is_acid_or_amide_or_ester = mol.adj[b].iter().any(|&(nb, bid2)| {
                    (mol.atoms[nb].atomic_number == 7 || mol.atoms[nb].atomic_number == 8)
                        && mol
                            .bonds
                            .get(bid2)
                            .is_some_and(|x| x.order == BondOrder::Single)
                });
                if !is_acid_or_amide_or_ester
                    && let Some(BondCleavage::TwoFragments(mut g_ar, ari, mut g_ac, aci)) =
                        cleave_bond(mol, bid, a, b)
                    && attach_cap(&mut g_ar, ari, CapGroup::NoneH)
                    && attach_cap(&mut g_ac, aci, CapGroup::Cl)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "FRIEDEL_CRAFTS",
                        "Friedel-Crafts Acylation (Ar-H + R-COCl)",
                        "Friedel-Crafts Acylation",
                        vec![g_ar, g_ac],
                        "AlCl3 (1.3 eq)",
                        "DCM or DCE",
                        "0 -> 40 °C, 4 h",
                        0.76,
                        "Pour onto ice/1M HCl, extract DCM, wash NaHCO3",
                    );
                }
            }
        }
    }

    // --- PART B: Heterocycle Click / Cycloaddition Disconnections ---
    for ring in &mol.rings {
        // 14. Tetrazole Cycloaddition: R - c1nnn[nH]1  <=  R - C#N + NaN3
        if ring.len() == 5 {
            let n_cnt = ring
                .iter()
                .filter(|&&i| mol.atoms[i].atomic_number == 7)
                .count();
            let c_cnt = ring
                .iter()
                .filter(|&&i| mol.atoms[i].atomic_number == 6)
                .count();
            if n_cnt == 4
                && c_cnt == 1
                && let Some(&ring_c) = ring.iter().find(|&&i| mol.atoms[i].atomic_number == 6)
                && let Some(&(ext_a, ext_bid)) =
                    mol.adj[ring_c].iter().find(|&&(nb, _)| !ring.contains(&nb))
                && let Some(BondCleavage::TwoFragments(mut g_ext, ei, _, _)) =
                    cleave_bond(mol, ext_bid, ext_a, ring_c)
                && attach_cap(&mut g_ext, ei, CapGroup::Cyano)
            {
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "TETRAZOLE_CLICK",
                    "Tetrazole Synthesis (Nitrile + NaN3 / ZnBr2)",
                    "Cycloaddition",
                    vec![g_ext],
                    "NaN3 (1.5 eq), ZnBr2 (1.0 eq)",
                    "i-PrOH / H2O (1:2)",
                    "100 °C, 18 h",
                    0.84,
                    "Acidify with 2M HCl to pH 2, extract EtOAc or filter solid",
                );
            }

            // 15. CuAAC Click 1,2,3-Triazole: R1-N3 + HC#C-R2
            if n_cnt == 3 && c_cnt == 2 {
                let mut n_ext = None;
                let mut c_ext = None;
                for &ri in ring {
                    for &(nb, bid) in &mol.adj[ri] {
                        if !ring.contains(&nb) {
                            if mol.atoms[ri].atomic_number == 7 {
                                n_ext = Some((nb, bid, ri));
                            } else if mol.atoms[ri].atomic_number == 6 {
                                c_ext = Some((nb, bid, ri));
                            }
                        }
                    }
                }
                if let (Some((na, nbid, rn)), Some((ca, cbid, rc))) = (n_ext, c_ext)
                    && let Some(BondCleavage::TwoFragments(mut g_az, azi, _, _)) =
                        cleave_bond(mol, nbid, na, rn)
                    && let Some(BondCleavage::TwoFragments(mut g_alk, alki, _, _)) =
                        cleave_bond(mol, cbid, ca, rc)
                {
                    // Attach -N3 to g_az and -C#CH to g_alk
                    let ok_az = attach_cap(&mut g_az, azi, CapGroup::Azide);
                    let ok_alk = {
                        if let (Ok(c1), Ok(c2)) = (
                            MolAtom::new(g_alk.atoms.len(), "C"),
                            MolAtom::new(g_alk.atoms.len() + 1, "C"),
                        ) {
                            let i1 = g_alk.add_atom(c1);
                            let i2 = g_alk.add_atom(c2);
                            g_alk
                                .add_bond(alki, i1, BondOrder::Single, BondStereo::None)
                                .is_ok()
                                && g_alk
                                    .add_bond(i1, i2, BondOrder::Triple, BondStereo::None)
                                    .is_ok()
                        } else {
                            false
                        }
                    };
                    if ok_az && ok_alk {
                        push_disconnection(
                            &mut out,
                            &mut seen,
                            &target_key,
                            "CUAAC_TRIAZOLE",
                            "CuAAC Click 1,2,3-Triazole (Azide + Alkyne)",
                            "Cycloaddition",
                            vec![g_az, g_alk],
                            "CuSO4·5H2O (5 mol%), Sodium Ascorbate (15 mol%)",
                            "t-BuOH / H2O (1:1)",
                            "25 °C, 8 h",
                            0.92,
                            "Dilute H2O/NH4OH, extract EtOAc, wash brine",
                        );
                    }
                }
            }
        }
    }

    // --- PART C: Strategic Functional Group Interconversions (FGIs) ---
    for (idx, atom) in mol.atoms.iter().enumerate() {
        // 16. Sulfone oxidation: R-SO2-R' <= R-S-R' (mCPBA)
        if is_sulfonyl_sulfur(mol, idx) {
            let mut has_only_c = true;
            let mut o_indices = Vec::new();
            for &(nb, bid) in &mol.adj[idx] {
                if mol.atoms[nb].atomic_number == 8 && mol.bonds[bid].order == BondOrder::Double {
                    o_indices.push(nb);
                } else if mol.atoms[nb].atomic_number != 6 {
                    has_only_c = false;
                }
            }
            if has_only_c && o_indices.len() == 2 {
                let keep: HashSet<usize> = (0..mol.atoms.len())
                    .filter(|i| !o_indices.contains(i))
                    .collect();
                let (g_thio, _) = extract_subgraph(mol, &keep, None);
                push_disconnection(
                    &mut out,
                    &mut seen,
                    &target_key,
                    "SULFONE_OX",
                    "Sulfide Oxidation to Sulfone (mCPBA / Oxone)",
                    "Oxidation (FGI)",
                    vec![g_thio],
                    "mCPBA (2.2 eq) or Oxone",
                    "DCM",
                    "0 -> 25 °C, 4 h",
                    0.91,
                    "Quench sat. Na2S2O3 & NaHCO3, extract DCM",
                );
            }
        }

        // 17. Aromatic primary amine reduction: Ar-NH2 <= Ar-NO2 (H2, Pd/C)
        if atom.atomic_number == 7
            && !atom.aromatic
            && atom.total_h() >= 2
            && mol.adj[idx].len() == 1
        {
            let (ar_i, bid) = mol.adj[idx][0];
            if mol.atoms[ar_i].aromatic && mol.atoms[ar_i].atomic_number == 6 {
                let keep: HashSet<usize> = (0..mol.atoms.len()).filter(|&i| i != idx).collect();
                let (mut g_no2, map) = extract_subgraph(mol, &keep, Some(bid));
                if let Some(new_ar) = map[ar_i]
                    && attach_cap(&mut g_no2, new_ar, CapGroup::Nitro)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "NITRO_REDUCTION",
                        "Nitro Reduction to Aniline (H2, Pd/C or Fe/NH4Cl)",
                        "Nitro Reduction",
                        vec![g_no2],
                        "10% Pd/C (5 wt%), H2 (1 atm) [or Fe / NH4Cl]",
                        "MeOH / EtOAc",
                        "25 °C, 4 h",
                        0.93,
                        "Filter through Celite under N2, concentrate in vacuo",
                    );
                }
            }
        }

        // 18. Secondary / primary alcohol reduction from ketone/aldehyde: R-CH(OH)-R' <= R-C(=O)-R' (NaBH4)
        if is_sp3_carbon(mol, idx)
            && atom.total_h() >= 1
            && let Some((o_idx, oh_bid)) = mol.adj[idx].iter().copied().find(|&(nb, bid)| {
                mol.atoms[nb].atomic_number == 8
                    && !mol.atoms[nb].aromatic
                    && mol.atoms[nb].total_h() == 1
                    && mol.adj[nb].len() == 1
                    && mol
                        .bonds
                        .get(bid)
                        .is_some_and(|b| b.order == BondOrder::Single)
            })
        {
            let mut g_ket = mol.clone();
            if let Some(bo) = g_ket.bonds.get_mut(oh_bid) {
                bo.order = BondOrder::Double;
            }
            if let Some(oa) = g_ket.atoms.get_mut(o_idx) {
                oa.from_bracket = false;
                oa.explicit_h = 0;
            }
            if let Some(ca) = g_ket.atoms.get_mut(idx) {
                ca.from_bracket = false;
                ca.explicit_h = 0;
            }
            push_disconnection(
                &mut out,
                &mut seen,
                &target_key,
                "KETONE_REDUCTION",
                "Ketone / Aldehyde Reduction (NaBH4 or H2/Raney Ni)",
                "Reduction (FGI)",
                vec![g_ket],
                "NaBH4 (1.2 eq)",
                "MeOH / THF",
                "0 -> 25 °C, 2 h",
                0.94,
                "Quench sat. NH4Cl at 0 °C, evaporate MeOH, extract EtOAc",
            );
        }

        // 19. Benzylic Carboxylic Acid via BHC Pd-Carbonylation & Nitrile Hydrolysis:
        //     Ar-CH(R)-C(=O)OH <= Ar-CH(OH)-R + CO (BHC Green Process) OR R-C#N (hydrolysis)
        if is_carbonyl_carbon(mol, idx) {
            let has_oh = mol.adj[idx].iter().any(|&(nb, bid)| {
                mol.atoms[nb].atomic_number == 8
                    && mol.atoms[nb].total_h() >= 1
                    && mol.adj[nb].len() == 1
                    && mol
                        .bonds
                        .get(bid)
                        .is_some_and(|b| b.order == BondOrder::Single)
            });
            if has_oh
                && let Some((c_adj, c_bid)) = mol.adj[idx].iter().copied().find(|&(nb, bid)| {
                    mol.atoms[nb].atomic_number == 6
                        && mol
                            .bonds
                            .get(bid)
                            .is_some_and(|b| b.order == BondOrder::Single && !b.in_ring)
                })
            {
                let is_benzylic = is_sp3_carbon(mol, c_adj)
                    && mol.adj[c_adj]
                        .iter()
                        .any(|&(nb, _)| mol.atoms[nb].atomic_number == 6 && mol.atoms[nb].aromatic);
                if is_benzylic
                    && let Some(BondCleavage::TwoFragments(_, _, mut g_alc, ai)) =
                        cleave_bond(mol, c_bid, idx, c_adj)
                    && attach_cap(&mut g_alc, ai, CapGroup::Oh)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "BENZYLIC_CARBONYLATION",
                        "Pd-Catalyzed Benzylic Carbonylation (BHC Green Process, R-CH(OH)Ar + CO)",
                        "Carbonylation",
                        vec![g_alc],
                        "CO (40 bar), PdCl2(PPh3)2 (0.5 mol%), HCl (cat.)",
                        "AcOH / H2O",
                        "130 °C, 4 h",
                        0.89,
                        "Vent CO safely, filter Pd catalyst, crystallize carboxylic acid from H2O/EtOH",
                    );
                }
                if (is_benzylic || mol.atoms[c_adj].aromatic)
                    && let Some(BondCleavage::TwoFragments(_, _, mut g_cn, ci)) =
                        cleave_bond(mol, c_bid, idx, c_adj)
                    && attach_cap(&mut g_cn, ci, CapGroup::Cyano)
                {
                    push_disconnection(
                        &mut out,
                        &mut seen,
                        &target_key,
                        "NITRILE_HYDROLYSIS",
                        "Nitrile Hydrolysis to Carboxylic Acid (HCl / H2O or NaOH)",
                        "Hydrolysis (FGI)",
                        vec![g_cn],
                        "6M HCl (aq) or 2M NaOH",
                        "EtOH / H2O (1:1)",
                        "90 °C, 8 h",
                        0.88,
                        "Adjust to pH 2, filter or extract EtOAc, dry Na2SO4",
                    );
                }
            }
        }
    }

    out.sort_by(|a, b| a.step_cost.total_cmp(&b.step_cost));
    out
}
