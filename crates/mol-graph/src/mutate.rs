use std::collections::BTreeSet;

use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, MolAtom, MolGraph};
use crate::smiles::parse_smiles;

/// Predefined medicinal chemistry bioisostere & substituent transforms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstituentKind {
    Methyl,
    Fluoro,
    Chloro,
    Bromo,
    Hydroxyl,
    Amino,
    Methoxy,
    Trifluoromethyl,
    Cyano,
    Carboxyl,
    Tetrazole,
}

impl SubstituentKind {
    pub fn from_name(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "methyl" | "me" | "ch3" => Some(Self::Methyl),
            "fluoro" | "f" => Some(Self::Fluoro),
            "chloro" | "cl" => Some(Self::Chloro),
            "bromo" | "br" => Some(Self::Bromo),
            "hydroxyl" | "oh" | "hydroxy" => Some(Self::Hydroxyl),
            "amino" | "nh2" => Some(Self::Amino),
            "methoxy" | "ome" | "och3" => Some(Self::Methoxy),
            "trifluoromethyl" | "cf3" => Some(Self::Trifluoromethyl),
            "cyano" | "cn" | "nitrile" => Some(Self::Cyano),
            "carboxyl" | "cooh" | "carboxylic_acid" => Some(Self::Carboxyl),
            "tetrazole" | "tetrazolyl" | "bioisostere_tetrazole" => Some(Self::Tetrazole),
            _ => None,
        }
    }

    pub fn fragment_smiles(self) -> &'static str {
        match self {
            Self::Methyl => "C",
            Self::Fluoro => "F",
            Self::Chloro => "Cl",
            Self::Bromo => "Br",
            Self::Hydroxyl => "O",
            Self::Amino => "N",
            Self::Methoxy => "OC",
            Self::Trifluoromethyl => "C(F)(F)F",
            Self::Cyano => "C#N",
            Self::Carboxyl => "C(=O)O",
            Self::Tetrazole => "c1nnn[nH]1",
        }
    }
}

/// Mutate a molecule by either:
/// - replacing a carboxylic acid (`-C(=O)OH` / `-C(=O)[O-]`) with a `1H-tetrazol-5-yl`
///   bioisostere (`group = "tetrazole"` / `replace_carboxyl = true`), OR
/// - changing an atom's element or attaching a named substituent onto `atom_index`.
pub fn mutate_molecule(
    mol: &MolGraph,
    atom_index: Option<usize>,
    group: &str,
) -> Result<MolGraph, MolError> {
    let kind = SubstituentKind::from_name(group).ok_or_else(|| {
        MolError::Mutation(format!("unsupported mutation/bioisostere group '{group}'"))
    })?;

    if kind == SubstituentKind::Tetrazole && atom_index.is_none() {
        return replace_carboxyl_with_tetrazole(mol);
    }

    // Check if user passed the carboxyl carbon index specifically for tetrazole bioisostere swap
    if kind == SubstituentKind::Tetrazole
        && let Some(idx) = atom_index
        && is_carboxyl_carbon(mol, idx).is_some()
    {
        return replace_specific_carboxyl_with_tetrazole(mol, idx);
    }

    let target = atom_index.unwrap_or(0);
    if target >= mol.atoms.len() {
        return Err(MolError::Mutation(format!(
            "atom_index {target} out of range (len={})",
            mol.atoms.len()
        )));
    }
    attach_fragment_at(mol, target, kind.fragment_smiles())
}

/// Replace the first carboxylic acid group (`R-C(=O)OH` or `R-C(=O)[O-]`) with `R-c1nnn[nH]1`.
pub fn replace_carboxyl_with_tetrazole(mol: &MolGraph) -> Result<MolGraph, MolError> {
    let Some(c_idx) = (0..mol.atoms.len()).find(|&idx| is_carboxyl_carbon(mol, idx).is_some())
    else {
        return Err(MolError::Mutation(
            "no carboxylic acid group (-C(=O)OH) found for tetrazole bioisostere replacement"
                .to_string(),
        ));
    };
    replace_specific_carboxyl_with_tetrazole(mol, c_idx)
}

fn is_carboxyl_carbon(mol: &MolGraph, c_idx: usize) -> Option<(usize, usize, Option<usize>)> {
    if c_idx >= mol.atoms.len() || mol.atoms[c_idx].symbol != "C" || mol.adj[c_idx].len() != 3 {
        return None;
    }
    let mut o_double: Option<usize> = None;
    let mut o_single: Option<usize> = None;
    let mut anchor: Option<usize> = None;

    for &(nbr, bid) in &mol.adj[c_idx] {
        if mol.atoms[nbr].symbol == "O" && mol.adj[nbr].len() == 1 {
            match mol.bonds[bid].order {
                BondOrder::Double => o_double = Some(nbr),
                BondOrder::Single => o_single = Some(nbr),
                _ => {}
            }
        } else {
            anchor = Some(nbr);
        }
    }

    match (o_double, o_single) {
        (Some(od), Some(os)) => Some((od, os, anchor)),
        _ => None,
    }
}

fn replace_specific_carboxyl_with_tetrazole(
    mol: &MolGraph,
    c_idx: usize,
) -> Result<MolGraph, MolError> {
    let (o_double, o_single, anchor_opt) = is_carboxyl_carbon(mol, c_idx).ok_or_else(|| {
        MolError::Mutation(format!("atom {c_idx} is not a carboxylic acid carbon"))
    })?;

    // Build filtered scaffold removing c_idx, o_double, o_single
    let remove_set: BTreeSet<usize> = [c_idx, o_double, o_single].into_iter().collect();
    let mut out = MolGraph::new();
    out.name = mol.name.clone();
    out.properties = mol.properties.clone();

    let mut remap: Vec<Option<usize>> = vec![None; mol.atoms.len()];
    for (old_id, atom) in mol.atoms.iter().enumerate() {
        if remove_set.contains(&old_id) {
            continue;
        }
        let mut cloned = atom.clone();
        cloned.coords_2d = None;
        let nid = out.add_atom(cloned);
        remap[old_id] = Some(nid);
    }
    for b in &mol.bonds {
        if let (Some(nu), Some(nv)) = (remap[b.u], remap[b.v]) {
            out.add_bond(nu, nv, b.order, b.stereo)?;
        }
    }

    let tet = parse_smiles("c1nnn[nH]1")?;
    let base = out.atoms.len();
    for a in &tet.atoms {
        let mut na = a.clone();
        na.coords_2d = None;
        out.add_atom(na);
    }
    for b in &tet.bonds {
        out.add_bond(base + b.u, base + b.v, b.order, b.stereo)?;
    }

    // Connect anchor (if any, e.g. phenyl ring) to carbon atom `0` of tetrazole (`c1nnn[nH]1`)
    if let Some(old_anchor) = anchor_opt
        && let Some(new_anchor) = remap[old_anchor]
    {
        let tet_c = (0..tet.atoms.len())
            .find(|&i| tet.atoms[i].symbol == "C")
            .unwrap_or(0);
        out.add_bond(
            new_anchor,
            base + tet_c,
            BondOrder::Single,
            BondStereo::None,
        )?;
    }

    out.perceive()?;
    out.input_smiles = Some(out.canonical_smiles());
    Ok(out)
}

/// Attach a SMILES fragment (connecting atom 0 of fragment to `target_atom` in `mol`).
pub fn attach_fragment_at(
    mol: &MolGraph,
    target_atom: usize,
    fragment_smiles: &str,
) -> Result<MolGraph, MolError> {
    if target_atom >= mol.atoms.len() {
        return Err(MolError::Mutation(format!(
            "target atom {target_atom} out of range"
        )));
    }
    let frag = parse_smiles(fragment_smiles)?;
    let mut out = mol.clone();
    for a in &mut out.atoms {
        a.coords_2d = None;
    }
    let offset = out.atoms.len();
    for a in &frag.atoms {
        let mut na = MolAtom::new(0, &a.symbol)?;
        na.aromatic = a.aromatic;
        na.formal_charge = a.formal_charge;
        na.explicit_h = a.explicit_h;
        na.from_bracket = a.from_bracket;
        out.add_atom(na);
    }
    for b in &frag.bonds {
        out.add_bond(offset + b.u, offset + b.v, b.order, b.stereo)?;
    }
    out.add_bond(target_atom, offset, BondOrder::Single, BondStereo::None)?;
    out.perceive()?;
    out.input_smiles = Some(out.canonical_smiles());
    Ok(out)
}
