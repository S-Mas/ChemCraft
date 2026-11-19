//! ChemDraw-equivalent Superatom / Functional Group / Protecting Group templates
//! (`Ph`, `Bn`, `Boc`, `Fmoc`, `Cbz`, `Ac`, `tBu`, `iPr`, `CF3`, `CO2H`, `CO2Me`, `SO2Me`, etc.)
//! that attach cleanly to any selected atom (with valence check & perception).

use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, MolGraph};
use crate::smiles::parse_smiles;

/// Metadata for a ChemDraw superatom functional group or protecting group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SuperatomGroup {
    pub abbr: &'static str,
    pub label: &'static str,
    pub category: &'static str,
    /// SMILES of the fragment whose atom index 0 is the attachment point.
    pub smiles: &'static str,
}

/// 26 standard organic & medicinal chemistry superatom groups matching ChemDraw hotkeys & nicknames.
pub const SUPERATOM_GROUPS: &[SuperatomGroup] = &[
    // Alkyl / Haloalkyl
    SuperatomGroup {
        abbr: "Me",
        label: "Methyl (-CH3)",
        category: "Alkyl",
        smiles: "C",
    },
    SuperatomGroup {
        abbr: "Et",
        label: "Ethyl (-CH2CH3)",
        category: "Alkyl",
        smiles: "CC",
    },
    SuperatomGroup {
        abbr: "iPr",
        label: "Isopropyl (-CH(CH3)2)",
        category: "Alkyl",
        smiles: "C(C)C",
    },
    SuperatomGroup {
        abbr: "tBu",
        label: "tert-Butyl (-C(CH3)3)",
        category: "Alkyl",
        smiles: "C(C)(C)C",
    },
    SuperatomGroup {
        abbr: "cPr",
        label: "Cyclopropyl",
        category: "Alkyl",
        smiles: "C1CC1",
    },
    SuperatomGroup {
        abbr: "CF3",
        label: "Trifluoromethyl (-CF3)",
        category: "Haloalkyl",
        smiles: "C(F)(F)F",
    },
    // Aryl / Heteroaryl
    SuperatomGroup {
        abbr: "Ph",
        label: "Phenyl (-C6H5)",
        category: "Aryl",
        smiles: "c1ccccc1",
    },
    SuperatomGroup {
        abbr: "Bn",
        label: "Benzyl (-CH2Ph)",
        category: "Aryl",
        smiles: "Cc1ccccc1",
    },
    SuperatomGroup {
        abbr: "Pyr",
        label: "4-Pyridyl",
        category: "Heteroaryl",
        smiles: "c1ccncc1",
    },
    SuperatomGroup {
        abbr: "Morph",
        label: "Morpholino",
        category: "Heterocycle",
        smiles: "N1CCOCC1",
    },
    SuperatomGroup {
        abbr: "Pip",
        label: "Piperazino",
        category: "Heterocycle",
        smiles: "N1CCNCC1",
    },
    // Oxygen / Carbonyl
    SuperatomGroup {
        abbr: "OMe",
        label: "Methoxy (-OCH3)",
        category: "Oxygen",
        smiles: "OC",
    },
    SuperatomGroup {
        abbr: "OEt",
        label: "Ethoxy (-OCH2CH3)",
        category: "Oxygen",
        smiles: "OCC",
    },
    SuperatomGroup {
        abbr: "Ac",
        label: "Acetyl (-C(=O)CH3)",
        category: "Carbonyl",
        smiles: "C(=O)C",
    },
    SuperatomGroup {
        abbr: "Bz",
        label: "Benzoyl (-C(=O)Ph)",
        category: "Carbonyl",
        smiles: "C(=O)c1ccccc1",
    },
    SuperatomGroup {
        abbr: "CO2H",
        label: "Carboxylic Acid (-COOH)",
        category: "Carbonyl",
        smiles: "C(=O)O",
    },
    SuperatomGroup {
        abbr: "CO2Me",
        label: "Methyl Ester (-COOMe)",
        category: "Carbonyl",
        smiles: "C(=O)OC",
    },
    SuperatomGroup {
        abbr: "CO2Et",
        label: "Ethyl Ester (-COOEt)",
        category: "Carbonyl",
        smiles: "C(=O)OCC",
    },
    SuperatomGroup {
        abbr: "CONH2",
        label: "Primary Amide (-CONH2)",
        category: "Carbonyl",
        smiles: "C(=O)N",
    },
    // Nitrogen / Sulfur / Bioisostere
    SuperatomGroup {
        abbr: "CN",
        label: "Cyano / Nitrile (-CN)",
        category: "Nitrogen",
        smiles: "C#N",
    },
    SuperatomGroup {
        abbr: "NO2",
        label: "Nitro (-NO2)",
        category: "Nitrogen",
        smiles: "[N+](=O)[O-]",
    },
    SuperatomGroup {
        abbr: "NMe2",
        label: "Dimethylamino (-NMe2)",
        category: "Nitrogen",
        smiles: "N(C)C",
    },
    SuperatomGroup {
        abbr: "SO2Me",
        label: "Mesyl (-SO2Me)",
        category: "Sulfur",
        smiles: "S(=O)(=O)C",
    },
    SuperatomGroup {
        abbr: "SO2NH2",
        label: "Sulfonamide (-SO2NH2)",
        category: "Sulfur",
        smiles: "S(=O)(=O)N",
    },
    SuperatomGroup {
        abbr: "OTf",
        label: "Triflate (-OSO2CF3)",
        category: "Leaving Group",
        smiles: "OS(=O)(=O)C(F)(F)F",
    },
    // Protecting Groups (ChemDraw peptide / synthesis staples)
    SuperatomGroup {
        abbr: "Boc",
        label: "tert-Butoxycarbonyl (Boc)",
        category: "Protecting",
        smiles: "C(=O)OC(C)(C)C",
    },
    SuperatomGroup {
        abbr: "Cbz",
        label: "Benzyloxycarbonyl (Cbz)",
        category: "Protecting",
        smiles: "C(=O)OCc1ccccc1",
    },
    SuperatomGroup {
        abbr: "Fmoc",
        label: "Fluorenylmethoxycarbonyl (Fmoc)",
        category: "Protecting",
        smiles: "C(=O)OCC1c2ccccc2-c2ccccc12",
    },
    SuperatomGroup {
        abbr: "Ts",
        label: "Tosyl (-SO2-C6H4-Me)",
        category: "Protecting",
        smiles: "S(=O)(=O)c1ccc(C)cc1",
    },
    SuperatomGroup {
        abbr: "TBS",
        label: "tert-Butyldimethylsilyl (TBS)",
        category: "Protecting",
        smiles: "[Si](C)(C)C(C)(C)C",
    },
];

/// Look up a superatom group by abbreviation (case-insensitive, e.g. `"Boc"`, `"co2me"`, `"Ph"`).
pub fn find_superatom_group(abbr: &str) -> Option<&'static SuperatomGroup> {
    let needle = abbr.trim();
    SUPERATOM_GROUPS
        .iter()
        .find(|g| g.abbr.eq_ignore_ascii_case(needle))
}

/// Expand a superatom group standalone as a capped molecule (`Me-<group>` or `<group>H`).
pub fn build_superatom_fragment(abbr: &str) -> Result<MolGraph, MolError> {
    let grp = find_superatom_group(abbr).ok_or_else(|| {
        MolError::Mutation(format!("Unknown superatom/functional group '{abbr}'"))
    })?;
    let mut frag = parse_smiles(grp.smiles)?;
    frag.name = Some(grp.abbr.to_string());
    Ok(frag)
}

/// Attach a ChemDraw superatom / protecting group directly onto `attach_atom` of `mol`.
/// If `mol` is empty, creates the standalone group molecule cleanly.
pub fn attach_superatom_group(
    mol: &MolGraph,
    attach_atom: Option<usize>,
    abbr: &str,
) -> Result<MolGraph, MolError> {
    let grp = find_superatom_group(abbr).ok_or_else(|| {
        MolError::Mutation(format!("Unknown superatom/functional group '{abbr}'"))
    })?;
    let frag = parse_smiles(grp.smiles)?;
    if mol.atoms.is_empty() {
        let mut out = frag;
        out.name = Some(grp.abbr.to_string());
        return Ok(out);
    }

    let anchor = attach_atom.unwrap_or_else(|| {
        // Default: last atom with at least one available hydrogen, else last atom
        mol.atoms
            .iter()
            .rev()
            .find(|a| a.total_h() > 0)
            .map_or(mol.atoms.len() - 1, |a| a.id)
    });
    if anchor >= mol.atoms.len() {
        return Err(MolError::Mutation(format!(
            "Attachment atom index {anchor} out of bounds (len={})",
            mol.atoms.len()
        )));
    }

    let mut out = mol.clone();
    // Clear stale coords so clean 2D/3D layout regenerates seamlessly
    for a in &mut out.atoms {
        a.coords_2d = None;
        a.coords_3d = None;
    }
    // If the anchor had explicit_h from bracket notation, decrement by 1 to accommodate the new bond
    if out.atoms[anchor].explicit_h > 0 {
        out.atoms[anchor].explicit_h -= 1;
    }
    let offset = out.atoms.len();
    for fa in &frag.atoms {
        let mut na = fa.clone();
        na.coords_2d = None;
        na.coords_3d = None;
        out.add_atom(na);
    }
    for fb in &frag.bonds {
        out.add_bond(offset + fb.u, offset + fb.v, fb.order, fb.stereo)?;
    }
    out.add_bond(anchor, offset, BondOrder::Single, BondStereo::None)?;
    out.perceive()?;
    Ok(out)
}
