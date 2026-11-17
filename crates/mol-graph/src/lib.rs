//! Layer 0 (`mol-graph`): Clean-room editable molecular graph, aromaticity/SSSR/stereo perception,
//! SMILES, SMARTS substructure matching, SDF V2000/V3000, DataWarrior `.dwar`, CSV, and bioisostere mutations.

#![forbid(unsafe_code)]

pub mod dwar;
pub mod element;
pub mod error;
pub mod mol;
pub mod mutate;
pub mod sdf;
pub mod smarts;
pub mod smiles;

pub use dwar::{TabularMolRecord, parse_csv_records, parse_dwar, write_dwar};
pub use element::{ElementMeta, element_by_number, element_by_symbol, max_neutral_valence};
pub use error::MolError;
pub use mol::{BondCisTrans, BondOrder, BondStereo, ChiralTag, MolAtom, MolBond, MolGraph};
pub use mutate::{
    SubstituentKind, attach_fragment_at, mutate_molecule, replace_carboxyl_with_tetrazole,
};
pub use sdf::{parse_mol_block, parse_sdf, write_sdf_v2000, write_sdf_v3000};
pub use smarts::SmartsQuery;
pub use smiles::{DEFAULT_MAX_ATOMS, parse_smiles, parse_smiles_line, parse_smiles_with_limit};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_aspirin_steroid_and_bioisostere_tetrazole() -> Result<(), MolError> {
        let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        assert_eq!(aspirin.formula(), "C9H8O4");
        assert!(aspirin.aromatic_atom_count() == 6);

        let tetrazole_analogue = replace_carboxyl_with_tetrazole(&aspirin)?;
        assert_eq!(tetrazole_analogue.formula(), "C9H8N4O2");
        assert_ne!(
            tetrazole_analogue.structural_key(),
            aspirin.structural_key()
        );

        let ibuprofen = parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")?;
        let ibu_tetrazole = mutate_molecule(&ibuprofen, None, "tetrazole")?;
        let smi_q = SmartsQuery::parse("c1nnn[nH]1")?;
        assert!(smi_q.matches(&ibu_tetrazole));
        Ok(())
    }

    #[test]
    fn sdf_v2000_v3000_and_dwar_roundtrip() -> Result<(), MolError> {
        let mut caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")?;
        caffeine.name = Some("Caffeine".to_string());
        caffeine
            .properties
            .insert("pIC50".to_string(), "6.42".to_string());

        let v2000 = write_sdf_v2000(&[&caffeine]);
        let parsed_v2 = parse_sdf(&v2000);
        assert_eq!(parsed_v2.len(), 1);
        let m2 = parsed_v2.into_iter().next().ok_or(MolError::EmptyInput)??;
        assert_eq!(m2.formula(), "C8H10N4O2");
        assert_eq!(m2.properties.get("pIC50").map(String::as_str), Some("6.42"));

        let v3000 = write_sdf_v3000(&caffeine);
        let parsed_v3 = parse_sdf(&v3000);
        assert_eq!(parsed_v3.len(), 1);
        let m3 = parsed_v3.into_iter().next().ok_or(MolError::EmptyInput)??;
        assert_eq!(m3.formula(), "C8H10N4O2");

        let dwar = write_dwar(&[(&caffeine, &caffeine.properties)]);
        let parsed_dw = parse_dwar(&dwar)?;
        assert_eq!(parsed_dw.len(), 1);
        let rec = parsed_dw.into_iter().next().ok_or(MolError::EmptyInput)??;
        assert_eq!(rec.name, "Caffeine");
        assert_eq!(rec.mol.formula(), "C8H10N4O2");
        Ok(())
    }

    #[test]
    fn valence_warning_flags_pentavalent_carbon() -> Result<(), MolError> {
        let mut mol = MolGraph::new();
        let c = mol.add_atom(MolAtom::new(0, "C")?);
        for _ in 0..5 {
            let f = mol.add_atom(MolAtom::new(0, "F")?);
            mol.add_bond(c, f, BondOrder::Single, BondStereo::None)?;
        }
        mol.perceive()?;
        assert!(mol.atoms[c].valence_exceeded);
        assert!(mol.warnings.iter().any(|w| w.contains("Valence error")));
        Ok(())
    }
}
