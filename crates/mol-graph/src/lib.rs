//! Layer 0 (`mol-graph`): Clean-room editable molecular graph, iterative Hückel aromaticity/SSSR,
//! CIP R/S & E/Z stereochemistry, SMILES/SMARTS, SDF V2000/V3000, DataWarrior `.dwar`, CSV, and mutations.

#![forbid(unsafe_code)]

pub mod dwar;
pub mod element;
pub mod error;
pub mod mol;
pub mod mutate;
pub mod sdf;
pub mod smarts;
pub mod smiles;

pub use dwar::{
    TabularMolRecord, parse_csv_records, parse_csv_records_with_limit, parse_dwar,
    parse_dwar_with_limit, write_dwar,
};
pub use element::{ElementMeta, element_by_number, element_by_symbol, max_neutral_valence};
pub use error::MolError;
pub use mol::{BondCisTrans, BondOrder, BondStereo, ChiralTag, MolAtom, MolBond, MolGraph};
pub use mutate::{
    SubstituentKind, attach_fragment_at, mutate_molecule, replace_carboxyl_with_tetrazole,
};
pub use sdf::{parse_mol_block, parse_sdf, parse_sdf_with_limit, write_sdf_v2000, write_sdf_v3000};
pub use smarts::SmartsQuery;
pub use smiles::{DEFAULT_MAX_ATOMS, parse_smiles, parse_smiles_line, parse_smiles_with_limit};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_aspirin_steroid_and_bioisostere_tetrazole() -> Result<(), MolError> {
        let aspirin = parse_smiles("CC(=O)Oc1ccccc1C(=O)O")?;
        assert_eq!(aspirin.formula(), "C9H8O4");
        assert_eq!(aspirin.aromatic_atom_count(), 6);

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
    fn kekule_equivalence_fused_valence_pyrrole_acetate_ferrocene_and_macrocycle()
    -> Result<(), MolError> {
        // 1. Naphthalene & Caffeine must have 0 valence errors and Kekule == Aromatic structural_key
        let naph_arom = parse_smiles("c1ccc2ccccc2c1")?;
        let naph_kek = parse_smiles("C1=CC=C2C=CC=CC2=C1")?;
        assert!(naph_arom.warnings.is_empty(), "{:?}", naph_arom.warnings);
        assert!(naph_kek.warnings.is_empty(), "{:?}", naph_kek.warnings);
        assert_eq!(naph_arom.structural_key(), naph_kek.structural_key());

        let caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")?;
        assert!(caffeine.warnings.is_empty(), "{:?}", caffeine.warnings);

        // 2. Pyrrole NH preserved in both aromatic and Kekule forms (C4H5N)
        let pyr_arom = parse_smiles("c1cc[nH]c1")?;
        let pyr_kek = parse_smiles("C1=CC=CN1")?;
        assert_eq!(pyr_arom.formula(), "C4H5N");
        assert_eq!(pyr_kek.formula(), "C4H5N");
        assert_eq!(pyr_arom.structural_key(), pyr_kek.structural_key());

        // 3. Carboxylate read back from SDF keeps exact 59.04 Da (no phantom H on O-)
        let acetate = parse_smiles("CC(=O)[O-]")?;
        assert!((acetate.molecular_weight() - 59.044).abs() < 0.05);
        let sdf = write_sdf_v2000(&[&acetate]);
        let back = parse_sdf(&sdf)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert!((back.molecular_weight() - 59.044).abs() < 0.05);
        assert_eq!(back.formula(), "C2H3O2");

        // 4. Ferrocene keeps Iron (Fe)
        let ferrocene = parse_smiles("[cH-]1cccc1.[cH-]1cccc1.[Fe+2]")?;
        assert_eq!(ferrocene.formula(), "C10H10Fe");

        // 5. Macrocycle >16 atoms detected by SSSR
        let c20 = parse_smiles("C1CCCCCCCCCCCCCCCCCCC1")?;
        assert_eq!(c20.rings.len(), 1);
        assert_eq!(c20.rings[0].len(), 20);
        Ok(())
    }

    #[test]
    fn cip_rs_and_ez_stereochemistry_and_canonical_distinction() -> Result<(), MolError> {
        let l_ala = parse_smiles("N[C@@H](C)C(=O)O")?;
        let l_ala_alt = parse_smiles("C[C@H](N)C(=O)O")?;
        let d_ala = parse_smiles("N[C@H](C)C(=O)O")?;

        assert_eq!(l_ala.atoms[1].cip_label, Some('S'));
        assert_eq!(l_ala_alt.atoms[1].cip_label, Some('S'));
        assert_eq!(d_ala.atoms[1].cip_label, Some('R'));
        assert_ne!(l_ala.structural_key(), d_ala.structural_key());

        let e_alkene = parse_smiles("F/C=C/Cl")?;
        let z_alkene = parse_smiles("F/C=C\\Cl")?;
        let db_e = e_alkene
            .bonds
            .iter()
            .find(|b| b.order == BondOrder::Double)
            .ok_or(MolError::EmptyInput)?;
        let db_z = z_alkene
            .bonds
            .iter()
            .find(|b| b.order == BondOrder::Double)
            .ok_or(MolError::EmptyInput)?;
        assert_eq!(db_e.cis_trans, BondCisTrans::E);
        assert_eq!(db_z.cis_trans, BondCisTrans::Z);
        assert_ne!(e_alkene.structural_key(), z_alkene.structural_key());
        Ok(())
    }

    #[test]
    fn smarts_rejects_min_charge_and_honors_precedence_and_r2() -> Result<(), MolError> {
        // Must return Err cleanly instead of panicking on -128
        assert!(matches!(
            SmartsQuery::parse("[-128]"),
            Err(MolError::InvalidSmarts { .. })
        ));

        // [C,N;H1] = (C OR N) AND H1 -> matches aniline NH2? no, matches methanol OH? no, matches dimethylamine [NH] (H1)!
        let q = SmartsQuery::parse("[C,N;H1]")?;
        let nh = parse_smiles("CNC")?;
        let pyridine = parse_smiles("c1ccncc1")?;
        assert!(q.matches(&nh));
        assert!(!q.matches(&pyridine));

        // [R2] matches naphthalene bridgehead carbons (2 SSSR rings) and NOT benzene (1 ring)
        let r2 = SmartsQuery::parse("[R2]")?;
        let naph = parse_smiles("c1ccc2ccccc2c1")?;
        let benz = parse_smiles("c1ccccc1")?;
        assert!(r2.matches(&naph));
        assert!(!r2.matches(&benz));
        Ok(())
    }

    #[test]
    fn sdf_v2000_v3000_dwar_and_headerless_rejection() -> Result<(), MolError> {
        let mut caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")?;
        caffeine.name = Some("Caffeine".to_string());
        caffeine
            .properties
            .insert("pIC50".to_string(), "6.42".to_string());

        let v2000 = write_sdf_v2000(&[&caffeine]);
        let m2 = parse_sdf(&v2000)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(m2.formula(), "C8H10N4O2");
        assert_eq!(m2.properties.get("pIC50").map(String::as_str), Some("6.42"));

        let v3000 = write_sdf_v3000(&caffeine);
        let m3 = parse_sdf(&v3000)
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(m3.formula(), "C8H10N4O2");

        let dwar = write_dwar(&[(&caffeine, &caffeine.properties)]);
        let rec = parse_dwar(&dwar)?
            .into_iter()
            .next()
            .ok_or(MolError::EmptyInput)??;
        assert_eq!(rec.name, "Caffeine");

        // Corrupt V2000 >999 columns rejected cleanly
        let corrupt = "bad\n  t\n\n1200  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        assert!(matches!(
            parse_sdf(corrupt).first(),
            Some(Err(MolError::InvalidSdf { .. }))
        ));

        // Headerless CSV without SMILES header rejected cleanly
        assert!(matches!(
            parse_csv_records("ID,Value\n1,42\n"),
            Err(MolError::InvalidCsv(_))
        ));
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
