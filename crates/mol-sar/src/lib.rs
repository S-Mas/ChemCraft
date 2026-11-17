//! Layer 2 (`mol-sar`): High-performance in-memory DataWarrior-style SAR engine:
//! ECFP4 Morgan fingerprints, parallel Tanimoto & SMARTS filters, Matched Molecular Pairs (MMP),
//! Activity Cliffs (SALI), Ligand Efficiency metrics, and AdmetCraft bridge enrichment.

#![forbid(unsafe_code)]

pub mod fp;
pub mod mmp;
pub mod table;

pub use fp::MorganFp;
pub use mmp::{
    ActivityCliff, MatchedMolecularPair, find_activity_cliffs, find_matched_molecular_pairs,
};
pub use table::{SarFilter, SarRow, SarTable};

#[cfg(test)]
mod tests {
    use super::*;
    use mol_admet_bridge::AdmetCraftProvider;
    use mol_graph::parse_smiles;

    #[test]
    fn ecfp4_tanimoto_substructure_mmp_and_cliffs_work_on_50_compounds() {
        let mut table = SarTable::demo_50_compounds();
        assert_eq!(table.rows.len(), 50);

        table.enrich_with_admet(&AdmetCraftProvider);
        assert!(table.rows.iter().all(|r| r.admet.is_some()));

        let filter = SarFilter {
            smarts_substructure: Some("c1ccccc1".to_string()),
            ..SarFilter::default()
        };
        let aromatic_hits = table.matching_indices(&filter);
        assert!(aromatic_hits.is_ok());
        let hits = aromatic_hits.unwrap_or_default();
        assert!(hits.len() >= 45);

        let cliffs = find_activity_cliffs(&table.rows, 0.35);
        assert!(!cliffs.is_empty());

        let mmps = find_matched_molecular_pairs(&table.rows);
        assert!(!mmps.is_empty());

        if let (Ok(m1), Ok(m2)) = (
            parse_smiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
            parse_smiles("CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1"),
        ) {
            let sim = MorganFp::ecfp4(&m1).tanimoto(&MorganFp::ecfp4(&m2));
            assert!(
                sim > 0.25 && sim < 1.0,
                "expected moderate ECFP4 similarity for bioisostere pair, got {sim}"
            );
            assert!((MorganFp::ecfp4(&m1).tanimoto(&MorganFp::ecfp4(&m1)) - 1.0).abs() < 1e-12);
        }
    }
}
