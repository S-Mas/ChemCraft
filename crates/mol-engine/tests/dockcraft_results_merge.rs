//! Folding DockCraft's results back into the SAR table.
//!
//! This is the return leg of the ChemCraft -> DockCraft -> ChemCraft loop. DockCraft writes a
//! versioned `*.dockcraft.json` sidecar rather than a `.chemcraft` project, because producing one
//! would mean fabricating `MolGraph` and `MorganFp` on the far side and the two schemas would
//! drift apart. The merge belongs here, where the document is owned.
//!
//! What the tests pin: the join really joins, nothing measured is overwritten, and a sidecar that
//! does not match is reported rather than silently half-applied.

use std::path::Path;

use mol_engine::{DockCraftResults, MergeOutcome, merge_dockcraft_results, read_dockcraft_results};
use mol_graph::parse_smiles;
use mol_sar::SarTable;

fn series() -> Result<SarTable, Box<dyn std::error::Error>> {
    let mut table = SarTable::new();
    for (name, smiles, pic50) in [
        ("Benzamidine", "NC(=N)c1ccccc1", 4.74),
        ("2-Naphthamidine", "NC(=N)c1ccc2ccccc2c1", 5.62),
        ("Tryptamine", "NCCc1c[nH]c2ccccc12", 3.10),
    ] {
        let mut mol = parse_smiles(smiles)?;
        mol.name = Some(name.to_string());
        table.add_molecule(mol, Some(pic50), Some("trypsin_S1".to_string()));
    }
    Ok(table)
}

fn sidecar_json(rows: &str) -> String {
    format!(
        r#"{{"format":"dockcraft-results","version":"1",
            "run":{{"receptor_name":"3PTB","scoring":"vina","exhaustiveness":8,
                   "box_center":[0.0,0.0,0.0],"box_size":[20.0,20.0,20.0],
                   "dockcraft_version":"0.1.0","elapsed_ms":1234}},
            "rows":[{rows}]}}"#
    )
}

#[test]
fn results_join_by_row_id_and_land_beside_the_measured_column()
-> Result<(), Box<dyn std::error::Error>> {
    let mut table = series()?;
    let rows = r#"
        {"canonical_smiles":"","source_row_id":0,"name":"Benzamidine","rank":3,
         "affinity_kcal":-5.672,"hyde_delta_g_kcal":-5.208,"estimated_ki_nm":151663.4,
         "ligand_efficiency":0.579,"mpo_score":0.537,"measured_pic50":4.74,
         "pose_sdf":"poses/3_Benzamidine.sdf"},
        {"canonical_smiles":"","source_row_id":1,"name":"2-Naphthamidine","rank":1,
         "affinity_kcal":-6.519,"hyde_delta_g_kcal":-5.770,"estimated_ki_nm":63451.1,
         "ligand_efficiency":0.438,"mpo_score":0.613,"measured_pic50":5.62,
         "pose_sdf":"poses/1_2-Naphthamidine.sdf"}"#;
    let results: DockCraftResults = serde_json::from_str(&sidecar_json(rows))?;

    let outcome = merge_dockcraft_results(&mut table, &results, Path::new("/tmp/results"));
    assert_eq!(
        outcome,
        MergeOutcome {
            merged: 2,
            unmatched: Vec::new()
        }
    );

    let benz = table
        .rows
        .iter()
        .find(|r| r.name == "Benzamidine")
        .ok_or("missing benzamidine")?;
    assert_eq!(
        benz.properties
            .get("DockCraft_affinity_kcal")
            .map(String::as_str),
        Some("-5.672")
    );
    assert_eq!(
        benz.properties.get("DockCraft_rank").map(String::as_str),
        Some("3")
    );
    assert_eq!(
        benz.properties
            .get("DockCraft_receptor")
            .map(String::as_str),
        Some("3PTB")
    );
    // Pose paths arrive relative to the sidecar and are stored resolved, so a row points at a
    // file rather than at a fragment.
    assert_eq!(
        benz.properties
            .get("DockCraft_pose_sdf")
            .map(String::as_str),
        Some("/tmp/results/poses/3_Benzamidine.sdf")
    );

    // The measured column is untouched. A computed score and an experimental one must stay
    // distinguishable, which is the entire reason the keys are prefixed.
    assert_eq!(benz.pic50, Some(4.74));
    assert_eq!(benz.series, "trypsin_S1");
    assert!(
        !benz.properties.contains_key("pIC50"),
        "the merge invented a measured column"
    );

    // A row with no results keeps its measured data and gains nothing.
    let tryp = table
        .rows
        .iter()
        .find(|r| r.name == "Tryptamine")
        .ok_or("missing tryptamine")?;
    assert_eq!(tryp.pic50, Some(3.10));
    assert!(!tryp.properties.contains_key("DockCraft_affinity_kcal"));
    Ok(())
}

#[test]
fn results_join_by_canonical_smiles_when_there_is_no_row_id()
-> Result<(), Box<dyn std::error::Error>> {
    let mut table = series()?;
    // Use whatever this program calls the structure, which is the only spelling the far side is
    // given to join on.
    let naph = table
        .rows
        .iter()
        .find(|r| r.name == "2-Naphthamidine")
        .ok_or("missing")?
        .canonical_smiles
        .clone();
    let rows = format!(
        r#"{{"canonical_smiles":"{naph}","source_row_id":null,"name":"whatever",
             "rank":1,"affinity_kcal":-6.5,"hyde_delta_g_kcal":-5.8,"estimated_ki_nm":1.0,
             "ligand_efficiency":0.44,"mpo_score":0.61,"measured_pic50":null,"pose_sdf":""}}"#
    );
    let results: DockCraftResults = serde_json::from_str(&sidecar_json(&rows))?;

    let outcome = merge_dockcraft_results(&mut table, &results, Path::new("/tmp/r"));
    assert_eq!(outcome.merged, 1);
    assert!(outcome.unmatched.is_empty());
    let naph_row = table
        .rows
        .iter()
        .find(|r| r.name == "2-Naphthamidine")
        .ok_or("missing")?;
    assert_eq!(
        naph_row
            .properties
            .get("DockCraft_affinity_kcal")
            .map(String::as_str),
        Some("-6.500")
    );
    // Joining does not rename the row: the SAR table's own naming wins.
    assert_eq!(naph_row.name, "2-Naphthamidine");
    Ok(())
}

#[test]
fn a_row_that_matches_nothing_is_reported_not_swallowed() -> Result<(), Box<dyn std::error::Error>>
{
    let mut table = series()?;
    let rows = r#"
        {"canonical_smiles":"CCO","source_row_id":99,"name":"Ethanol","rank":1,
         "affinity_kcal":-2.0,"hyde_delta_g_kcal":-1.0,"estimated_ki_nm":1.0,
         "ligand_efficiency":0.1,"mpo_score":0.1,"measured_pic50":null,"pose_sdf":""}"#;
    let results: DockCraftResults = serde_json::from_str(&sidecar_json(rows))?;

    let outcome = merge_dockcraft_results(&mut table, &results, Path::new("/tmp/r"));
    assert_eq!(outcome.merged, 0);
    assert_eq!(outcome.unmatched, vec!["Ethanol".to_string()]);
    assert!(
        table
            .rows
            .iter()
            .all(|r| !r.properties.contains_key("DockCraft_affinity_kcal")),
        "an unmatched row was applied somewhere"
    );
    Ok(())
}

#[test]
fn a_file_of_the_wrong_format_is_refused() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("chemcraft_sidecar_{}", std::process::id()));
    std::fs::create_dir_all(&dir)?;

    let wrong = dir.join("wrong.json");
    std::fs::write(
        &wrong,
        br#"{"format":"something-else","version":"1","rows":[]}"#,
    )?;
    assert!(
        read_dockcraft_results(&wrong).is_err(),
        "a sidecar of another format was accepted"
    );

    let good = dir.join("good.dockcraft.json");
    std::fs::write(&good, sidecar_json("").as_bytes())?;
    let parsed = read_dockcraft_results(&good)?;
    assert_eq!(parsed.run.receptor_name, "3PTB");
    assert!(parsed.rows.is_empty());

    // A field DockCraft adds later must not stop an older ChemCraft reading the rest.
    let future = dir.join("future.dockcraft.json");
    std::fs::write(
        &future,
        br#"{"format":"dockcraft-results","version":"2","run":{"receptor_name":"X"},
             "rows":[{"name":"A","affinity_kcal":-1.5,"something_new":{"a":1}}],
             "another_new_block":[1,2,3]}"#,
    )?;
    let forward = read_dockcraft_results(&future)?;
    assert_eq!(forward.rows.len(), 1);
    assert!((forward.rows[0].affinity_kcal + 1.5).abs() < 1e-9);

    let _ = std::fs::remove_dir_all(&dir);
    Ok(())
}
