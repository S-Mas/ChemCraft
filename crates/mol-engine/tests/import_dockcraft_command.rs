//! The import as a `Command`, which is what the GUI button, the palette and any other caller run.
//!
//! Going through the command registry rather than mutating the table from the UI is what gives
//! the import undo for free, puts it in the command palette, and keeps the GUI and the headless
//! paths doing the same thing.

use mol_engine::{ChemSession, Command, CommandRegistry, default_dockcraft_results_dir};
use mol_graph::parse_smiles;

fn sidecar(dir: &std::path::Path, rows: &str) -> std::io::Result<std::path::PathBuf> {
    let path = dir.join("series.dockcraft.json");
    std::fs::write(
        &path,
        format!(
            r#"{{"format":"dockcraft-results","version":"1",
                "run":{{"receptor_name":"3PTB","scoring":"vina","exhaustiveness":8,
                       "box_center":[0.0,0.0,0.0],"box_size":[20.0,20.0,20.0],
                       "dockcraft_version":"0.1.0","elapsed_ms":10}},
                "rows":[{rows}]}}"#
        )
        .as_bytes(),
    )?;
    Ok(path)
}

fn session_with_series() -> Result<ChemSession, Box<dyn std::error::Error>> {
    let mut session = ChemSession::new();
    for (name, smiles, pic50) in [
        ("Benzamidine", "NC(=N)c1ccccc1", 4.74),
        ("2-Naphthamidine", "NC(=N)c1ccc2ccccc2c1", 5.62),
    ] {
        let mut mol = parse_smiles(smiles)?;
        mol.name = Some(name.to_string());
        let mut table = (*session.sar_table).clone();
        table.add_molecule(mol, Some(pic50), Some("trypsin_S1".to_string()));
        session.sar_table = std::sync::Arc::new(table);
    }
    Ok(session)
}

#[test]
fn importing_results_adds_columns_and_can_be_undone() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("chemcraft_import_cmd_{}", std::process::id()));
    std::fs::create_dir_all(&dir)?;
    let path = sidecar(
        &dir,
        r#"{"canonical_smiles":"","source_row_id":0,"name":"Benzamidine","rank":2,
            "affinity_kcal":-5.672,"hyde_delta_g_kcal":-5.208,"estimated_ki_nm":151663.4,
            "ligand_efficiency":0.579,"mpo_score":0.537,"measured_pic50":4.74,
            "pose_sdf":"poses/2_Benzamidine.sdf"},
           {"canonical_smiles":"","source_row_id":1,"name":"2-Naphthamidine","rank":1,
            "affinity_kcal":-6.519,"hyde_delta_g_kcal":-5.770,"estimated_ki_nm":63451.1,
            "ligand_efficiency":0.438,"mpo_score":0.613,"measured_pic50":5.62,
            "pose_sdf":"poses/1_2-Naphthamidine.sdf"}"#,
    )?;

    let mut session = session_with_series()?;
    let outcome = CommandRegistry::execute(
        &mut session,
        Command::ImportDockCraftResults {
            path: Some(path.display().to_string()),
        },
    )?;
    assert!(
        outcome.message.contains("2 of 2"),
        "unexpected message: {}",
        outcome.message
    );
    assert!(
        outcome.message.contains("3PTB") && outcome.message.contains("vina"),
        "the message should say what the compounds were docked against: {}",
        outcome.message
    );

    let benz = session
        .sar_table
        .rows
        .iter()
        .find(|r| r.name == "Benzamidine")
        .ok_or("missing row")?;
    assert_eq!(
        benz.properties
            .get("DockCraft_affinity_kcal")
            .map(String::as_str),
        Some("-5.672")
    );
    assert_eq!(benz.pic50, Some(4.74), "the measured column was disturbed");

    // Undo has to restore the table, which is the reason this goes through the registry at all.
    CommandRegistry::execute(&mut session, Command::Undo)?;
    let after_undo = session
        .sar_table
        .rows
        .iter()
        .find(|r| r.name == "Benzamidine")
        .ok_or("missing row")?;
    assert!(
        !after_undo
            .properties
            .contains_key("DockCraft_affinity_kcal"),
        "undo did not remove the imported columns"
    );
    assert_eq!(after_undo.pic50, Some(4.74));

    let _ = std::fs::remove_dir_all(&dir);
    Ok(())
}

#[test]
fn importing_into_an_empty_table_says_so_instead_of_doing_nothing()
-> Result<(), Box<dyn std::error::Error>> {
    let mut session = ChemSession::new();
    let result =
        CommandRegistry::execute(&mut session, Command::ImportDockCraftResults { path: None });
    let err = match result {
        Err(e) => e,
        Ok(outcome) => {
            return Err(format!(
                "importing into an empty SAR table should fail loudly, got: {}",
                outcome.message
            )
            .into());
        }
    };
    assert!(err.to_string().contains("empty"), "unhelpful error: {err}");
    Ok(())
}

#[test]
fn a_missing_sidecar_names_the_directory_it_looked_in() -> Result<(), Box<dyn std::error::Error>> {
    let mut session = session_with_series()?;
    let result = CommandRegistry::execute(
        &mut session,
        Command::ImportDockCraftResults {
            path: Some("/nonexistent/nope.dockcraft.json".to_string()),
        },
    );
    let err = match result {
        Err(e) => e,
        Ok(outcome) => {
            return Err(format!("a missing file should fail, got: {}", outcome.message).into());
        }
    };
    assert!(err.to_string().contains("nope.dockcraft.json"), "{err}");

    // And the default directory is the mirror of the handoff one, so the pair reads symmetrically.
    let dir = default_dockcraft_results_dir().display().to_string();
    assert!(
        dir.contains("ChemCraft") && dir.ends_with("dockcraft_results"),
        "unexpected default results directory {dir}"
    );
    Ok(())
}
