//! Layer 3 (`mol-engine`): Session state, Copy-on-Write Undo/Redo, `CommandRegistry`,
//! atomic `.chemcraft` ZIP+JSON project persistence, and hardened `PathSandbox`.

#![forbid(unsafe_code)]

pub mod chimerax;
pub mod commands;
pub mod project;
pub mod sandbox;

pub use chimerax::{
    ChimeraXBundle, ChimeraXError, default_chimerax_export_dir, is_chimerax_available,
    launch_chimerax_for_mol, write_chimerax_bundle,
};
pub use commands::{
    ChemSession, Command, CommandError, CommandInfo, CommandOutcome, CommandRegistry,
};
pub use project::{
    ChemCraftDocument, ProjectError, decode_chemcraft_zip, encode_chemcraft_zip,
    save_chemcraft_atomic,
};
pub use sandbox::{MAX_READ_BYTES, PathSandbox, SandboxError};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn undo_redo_and_chemcraft_zip_roundtrip_work() -> Result<(), Box<dyn std::error::Error>> {
        let mut session = ChemSession::new();
        CommandRegistry::execute(&mut session, Command::LoadDemo50Sar)?;
        assert_eq!(session.sar_table.rows.len(), 50);

        // Mutate carboxylic acid -> tetrazole on Ibuprofen and check ADMET delta
        CommandRegistry::execute(
            &mut session,
            Command::MutateBioisostere {
                atom_index: None,
                group: "tetrazole".to_string(),
            },
        )?;
        assert!(session.last_mutation_delta.is_some());
        assert!(session.can_undo());

        CommandRegistry::execute(&mut session, Command::Undo)?;
        assert!(session.can_redo());

        let doc = ChemCraftDocument {
            format: "chemcraft".to_string(),
            version: "0.1.0".to_string(),
            canvas_mol: (*session.canvas_mol).clone(),
            reaction_arrows: (*session.reaction_arrows).clone(),
            sar_table: (*session.sar_table).clone(),
            sar_filter: session.sar_filter.clone(),
        };
        let bytes = encode_chemcraft_zip(&doc)?;
        let decoded = decode_chemcraft_zip(&bytes)?;
        assert_eq!(decoded.sar_table.rows.len(), 50);
        assert_eq!(
            decoded.canvas_mol.canonical_smiles(),
            session.canvas_mol.canonical_smiles()
        );
        Ok(())
    }

    #[test]
    fn sandbox_blocks_symlinks_and_path_traversal() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_sandbox_test");
        let _ = std::fs::remove_dir_all(&base);
        let r = base.join("read");
        let w = base.join("write");
        let outside = base.join("outside.txt");
        std::fs::create_dir_all(&r)?;
        std::fs::create_dir_all(&w)?;
        std::fs::write(&outside, "secret")?;
        #[cfg(unix)]
        std::os::unix::fs::symlink(&outside, r.join("link.smi"))?;

        let sb = PathSandbox::new(Some(&r), Some(&w))?;
        assert!(
            sb.read_to_string_limited(std::path::Path::new("../outside.txt"), 1024)
                .is_err()
        );
        #[cfg(unix)]
        assert!(
            sb.read_to_string_limited(std::path::Path::new("link.smi"), 1024)
                .is_err()
        );
        assert!(sb.open_write(std::path::Path::new("../pwned.svg")).is_err());

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn zip_decompression_bomb_and_4000_openat2_dir_races_are_blocked()
    -> Result<(), Box<dyn std::error::Error>> {
        use flate2::{Compression, write::DeflateEncoder};
        use std::io::Write;

        // 1. Decompression bomb test: ~40 MB zeros compressed to tiny archive declaring uncomp=42 bytes
        let huge = vec![b' '; 40 * 1024 * 1024];
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::fast());
        enc.write_all(&huge)?;
        let comp = enc.finish()?;
        let mut zip = Vec::new();
        zip.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        zip.extend_from_slice(&20u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(&8u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(&0u32.to_le_bytes());
        zip.extend_from_slice(&(comp.len() as u32).to_le_bytes());
        zip.extend_from_slice(&42u32.to_le_bytes()); // spoofed 42-byte header!
        zip.extend_from_slice(&13u16.to_le_bytes());
        zip.extend_from_slice(&0u16.to_le_bytes());
        zip.extend_from_slice(b"manifest.json");
        zip.extend_from_slice(&comp);
        assert!(matches!(
            decode_chemcraft_zip(&zip),
            Err(ProjectError::InvalidArchive(_))
        ));

        // 2. 4,000 intermediate-directory symlink swap race test
        #[cfg(unix)]
        {
            use std::os::unix::fs::symlink;
            use std::sync::{
                Arc,
                atomic::{AtomicBool, Ordering},
            };
            let base = std::env::temp_dir().join("chemcraft_sb_race_4k");
            let _ = std::fs::remove_dir_all(&base);
            let r = base.join("read");
            let w = base.join("write");
            let outside = base.join("outside");
            std::fs::create_dir_all(r.join("sub"))?;
            std::fs::create_dir_all(w.join("sub"))?;
            std::fs::create_dir_all(&outside)?;
            std::fs::write(outside.join("leak.smi"), "SECRET_OUTSIDE")?;
            std::fs::write(r.join("sub/leak.smi"), "SAFE_INSIDE")?;

            let sb = PathSandbox::new(Some(&r), Some(&w))?;
            let stop = Arc::new(AtomicBool::new(false));
            let stop_r = Arc::clone(&stop);
            let r_link = r.join("sub");
            let w_link = w.join("sub");

            let racer = std::thread::spawn(move || {
                let mut t = false;
                while !stop_r.load(Ordering::Relaxed) {
                    t = !t;
                    let _ = std::fs::remove_file(&r_link);
                    let _ = std::fs::remove_dir_all(&r_link);
                    let _ = std::fs::remove_file(&w_link);
                    let _ = std::fs::remove_dir_all(&w_link);
                    if t {
                        let _ = symlink(&outside, &r_link);
                        let _ = symlink(&outside, &w_link);
                    } else {
                        let _ = std::fs::create_dir(&r_link);
                        let _ = std::fs::write(r_link.join("leak.smi"), "SAFE_INSIDE");
                        let _ = std::fs::create_dir(&w_link);
                    }
                }
            });

            for _ in 0..4_000 {
                if let Ok((_, txt)) =
                    sb.read_to_string_limited(std::path::Path::new("sub/leak.smi"), 1024)
                {
                    assert_ne!(txt, "SECRET_OUTSIDE", "leaked outside read-root!");
                }
                if let Ok((_, mut f)) = sb.open_write(std::path::Path::new("sub/pwn.svg")) {
                    let _ = f.write_all(b"x");
                }
            }
            stop.store(true, Ordering::Relaxed);
            let _ = racer.join();
            assert!(
                !base.join("outside/pwn.svg").exists(),
                "wrote outside write-root!"
            );
            let _ = std::fs::remove_dir_all(&base);
        }
        Ok(())
    }

    #[test]
    fn peptide_command_and_chimerax_bundle_export_work() -> Result<(), Box<dyn std::error::Error>> {
        let mut session = ChemSession::new();
        CommandRegistry::execute(
            &mut session,
            Command::BuildPeptide {
                sequence: "Ac-Ala-Gly-Tyr-Ser-NH2".to_string(),
                secondary_structure: mol_graph::SecondaryStructure::AlphaHelix,
            },
        )?;
        assert!(!session.canvas_mol.atoms.is_empty());

        let out_dir = std::env::temp_dir().join("chemcraft_chimerax_unit_test");
        let _ = std::fs::remove_dir_all(&out_dir);
        let bundle = write_chimerax_bundle(&session.canvas_mol, &out_dir, Some("agys_helix"))?;
        assert!(bundle.is_peptide);
        assert!(bundle.script_path.exists());
        assert!(bundle.structure_path.exists());
        let cxc = std::fs::read_to_string(&bundle.script_path)?;
        assert!(cxc.contains("open agys_helix.pdb"));
        assert!(cxc.contains("cartoon"));
        assert!(cxc.contains("hbonds reveal true"));
        let _ = std::fs::remove_dir_all(&out_dir);
        Ok(())
    }
}
