//! Layer 3 (`mol-engine`): Session state, Copy-on-Write Undo/Redo, `CommandRegistry`,
//! atomic `.chemcraft` ZIP+JSON project persistence, and hardened `PathSandbox`.

#![forbid(unsafe_code)]

pub mod commands;
pub mod project;
pub mod sandbox;

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
}
