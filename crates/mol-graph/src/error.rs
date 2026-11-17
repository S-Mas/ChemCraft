use thiserror::Error;

/// Errors returned by `mol-graph` parsing, perception, and editing operations.
#[derive(Debug, Error)]
pub enum MolError {
    #[error("empty chemical input")]
    EmptyInput,

    #[error("invalid SMILES at byte {pos}: {msg}")]
    InvalidSmiles { pos: usize, msg: String },

    #[error("invalid SMARTS at byte {pos}: {msg}")]
    InvalidSmarts { pos: usize, msg: String },

    #[error("invalid SDF/MOL at line {line}: {msg}")]
    InvalidSdf { line: usize, msg: String },

    #[error("invalid DataWarrior .dwar file at line {line}: {msg}")]
    InvalidDwar { line: usize, msg: String },

    #[error("invalid CSV input: {0}")]
    InvalidCsv(String),

    #[error("molecule exceeds atom safety limit ({count} > {limit})")]
    TooManyAtoms { count: usize, limit: usize },

    #[error("canonical SMILES exceeded 99 simultaneous open ring closures")]
    TooManyOpenRings,

    #[error("molecular graph error: {0}")]
    Graph(String),

    #[error("mutation error: {0}")]
    Mutation(String),
}
