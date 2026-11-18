//! Clean-room 2D + 3D peptide builder for canonical L-amino acids.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::coords3d::{covalent_radius, generate_3d_coords, place_nerf_atom};
use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, ChiralTag, MolAtom, MolGraph, PdbAtomMeta};

/// Maximum residues permitted per peptide build for interactive/MCP safety.
pub const MAX_PEPTIDE_RESIDUES: usize = 50;

/// Secondary structure backbone torsion preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum SecondaryStructure {
    #[default]
    BetaStrand,
    AlphaHelix,
    Extended,
}

impl SecondaryStructure {
    pub fn label(self) -> &'static str {
        match self {
            Self::BetaStrand => "Beta-Strand",
            Self::AlphaHelix => "Alpha-Helix",
            Self::Extended => "Extended",
        }
    }

    /// Standard Ramachandran `(phi, psi, omega)` torsion angles in degrees.
    pub fn phi_psi_omega_deg(self) -> (f64, f64, f64) {
        match self {
            Self::BetaStrand => (-139.0, 135.0, 180.0),
            Self::AlphaHelix => (-57.0, -47.0, 180.0),
            Self::Extended => (-180.0, 180.0, 180.0),
        }
    }

    pub fn from_str_lossy(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "beta" | "beta-strand" | "beta_strand" | "strand" | "sheet" => Some(Self::BetaStrand),
            "alpha" | "alpha-helix" | "alpha_helix" | "helix" => Some(Self::AlphaHelix),
            "extended" | "ext" | "linear" => Some(Self::Extended),
            _ => None,
        }
    }
}

/// Sidechain internal atom descriptor (`sc_idx = 0` denotes backbone `CA`;
/// `sc_idx = 1..` denotes each sidechain heavy atom in build order).
struct ScAtom {
    pdb_name: &'static str,
    symbol: &'static str,
    parent_sc_idx: usize,
    order: BondOrder,
    aromatic: bool,
    /// Relative 2D offset from `CA` in standard bond-length units `(dx_right, dy_outward)`.
    pos_2d: [f64; 2],
    /// NeRF build parameters relative to `(gp_ref, parent)`:
    /// `ref_a` (`None` = use backbone `N`), `ref_b` (`None` = use `CA`), `angle_deg`, `dihedral_deg`.
    nerf_angle_deg: f64,
    nerf_dihedral_deg: f64,
    nerf_gp_idx: Option<usize>,
}

#[allow(clippy::too_many_arguments)]
const fn sc(
    pdb_name: &'static str,
    symbol: &'static str,
    parent_sc_idx: usize,
    order: BondOrder,
    aromatic: bool,
    pos_2d: [f64; 2],
    nerf_angle_deg: f64,
    nerf_dihedral_deg: f64,
    nerf_gp_idx: Option<usize>,
) -> ScAtom {
    ScAtom {
        pdb_name,
        symbol,
        parent_sc_idx,
        order,
        aromatic,
        pos_2d,
        nerf_angle_deg,
        nerf_dihedral_deg,
        nerf_gp_idx,
    }
}

struct ResidueSpec {
    code1: char,
    code3: &'static str,
    name: &'static str,
    charge_ph74: f64,
    sidechain: &'static [ScAtom],
    /// Extra bond(s) inside sidechain `(sc_idx_a, sc_idx_b, order)` for rings.
    ring_closures: &'static [(usize, usize, BondOrder)],
    /// For Proline: sidechain `sc_idx` bonded back to backbone `N`.
    close_to_backbone_n: Option<usize>,
    /// Optional explicit sidechain stereocenter `(sc_idx, ChiralTag, [a, b, c, d])` using `usize::MAX` for H.
    extra_chiral: Option<(usize, ChiralTag, [usize; 4])>,
}

const AMINO_ACIDS: &[ResidueSpec] = &[
    ResidueSpec {
        code1: 'G',
        code3: "GLY",
        name: "Glycine",
        charge_ph74: 0.0,
        sidechain: &[],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'A',
        code3: "ALA",
        name: "Alanine",
        charge_ph74: 0.0,
        sidechain: &[sc(
            "CB",
            "C",
            0,
            BondOrder::Single,
            false,
            [0.0, 1.0],
            110.5,
            -122.0,
            None,
        )],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'V',
        code3: "VAL",
        name: "Valine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG1",
                "C",
                1,
                BondOrder::Single,
                false,
                [-0.87, 1.5],
                111.0,
                180.0,
                None,
            ),
            sc(
                "CG2",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                111.0,
                -60.0,
                None,
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'L',
        code3: "LEU",
        name: "Leucine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                116.0,
                -65.0,
                None,
            ),
            sc(
                "CD1",
                "C",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                111.0,
                175.0,
                Some(0),
            ),
            sc(
                "CD2",
                "C",
                2,
                BondOrder::Single,
                false,
                [1.73, 1.0],
                111.0,
                55.0,
                Some(0),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'I',
        code3: "ILE",
        name: "Isoleucine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG1",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                111.0,
                -60.0,
                None,
            ),
            sc(
                "CG2",
                "C",
                1,
                BondOrder::Single,
                false,
                [-0.87, 1.5],
                111.0,
                60.0,
                None,
            ),
            sc(
                "CD1",
                "C",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                113.0,
                170.0,
                Some(0),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        // L-Ile CB is (3S): neighbors in seq: [CA(0), H(MAX), CG2(3), CG1(2)] -> CounterClockwise?
        // Using open-smiles orientation around CA(0): [0, MAX, 3, 2]
        extra_chiral: Some((1, ChiralTag::CounterClockwise, [0, usize::MAX, 3, 2])),
    },
    ResidueSpec {
        code1: 'P',
        code3: "PRO",
        name: "Proline",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.25, 0.95],
                103.5,
                -118.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [-0.45, 1.58],
                105.0,
                28.0,
                None,
            ),
            sc(
                "CD",
                "C",
                2,
                BondOrder::Single,
                false,
                [-1.10, 0.90],
                105.0,
                -28.0,
                Some(0),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: Some(3),
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'F',
        code3: "PHE",
        name: "Phenylalanine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                true,
                [0.0, 2.0],
                113.8,
                -65.0,
                None,
            ),
            sc(
                "CD1",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [-0.87, 2.5],
                120.0,
                90.0,
                Some(0),
            ),
            sc(
                "CD2",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [0.87, 2.5],
                120.0,
                -90.0,
                Some(0),
            ),
            sc(
                "CE1",
                "C",
                3,
                BondOrder::Aromatic,
                true,
                [-0.87, 3.5],
                120.0,
                180.0,
                Some(1),
            ),
            sc(
                "CE2",
                "C",
                4,
                BondOrder::Aromatic,
                true,
                [0.87, 3.5],
                120.0,
                180.0,
                Some(1),
            ),
            sc(
                "CZ",
                "C",
                5,
                BondOrder::Aromatic,
                true,
                [0.0, 4.0],
                120.0,
                0.0,
                Some(2),
            ),
        ],
        ring_closures: &[(6, 7, BondOrder::Aromatic)],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'Y',
        code3: "TYR",
        name: "Tyrosine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                true,
                [0.0, 2.0],
                113.8,
                -65.0,
                None,
            ),
            sc(
                "CD1",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [-0.87, 2.5],
                120.0,
                90.0,
                Some(0),
            ),
            sc(
                "CD2",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [0.87, 2.5],
                120.0,
                -90.0,
                Some(0),
            ),
            sc(
                "CE1",
                "C",
                3,
                BondOrder::Aromatic,
                true,
                [-0.87, 3.5],
                120.0,
                180.0,
                Some(1),
            ),
            sc(
                "CE2",
                "C",
                4,
                BondOrder::Aromatic,
                true,
                [0.87, 3.5],
                120.0,
                180.0,
                Some(1),
            ),
            sc(
                "CZ",
                "C",
                5,
                BondOrder::Aromatic,
                true,
                [0.0, 4.0],
                120.0,
                0.0,
                Some(2),
            ),
            sc(
                "OH",
                "O",
                7,
                BondOrder::Single,
                false,
                [0.0, 5.0],
                120.0,
                180.0,
                Some(5),
            ),
        ],
        ring_closures: &[(6, 7, BondOrder::Aromatic)],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'W',
        code3: "TRP",
        name: "Tryptophan",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                true,
                [0.0, 2.0],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "CD1",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [-0.88, 2.58],
                126.0,
                90.0,
                Some(0),
            ),
            sc(
                "CD2",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [0.88, 2.58],
                126.0,
                -90.0,
                Some(0),
            ),
            sc(
                "NE1",
                "N",
                3,
                BondOrder::Aromatic,
                true,
                [-0.54, 3.52],
                108.0,
                180.0,
                Some(1),
            ),
            sc(
                "CE2",
                "C",
                4,
                BondOrder::Aromatic,
                true,
                [0.54, 3.52],
                108.0,
                180.0,
                Some(1),
            ),
            sc(
                "CE3",
                "C",
                4,
                BondOrder::Aromatic,
                true,
                [1.80, 2.25],
                132.0,
                0.0,
                Some(1),
            ),
            sc(
                "CZ2",
                "C",
                6,
                BondOrder::Aromatic,
                true,
                [1.20, 4.32],
                120.0,
                180.0,
                Some(4),
            ),
            sc(
                "CZ3",
                "C",
                7,
                BondOrder::Aromatic,
                true,
                [2.42, 3.05],
                120.0,
                180.0,
                Some(4),
            ),
            sc(
                "CH2",
                "C",
                8,
                BondOrder::Aromatic,
                true,
                [2.12, 4.05],
                120.0,
                0.0,
                Some(6),
            ),
        ],
        ring_closures: &[(5, 6, BondOrder::Aromatic), (9, 10, BondOrder::Aromatic)],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'S',
        code3: "SER",
        name: "Serine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "OG",
                "O",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                111.0,
                65.0,
                None,
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'T',
        code3: "THR",
        name: "Threonine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "OG1",
                "O",
                1,
                BondOrder::Single,
                false,
                [-0.87, 1.5],
                109.5,
                60.0,
                None,
            ),
            sc(
                "CG2",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                111.0,
                -60.0,
                None,
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: Some((1, ChiralTag::Clockwise, [0, usize::MAX, 2, 3])),
    },
    ResidueSpec {
        code1: 'C',
        code3: "CYS",
        name: "Cysteine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "SG",
                "S",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                114.0,
                -65.0,
                None,
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'M',
        code3: "MET",
        name: "Methionine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "SD",
                "S",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                112.5,
                180.0,
                Some(0),
            ),
            sc(
                "CE",
                "C",
                3,
                BondOrder::Single,
                false,
                [1.73, 3.0],
                100.0,
                180.0,
                Some(1),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'N',
        code3: "ASN",
        name: "Asparagine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.0, 2.0],
                112.5,
                -65.0,
                None,
            ),
            sc(
                "OD1",
                "O",
                2,
                BondOrder::Double,
                false,
                [-0.87, 2.5],
                121.0,
                15.0,
                Some(0),
            ),
            sc(
                "ND2",
                "N",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                117.0,
                -165.0,
                Some(0),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'Q',
        code3: "GLN",
        name: "Glutamine",
        charge_ph74: 0.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "CD",
                "C",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                112.5,
                180.0,
                Some(0),
            ),
            sc(
                "OE1",
                "O",
                3,
                BondOrder::Double,
                false,
                [0.0, 3.0],
                121.0,
                15.0,
                Some(1),
            ),
            sc(
                "NE2",
                "N",
                3,
                BondOrder::Single,
                false,
                [1.73, 3.0],
                117.0,
                -165.0,
                Some(1),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'D',
        code3: "ASP",
        name: "Aspartate",
        charge_ph74: -1.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.0, 2.0],
                112.5,
                -65.0,
                None,
            ),
            sc(
                "OD1",
                "O",
                2,
                BondOrder::Double,
                false,
                [-0.87, 2.5],
                122.0,
                15.0,
                Some(0),
            ),
            sc(
                "OD2",
                "O",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                115.0,
                -165.0,
                Some(0),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'E',
        code3: "GLU",
        name: "Glutamate",
        charge_ph74: -1.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "CD",
                "C",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                112.5,
                180.0,
                Some(0),
            ),
            sc(
                "OE1",
                "O",
                3,
                BondOrder::Double,
                false,
                [0.0, 3.0],
                122.0,
                15.0,
                Some(1),
            ),
            sc(
                "OE2",
                "O",
                3,
                BondOrder::Single,
                false,
                [1.73, 3.0],
                115.0,
                -165.0,
                Some(1),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'K',
        code3: "LYS",
        name: "Lysine",
        charge_ph74: 1.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "CD",
                "C",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                112.5,
                180.0,
                Some(0),
            ),
            sc(
                "CE",
                "C",
                3,
                BondOrder::Single,
                false,
                [1.73, 3.0],
                112.5,
                180.0,
                Some(1),
            ),
            sc(
                "NZ",
                "N",
                4,
                BondOrder::Single,
                false,
                [1.73, 4.0],
                112.0,
                180.0,
                Some(2),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'R',
        code3: "ARG",
        name: "Arginine",
        charge_ph74: 1.0,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                false,
                [0.87, 1.5],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "CD",
                "C",
                2,
                BondOrder::Single,
                false,
                [0.87, 2.5],
                112.0,
                180.0,
                Some(0),
            ),
            sc(
                "NE",
                "N",
                3,
                BondOrder::Single,
                false,
                [1.73, 3.0],
                112.0,
                180.0,
                Some(1),
            ),
            sc(
                "CZ",
                "C",
                4,
                BondOrder::Single,
                false,
                [1.73, 4.0],
                124.0,
                180.0,
                Some(2),
            ),
            sc(
                "NH1",
                "N",
                5,
                BondOrder::Double,
                false,
                [0.87, 4.5],
                120.0,
                0.0,
                Some(3),
            ),
            sc(
                "NH2",
                "N",
                5,
                BondOrder::Single,
                false,
                [2.60, 4.5],
                120.0,
                180.0,
                Some(3),
            ),
        ],
        ring_closures: &[],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
    ResidueSpec {
        code1: 'H',
        code3: "HIS",
        name: "Histidine",
        charge_ph74: 0.1,
        sidechain: &[
            sc(
                "CB",
                "C",
                0,
                BondOrder::Single,
                false,
                [0.0, 1.0],
                110.5,
                -122.0,
                None,
            ),
            sc(
                "CG",
                "C",
                1,
                BondOrder::Single,
                true,
                [0.0, 2.0],
                114.0,
                -65.0,
                None,
            ),
            sc(
                "ND1",
                "N",
                2,
                BondOrder::Aromatic,
                true,
                [-0.88, 2.58],
                126.0,
                90.0,
                Some(0),
            ),
            sc(
                "CD2",
                "C",
                2,
                BondOrder::Aromatic,
                true,
                [0.88, 2.58],
                126.0,
                -90.0,
                Some(0),
            ),
            sc(
                "CE1",
                "C",
                3,
                BondOrder::Aromatic,
                true,
                [-0.54, 3.52],
                108.0,
                180.0,
                Some(1),
            ),
            sc(
                "NE2",
                "N",
                4,
                BondOrder::Aromatic,
                true,
                [0.54, 3.52],
                108.0,
                180.0,
                Some(1),
            ),
        ],
        ring_closures: &[(5, 6, BondOrder::Aromatic)],
        close_to_backbone_n: None,
        extra_chiral: None,
    },
];

fn lookup_1(c: char) -> Option<&'static ResidueSpec> {
    let cu = c.to_ascii_uppercase();
    AMINO_ACIDS.iter().find(|r| r.code1 == cu)
}

fn lookup_3(code: &str) -> Option<&'static ResidueSpec> {
    let cu = code.to_ascii_uppercase();
    AMINO_ACIDS.iter().find(|r| r.code3 == cu)
}

/// Parsed sequence containing terminal cap flags and validated residues.
#[derive(Debug, Clone)]
pub struct ParsedPeptide {
    pub residues_idx: Vec<usize>,
    pub n_acetyl: bool,
    pub c_amide: bool,
    pub one_letter: String,
    pub three_letter: String,
    pub est_charge_ph74: f64,
}

/// Parse a 1-letter (`AGYS`, `Ac-AGYS-NH2`) or dash/space-separated 3-letter (`Ala-Gly-Ser`) sequence.
pub fn parse_peptide_sequence(raw: &str) -> Result<ParsedPeptide, MolError> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(MolError::EmptyInput);
    }

    let mut s = trimmed.to_string();
    let mut n_acetyl = false;
    let mut c_amide = false;

    if let Some(rest) = s
        .strip_prefix("Ac-")
        .or_else(|| s.strip_prefix("AC-"))
        .or_else(|| s.strip_prefix("ac-"))
    {
        n_acetyl = true;
        s = rest.to_string();
    } else if let Some(rest) = s.strip_prefix("H-").or_else(|| s.strip_prefix("h-")) {
        s = rest.to_string();
    }

    if let Some(rest) = s
        .strip_suffix("-NH2")
        .or_else(|| s.strip_suffix("-nh2"))
        .or_else(|| s.strip_suffix("-NH₂"))
    {
        c_amide = true;
        s = rest.to_string();
    } else if let Some(rest) = s.strip_suffix("-OH").or_else(|| s.strip_suffix("-oh")) {
        s = rest.to_string();
    }

    let s = s.trim();
    if s.is_empty() {
        return Err(MolError::EmptyInput);
    }

    let mut residues_idx = Vec::new();
    let has_delims = s.contains('-') || s.contains(' ') || s.contains(',') || s.contains(';');
    if has_delims {
        for tok in s
            .split(['-', ' ', ',', ';'])
            .map(str::trim)
            .filter(|t| !t.is_empty())
        {
            if let Some(pos) = AMINO_ACIDS
                .iter()
                .position(|r| r.code3.eq_ignore_ascii_case(tok))
            {
                residues_idx.push(pos);
            } else if tok.len() == 1
                && let Some(ch) = tok.chars().next()
                && let Some(pos) = AMINO_ACIDS
                    .iter()
                    .position(|r| r.code1 == ch.to_ascii_uppercase())
            {
                residues_idx.push(pos);
            } else {
                return Err(MolError::Graph(format!(
                    "unknown amino acid token '{tok}' in peptide sequence"
                )));
            }
        }
    } else if s.len() == 3
        && let Some(pos) = AMINO_ACIDS
            .iter()
            .position(|r| r.code3.eq_ignore_ascii_case(s))
        && s.chars().skip(1).any(|c| c.is_ascii_lowercase())
    {
        // Single title-cased 3-letter residue like "Ala" or "Trp"
        residues_idx.push(pos);
    } else {
        for ch in s.chars() {
            if let Some(pos) = AMINO_ACIDS
                .iter()
                .position(|r| r.code1 == ch.to_ascii_uppercase())
            {
                residues_idx.push(pos);
            } else {
                return Err(MolError::Graph(format!(
                    "invalid 1-letter amino acid code '{ch}' (expected A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)"
                )));
            }
        }
    }

    if residues_idx.is_empty() {
        return Err(MolError::EmptyInput);
    }
    if residues_idx.len() > MAX_PEPTIDE_RESIDUES {
        return Err(MolError::Graph(format!(
            "peptide length {} exceeds interactive limit of {MAX_PEPTIDE_RESIDUES} residues",
            residues_idx.len()
        )));
    }

    let one_letter: String = residues_idx.iter().map(|&i| AMINO_ACIDS[i].code1).collect();
    let three_vec: Vec<&str> = residues_idx.iter().map(|&i| AMINO_ACIDS[i].code3).collect();
    let mut est_charge_ph74: f64 = residues_idx
        .iter()
        .map(|&i| AMINO_ACIDS[i].charge_ph74)
        .sum();
    if !n_acetyl {
        est_charge_ph74 += 0.9; // N-terminal NH3+
    }
    if !c_amide {
        est_charge_ph74 -= 1.0; // C-terminal COO-
    }

    let prefix = if n_acetyl { "Ac-" } else { "" };
    let suffix = if c_amide { "-NH2" } else { "" };
    let three_letter = format!("{prefix}{}{suffix}", three_vec.join("-"));

    Ok(ParsedPeptide {
        residues_idx,
        n_acetyl,
        c_amide,
        one_letter,
        three_letter,
        est_charge_ph74,
    })
}

/// Build a peptide `MolGraph` with clean 2D zig-zag depiction coordinates, authentic 3D
/// secondary-structure (`BetaStrand`, `AlphaHelix`, `Extended`) coordinates, and PDB atom metadata.
pub fn build_peptide(sequence: &str, ss: SecondaryStructure) -> Result<MolGraph, MolError> {
    let parsed = parse_peptide_sequence(sequence)?;
    let mut mol = MolGraph::new();
    let prefix = if parsed.n_acetyl { "Ac-" } else { "" };
    let suffix = if parsed.c_amide { "-NH2" } else { "" };
    let disp_seq = format!("{prefix}{}{suffix}", parsed.one_letter);
    mol.name = Some(format!("Peptide {disp_seq} ({})", ss.label()));
    mol.properties
        .insert("PEPTIDE_SEQUENCE".to_string(), parsed.three_letter.clone());
    mol.properties
        .insert("PEPTIDE_1LETTER".to_string(), disp_seq);
    mol.properties.insert(
        "PEPTIDE_RESIDUES".to_string(),
        parsed.residues_idx.len().to_string(),
    );
    mol.properties
        .insert("SECONDARY_STRUCTURE".to_string(), ss.label().to_string());
    mol.properties.insert(
        "EST_CHARGE_PH74".to_string(),
        format!("{:+.1}", parsed.est_charge_ph74),
    );

    let (phi_deg, psi_deg, omega_deg) = ss.phi_psi_omega_deg();
    let phi = phi_deg.to_radians();
    let psi = psi_deg.to_radians();
    let omega = omega_deg.to_radians();

    // Engh & Huber peptide backbone bond lengths (A) & angles (deg -> rad)
    let d_n_ca = 1.458;
    let d_ca_c = 1.525;
    let d_c_n = 1.329;
    let d_c_o = 1.231;
    let ang_c_n_ca = 121.7_f64.to_radians();
    let ang_n_ca_c = 111.2_f64.to_radians();
    let ang_ca_c_n = 116.2_f64.to_radians();
    let ang_ca_c_o = 120.8_f64.to_radians();

    // 2D horizontal zig-zag parameters (~1.5 A bond length)
    let bl2 = 1.50;
    let dx_step = bl2 * 0.866;
    let dy_step = bl2 * 0.5;
    let mut x2 = 0.0_f64;
    let mut backbone_2d_idx: usize = 0;

    // Seed NeRF history before first N: imaginary frame behind residue 1
    let mut prev_c_3d = [-1.329, 0.0, 0.0];
    let mut prev_ca_3d = [-2.25, 1.15, 0.0];
    let mut prev_n_3d = [-3.60, 0.80, 0.0];
    let mut prev_c_atom: Option<usize> = None;

    // Optional N-terminal Acetyl cap (ACE)
    if parsed.n_acetyl {
        let mut ch3 = MolAtom::new(0, "C")?;
        ch3.coords_2d = Some([-2.25, -0.5]);
        ch3.coords_3d = Some(prev_ca_3d);
        ch3.pdb_meta = Some(PdbAtomMeta {
            atom_name: "CH3".to_string(),
            res_name: "ACE".to_string(),
            res_seq: 0,
            chain: 'A',
            is_hetatm: false,
        });
        let ch3_id = mol.add_atom(ch3);

        let mut c_cap = MolAtom::new(0, "C")?;
        c_cap.coords_2d = Some([-1.05, 0.2]);
        c_cap.coords_3d = Some(prev_c_3d);
        c_cap.pdb_meta = Some(PdbAtomMeta {
            atom_name: "C".to_string(),
            res_name: "ACE".to_string(),
            res_seq: 0,
            chain: 'A',
            is_hetatm: false,
        });
        let c_id = mol.add_atom(c_cap);

        let o_cap_3d = place_nerf_atom(prev_n_3d, prev_ca_3d, prev_c_3d, d_c_o, ang_ca_c_o, 0.0);
        let mut o_cap = MolAtom::new(0, "O")?;
        o_cap.coords_2d = Some([-1.05, 1.55]);
        o_cap.coords_3d = Some(o_cap_3d);
        o_cap.pdb_meta = Some(PdbAtomMeta {
            atom_name: "O".to_string(),
            res_name: "ACE".to_string(),
            res_seq: 0,
            chain: 'A',
            is_hetatm: false,
        });
        let o_id = mol.add_atom(o_cap);

        mol.add_bond(ch3_id, c_id, BondOrder::Single, BondStereo::None)?;
        mol.add_bond(c_id, o_id, BondOrder::Double, BondStereo::None)?;
        prev_c_atom = Some(c_id);
    }

    let total_res = parsed.residues_idx.len();
    let mut last_c_id = 0usize;
    let mut last_c_3d = [0.0; 3];
    let mut last_ca_3d = [0.0; 3];
    let mut last_n_3d = [0.0; 3];
    let mut last_c_2d = [0.0; 2];

    for (res_pos, &aa_idx) in parsed.residues_idx.iter().enumerate() {
        let spec = &AMINO_ACIDS[aa_idx];
        let res_seq = (res_pos + 1) as u16;
        let side_sign = if res_pos % 2 == 0 { 1.0 } else { -1.0 };

        // 1. Backbone N
        let n_3d = if res_pos == 0 {
            [0.0, 0.0, 0.0]
        } else {
            place_nerf_atom(prev_n_3d, prev_ca_3d, prev_c_3d, d_c_n, ang_ca_c_n, psi)
        };
        let n_2d = [
            x2,
            if backbone_2d_idx.is_multiple_of(2) {
                -0.42
            } else {
                0.42
            },
        ];
        x2 += dx_step;
        backbone_2d_idx += 1;

        let mut a_n = MolAtom::new(0, "N")?;
        a_n.coords_2d = Some(n_2d);
        a_n.coords_3d = Some(n_3d);
        a_n.pdb_meta = Some(PdbAtomMeta {
            atom_name: "N".to_string(),
            res_name: spec.code3.to_string(),
            res_seq,
            chain: 'A',
            is_hetatm: false,
        });
        let n_id = mol.add_atom(a_n);
        if let Some(pc) = prev_c_atom {
            mol.add_bond(pc, n_id, BondOrder::Single, BondStereo::None)?;
        }

        // 2. Backbone CA
        let ca_3d = if res_pos == 0 {
            let rad = (PI - ang_c_n_ca).max(0.2);
            [d_n_ca * rad.cos(), d_n_ca * rad.sin(), 0.0]
        } else {
            place_nerf_atom(prev_ca_3d, prev_c_3d, n_3d, d_n_ca, ang_c_n_ca, omega)
        };
        let ca_2d = [x2, side_sign * dy_step * 0.45];
        x2 += dx_step;
        backbone_2d_idx += 1;

        let mut a_ca = MolAtom::new(0, "C")?;
        a_ca.coords_2d = Some(ca_2d);
        a_ca.coords_3d = Some(ca_3d);
        a_ca.pdb_meta = Some(PdbAtomMeta {
            atom_name: "CA".to_string(),
            res_name: spec.code3.to_string(),
            res_seq,
            chain: 'A',
            is_hetatm: false,
        });
        let ca_id = mol.add_atom(a_ca);
        mol.add_bond(n_id, ca_id, BondOrder::Single, BondStereo::None)?;

        // 3. Backbone C
        let c_3d = place_nerf_atom(prev_c_3d, n_3d, ca_3d, d_ca_c, ang_n_ca_c, phi);
        let c_2d = [x2, -side_sign * dy_step * 0.45];
        x2 += dx_step;
        backbone_2d_idx += 1;

        let mut a_c = MolAtom::new(0, "C")?;
        a_c.coords_2d = Some(c_2d);
        a_c.coords_3d = Some(c_3d);
        a_c.pdb_meta = Some(PdbAtomMeta {
            atom_name: "C".to_string(),
            res_name: spec.code3.to_string(),
            res_seq,
            chain: 'A',
            is_hetatm: false,
        });
        let c_id = mol.add_atom(a_c);

        // Sidechain atoms (`sc_ids[0] = ca_id`, `sc_ids[1..] = sidechain atoms`)
        let mut sc_ids = Vec::with_capacity(1 + spec.sidechain.len());
        let mut sc_3d = Vec::with_capacity(1 + spec.sidechain.len());
        sc_ids.push(ca_id);
        sc_3d.push(ca_3d);

        // If chiral non-Gly residue, record L-amino-acid stereochemistry BEFORE adding C bond?
        // Wait: bond insertion order decides `adj[ca_id]` unless `chiral_seq` is set!
        // We set `chiral_seq` explicitly so order in `adj` does not limit chirality.
        for (sc_pos, sc_spec) in spec.sidechain.iter().enumerate() {
            let parent_id = sc_ids[sc_spec.parent_sc_idx];
            let p_pos_3d = sc_3d[sc_spec.parent_sc_idx];

            let bl = covalent_radius(&mol.atoms[parent_id].symbol)
                + covalent_radius(sc_spec.symbol)
                    * match sc_spec.order {
                        BondOrder::Aromatic => 0.91,
                        BondOrder::Double => 0.86,
                        _ => 1.0,
                    };
            let bl = bl.clamp(1.22, 1.85);

            let pos_3d = if sc_pos == 0 {
                // Place CB off (C, N, CA) to guarantee L-amino acid CORN chirality!
                place_nerf_atom(
                    c_3d,
                    n_3d,
                    ca_3d,
                    bl,
                    sc_spec.nerf_angle_deg.to_radians(),
                    sc_spec.nerf_dihedral_deg.to_radians(),
                )
            } else {
                let gp_3d = match sc_spec.nerf_gp_idx {
                    Some(g) => sc_3d[g],
                    None => n_3d,
                };
                place_nerf_atom(
                    gp_3d,
                    if sc_spec.parent_sc_idx == 1 && sc_spec.nerf_gp_idx.is_none() {
                        ca_3d
                    } else {
                        let parent_of_parent =
                            spec.sidechain[sc_spec.parent_sc_idx - 1].parent_sc_idx;
                        if sc_spec.nerf_gp_idx.is_some() {
                            sc_3d[parent_of_parent]
                        } else {
                            ca_3d
                        }
                    },
                    p_pos_3d,
                    bl,
                    sc_spec.nerf_angle_deg.to_radians(),
                    sc_spec.nerf_dihedral_deg.to_radians(),
                )
            };

            let pos_2d = [
                ca_2d[0] + sc_spec.pos_2d[0] * bl2 * 0.92,
                ca_2d[1] + side_sign * sc_spec.pos_2d[1] * bl2 * 0.92,
            ];

            let mut sa = MolAtom::new(0, sc_spec.symbol)?;
            sa.aromatic = sc_spec.aromatic;
            sa.coords_2d = Some(pos_2d);
            sa.coords_3d = Some(pos_3d);
            sa.pdb_meta = Some(PdbAtomMeta {
                atom_name: sc_spec.pdb_name.to_string(),
                res_name: spec.code3.to_string(),
                res_seq,
                chain: 'A',
                is_hetatm: false,
            });
            let sid = mol.add_atom(sa);
            let stereo = if sc_pos == 0 {
                if side_sign > 0.0 {
                    BondStereo::Up
                } else {
                    BondStereo::Down
                }
            } else {
                BondStereo::None
            };
            mol.add_bond(parent_id, sid, sc_spec.order, stereo)?;
            sc_ids.push(sid);
            sc_3d.push(pos_3d);
        }

        for &(a_sc, b_sc, ord) in spec.ring_closures {
            mol.add_bond(sc_ids[a_sc], sc_ids[b_sc], ord, BondStereo::None)?;
        }
        if let Some(cd_sc) = spec.close_to_backbone_n {
            mol.add_bond(sc_ids[cd_sc], n_id, BondOrder::Single, BondStereo::None)?;
        }

        // Now connect CA -> C
        mol.add_bond(ca_id, c_id, BondOrder::Single, BondStereo::None)?;

        // Record authentic L-alpha tetrahedral stereochemistry when sidechain CB exists
        if let Some(&cb_id) = sc_ids.get(1) {
            mol.atoms[ca_id].chiral = ChiralTag::Clockwise;
            mol.atoms[ca_id].chiral_h_slot = Some(1);
            mol.atoms[ca_id].chiral_seq = vec![n_id, usize::MAX, cb_id, c_id];
        }
        if let Some((sc_idx, tag, seq_local)) = spec.extra_chiral {
            let target_atom = sc_ids[sc_idx];
            let mapped = seq_local.map(|k| {
                if k == usize::MAX {
                    usize::MAX
                } else {
                    sc_ids[k]
                }
            });
            mol.atoms[target_atom].chiral = tag;
            mol.atoms[target_atom].chiral_h_slot = Some(1);
            mol.atoms[target_atom].chiral_seq = mapped.to_vec();
        }

        // 4. Backbone carbonyl =O
        let o_dihedral = if res_pos + 1 < total_res {
            psi + PI
        } else {
            0.0
        };
        let o_3d = place_nerf_atom(n_3d, ca_3d, c_3d, d_c_o, ang_ca_c_o, o_dihedral);
        let o_2d = [c_2d[0], c_2d[1] - side_sign * bl2 * 0.90];
        let mut a_o = MolAtom::new(0, "O")?;
        a_o.coords_2d = Some(o_2d);
        a_o.coords_3d = Some(o_3d);
        a_o.pdb_meta = Some(PdbAtomMeta {
            atom_name: "O".to_string(),
            res_name: spec.code3.to_string(),
            res_seq,
            chain: 'A',
            is_hetatm: false,
        });
        let o_id = mol.add_atom(a_o);
        mol.add_bond(c_id, o_id, BondOrder::Double, BondStereo::None)?;

        prev_n_3d = n_3d;
        prev_ca_3d = ca_3d;
        prev_c_3d = c_3d;
        prev_c_atom = Some(c_id);
        last_c_id = c_id;
        last_c_3d = c_3d;
        last_ca_3d = ca_3d;
        last_n_3d = n_3d;
        last_c_2d = c_2d;
    }

    // Terminal cap: C-terminal amide -NH2 or carboxyl -OH (OXT)
    let term_3d = place_nerf_atom(last_n_3d, last_ca_3d, last_c_3d, 1.33, ang_ca_c_n, PI);
    let term_2d = [last_c_2d[0] + dx_step, -last_c_2d[1] * 0.6 + 0.2];
    if parsed.c_amide {
        let mut nh2 = MolAtom::new(0, "N")?;
        nh2.coords_2d = Some(term_2d);
        nh2.coords_3d = Some(term_3d);
        nh2.pdb_meta = Some(PdbAtomMeta {
            atom_name: "N".to_string(),
            res_name: "NH2".to_string(),
            res_seq: (total_res + 1) as u16,
            chain: 'A',
            is_hetatm: false,
        });
        let tid = mol.add_atom(nh2);
        mol.add_bond(last_c_id, tid, BondOrder::Single, BondStereo::None)?;
    } else {
        let last_spec = &AMINO_ACIDS[parsed.residues_idx[total_res - 1]];
        let mut oxt = MolAtom::new(0, "O")?;
        oxt.coords_2d = Some(term_2d);
        oxt.coords_3d = Some(term_3d);
        oxt.pdb_meta = Some(PdbAtomMeta {
            atom_name: "OXT".to_string(),
            res_name: last_spec.code3.to_string(),
            res_seq: total_res as u16,
            chain: 'A',
            is_hetatm: false,
        });
        let tid = mol.add_atom(oxt);
        mol.add_bond(last_c_id, tid, BondOrder::Single, BondStereo::None)?;
    }

    // Center 2D and 3D centroids cleanly
    center_2d_and_3d(&mut mol);
    // Run perception (SSSR, implicit H, CIP R/S) without clobbering 2D/3D coordinates!
    mol.perceive()?;

    // Verify 3D finite sanity; fallback to conformer relaxation only if any coordinate failed
    if mol
        .atoms
        .iter()
        .any(|a| a.coords_3d.is_none_or(|p| p.iter().any(|x| !x.is_finite())))
    {
        generate_3d_coords(&mut mol);
    }

    Ok(mol)
}

fn center_2d_and_3d(mol: &mut MolGraph) {
    if mol.atoms.is_empty() {
        return;
    }
    let mut cx2 = 0.0;
    let mut cy2 = 0.0;
    let mut c3 = [0.0_f64; 3];
    let n = mol.atoms.len() as f64;
    for a in &mol.atoms {
        if let Some([x, y]) = a.coords_2d {
            cx2 += x;
            cy2 += y;
        }
        if let Some([x, y, z]) = a.coords_3d {
            c3[0] += x;
            c3[1] += y;
            c3[2] += z;
        }
    }
    for a in &mut mol.atoms {
        if let Some([x, y]) = &mut a.coords_2d {
            *x -= cx2 / n;
            *y -= cy2 / n;
        }
        if let Some([x, y, z]) = &mut a.coords_3d {
            *x -= c3[0] / n;
            *y -= c3[1] / n;
            *z -= c3[2] / n;
        }
    }
}

/// Expose lookup metadata list for GUI peptide amino-acid palette buttons.
pub fn canonical_amino_acids() -> Vec<(char, &'static str, &'static str)> {
    AMINO_ACIDS
        .iter()
        .map(|r| (r.code1, r.code3, r.name))
        .collect()
}

#[allow(dead_code)]
fn _keep_lookups() -> bool {
    lookup_1('A').is_some() && lookup_3("ALA").is_some()
}
