use serde::{Deserialize, Serialize};

use mol_graph::{BondOrder, MolGraph};

/// 2048-bit folded Morgan / ECFP4 circular fingerprint (`radius = 2`).
#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct MorganFp {
    pub bits: [u64; 32],
}

impl MorganFp {
    /// Compute ECFP4 (radius 2) circular fingerprint for `mol`.
    pub fn ecfp4(mol: &MolGraph) -> Self {
        Self::morgan(mol, 2)
    }

    /// Compute folded Morgan circular fingerprint up to `radius`.
    pub fn morgan(mol: &MolGraph, radius: usize) -> Self {
        let mut fp = Self::default();
        let n = mol.atoms.len();
        if n == 0 {
            return fp;
        }

        // Round 0 atom invariants (atomic number, heavy degree, total H, charge, aromatic, ring)
        let mut inv: Vec<u64> = (0..n)
            .map(|i| {
                let a = &mol.atoms[i];
                let deg = mol.adj[i].len() as u64;
                let cip = match a.cip_label {
                    Some('R') => 1u64,
                    Some('S') => 2u64,
                    _ => 0u64,
                };
                let mut h = mix64(0x9e37_79b9_7f4a_7c15 ^ u64::from(a.atomic_number));
                h = mix64(h ^ (deg << 8));
                h = mix64(h ^ (u64::from(a.total_h()) << 16));
                h = mix64(h ^ ((a.formal_charge as i64 as u64) << 24));
                h = mix64(
                    h ^ (u64::from(a.aromatic) << 32) ^ (u64::from(a.in_ring) << 33) ^ (cip << 36),
                );
                fp.set_hash(h);
                h
            })
            .collect();

        for _r in 1..=radius {
            let mut next = vec![0u64; n];
            for i in 0..n {
                let mut nbr_tags: Vec<(u64, u64)> = mol.adj[i]
                    .iter()
                    .map(|&(nbr, bid)| {
                        let bo = match mol.bonds[bid].order {
                            BondOrder::Single => 1u64,
                            BondOrder::Double => 2,
                            BondOrder::Triple => 3,
                            BondOrder::Aromatic => 4,
                            BondOrder::Dative => 5,
                        };
                        let ez = match mol.bonds[bid].cis_trans {
                            mol_graph::BondCisTrans::None => 0u64,
                            mol_graph::BondCisTrans::E => 1,
                            mol_graph::BondCisTrans::Z => 2,
                        };
                        ((bo << 2) | ez, inv[nbr])
                    })
                    .collect();
                nbr_tags.sort_unstable();
                let mut acc = mix64(inv[i].wrapping_mul(0xbf58_476d_1ce4_e5b9));
                for (bo, ninv) in nbr_tags {
                    acc = mix64(acc ^ bo.wrapping_mul(0x94d0_49bb_1331_11eb) ^ ninv);
                }
                fp.set_hash(acc);
                next[i] = acc;
            }
            inv = next;
        }
        fp
    }

    #[inline]
    fn set_hash(&mut self, h: u64) {
        let bit = (h as usize) & 2047;
        self.bits[bit >> 6] |= 1u64 << (bit & 63);
    }

    pub fn bit_count(&self) -> u32 {
        self.bits.iter().map(|w| w.count_ones()).sum()
    }

    /// Tanimoto similarity in `[0.0, 1.0]`.
    pub fn tanimoto(&self, other: &Self) -> f64 {
        let mut inter = 0u32;
        let mut union = 0u32;
        for (a, b) in self.bits.iter().zip(other.bits.iter()) {
            inter += (a & b).count_ones();
            union += (a | b).count_ones();
        }
        if union == 0 {
            0.0
        } else {
            f64::from(inter) / f64::from(union)
        }
    }
}

#[inline]
fn mix64(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}
