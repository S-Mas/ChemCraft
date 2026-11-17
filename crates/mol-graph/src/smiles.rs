use std::collections::HashMap;

use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, ChiralTag, MolAtom, MolGraph};

/// Default safety limit on atoms in a single SMILES input.
pub const DEFAULT_MAX_ATOMS: usize = 2_000;

/// Parse a SMILES string (supporting optional whitespace-separated name when using `parse_smiles_line`).
pub fn parse_smiles(input: &str) -> Result<MolGraph, MolError> {
    parse_smiles_with_limit(input, DEFAULT_MAX_ATOMS)
}

/// Parse a SMILES string with an explicit maximum atom limit.
pub fn parse_smiles_with_limit(input: &str, max_atoms: usize) -> Result<MolGraph, MolError> {
    let smi = input.trim();
    if smi.is_empty() {
        return Err(MolError::EmptyInput);
    }
    let bytes = smi.as_bytes();
    let mut i = 0usize;
    let mut mol = MolGraph::new();
    mol.input_smiles = Some(smi.to_string());

    let mut stack: Vec<usize> = Vec::new();
    let mut current: Option<usize> = None;
    let mut pending_bond: Option<BondOrder> = None;
    let mut pending_stereo = BondStereo::None;
    let mut rings: HashMap<u16, (usize, Option<BondOrder>, BondStereo, usize)> = HashMap::new();
    #[derive(Clone, Copy)]
    enum ChiralSlot {
        Nbr(usize),
        H,
        RingOpen(u16),
    }
    let mut nbr_slots: Vec<Vec<ChiralSlot>> = Vec::new();

    while i < bytes.len() {
        let ch = bytes[i];
        match ch {
            b'.' => {
                if pending_bond.is_some() {
                    return Err(MolError::InvalidSmiles {
                        pos: i,
                        msg: "dangling bond before '.'".to_string(),
                    });
                }
                current = None;
                pending_stereo = BondStereo::None;
                i += 1;
            }
            b'(' => {
                let Some(cur) = current else {
                    return Err(MolError::InvalidSmiles {
                        pos: i,
                        msg: "branch '(' without preceding atom".to_string(),
                    });
                };
                stack.push(cur);
                i += 1;
            }
            b')' => {
                current = Some(stack.pop().ok_or(MolError::InvalidSmiles {
                    pos: i,
                    msg: "unmatched ')'".to_string(),
                })?);
                i += 1;
            }
            b'-' => {
                pending_bond = Some(BondOrder::Single);
                i += 1;
            }
            b'=' => {
                pending_bond = Some(BondOrder::Double);
                i += 1;
            }
            b'#' => {
                pending_bond = Some(BondOrder::Triple);
                i += 1;
            }
            b':' => {
                pending_bond = Some(BondOrder::Aromatic);
                i += 1;
            }
            b'/' => {
                pending_stereo = BondStereo::CisTransUp;
                i += 1;
            }
            b'\\' => {
                pending_stereo = BondStereo::CisTransDown;
                i += 1;
            }
            b'%' | b'0'..=b'9' => {
                let (rnum, next_i) = if ch == b'%' {
                    if i + 2 >= bytes.len()
                        || !bytes[i + 1].is_ascii_digit()
                        || !bytes[i + 2].is_ascii_digit()
                    {
                        return Err(MolError::InvalidSmiles {
                            pos: i,
                            msg: "invalid '%' ring closure".to_string(),
                        });
                    }
                    let d = ((bytes[i + 1] - b'0') * 10 + (bytes[i + 2] - b'0')) as u16;
                    (d, i + 3)
                } else {
                    (u16::from(ch - b'0'), i + 1)
                };
                let Some(cur) = current else {
                    return Err(MolError::InvalidSmiles {
                        pos: i,
                        msg: "ring digit without preceding atom".to_string(),
                    });
                };
                if let Some((other, bo1, st1, _)) = rings.remove(&rnum) {
                    if other == cur {
                        return Err(MolError::InvalidSmiles {
                            pos: i,
                            msg: "self-loop ring closure".to_string(),
                        });
                    }
                    let order = resolve_bond_order(
                        bo1.or(pending_bond),
                        mol.atoms[other].aromatic,
                        mol.atoms[cur].aromatic,
                    );
                    let stereo = if st1 != BondStereo::None {
                        st1
                    } else {
                        pending_stereo
                    };
                    mol.add_bond(other, cur, order, stereo)?;
                    if let Some(slot) = nbr_slots[other]
                        .iter_mut()
                        .rev()
                        .find(|s| matches!(s, ChiralSlot::RingOpen(r) if *r == rnum))
                    {
                        *slot = ChiralSlot::Nbr(cur);
                    }
                    nbr_slots[cur].push(ChiralSlot::Nbr(other));
                } else {
                    rings.insert(rnum, (cur, pending_bond, pending_stereo, i));
                    nbr_slots[cur].push(ChiralSlot::RingOpen(rnum));
                }
                pending_bond = None;
                pending_stereo = BondStereo::None;
                i = next_i;
            }
            b'[' => {
                let (mut atom, next_i) = parse_bracket_atom(bytes, i, mol.atoms.len())?;
                let has_bracket_chiral_h = atom.chiral != ChiralTag::None && atom.explicit_h > 0;
                if has_bracket_chiral_h {
                    atom.chiral_h_slot = Some(if current.is_some() { 1 } else { 0 });
                }
                if mol.atoms.len() >= max_atoms {
                    return Err(MolError::TooManyAtoms {
                        count: mol.atoms.len() + 1,
                        limit: max_atoms,
                    });
                }
                let nid = mol.add_atom(atom);
                nbr_slots.push(Vec::new());
                if let Some(prev) = current {
                    let order = resolve_bond_order(
                        pending_bond,
                        mol.atoms[prev].aromatic,
                        mol.atoms[nid].aromatic,
                    );
                    mol.add_bond(prev, nid, order, pending_stereo)?;
                    nbr_slots[prev].push(ChiralSlot::Nbr(nid));
                    nbr_slots[nid].push(ChiralSlot::Nbr(prev));
                }
                if has_bracket_chiral_h {
                    nbr_slots[nid].push(ChiralSlot::H);
                }
                current = Some(nid);
                pending_bond = None;
                pending_stereo = BondStereo::None;
                i = next_i;
            }
            _ => {
                let (atom, next_i) = parse_organic_atom(bytes, i, mol.atoms.len())?;
                if mol.atoms.len() >= max_atoms {
                    return Err(MolError::TooManyAtoms {
                        count: mol.atoms.len() + 1,
                        limit: max_atoms,
                    });
                }
                let nid = mol.add_atom(atom);
                nbr_slots.push(Vec::new());
                if let Some(prev) = current {
                    let order = resolve_bond_order(
                        pending_bond,
                        mol.atoms[prev].aromatic,
                        mol.atoms[nid].aromatic,
                    );
                    mol.add_bond(prev, nid, order, pending_stereo)?;
                    nbr_slots[prev].push(ChiralSlot::Nbr(nid));
                    nbr_slots[nid].push(ChiralSlot::Nbr(prev));
                }
                current = Some(nid);
                pending_bond = None;
                pending_stereo = BondStereo::None;
                i = next_i;
            }
        }
    }

    if !stack.is_empty() {
        return Err(MolError::InvalidSmiles {
            pos: bytes.len(),
            msg: "unclosed branch '('".to_string(),
        });
    }
    if let Some((rnum, (_, _, _, pos))) = rings.into_iter().next() {
        return Err(MolError::InvalidSmiles {
            pos,
            msg: format!("unclosed ring closure {rnum}"),
        });
    }
    if mol.atoms.is_empty() {
        return Err(MolError::EmptyInput);
    }

    for (idx, slots) in nbr_slots.into_iter().enumerate() {
        if mol.atoms[idx].chiral != ChiralTag::None {
            mol.atoms[idx].chiral_seq = slots
                .into_iter()
                .filter_map(|s| match s {
                    ChiralSlot::Nbr(n) => Some(n),
                    ChiralSlot::H => Some(usize::MAX),
                    ChiralSlot::RingOpen(_) => None,
                })
                .collect();
        }
    }
    mol.suppress_explicit_hydrogens();
    mol.desalt();
    mol.perceive()?;
    Ok(mol)
}

/// Parse a batch line `SMILES [Optional Name]`.
pub fn parse_smiles_line(line: &str, default_name: &str) -> Result<MolGraph, MolError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Err(MolError::EmptyInput);
    }
    let mut parts = trimmed.splitn(2, char::is_whitespace);
    let smi = parts.next().unwrap_or("");
    let name = parts
        .next()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or(default_name);
    let mut mol = parse_smiles(smi)?;
    mol.name = Some(name.to_string());
    Ok(mol)
}

fn resolve_bond_order(explicit: Option<BondOrder>, a1: bool, a2: bool) -> BondOrder {
    match explicit {
        Some(bo) => bo,
        None => {
            if a1 && a2 {
                BondOrder::Aromatic
            } else {
                BondOrder::Single
            }
        }
    }
}

fn parse_organic_atom(bytes: &[u8], pos: usize, id: usize) -> Result<(MolAtom, usize), MolError> {
    let b = bytes[pos];
    let (sym, arom, len) = if b == b'C' && bytes.get(pos + 1) == Some(&b'l') {
        ("Cl", false, 2)
    } else if b == b'B' && bytes.get(pos + 1) == Some(&b'r') {
        ("Br", false, 2)
    } else {
        let s = match b {
            b'*' => ("*", false),
            b'B' => ("B", false),
            b'C' => ("C", false),
            b'N' => ("N", false),
            b'O' => ("O", false),
            b'P' => ("P", false),
            b'S' => ("S", false),
            b'F' => ("F", false),
            b'I' => ("I", false),
            b'b' => ("B", true),
            b'c' => ("C", true),
            b'n' => ("N", true),
            b'o' => ("O", true),
            b'p' => ("P", true),
            b's' => ("S", true),
            _ => {
                return Err(MolError::InvalidSmiles {
                    pos,
                    msg: format!("unexpected SMILES token '{}'", b as char),
                });
            }
        };
        (s.0, s.1, 1)
    };
    let mut atom = MolAtom::new(id, sym)?;
    atom.aromatic = arom;
    Ok((atom, pos + len))
}

fn parse_bracket_atom(bytes: &[u8], start: usize, id: usize) -> Result<(MolAtom, usize), MolError> {
    let mut i = start + 1;
    let mut isotope: Option<u16> = None;
    let iso_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i > iso_start {
        let txt =
            std::str::from_utf8(&bytes[iso_start..i]).map_err(|_| MolError::InvalidSmiles {
                pos: iso_start,
                msg: "invalid UTF-8 isotope".to_string(),
            })?;
        isotope = txt.parse::<u16>().ok();
    }
    if i >= bytes.len() {
        return Err(MolError::InvalidSmiles {
            pos: start,
            msg: "unclosed '['".to_string(),
        });
    }

    let first = bytes[i];
    if first == b'*' {
        i += 1;
        let mut atom = MolAtom::new(id, "*")?;
        atom.from_bracket = true;
        if bytes.get(i) == Some(&b':') {
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
        }
        if bytes.get(i) != Some(&b']') {
            return Err(MolError::InvalidSmiles {
                pos: i,
                msg: "expected ']' closing '*' bracket atom".to_string(),
            });
        }
        return Ok((atom, i + 1));
    }
    if !first.is_ascii_alphabetic() {
        return Err(MolError::InvalidSmiles {
            pos: i,
            msg: "missing element symbol inside '[]'".to_string(),
        });
    }
    let arom = first.is_ascii_lowercase();
    let sym = if let Some(&second) = bytes.get(i + 1)
        && second.is_ascii_lowercase()
        && crate::element::element_by_symbol(&format!(
            "{}{}",
            first.to_ascii_uppercase() as char,
            second as char
        ))
        .is_some()
    {
        i += 2;
        format!("{}{}", first.to_ascii_uppercase() as char, second as char)
    } else {
        i += 1;
        (first.to_ascii_uppercase() as char).to_string()
    };

    let mut atom = MolAtom::new(id, &sym)?;
    atom.aromatic = arom;
    atom.isotope = isotope;
    atom.from_bracket = true;

    if bytes.get(i) == Some(&b'@') {
        if bytes.get(i + 1) == Some(&b'@') {
            atom.chiral = ChiralTag::Clockwise;
            i += 2;
        } else {
            atom.chiral = ChiralTag::CounterClockwise;
            i += 1;
        }
    }

    if bytes.get(i) == Some(&b'H') {
        i += 1;
        let d_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i > d_start {
            let txt = std::str::from_utf8(&bytes[d_start..i]).unwrap_or("1");
            atom.explicit_h = txt.parse::<u8>().unwrap_or(1);
        } else {
            atom.explicit_h = 1;
        }
    }

    if let Some(&sign) = bytes.get(i)
        && (sign == b'+' || sign == b'-')
    {
        let mult: i16 = if sign == b'+' { 1 } else { -1 };
        i += 1;
        let d_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let charge_val: i16 = if i > d_start {
            let txt = std::str::from_utf8(&bytes[d_start..i]).unwrap_or("1");
            let mag = txt.parse::<i16>().map_err(|_| MolError::InvalidSmiles {
                pos: d_start,
                msg: "formal charge out of range".to_string(),
            })?;
            mult.checked_mul(mag).ok_or(MolError::InvalidSmiles {
                pos: d_start,
                msg: "formal charge overflow".to_string(),
            })?
        } else {
            let mut count: i16 = 1;
            while bytes.get(i) == Some(&sign) {
                count = count.checked_add(1).ok_or(MolError::InvalidSmiles {
                    pos: i,
                    msg: "formal charge count overflow".to_string(),
                })?;
                if count > 15 {
                    return Err(MolError::InvalidSmiles {
                        pos: i,
                        msg: "formal charge magnitude exceeds [-15, +15]".to_string(),
                    });
                }
                i += 1;
            }
            mult * count
        };
        if !(-15..=15).contains(&charge_val) {
            return Err(MolError::InvalidSmiles {
                pos: start,
                msg: format!("formal charge {charge_val} outside [-15, +15]"),
            });
        }
        atom.formal_charge = charge_val as i8;
    }

    // Optional atom map `:N`
    if bytes.get(i) == Some(&b':') {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
    }

    if bytes.get(i) != Some(&b']') {
        return Err(MolError::InvalidSmiles {
            pos: i,
            msg: "expected ']' closing bracket atom".to_string(),
        });
    }
    Ok((atom, i + 1))
}
