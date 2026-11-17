use std::fmt::Write;

use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, MolAtom, MolGraph};
use crate::smiles::DEFAULT_MAX_ATOMS;

/// Parse multi-record SDF text (supports both CTAB **V2000** and **V3000** blocks plus SD tags).
pub fn parse_sdf(input: &str) -> Vec<Result<MolGraph, MolError>> {
    if input.trim().is_empty() {
        return vec![Err(MolError::EmptyInput)];
    }
    let normalized = input.replace("\r\n", "\n");
    normalized
        .split("$$$$")
        .filter_map(|rec| {
            let block = rec.strip_prefix('\n').unwrap_or(rec);
            if block.trim().is_empty() {
                None
            } else {
                Some(parse_mol_block(block))
            }
        })
        .collect()
}

/// Parse a single CTAB MOL block (V2000 or V3000) + optional `> <TAG>` SD properties.
pub fn parse_mol_block(block: &str) -> Result<MolGraph, MolError> {
    let lines: Vec<&str> = block.lines().collect();
    if lines.len() < 4 {
        return Err(MolError::InvalidSdf {
            line: 1,
            msg: "MOL block shorter than 4 header lines".to_string(),
        });
    }
    let title = lines[0].trim();
    let counts_line = lines[3];

    if counts_line.contains("V3000") {
        return parse_v3000_block(&lines, title);
    }
    parse_v2000_block(&lines, title)
}

fn parse_v2000_block(lines: &[&str], title: &str) -> Result<MolGraph, MolError> {
    let counts_line = lines[3];
    let (atom_count, bond_count) = if counts_line.len() >= 6 {
        let a = counts_line
            .get(0..3)
            .unwrap_or("")
            .trim()
            .parse::<usize>()
            .map_err(|_| MolError::InvalidSdf {
                line: 4,
                msg: "invalid V2000 atom count".to_string(),
            })?;
        let b = counts_line
            .get(3..6)
            .unwrap_or("")
            .trim()
            .parse::<usize>()
            .map_err(|_| MolError::InvalidSdf {
                line: 4,
                msg: "invalid V2000 bond count".to_string(),
            })?;
        (a, b)
    } else {
        return Err(MolError::InvalidSdf {
            line: 4,
            msg: "truncated V2000 counts line".to_string(),
        });
    };

    if atom_count > DEFAULT_MAX_ATOMS {
        return Err(MolError::InvalidSdf {
            line: 4,
            msg: format!("atom count {atom_count} exceeds limit {DEFAULT_MAX_ATOMS}"),
        });
    }
    if bond_count > DEFAULT_MAX_ATOMS * 4 {
        return Err(MolError::InvalidSdf {
            line: 4,
            msg: format!("bond count {bond_count} exceeds safety limit"),
        });
    }
    let min_lines = 4usize.saturating_add(atom_count).saturating_add(bond_count);
    if lines.len() < min_lines {
        return Err(MolError::InvalidSdf {
            line: lines.len(),
            msg: "MOL block truncated before atom/bond table ends".to_string(),
        });
    }

    let mut mol = MolGraph::new();
    if !title.is_empty() {
        mol.name = Some(title.to_string());
    }

    for idx in 0..atom_count {
        let ln = 4 + idx;
        let line = lines[ln];
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 4 {
            return Err(MolError::InvalidSdf {
                line: ln + 1,
                msg: "atom line has fewer than 4 fields".to_string(),
            });
        }
        let x: f64 = toks[0].parse().map_err(|_| MolError::InvalidSdf {
            line: ln + 1,
            msg: "invalid X coordinate".to_string(),
        })?;
        let y: f64 = toks[1].parse().map_err(|_| MolError::InvalidSdf {
            line: ln + 1,
            msg: "invalid Y coordinate".to_string(),
        })?;
        let sym = toks[3];
        let mut atom = MolAtom::new(idx, sym)?;
        if x.is_finite() && y.is_finite() {
            atom.coords_2d = Some([x, y]);
        }
        if let Some(chg_s) = toks.get(5)
            && let Ok(code) = chg_s.parse::<u8>()
        {
            atom.formal_charge = match code {
                1 => 3,
                2 => 2,
                3 => 1,
                5 => -1,
                6 => -2,
                7 => -3,
                _ => 0,
            };
        }
        mol.add_atom(atom);
    }

    for idx in 0..bond_count {
        let ln = 4 + atom_count + idx;
        let line = lines[ln];
        let (u_raw, v_raw, t_raw, s_raw) = if line.len() >= 9 {
            let s1 = line.get(0..3).unwrap_or("").trim();
            let s2 = line.get(3..6).unwrap_or("").trim();
            let s3 = line.get(6..9).unwrap_or("").trim();
            let s4 = line.get(9..12).unwrap_or("0").trim();
            if let (Ok(u), Ok(v), Ok(t)) =
                (s1.parse::<usize>(), s2.parse::<usize>(), s3.parse::<u8>())
            {
                (u, v, t, s4.parse::<u8>().unwrap_or(0))
            } else {
                parse_v2000_bond_tokens(line, ln + 1)?
            }
        } else {
            parse_v2000_bond_tokens(line, ln + 1)?
        };
        if u_raw == 0 || v_raw == 0 || u_raw > atom_count || v_raw > atom_count {
            return Err(MolError::InvalidSdf {
                line: ln + 1,
                msg: "bond atom index out of range".to_string(),
            });
        }
        let order = v2000_bond_order(t_raw);
        let stereo = match s_raw {
            1 => BondStereo::Up,
            6 => BondStereo::Down,
            4 => BondStereo::Wavy,
            _ => BondStereo::None,
        };
        mol.add_bond(u_raw - 1, v_raw - 1, order, stereo)?;
    }

    parse_props_and_sd_tags(&mut mol, lines, 4 + atom_count + bond_count, title);
    finish_sdf_mol(mol)
}

fn parse_v3000_block(lines: &[&str], title: &str) -> Result<MolGraph, MolError> {
    let mut mol = MolGraph::new();
    if !title.is_empty() {
        mol.name = Some(title.to_string());
    }
    let mut atom_id_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut in_atom = false;
    let mut in_bond = false;
    let mut end_idx = lines.len();

    for (ln, raw) in lines.iter().enumerate().skip(4) {
        let line = raw.trim();
        if line == "M  END" {
            end_idx = ln;
            break;
        }
        let Some(rest) = line.strip_prefix("M  V30 ") else {
            continue;
        };
        let rest = rest.trim();
        if rest == "BEGIN ATOM" {
            in_atom = true;
            in_bond = false;
            continue;
        } else if rest == "END ATOM" {
            in_atom = false;
            continue;
        } else if rest == "BEGIN BOND" {
            in_bond = true;
            in_atom = false;
            continue;
        } else if rest == "END BOND" {
            in_bond = false;
            continue;
        }

        if in_atom {
            let toks: Vec<&str> = rest.split_whitespace().collect();
            if toks.len() < 5 {
                return Err(MolError::InvalidSdf {
                    line: ln + 1,
                    msg: "malformed V3000 ATOM record".to_string(),
                });
            }
            let ext_id: usize = toks[0].parse().map_err(|_| MolError::InvalidSdf {
                line: ln + 1,
                msg: "invalid V3000 atom id".to_string(),
            })?;
            if mol.atoms.len() >= DEFAULT_MAX_ATOMS {
                return Err(MolError::TooManyAtoms {
                    count: mol.atoms.len() + 1,
                    limit: DEFAULT_MAX_ATOMS,
                });
            }
            let sym = toks[1];
            let x: f64 = toks[2].parse().unwrap_or(0.0);
            let y: f64 = toks[3].parse().unwrap_or(0.0);
            let mut atom = MolAtom::new(mol.atoms.len(), sym)?;
            if x.is_finite() && y.is_finite() {
                atom.coords_2d = Some([x, y]);
            }
            for kv in &toks[5..] {
                if let Some(chg_s) = kv.strip_prefix("CHG=")
                    && let Ok(c) = chg_s.parse::<i32>()
                {
                    atom.formal_charge = c.clamp(-15, 15) as i8;
                }
            }
            let aid = mol.add_atom(atom);
            atom_id_map.insert(ext_id, aid);
        } else if in_bond {
            let toks: Vec<&str> = rest.split_whitespace().collect();
            if toks.len() < 4 {
                return Err(MolError::InvalidSdf {
                    line: ln + 1,
                    msg: "malformed V3000 BOND record".to_string(),
                });
            }
            let bo_num: u8 = toks[1].parse().unwrap_or(1);
            let u_ext: usize = toks[2].parse().unwrap_or(0);
            let v_ext: usize = toks[3].parse().unwrap_or(0);
            let (Some(&u), Some(&v)) = (atom_id_map.get(&u_ext), atom_id_map.get(&v_ext)) else {
                return Err(MolError::InvalidSdf {
                    line: ln + 1,
                    msg: "V3000 bond references unknown atom id".to_string(),
                });
            };
            let mut stereo = BondStereo::None;
            for kv in &toks[4..] {
                if let Some(cfg) = kv.strip_prefix("CFG=") {
                    stereo = match cfg {
                        "1" => BondStereo::Up,
                        "3" => BondStereo::Down,
                        "2" => BondStereo::Wavy,
                        _ => BondStereo::None,
                    };
                }
            }
            mol.add_bond(u, v, v2000_bond_order(bo_num), stereo)?;
        }
    }

    parse_props_and_sd_tags(&mut mol, lines, end_idx, title);
    finish_sdf_mol(mol)
}

fn finish_sdf_mol(mut mol: MolGraph) -> Result<MolGraph, MolError> {
    if mol.atoms.is_empty() {
        return Err(MolError::EmptyInput);
    }
    mol.suppress_explicit_hydrogens();
    mol.desalt();
    mol.perceive()?;
    Ok(mol)
}

fn parse_props_and_sd_tags(mol: &mut MolGraph, lines: &[&str], start_idx: usize, title: &str) {
    let mut i = start_idx;
    while i < lines.len() {
        let line = lines[i];
        if let Some(rest) = line.strip_prefix("M  CHG") {
            let nums: Vec<i32> = rest
                .split_whitespace()
                .filter_map(|t| t.parse::<i32>().ok())
                .collect();
            if let Some(&count) = nums.first()
                && count > 0
            {
                let safe_count = usize::try_from(count).unwrap_or(0).min(16);
                for pair in nums[1..].as_chunks::<2>().0.iter().take(safe_count) {
                    if let Some(one_based) = pair[0].checked_sub(1)
                        && let Ok(a_idx) = usize::try_from(one_based)
                        && let Some(atom) = mol.atoms.get_mut(a_idx)
                    {
                        atom.formal_charge = pair[1].clamp(-15, 15) as i8;
                    }
                }
            }
        } else if line.trim_start().starts_with("> ")
            && let Some(tag_start) = line.find('<')
            && let Some(tag_end) = line[tag_start + 1..].find('>')
        {
            let tag = line[tag_start + 1..tag_start + 1 + tag_end].trim();
            if let Some(val_line) = lines.get(i + 1) {
                let val = val_line.trim();
                if !val.is_empty() {
                    if matches!(
                        tag.to_ascii_uppercase().as_str(),
                        "NAME" | "ID" | "COMMON_NAME"
                    ) && title.is_empty()
                        && mol.name.is_none()
                    {
                        mol.name = Some(val.to_string());
                    } else {
                        mol.properties.insert(tag.to_string(), val.to_string());
                    }
                }
                i += 1;
            }
        }
        i += 1;
    }
}

fn parse_v2000_bond_tokens(
    line: &str,
    line_num: usize,
) -> Result<(usize, usize, u8, u8), MolError> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.len() < 3 {
        return Err(MolError::InvalidSdf {
            line: line_num,
            msg: "bond line has fewer than 3 fields".to_string(),
        });
    }
    let u = toks[0].parse().map_err(|_| MolError::InvalidSdf {
        line: line_num,
        msg: "invalid bond u".to_string(),
    })?;
    let v = toks[1].parse().map_err(|_| MolError::InvalidSdf {
        line: line_num,
        msg: "invalid bond v".to_string(),
    })?;
    let t = toks[2].parse().unwrap_or(1);
    let s = toks.get(3).and_then(|x| x.parse().ok()).unwrap_or(0);
    Ok((u, v, t, s))
}

fn v2000_bond_order(t: u8) -> BondOrder {
    match t {
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        9 => BondOrder::Dative,
        _ => BondOrder::Single,
    }
}

/// Write a list of `MolGraph` structures as an SDF V2000 multi-record file.
pub fn write_sdf_v2000(mols: &[&MolGraph]) -> String {
    let mut out = String::new();
    for mol in mols {
        let name = mol.name.as_deref().unwrap_or("ChemCraft");
        let _ = writeln!(out, "{name}");
        let _ = writeln!(out, "  ChemCraft-V2000");
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
            mol.atoms.len(),
            mol.bonds.len()
        );
        for a in &mol.atoms {
            let [x, y] = a.coords_2d.unwrap_or([0.0, 0.0]);
            let _ = writeln!(
                out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                x, y, 0.0, a.symbol
            );
        }
        for b in &mol.bonds {
            let t = match b.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Dative => 9,
            };
            let s = match b.stereo {
                BondStereo::Up => 1,
                BondStereo::Down => 6,
                BondStereo::Wavy => 4,
                _ => 0,
            };
            let _ = writeln!(out, "{:>3}{:>3}{:>3}{:>3}  0  0  0", b.u + 1, b.v + 1, t, s);
        }
        let charged: Vec<(usize, i8)> = mol
            .atoms
            .iter()
            .filter(|a| a.formal_charge != 0)
            .map(|a| (a.id + 1, a.formal_charge))
            .collect();
        for chunk in charged.chunks(8) {
            let _ = write!(out, "M  CHG{:>3}", chunk.len());
            for (aid, chg) in chunk {
                let _ = write!(out, "{aid:>4}{chg:>4}");
            }
            let _ = writeln!(out);
        }
        let _ = writeln!(out, "M  END");
        for (k, v) in &mol.properties {
            let _ = writeln!(out, "> <{k}>\n{v}\n");
        }
        let _ = writeln!(out, "$$$$");
    }
    out
}

/// Write a `MolGraph` as a CTAB V3000 block.
pub fn write_sdf_v3000(mol: &MolGraph) -> String {
    let mut out = String::new();
    let name = mol.name.as_deref().unwrap_or("ChemCraft");
    let _ = writeln!(out, "{name}");
    let _ = writeln!(out, "  ChemCraft-V3000");
    let _ = writeln!(out);
    let _ = writeln!(out, "  0  0  0     0  0            999 V3000");
    let _ = writeln!(out, "M  V30 BEGIN CTAB");
    let _ = writeln!(
        out,
        "M  V30 COUNTS {} {} 0 0 0",
        mol.atoms.len(),
        mol.bonds.len()
    );
    let _ = writeln!(out, "M  V30 BEGIN ATOM");
    for a in &mol.atoms {
        let [x, y] = a.coords_2d.unwrap_or([0.0, 0.0]);
        if a.formal_charge != 0 {
            let _ = writeln!(
                out,
                "M  V30 {} {} {:.4} {:.4} 0 0 CHG={}",
                a.id + 1,
                a.symbol,
                x,
                y,
                a.formal_charge
            );
        } else {
            let _ = writeln!(
                out,
                "M  V30 {} {} {:.4} {:.4} 0 0",
                a.id + 1,
                a.symbol,
                x,
                y
            );
        }
    }
    let _ = writeln!(out, "M  V30 END ATOM");
    let _ = writeln!(out, "M  V30 BEGIN BOND");
    for (idx, b) in mol.bonds.iter().enumerate() {
        let t = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::Dative => 9,
        };
        let cfg = match b.stereo {
            BondStereo::Up => " CFG=1",
            BondStereo::Down => " CFG=3",
            BondStereo::Wavy => " CFG=2",
            _ => "",
        };
        let _ = writeln!(
            out,
            "M  V30 {} {} {} {}{}",
            idx + 1,
            t,
            b.u + 1,
            b.v + 1,
            cfg
        );
    }
    let _ = writeln!(out, "M  V30 END BOND");
    let _ = writeln!(out, "M  V30 END CTAB");
    let _ = writeln!(out, "M  END");
    let _ = writeln!(out, "$$$$");
    out
}
