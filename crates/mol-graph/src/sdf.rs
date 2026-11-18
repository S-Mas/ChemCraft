use std::collections::VecDeque;
use std::f64::consts::PI;
use std::fmt::Write;

use crate::error::MolError;
use crate::mol::{
    BondCisTrans, BondOrder, BondStereo, ChiralTag, MolAtom, MolGraph, highest_cip_end_nbr,
};
use crate::smiles::DEFAULT_MAX_ATOMS;

/// Parse multi-record SDF text (supports both CTAB **V2000** and **V3000** blocks plus SD tags).
pub fn parse_sdf(input: &str) -> Vec<Result<MolGraph, MolError>> {
    parse_sdf_with_limit(input, usize::MAX)
}

/// Parse at most `max_records` SDF records lazily without allocating beyond the limit.
pub fn parse_sdf_with_limit(input: &str, max_records: usize) -> Vec<Result<MolGraph, MolError>> {
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
        .take(max_records)
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
    if let Some(first_tok) = counts_line.split_whitespace().next()
        && first_tok.len() > 3
        && counts_line
            .get(0..3)
            .is_some_and(|s| s.trim().chars().all(|c| c.is_ascii_digit()))
        && counts_line
            .as_bytes()
            .get(3)
            .is_some_and(u8::is_ascii_digit)
        && counts_line
            .as_bytes()
            .get(6)
            .is_some_and(u8::is_ascii_digit)
    {
        return Err(MolError::InvalidSdf {
            line: 4,
            msg: "corrupt V2000 counts line (>999 columns overflow); use V3000".to_string(),
        });
    }
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
        if counts_line.as_bytes().get(3).is_some_and(|&b| b == b' ')
            && let Some(tok0) = counts_line.split_whitespace().next()
            && tok0.parse::<usize>().is_ok_and(|n| n > 999)
        {
            return Err(MolError::InvalidSdf {
                line: 4,
                msg: "V2000 atom count exceeds 999; use V3000".to_string(),
            });
        }
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
        let z: f64 = toks[2].parse().unwrap_or(0.0);
        let sym = toks[3];
        let mut atom = MolAtom::new(idx, sym)?;
        if x.is_finite() && y.is_finite() {
            atom.coords_2d = Some([x, y]);
            if z.is_finite() && z.abs() > 1e-4 {
                atom.coords_3d = Some([x, y, z]);
            }
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
            3 | 4 => BondStereo::Wavy,
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
            let z: f64 = toks[4].parse().unwrap_or(0.0);
            let mut atom = MolAtom::new(mol.atoms.len(), sym)?;
            if x.is_finite() && y.is_finite() {
                atom.coords_2d = Some([x, y]);
                if z.is_finite() && z.abs() > 1e-4 {
                    atom.coords_3d = Some([x, y, z]);
                }
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
    perceive_sdf_wedges(&mut mol);
    mol.suppress_explicit_hydrogens();
    mol.desalt();
    mol.perceive()?;
    Ok(mol)
}

fn det3_vol(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn stereo_signed_vol_for_seq(u_pos: [f64; 2], nbr_xy: &[[f64; 2]; 4], z: &[f64; 4]) -> f64 {
    let p0 = [nbr_xy[0][0] - u_pos[0], nbr_xy[0][1] - u_pos[1], z[0]];
    let p1 = [nbr_xy[1][0] - u_pos[0], nbr_xy[1][1] - u_pos[1], z[1]];
    let p2 = [nbr_xy[2][0] - u_pos[0], nbr_xy[2][1] - u_pos[1], z[2]];
    let p3 = [nbr_xy[3][0] - u_pos[0], nbr_xy[3][1] - u_pos[1], z[3]];
    let v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let v3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
    det3_vol(v1, v2, v3)
}

fn away_2d_dir(u: usize, mol: &MolGraph, coords: &[[f64; 2]]) -> [f64; 2] {
    let pu = coords[u];
    let mut sx = 0.0;
    let mut sy = 0.0;
    for &(nbr, _) in &mol.adj[u] {
        let pn = coords[nbr];
        let dx = pn[0] - pu[0];
        let dy = pn[1] - pu[1];
        let d = dx.hypot(dy).max(1e-6);
        sx += dx / d;
        sy += dy / d;
    }
    let len = sx.hypot(sy);
    if len < 1e-4 {
        [pu[0] + 1.0, pu[1]]
    } else {
        [pu[0] - 1.2 * sx / len, pu[1] - 1.2 * sy / len]
    }
}

/// Reconstruct tetrahedral SMILES parity (`ChiralTag` + `chiral_seq`) from 2D CTAB coordinates
/// and Up/Down stereo wedges (`b.u` is the narrow end at the chiral center).
fn perceive_sdf_wedges(mol: &mut MolGraph) {
    let coords: Vec<[f64; 2]> = mol
        .atoms
        .iter()
        .map(|a| a.coords_2d.unwrap_or([0.0, 0.0]))
        .collect();
    for u in 0..mol.atoms.len() {
        let mut wedged: Vec<(usize, f64)> = Vec::new();
        for &(nbr, bid) in &mol.adj[u] {
            let b = &mol.bonds[bid];
            if b.u == u && matches!(b.stereo, BondStereo::Up | BondStereo::Down) {
                let zv = if b.stereo == BondStereo::Up {
                    1.0
                } else {
                    -1.0
                };
                wedged.push((nbr, zv));
            }
        }
        if wedged.is_empty() {
            continue;
        }
        let mut seq: Vec<usize> = mol.adj[u].iter().map(|&(nb, _)| nb).collect();
        if seq.len() == 3 {
            seq.push(usize::MAX);
        }
        if seq.len() != 4 {
            continue;
        }
        let mut xy = [[0.0; 2]; 4];
        let mut zs = [0.0; 4];
        for (i, &nid) in seq.iter().enumerate() {
            if nid == usize::MAX {
                xy[i] = away_2d_dir(u, mol, &coords);
            } else {
                xy[i] = coords[nid];
                if let Some(&(_, zv)) = wedged.iter().find(|&&(wn, _)| wn == nid) {
                    zs[i] = zv;
                }
            }
        }
        let vol = stereo_signed_vol_for_seq(coords[u], &xy, &zs);
        if vol.abs() > 1e-5 {
            mol.atoms[u].chiral = if vol > 0.0 {
                ChiralTag::Clockwise
            } else {
                ChiralTag::CounterClockwise
            };
            if seq[3] == usize::MAX {
                mol.atoms[u].chiral_h_slot = Some(3);
            }
            mol.atoms[u].chiral_seq = seq;
        }
    }
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

/// Collect explicit H atoms needed to disambiguate aromatic `[nH]` nitrogens in SDF CTABs.
fn aromatic_nh_hydrogens(mol: &MolGraph, coords: &[[f64; 2]]) -> Vec<(usize, [f64; 2])> {
    let mut extras = Vec::new();
    for a in &mol.atoms {
        if a.aromatic && a.symbol == "N" && a.total_h() > 0 {
            let [x, y] = coords[a.id];
            for k in 0..a.total_h() {
                extras.push((a.id, [x + 0.8 * f64::from(k + 1), y + 0.4]));
            }
        }
    }
    extras
}

/// Generate lightweight trigonal 2D coordinates in pure `mol-graph` when exporting SDF from
/// coordinate-less SMILES so E/Z double-bond geometry and chiral wedges encode faithfully.
fn ensure_export_2d_coords(mol: &MolGraph) -> Vec<[f64; 2]> {
    if mol.atoms.iter().all(|a| a.coords_2d.is_some()) {
        return mol
            .atoms
            .iter()
            .map(|a| a.coords_2d.unwrap_or([0.0, 0.0]))
            .collect();
    }
    let n = mol.atoms.len();
    let bond_len = 1.5_f64;
    let mut pos: Vec<Option<[f64; 2]>> = vec![None; n];
    let mut cursor_x = 0.0_f64;

    // 1. Place rings first (monocyclic regular + fused shared-edge rings)
    if !mol.rings.is_empty() {
        let mut ring_done = vec![false; mol.rings.len()];
        loop {
            let mut progressed = false;
            for (ri, r) in mol.rings.iter().enumerate() {
                if ring_done[ri] || r.is_empty() {
                    continue;
                }
                let placed: Vec<usize> = r.iter().copied().filter(|&v| pos[v].is_some()).collect();
                if placed.is_empty() {
                    if ring_done.iter().all(|&d| !d)
                        || mol
                            .rings
                            .iter()
                            .enumerate()
                            .all(|(rj, _)| ring_done[rj] || rj >= ri)
                    {
                        let k = r.len() as f64;
                        let rad = bond_len / (2.0 * (PI / k).sin());
                        for (i, &v) in r.iter().enumerate() {
                            let ang = 2.0 * PI * (i as f64) / k;
                            pos[v] = Some([cursor_x + rad * ang.cos(), rad * ang.sin()]);
                        }
                        ring_done[ri] = true;
                        progressed = true;
                    }
                    continue;
                }
                if let Some((k1, k2)) = (0..r.len()).find_map(|i| {
                    let j = (i + 1) % r.len();
                    (pos[r[i]].is_some() && pos[r[j]].is_some()).then_some((i, j))
                }) {
                    let p1 = pos[r[k1]].unwrap_or([0.0, 0.0]);
                    let p2 = pos[r[k2]].unwrap_or([bond_len, 0.0]);
                    let mid = [(p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5];
                    let dx = p2[0] - p1[0];
                    let dy = p2[1] - p1[1];
                    let elen = dx.hypot(dy).max(1e-6);
                    let n1 = [-dy / elen, dx / elen];
                    let n2 = [dy / elen, -dx / elen];
                    let mut sx = 0.0;
                    let mut sy = 0.0;
                    let mut cnt = 0.0;
                    for p in pos.iter().flatten() {
                        sx += p[0];
                        sy += p[1];
                        cnt += 1.0;
                    }
                    let c = if cnt > 2.0 { [sx / cnt, sy / cnt] } else { mid };
                    let normal = if (mid[0] + n1[0] - c[0]).hypot(mid[1] + n1[1] - c[1])
                        >= (mid[0] + n2[0] - c[0]).hypot(mid[1] + n2[1] - c[1])
                    {
                        n1
                    } else {
                        n2
                    };
                    let sides = r.len() as f64;
                    let rad = elen / (2.0 * (PI / sides).sin());
                    let apo = (rad * rad - (elen * 0.5).powi(2)).max(0.0).sqrt();
                    let center = [mid[0] + normal[0] * apo, mid[1] + normal[1] * apo];
                    let a1 = (p1[1] - center[1]).atan2(p1[0] - center[0]);
                    let a2 = (p2[1] - center[1]).atan2(p2[0] - center[0]);
                    let mut diff = a2 - a1;
                    while diff > PI {
                        diff -= 2.0 * PI;
                    }
                    while diff < -PI {
                        diff += 2.0 * PI;
                    }
                    let step = if diff < 0.0 {
                        -2.0 * PI / sides
                    } else {
                        2.0 * PI / sides
                    };
                    for (idx, &v) in r.iter().enumerate() {
                        if pos[v].is_none() {
                            let d =
                                (idx as isize - k1 as isize).rem_euclid(r.len() as isize) as f64;
                            let ang = a1 + d * step;
                            pos[v] =
                                Some([center[0] + rad * ang.cos(), center[1] + rad * ang.sin()]);
                        }
                    }
                    ring_done[ri] = true;
                    progressed = true;
                }
            }
            if !progressed {
                if let Some(ri) = ring_done.iter().position(|&d| !d) {
                    let r = &mol.rings[ri];
                    let k = r.len().max(3) as f64;
                    let rad = bond_len / (2.0 * (PI / k).sin());
                    let max_x = pos.iter().flatten().map(|p| p[0]).fold(cursor_x, f64::max);
                    let cx = max_x + bond_len * 2.5;
                    for (i, &v) in r.iter().enumerate() {
                        if pos[v].is_none() {
                            let ang = 2.0 * PI * (i as f64) / k;
                            pos[v] = Some([cx + rad * ang.cos(), rad * ang.sin()]);
                        }
                    }
                    ring_done[ri] = true;
                } else {
                    break;
                }
            }
        }
    }

    // 2. BFS place acyclic components & substituents respecting E/Z double bonds!
    for s in 0..n {
        let mut q = VecDeque::new();
        if pos[s].is_none() {
            let max_x = pos.iter().flatten().map(|p| p[0]).fold(cursor_x, f64::max);
            cursor_x = if pos.iter().any(Option::is_some) {
                max_x + bond_len * 2.2
            } else {
                0.0
            };
            pos[s] = Some([cursor_x, 0.0]);
            q.push_back(s);
        } else {
            // Drain neighbors of already placed ring atoms too
            for (v, pv) in pos.iter().enumerate() {
                if pv.is_some() {
                    q.push_back(v);
                }
            }
        }
        while let Some(u) = q.pop_front() {
            let pu = pos[u].unwrap_or([0.0, 0.0]);
            let placed_nbrs: Vec<usize> = mol.adj[u]
                .iter()
                .map(|&(v, _)| v)
                .filter(|&v| pos[v].is_some())
                .collect();
            let unplaced: Vec<(usize, usize)> = mol.adj[u]
                .iter()
                .copied()
                .filter(|&(v, _)| pos[v].is_none())
                .collect();
            if unplaced.is_empty() {
                continue;
            }

            // Check if `u` is one end of an acyclic E/Z double bond with its other end ALREADY placed (`parent = v_dbl`)
            let dbl_parent = placed_nbrs.iter().copied().find_map(|pn| {
                mol.adj[u].iter().find(|&&(nb, bid)| {
                    nb == pn
                        && mol.bonds[bid].order == BondOrder::Double
                        && matches!(mol.bonds[bid].cis_trans, BondCisTrans::E | BondCisTrans::Z)
                })
            });

            if let Some(&(v_dbl, dbl_bid)) = dbl_parent {
                let pv = pos[v_dbl].unwrap_or([0.0, 0.0]);
                let dx = pu[0] - pv[0];
                let dy = pu[1] - pv[1];
                let alpha = dy.atan2(dx);
                let len = dx.hypot(dy).max(1e-6);
                let nx = -dy / len;
                let ny = dx / len;
                let hi_other = highest_cip_end_nbr(mol, v_dbl, u);
                let hi_self = highest_cip_end_nbr(mol, u, v_dbl);
                let other_side = hi_other
                    .and_then(|ho| pos[ho])
                    .map_or(1.0, |pho| (pho[0] - pv[0]) * nx + (pho[1] - pv[1]) * ny);
                let want_self_hi_sign = if mol.bonds[dbl_bid].cis_trans == BondCisTrans::Z {
                    if other_side >= 0.0 { 1.0 } else { -1.0 }
                } else if other_side >= 0.0 {
                    -1.0
                } else {
                    1.0
                };
                for (i, (w, _)) in unplaced.into_iter().enumerate() {
                    let is_hi = hi_self == Some(w) || (hi_self.is_none() && i == 0);
                    let target_sign = if is_hi {
                        want_self_hi_sign
                    } else {
                        -want_self_hi_sign
                    };
                    // `alpha + PI/3` points in +n direction (`dot = +sin(PI/3)`), `alpha - PI/3` points in -n
                    let ang = if target_sign > 0.0 {
                        alpha + PI / 3.0
                    } else {
                        alpha - PI / 3.0
                    };
                    pos[w] = Some([pu[0] + bond_len * ang.cos(), pu[1] + bond_len * ang.sin()]);
                    q.push_back(w);
                }
                continue;
            }

            let base_ang = if placed_nbrs.is_empty() {
                PI / 6.0
            } else {
                let mut sx = 0.0;
                let mut sy = 0.0;
                for &pn in &placed_nbrs {
                    let pp = pos[pn].unwrap_or(pu);
                    let dx = pp[0] - pu[0];
                    let dy = pp[1] - pu[1];
                    let d = dx.hypot(dy).max(1e-6);
                    sx += dx / d;
                    sy += dy / d;
                }
                (-sy).atan2(-sx)
            };
            let k = unplaced.len();
            for (i, (v, _)) in unplaced.into_iter().enumerate() {
                let offset = if k == 1 {
                    // Tilt slightly on entering/exiting double bonds so E/Z vector never stays collinear
                    let is_dbl_adj = mol.adj[v]
                        .iter()
                        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Double);
                    if is_dbl_adj && placed_nbrs.len() == 1 {
                        PI / 6.0
                    } else {
                        0.0
                    }
                } else {
                    let spread = if k == 2 { PI * 2.0 / 3.0 } else { PI * 0.55 };
                    (i as f64 - (k as f64 - 1.0) * 0.5) * spread
                };
                let ang = base_ang + offset;
                pos[v] = Some([pu[0] + bond_len * ang.cos(), pu[1] + bond_len * ang.sin()]);
                q.push_back(v);
            }
        }
    }

    pos.into_iter().map(|p| p.unwrap_or([0.0, 0.0])).collect()
}

/// Compute per-bond oriented endpoints `(u_out, v_out, stereo_out)` carrying CTAB stereo wedges
/// for every tetrahedral stereocenter in `mol`.
fn sdf_bond_records(
    mol: &MolGraph,
    coords: &[[f64; 2]],
) -> Vec<(usize, usize, BondOrder, BondStereo)> {
    let mut records: Vec<(usize, usize, BondOrder, BondStereo)> = mol
        .bonds
        .iter()
        .map(|b| (b.u, b.v, b.order, b.stereo))
        .collect();

    for u in 0..mol.atoms.len() {
        let a = &mol.atoms[u];
        if !matches!(a.chiral, ChiralTag::Clockwise | ChiralTag::CounterClockwise) {
            continue;
        }
        // Keep existing explicit wedge if already present with narrow tip at `u`
        if records
            .iter()
            .any(|&(bu, _, _, st)| bu == u && matches!(st, BondStereo::Up | BondStereo::Down))
        {
            continue;
        }
        let Some(seq) = mol.effective_chiral_seq(u) else {
            continue;
        };
        let orig_sign: i32 = if a.chiral == ChiralTag::Clockwise {
            1
        } else {
            -1
        };
        // Pick a non-H neighbor `w` to carry the wedge so |det3| is non-degenerate
        let mut candidates: Vec<usize> = seq
            .iter()
            .copied()
            .filter(|&nid| nid != usize::MAX)
            .collect();
        // Prefer degree-1 or non-ring neighbor for clean 2D depiction
        candidates.sort_by_key(|&w| {
            let in_ring = mol.adj[u]
                .iter()
                .find(|&&(nb, _)| nb == w)
                .is_some_and(|&(_, bid)| mol.bonds[bid].in_ring);
            let used_as_wedge = records.iter().any(|&(bu, bv, _, st)| {
                ((bu == u && bv == w) || (bv == u && bu == w))
                    && matches!(st, BondStereo::Up | BondStereo::Down)
            });
            (used_as_wedge, in_ring, mol.adj[w].len())
        });

        let mut chosen: Option<(usize, BondStereo)> = None;
        for w in candidates {
            let mut xy = [[0.0; 2]; 4];
            let mut zs_up = [0.0; 4];
            for (i, &nid) in seq.iter().enumerate() {
                if nid == usize::MAX {
                    xy[i] = away_2d_dir(u, mol, coords);
                } else {
                    xy[i] = coords[nid];
                    if nid == w {
                        zs_up[i] = 1.0;
                    }
                }
            }
            let vol_up = stereo_signed_vol_for_seq(coords[u], &xy, &zs_up);
            if vol_up.abs() > 1e-5 {
                let st = if (vol_up > 0.0) == (orig_sign > 0) {
                    BondStereo::Up
                } else {
                    BondStereo::Down
                };
                chosen = Some((w, st));
                break;
            }
        }
        if let Some((w, st)) = chosen
            && let Some(rec) = records
                .iter_mut()
                .find(|(bu, bv, _, _)| (*bu == u && *bv == w) || (*bu == w && *bv == u))
        {
            *rec = (u, w, rec.2, st);
        }
    }
    records
}

/// Write a list of `MolGraph` structures as an SDF V2000 multi-record file.
pub fn write_sdf_v2000(mols: &[&MolGraph]) -> String {
    let mut out = String::new();
    for mol in mols {
        let coords = ensure_export_2d_coords(mol);
        let bond_recs = sdf_bond_records(mol, &coords);
        let extra_hs = aromatic_nh_hydrogens(mol, &coords);
        let name = mol.name.as_deref().unwrap_or("ChemCraft");
        let _ = writeln!(out, "{name}");
        let _ = writeln!(out, "  ChemCraft-V2000");
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
            mol.atoms.len() + extra_hs.len(),
            bond_recs.len() + extra_hs.len()
        );
        for (idx, a) in mol.atoms.iter().enumerate() {
            let [x, y] = coords[idx];
            let _ = writeln!(
                out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                x, y, 0.0, a.symbol
            );
        }
        for &(_, [hx, hy]) in &extra_hs {
            let _ = writeln!(
                out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                hx, hy, 0.0, "H"
            );
        }
        for (bu, bv, order, stereo) in bond_recs {
            let t = match order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Dative => 9,
            };
            let s = match stereo {
                BondStereo::Up => 1,
                BondStereo::Down => 6,
                BondStereo::Wavy => 4,
                _ => 0,
            };
            let _ = writeln!(out, "{:>3}{:>3}{:>3}{:>3}  0  0  0", bu + 1, bv + 1, t, s);
        }
        for (h_i, &(parent, _)) in extra_hs.iter().enumerate() {
            let h_id = mol.atoms.len() + h_i + 1;
            let _ = writeln!(out, "{:>3}{:>3}{:>3}{:>3}  0  0  0", parent + 1, h_id, 1, 0);
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
    let coords = ensure_export_2d_coords(mol);
    let bond_recs = sdf_bond_records(mol, &coords);
    let extra_hs = aromatic_nh_hydrogens(mol, &coords);
    let name = mol.name.as_deref().unwrap_or("ChemCraft");
    let _ = writeln!(out, "{name}");
    let _ = writeln!(out, "  ChemCraft-V3000");
    let _ = writeln!(out);
    let _ = writeln!(out, "  0  0  0     0  0            999 V3000");
    let _ = writeln!(out, "M  V30 BEGIN CTAB");
    let _ = writeln!(
        out,
        "M  V30 COUNTS {} {} 0 0 0",
        mol.atoms.len() + extra_hs.len(),
        bond_recs.len() + extra_hs.len()
    );
    let _ = writeln!(out, "M  V30 BEGIN ATOM");
    for (idx, a) in mol.atoms.iter().enumerate() {
        let [x, y] = coords[idx];
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
    for (h_i, &(_, [hx, hy])) in extra_hs.iter().enumerate() {
        let _ = writeln!(
            out,
            "M  V30 {} H {:.4} {:.4} 0 0",
            mol.atoms.len() + h_i + 1,
            hx,
            hy
        );
    }
    let _ = writeln!(out, "M  V30 END ATOM");
    let _ = writeln!(out, "M  V30 BEGIN BOND");
    for (idx, (bu, bv, order, stereo)) in bond_recs.into_iter().enumerate() {
        let t = match order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::Dative => 9,
        };
        let cfg = match stereo {
            BondStereo::Up => " CFG=1",
            BondStereo::Down => " CFG=3",
            BondStereo::Wavy => " CFG=2",
            _ => "",
        };
        let _ = writeln!(out, "M  V30 {} {} {} {}{}", idx + 1, t, bu + 1, bv + 1, cfg);
    }
    for (h_i, &(parent, _)) in extra_hs.iter().enumerate() {
        let _ = writeln!(
            out,
            "M  V30 {} 1 {} {}",
            mol.bonds.len() + h_i + 1,
            parent + 1,
            mol.atoms.len() + h_i + 1
        );
    }
    let _ = writeln!(out, "M  V30 END BOND");
    let _ = writeln!(out, "M  V30 END CTAB");
    let _ = writeln!(out, "M  END");
    let _ = writeln!(out, "$$$$");
    out
}
