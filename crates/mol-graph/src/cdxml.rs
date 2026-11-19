//! Clean-room ChemDraw `.cdxml` XML reader and writer (interoperable with ChemDraw & EPAM Indigo/Ketcher).
//!
//! Coordinate convention:
//! - ChemCraft stores 2D coordinates in Cartesian Y-up units (`bond_len ≈ 1.5`).
//! - CDXML `p="x y"` uses screen Y-down points scaled to ACS Document 1996 `BondLength="14.400000"`.
//! - Negating Y on export/import preserves exact 2D page appearance and tetrahedral wedge parity.

use std::collections::HashMap;
use std::fmt::Write;

use crate::element::element_by_number;
use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, MolAtom, MolGraph};
use crate::sdf::{ensure_export_2d_coords, perceive_sdf_wedges, sdf_bond_records};

/// ACS Document 1996 default bond length in CDXML points (0.20 in = 14.4 pt).
const ACS_CDXML_BOND_LEN: f64 = 14.4;

/// Serialize `mol` into a standards-compliant ChemDraw `.cdxml` XML document (ACS Document 1996 settings).
pub fn write_cdxml(mol: &MolGraph) -> String {
    let coords = ensure_export_2d_coords(mol);
    let bond_recs = sdf_bond_records(mol, &coords);

    // Estimate current mean bond length so we scale cleanly to ACS 14.4 pt
    let mut bl_sum = 0.0_f64;
    let mut bl_cnt = 0usize;
    for b in &mol.bonds {
        let p0 = coords[b.u];
        let p1 = coords[b.v];
        let d = (p1[0] - p0[0]).hypot(p1[1] - p0[1]);
        if d > 1e-4 {
            bl_sum += d;
            bl_cnt += 1;
        }
    }
    let avg_bl = if bl_cnt > 0 {
        bl_sum / (bl_cnt as f64)
    } else {
        1.5
    };
    let scale = ACS_CDXML_BOND_LEN / avg_bl.max(1e-3);

    // Transform Cartesian Y-up to CDXML Y-down points with generous top-left page margin (72 pt)
    let mut pts_cd: Vec<[f64; 2]> = coords
        .iter()
        .map(|&[x, y]| [x * scale, -y * scale])
        .collect();
    let (mut min_x, mut min_y, mut max_x, mut max_y) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for &[x, y] in &pts_cd {
        min_x = min_x.min(x);
        min_y = min_y.min(y);
        max_x = max_x.max(x);
        max_y = max_y.max(y);
    }
    if min_x.is_finite() && min_y.is_finite() {
        let ox = 72.0 - min_x;
        let oy = 72.0 - min_y;
        for p in &mut pts_cd {
            p[0] += ox;
            p[1] += oy;
        }
        max_x += ox;
        max_y += oy;
        min_x = 72.0;
        min_y = 72.0;
    } else {
        min_x = 72.0;
        min_y = 72.0;
        max_x = 144.0;
        max_y = 144.0;
    }

    let mut out = String::new();
    let _ = writeln!(out, r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    let _ = writeln!(
        out,
        r#"<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">"#
    );
    let _ = writeln!(
        out,
        r#"<CDXML CreationProgram="ChemCraft" BondLength="{ACS_CDXML_BOND_LEN:.6}" HashSpacing="2.500000" LineWidth="0.600000" BoldWidth="2.000000" MarginWidth="1.600000" LabelFont="3" LabelSize="10" CaptionFont="3" CaptionSize="10">"#
    );
    let _ = writeln!(
        out,
        r#"  <colortable>
    <color r="1" g="1" b="1"/>
    <color r="0" g="0" b="0"/>
    <color r="1" g="0" b="0"/>
    <color r="1" g="1" b="0"/>
    <color r="0" g="1" b="0"/>
    <color r="0" g="1" b="1"/>
    <color r="0" g="0" b="1"/>
    <color r="1" g="0" b="1"/>
  </colortable>"#
    );
    let _ = writeln!(
        out,
        r#"  <fonttable>
    <font id="3" charset="iso-8859-1" name="Arial"/>
    <font id="4" charset="iso-8859-1" name="Times New Roman"/>
  </fonttable>"#
    );
    let _ = writeln!(
        out,
        r#"  <page id="1" BoundingBox="{min_x:.4} {min_y:.4} {max_x:.4} {max_y:.4}" WidthPages="1" HeightPages="1">"#
    );
    let _ = writeln!(out, r#"    <fragment id="2">"#);

    let mut next_id = 10usize;
    let mut atom_xml_ids = Vec::with_capacity(mol.atoms.len());
    for (idx, a) in mol.atoms.iter().enumerate() {
        let xml_id = next_id;
        next_id += 1;
        atom_xml_ids.push(xml_id);
        let [px, py] = pts_cd[idx];
        let _ = write!(out, r#"      <n id="{xml_id}" p="{px:.4} {py:.4}""#);
        if a.atomic_number != 6 {
            let _ = write!(out, r#" Element="{}""#, a.atomic_number);
        }
        if a.formal_charge != 0 {
            let _ = write!(out, r#" Charge="{}""#, a.formal_charge);
        }
        if let Some(iso) = a.isotope {
            let _ = write!(out, r#" Isotope="{iso}""#);
        }
        // Emit explicit NumHydrogens for heteroatoms & aromatic [nH] so ChemDraw/Indigo keep exact valences
        let need_nh = a.atomic_number != 6 || a.from_bracket || a.formal_charge != 0;
        if need_nh {
            let _ = write!(out, r#" NumHydrogens="{}""#, a.total_h());
        }
        if a.atomic_number != 6 || a.formal_charge != 0 {
            let label = format_cdxml_label(a);
            let _ = writeln!(out, ">");
            let _ = writeln!(
                out,
                r#"        <t p="{px:.4} {py:.4}" LabelJustification="Left"><s font="3" size="10">{label}</s></t>"#
            );
            let _ = writeln!(out, "      </n>");
        } else {
            let _ = writeln!(out, "/>");
        }
    }

    for (bu, bv, order, stereo) in bond_recs {
        let bid = next_id;
        next_id += 1;
        let b_xml = atom_xml_ids[bu];
        let e_xml = atom_xml_ids[bv];
        let _ = write!(out, r#"      <b id="{bid}" B="{b_xml}" E="{e_xml}""#);
        match order {
            BondOrder::Single => {}
            BondOrder::Double => {
                let _ = write!(out, r#" Order="2""#);
            }
            BondOrder::Triple => {
                let _ = write!(out, r#" Order="3""#);
            }
            BondOrder::Aromatic => {
                let _ = write!(out, r#" Order="1.5""#);
            }
            BondOrder::Dative => {
                let _ = write!(out, r#" Order="dative""#);
            }
        }
        match stereo {
            BondStereo::Up => {
                let _ = write!(out, r#" Display="WedgeBegin""#);
            }
            BondStereo::Down => {
                let _ = write!(out, r#" Display="WedgedHashBegin""#);
            }
            BondStereo::Wavy => {
                let _ = write!(out, r#" Display="Wavy""#);
            }
            _ => {}
        }
        let _ = writeln!(out, "/>");
    }

    let _ = writeln!(out, "    </fragment>");
    let _ = writeln!(out, "  </page>");
    let _ = writeln!(out, "</CDXML>");
    out
}

fn format_cdxml_label(a: &MolAtom) -> String {
    let mut s = a.symbol.clone();
    let h = a.total_h();
    if h == 1 {
        s.push('H');
    } else if h > 1 {
        let _ = write!(s, "H{h}");
    }
    if a.formal_charge == 1 {
        s.push('+');
    } else if a.formal_charge == -1 {
        s.push('-');
    } else if a.formal_charge > 1 {
        let _ = write!(s, "{}+", a.formal_charge);
    } else if a.formal_charge < -1 {
        let _ = write!(s, "{}-", -a.formal_charge);
    }
    s
}

/// Parse a ChemDraw `.cdxml` XML document into a perceived `MolGraph` with 2D coords and stereochemistry.
pub fn parse_cdxml(xml: &str) -> Result<MolGraph, MolError> {
    if !xml.contains("<CDXML") && !xml.contains("<fragment") {
        return Err(MolError::InvalidFormat(
            "Input does not look like a CDXML document (missing <CDXML> / <fragment>)".to_string(),
        ));
    }
    let mut mol = MolGraph::new();
    let mut bond_length_pt = ACS_CDXML_BOND_LEN;
    let mut id_to_atom: HashMap<String, usize> = HashMap::new();
    let mut pending_bonds: Vec<(String, String, BondOrder, BondStereo)> = Vec::new();

    // Lightweight, zero-dependency XML tag scanner safe against arbitrary whitespace/newlines
    let bytes = xml.as_bytes();
    let mut pos = 0usize;
    while pos < bytes.len() {
        let Some(lt_rel) = xml[pos..].find('<') else {
            break;
        };
        let tag_start = pos + lt_rel;
        if xml[tag_start..].starts_with("<!--") {
            if let Some(end_rel) = xml[tag_start + 4..].find("-->") {
                pos = tag_start + 4 + end_rel + 3;
                continue;
            }
            break;
        }
        let Some(gt_rel) = xml[tag_start..].find('>') else {
            break;
        };
        let tag_body = xml[tag_start + 1..tag_start + gt_rel].trim();
        pos = tag_start + gt_rel + 1;

        if tag_body.starts_with('/') || tag_body.starts_with('?') || tag_body.starts_with('!') {
            continue;
        }
        let tag_body = tag_body.trim_end_matches('/').trim();
        let mut split = tag_body.splitn(2, |c: char| c.is_ascii_whitespace());
        let tag_name = split.next().unwrap_or("");
        let attr_str = split.next().unwrap_or("");
        let attrs = parse_xml_attrs(attr_str);

        match tag_name {
            "CDXML" => {
                if let Some(bl) = attrs.get("BondLength").and_then(|v| v.parse::<f64>().ok())
                    && bl > 1.0
                {
                    bond_length_pt = bl;
                }
            }
            "n" => {
                let Some(xml_id) = attrs.get("id") else {
                    continue;
                };
                // Skip non-atom graphic/nickname wrapper nodes if NodeType is ExternalConnectionPoint
                if attrs
                    .get("NodeType")
                    .is_some_and(|nt| nt == "ExternalConnectionPoint")
                {
                    continue;
                }
                let atomic_num = attrs
                    .get("Element")
                    .and_then(|v| v.parse::<u8>().ok())
                    .unwrap_or(6);
                let sym = element_by_number(atomic_num).map_or("C", |e| e.symbol);
                let mut atom = MolAtom::new(0, sym)?;
                if let Some(chg) = attrs.get("Charge").and_then(|v| v.parse::<i8>().ok()) {
                    atom.formal_charge = chg;
                }
                if let Some(iso) = attrs.get("Isotope").and_then(|v| v.parse::<u16>().ok()) {
                    atom.isotope = Some(iso);
                }
                if let Some(nh) = attrs.get("NumHydrogens").and_then(|v| v.parse::<u8>().ok()) {
                    atom.explicit_h = nh;
                    atom.from_bracket = true;
                }
                if let Some(p_str) = attrs.get("p") {
                    let nums: Vec<f64> = p_str
                        .split_whitespace()
                        .filter_map(|t| t.parse::<f64>().ok())
                        .collect();
                    if nums.len() >= 2 && nums[0].is_finite() && nums[1].is_finite() {
                        // Convert CDXML Y-down points to Cartesian Y-up coordinates (1.5 Å scale)
                        let s = 1.5 / bond_length_pt.max(1.0);
                        atom.coords_2d = Some([nums[0] * s, -nums[1] * s]);
                    }
                }
                let aid = mol.add_atom(atom);
                id_to_atom.insert(xml_id.clone(), aid);
            }
            "b" => {
                let (Some(b_id), Some(e_id)) = (attrs.get("B"), attrs.get("E")) else {
                    continue;
                };
                let order = match attrs.get("Order").map(String::as_str) {
                    Some("2") => BondOrder::Double,
                    Some("3") => BondOrder::Triple,
                    Some("1.5") | Some("aromatic") => BondOrder::Aromatic,
                    Some("dative") => BondOrder::Dative,
                    _ => BondOrder::Single,
                };
                let (u_id, v_id, stereo) = match attrs.get("Display").map(String::as_str) {
                    Some("WedgeBegin") => (b_id.clone(), e_id.clone(), BondStereo::Up),
                    Some("WedgeEnd") => (e_id.clone(), b_id.clone(), BondStereo::Up),
                    Some("WedgedHashBegin") | Some("HashBegin") => {
                        (b_id.clone(), e_id.clone(), BondStereo::Down)
                    }
                    Some("WedgedHashEnd") | Some("HashEnd") => {
                        (e_id.clone(), b_id.clone(), BondStereo::Down)
                    }
                    Some("Wavy") => (b_id.clone(), e_id.clone(), BondStereo::Wavy),
                    _ => (b_id.clone(), e_id.clone(), BondStereo::None),
                };
                pending_bonds.push((u_id, v_id, order, stereo));
            }
            _ => {}
        }
    }

    if mol.atoms.is_empty() {
        return Err(MolError::EmptyInput);
    }

    for (u_xml, v_xml, order, stereo) in pending_bonds {
        if let (Some(&u), Some(&v)) = (id_to_atom.get(&u_xml), id_to_atom.get(&v_xml))
            && u != v
        {
            let _ = mol.add_bond(u, v, order, stereo);
        }
    }

    // Center 2D coordinates neatly around origin
    let mut sx = 0.0_f64;
    let mut sy = 0.0_f64;
    let mut cnt = 0usize;
    for a in &mol.atoms {
        if let Some([x, y]) = a.coords_2d {
            sx += x;
            sy += y;
            cnt += 1;
        }
    }
    if cnt > 0 {
        let cx = sx / (cnt as f64);
        let cy = sy / (cnt as f64);
        for a in &mut mol.atoms {
            if let Some([x, y]) = &mut a.coords_2d {
                *x -= cx;
                *y -= cy;
            }
        }
    }

    // Reconstruct tetrahedral chirality from CDXML wedges + 2D geometry, then perceive rings/aromaticity/CIP
    perceive_sdf_wedges(&mut mol);
    mol.suppress_explicit_hydrogens();
    mol.perceive()?;
    Ok(mol)
}

fn parse_xml_attrs(s: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let bytes = s.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        let k_start = i;
        while i < bytes.len() && bytes[i] != b'=' && !bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        let key = &s[k_start..i];
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'=' {
            continue;
        }
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        let quote = bytes[i];
        if quote != b'"' && quote != b'\'' {
            continue;
        }
        i += 1;
        let v_start = i;
        while i < bytes.len() && bytes[i] != quote {
            i += 1;
        }
        let val = &s[v_start..i];
        if i < bytes.len() {
            i += 1;
        }
        if !key.is_empty() {
            map.insert(key.to_string(), val.to_string());
        }
    }
    map
}
