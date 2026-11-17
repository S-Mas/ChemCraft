//! Layer 1 (`mol-vector`): Clean-room 2D coordinate layout, interactive ring snapping,
//! `Clean Structure 2D`, VectorCraft-style vector scene IR, and SVG/EPS/PDF/PNG export.

#![forbid(unsafe_code)]

pub mod export;
pub mod layout;
pub mod scene;

pub use export::{
    VectorError, export_sar_grid_eps, export_sar_grid_pdf_bytes, export_sar_grid_png_bytes,
    export_sar_grid_svg, scene_to_eps, scene_to_pdf_bytes, scene_to_png_bytes, scene_to_svg,
};
pub use layout::{
    RingTemplate, STANDARD_BOND_LEN, SnapTarget, clean_structure_2d, ensure_2d_coords,
    generate_2d_coords, stamp_ring_template,
};
pub use scene::{
    ReactionArrow, Rect2, Rgba, SceneStyle, VectorElement, VectorScene, build_mol_scene,
};

#[cfg(test)]
mod tests {
    use super::*;
    use mol_graph::{MolError, parse_smiles};

    #[test]
    fn renders_10_complex_molecules_to_svg_eps_pdf_png_without_text_overlaps()
    -> Result<(), Box<dyn std::error::Error>> {
        let complex = [
            ("Testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
            ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
            (
                "Sildenafil",
                "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
            ),
            ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
            ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
            ("Morphine", "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"),
            ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
            ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
            ("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
            (
                "Ciprofloxacin",
                "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
            ),
        ];

        let style = SceneStyle::default();
        let mut grid_mols = Vec::new();
        for (name, smi) in complex {
            let mut mol = parse_smiles(smi)?;
            mol.name = Some(name.to_string());
            generate_2d_coords(&mut mol);
            let scene = build_mol_scene(&mol, &style, &[]);
            assert_eq!(
                scene.label_overlap_count(),
                0,
                "{name} had overlapping atom labels"
            );
            let svg = scene_to_svg(&scene);
            assert!(svg.contains("<svg") && svg.contains("</svg>"));
            let eps = scene_to_eps(&scene);
            assert!(eps.starts_with("%!PS-Adobe-3.0 EPSF-3.0"));
            let pdf = scene_to_pdf_bytes(&scene)?;
            assert!(pdf.starts_with(b"%PDF-1.4") && pdf.contains(&b'/'));
            let png = scene_to_png_bytes(&scene, 1.0)?;
            assert!(png.starts_with(&[137, 80, 78, 71]));
            grid_mols.push((name, mol));
        }

        let refs: Vec<(&str, &mol_graph::MolGraph)> =
            grid_mols.iter().map(|(n, m)| (*n, m)).collect();
        let grid_svg = export_sar_grid_svg(&refs, 4, 240.0, 200.0, false);
        assert!(grid_svg.contains("Testosterone") && grid_svg.contains("Sildenafil"));
        let grid_pdf = export_sar_grid_pdf_bytes(&refs, 4, 240.0, 200.0)?;
        assert!(grid_pdf.starts_with(b"%PDF-1.4"));
        let grid_eps = export_sar_grid_eps(&refs, 4, 240.0, 200.0, false);
        assert!(grid_eps.contains("Testosterone") && grid_eps.contains("Ciprofloxacin"));
        let grid_png = export_sar_grid_png_bytes(&refs, 4, 240.0, 200.0, false, 1.0)?;
        assert!(grid_png.starts_with(&[137, 80, 78, 71]));
        Ok(())
    }

    #[test]
    fn pdf_renders_circles_winansi_accents_and_screen_to_world_roundtrips()
    -> Result<(), Box<dyn std::error::Error>> {
        let style = SceneStyle::default();
        let mut mol = parse_smiles("c1ccccc1")?;
        generate_2d_coords(&mut mol);
        let mut scene = build_mol_scene(&mol, &style, &[]);
        // Inverse screen -> world roundtrip check
        let sp = scene.atom_screen_pos[0];
        let wp = scene.screen_to_world(sp[0], sp[1]);
        let orig = mol.atoms[0].coords_2d.unwrap_or([0.0, 0.0]);
        assert!((wp[0] - orig[0]).abs() < 1e-6);
        assert!((wp[1] - orig[1]).abs() < 1e-6);

        // Add explicit Circle + accented Label
        scene.elements.push(VectorElement::Circle {
            cx: 120.0,
            cy: 90.0,
            r: 14.0,
            fill: Some(Rgba(225, 45, 45, 55)),
            stroke: Some((Rgba(220, 38, 38, 220), 1.5)),
        });
        scene.elements.push(VectorElement::Label {
            text: "Cafeína".to_string(),
            cx: 200.0,
            cy: 40.0,
            font_size: 12.0,
            color: Rgba(20, 20, 20, 255),
            bbox: Rect2 {
                x: 170.0,
                y: 34.0,
                w: 60.0,
                h: 12.0,
            },
        });
        let eps = scene_to_eps(&scene);
        assert!(eps.contains("Cafe\\355na"));
        assert!(eps.contains("arc fill") && eps.contains("arc stroke"));

        // Decompress PDF content stream to verify circle Bezier curves + WinAnsi octal!
        let pdf = scene_to_pdf_bytes(&scene)?;
        let start = pdf
            .windows(7)
            .position(|w| w == b"stream\n")
            .ok_or("missing pdf stream")?
            + 7;
        let end = pdf
            .windows(10)
            .position(|w| w == b"\nendstream")
            .ok_or("missing pdf endstream")?;
        let mut dec = flate2::read::ZlibDecoder::new(&pdf[start..end]);
        let mut text = String::new();
        std::io::Read::read_to_string(&mut dec, &mut text)?;
        assert!(text.contains("Cafe\\355na"));
        assert!(text.contains(" c f") && text.contains(" c S"));
        Ok(())
    }

    #[test]
    fn ring_snapping_fuses_benzene_to_naphthalene() -> Result<(), MolError> {
        let mut mol = parse_smiles("c1ccccc1")?;
        generate_2d_coords(&mut mol);
        stamp_ring_template(
            &mut mol,
            RingTemplate::Benzene,
            SnapTarget::Bond {
                bond_idx: 0,
                cursor_x: 5.0,
                cursor_y: 0.0,
            },
        )?;
        assert_eq!(mol.atoms.len(), 10);
        assert_eq!(mol.formula(), "C10H8");
        Ok(())
    }
}
