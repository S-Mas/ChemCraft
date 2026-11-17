//! Layer 1 (`mol-vector`): Clean-room 2D coordinate layout, interactive ring snapping,
//! `Clean Structure 2D`, VectorCraft-style vector scene IR, and SVG/EPS/PDF/PNG export.

#![forbid(unsafe_code)]

pub mod export;
pub mod layout;
pub mod scene;

pub use export::{
    VectorError, export_sar_grid_pdf_bytes, export_sar_grid_svg, scene_to_eps, scene_to_pdf_bytes,
    scene_to_png_bytes, scene_to_svg,
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
