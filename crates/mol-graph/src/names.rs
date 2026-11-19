//! Clean-room ChemDraw-equivalent Name-to-Structure converter:
//! 1. Curated registry of 110+ drugs, reagents, solvents, amino acids, heterocycles, and building blocks.
//! 2. Algorithmic IUPAC linear organic nomenclature parser supporting:
//!    - Carbon chains (`meth`..`dodeca`)
//!    - Unsaturation (`-ane`, `-ene`, `-yne`, positional locants like `2-butene`, `1-pentyne`)
//!    - Functional suffixes (`-ol`, `-al`, `-one`, `-oic acid`, `-amine`, `-nitrile`, `-thiol`, `-amide`)
//!    - Substituent prefixes with locants & multipliers (`2-methyl`, `2,2-dimethyl`, `3-chloro`, `4-fluoro`, `2-bromo`, `3-iodo`, `2-hydroxy`, `3-methoxy`, `2-amino`, `4-nitro`, `2-phenyl`, `3-ethyl`, `2-propyl`, `cyclo-` rings, `o-/m-/p-` disubstituted benzenes).

use crate::error::MolError;
use crate::mol::{BondOrder, BondStereo, MolAtom, MolGraph};
use crate::smiles::parse_smiles;

/// Curated dictionary of common medicinal chemistry, organic synthesis, solvent, and biochemical names.
const COMMON_NAMES: &[(&str, &str)] = &[
    // Reference Drugs & Natural Products
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("acetylsalicylic acid", "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("(s)-ibuprofen", "CC(C)Cc1ccc([C@H](C)C(=O)O)cc1"),
    ("(r)-ibuprofen", "CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
    ("naproxen", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
    ("diclofenac", "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    (
        "celecoxib",
        "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
    ),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("theophylline", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"),
    ("theobromine", "Cn1cnc2c1c(=O)[nH]c(=O)n2C"),
    ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
    ("dopamine", "NCCc1ccc(O)c(O)c1"),
    ("serotonin", "NCCc1c[nH]c2ccc(O)cc12"),
    ("histamine", "NCCc1c[nH]cn1"),
    ("adrenaline", "CNC[C@H](O)c1ccc(O)c(O)c1"),
    ("epinephrine", "CNC[C@H](O)c1ccc(O)c(O)c1"),
    ("melatonin", "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"),
    ("diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
    (
        "sildenafil",
        "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
    ),
    (
        "atorvastatin",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    ),
    (
        "imatinib",
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    ),
    ("erlotinib", "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"),
    (
        "ciprofloxacin",
        "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    ),
    (
        "penicillin g",
        "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    ),
    (
        "amoxicillin",
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O",
    ),
    ("metformin", "CN(C)C(=N)NC(=N)N"),
    ("omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
    (
        "losartan",
        "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1",
    ),
    ("propranolol", "CC(C)NCC(O)COc1cccc2ccccc12"),
    ("lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    (
        "morphine",
        "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
    ),
    (
        "testosterone",
        "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
    ),
    (
        "estradiol",
        "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O",
    ),
    (
        "cholesterol",
        "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(C)CC[C@H]3[C@H]2CC=C4C[C@@H](O)CC[C@]34C",
    ),
    (
        "dexamethasone",
        "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO",
    ),
    (
        "artemisinin",
        "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]23OO4",
    ),
    ("warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    (
        "tamoxifen",
        "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
    ),
    ("sulfamethoxazole", "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"),
    ("fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"),
    (
        "haloperidol",
        "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
    ),
    // Common Solvents & Reagents
    ("water", "O"),
    ("ammonia", "N"),
    ("methanol", "CO"),
    ("ethanol", "CCO"),
    ("isopropanol", "CC(O)C"),
    ("ipa", "CC(O)C"),
    ("acetone", "CC(=O)C"),
    ("acetonitrile", "CC#N"),
    ("mecn", "CC#N"),
    ("dmso", "CS(=O)C"),
    ("dimethyl sulfoxide", "CS(=O)C"),
    ("dmf", "CN(C)C=O"),
    ("dimethylformamide", "CN(C)C=O"),
    ("thf", "C1CCOC1"),
    ("tetrahydrofuran", "C1CCOC1"),
    ("dioxane", "C1COCCO1"),
    ("1,4-dioxane", "C1COCCO1"),
    ("dcm", "ClCCl"),
    ("dichloromethane", "ClCCl"),
    ("chloroform", "ClC(Cl)Cl"),
    ("carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
    ("diethyl ether", "CCOCC"),
    ("ether", "CCOCC"),
    ("ethyl acetate", "CCOC(=O)C"),
    ("etoac", "CCOC(=O)C"),
    ("acetic acid", "CC(=O)O"),
    ("formic acid", "O=CO"),
    ("trifluoroacetic acid", "O=C(O)C(F)(F)F"),
    ("tfa", "O=C(O)C(F)(F)F"),
    ("triethylamine", "CCN(CC)CC"),
    ("tea", "CCN(CC)CC"),
    ("dipea", "CCN(C(C)C)C(C)C"),
    ("hunig's base", "CCN(C(C)C)C(C)C"),
    ("dmap", "CN(C)c1ccncc1"),
    ("urea", "NC(=O)N"),
    ("guanidine", "NC(=N)N"),
    ("formaldehyde", "C=O"),
    ("acetaldehyde", "CC=O"),
    ("benzaldehyde", "O=Cc1ccccc1"),
    ("benzoic acid", "O=C(O)c1ccccc1"),
    ("salicylic acid", "O=C(O)c1ccccc1O"),
    ("cinnamic acid", "O=C(O)/C=C/c1ccccc1"),
    ("citric acid", "O=C(O)CC(O)(CC(=O)O)C(=O)O"),
    ("lactic acid", "CC(O)C(=O)O"),
    ("pyruvic acid", "CC(=O)C(=O)O"),
    ("oxalic acid", "O=C(O)C(=O)O"),
    ("malonic acid", "O=C(O)CC(=O)O"),
    ("succinic acid", "O=C(O)CCC(=O)O"),
    // Aromatics, Heterocycles & Cages
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("aniline", "Nc1ccccc1"),
    ("anisole", "COc1ccccc1"),
    ("styrene", "C=Cc1ccccc1"),
    ("nitrobenzene", "O=[N+]([O-])c1ccccc1"),
    ("chlorobenzene", "Clc1ccccc1"),
    ("bromobenzene", "Brc1ccccc1"),
    ("fluorobenzene", "Fc1ccccc1"),
    ("benzonitrile", "N#Cc1ccccc1"),
    ("acetophenone", "CC(=O)c1ccccc1"),
    ("benzophenone", "O=C(c1ccccc1)c1ccccc1"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("phenanthrene", "c1ccc2c(c1)ccc1ccccc12"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("pyridine", "c1ccncc1"),
    ("pyrimidine", "c1cncnc1"),
    ("pyrazine", "c1cnccn1"),
    ("pyridazine", "c1ccnnc1"),
    ("pyrrole", "c1cc[nH]c1"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("imidazole", "c1c[nH]cn1"),
    ("pyrazole", "c1cn[nH]c1"),
    ("oxazole", "c1cocn1"),
    ("isoxazole", "c1ccon1"),
    ("thiazole", "c1cscn1"),
    ("tetrazole", "c1nnn[nH]1"),
    ("triazole", "c1cn[nH]n1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("benzimidazole", "c1ccc2[nH]cnc2c1"),
    ("quinoline", "c1ccc2ncccc2c1"),
    ("isoquinoline", "c1ccc2cnccc2c1"),
    ("purine", "c1ncc2[nH]cnc2n1"),
    ("piperidine", "C1CCNCC1"),
    ("piperazine", "C1CNCCN1"),
    ("morpholine", "C1COCCN1"),
    ("pyrrolidine", "C1CCNC1"),
    ("adamantane", "C1C2CC3CC1CC(C2)C3"),
    ("cubane", "C12C3C4C1C5C2C3C45"),
    ("norbornane", "C1CC2CCC1C2"),
    // Nucleobases, Sugars & Amino Acids
    ("adenine", "Nc1ncnc2[nH]cnc12"),
    ("guanine", "Nc1nc2[nH]cnc2c(=O)[nH]1"),
    ("cytosine", "Nc1cc[nH]c(=O)n1"),
    ("thymine", "Cc1c[nH]c(=O)[nH]c1=O"),
    ("uracil", "O=c1cc[nH]c(=O)[nH]1"),
    ("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("d-glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (
        "atp",
        "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    ),
    ("glycine", "NCC(=O)O"),
    ("alanine", "C[C@H](N)C(=O)O"),
    ("l-alanine", "C[C@H](N)C(=O)O"),
    ("valine", "CC(C)[C@H](N)C(=O)O"),
    ("leucine", "CC(C)C[C@H](N)C(=O)O"),
    ("isoleucine", "CC[C@H](C)[C@H](N)C(=O)O"),
    ("proline", "O=C(O)[C@@H]1CCCN1"),
    ("phenylalanine", "N[C@@H](Cc1ccccc1)C(=O)O"),
    ("tyrosine", "N[C@@H](Cc1ccc(O)cc1)C(=O)O"),
    ("tryptophan", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    ("serine", "N[C@@H](CO)C(=O)O"),
    ("threonine", "C[C@@H](O)[C@H](N)C(=O)O"),
    ("cysteine", "N[C@@H](CS)C(=O)O"),
    ("methionine", "CSCC[C@H](N)C(=O)O"),
    ("asparagine", "NC(=O)C[C@H](N)C(=O)O"),
    ("glutamine", "NC(=O)CC[C@H](N)C(=O)O"),
    ("aspartic acid", "N[C@@H](CC(=O)O)C(=O)O"),
    ("glutamic acid", "N[C@@H](CCC(=O)O)C(=O)O"),
    ("lysine", "NCCCC[C@H](N)C(=O)O"),
    ("arginine", "NC(=N)NCCC[C@H](N)C(=O)O"),
    ("histidine", "N[C@@H](Cc1c[nH]cn1)C(=O)O"),
];

/// Convert a common drug/reagent name or systematic IUPAC name into a perceived `MolGraph`.
pub fn resolve_chemical_name(name: &str) -> Result<MolGraph, MolError> {
    let raw = name.trim();
    if raw.is_empty() {
        return Err(MolError::EmptyInput);
    }
    let lower = raw.to_ascii_lowercase();
    let norm: String = lower.split_whitespace().collect::<Vec<_>>().join(" ");

    // 1. Exact curated dictionary match
    if let Some(&(_, smi)) = COMMON_NAMES.iter().find(|&&(k, _)| k == norm) {
        let mut mol = parse_smiles(smi)?;
        mol.name = Some(raw.to_string());
        return Ok(mol);
    }

    // 2. Algorithmic IUPAC systematic name parser (o/m/p-disubstituted benzenes + linear/cyclic alkanes/alkenes/alkynes/functional groups)
    if let Some(mut mol) = parse_substituted_benzene(&norm)? {
        mol.name = Some(raw.to_string());
        return Ok(mol);
    }
    if let Some(mut mol) = parse_iupac_aliphatic(&norm)? {
        mol.name = Some(raw.to_string());
        return Ok(mol);
    }

    Err(MolError::InvalidFormat(format!(
        "Could not resolve chemical name '{raw}' (not in common registry or supported IUPAC subset)"
    )))
}

/// Return list of built-in dictionary names for autocomplete / CLI help.
pub fn known_chemical_names() -> Vec<&'static str> {
    COMMON_NAMES.iter().map(|&(k, _)| k).collect()
}

fn subst_smiles(prefix: &str) -> Option<&'static str> {
    match prefix {
        "methyl" => Some("C"),
        "ethyl" => Some("CC"),
        "propyl" => Some("CCC"),
        "isopropyl" => Some("C(C)C"),
        "butyl" => Some("CCCC"),
        "tert-butyl" | "t-butyl" => Some("C(C)(C)C"),
        "fluoro" => Some("F"),
        "chloro" => Some("Cl"),
        "bromo" => Some("Br"),
        "iodo" => Some("I"),
        "hydroxy" => Some("O"),
        "methoxy" => Some("OC"),
        "ethoxy" => Some("OCC"),
        "amino" => Some("N"),
        "dimethylamino" => Some("N(C)C"),
        "nitro" => Some("[N+](=O)[O-]"),
        "cyano" => Some("C#N"),
        "trifluoromethyl" => Some("C(F)(F)F"),
        "phenyl" => Some("c1ccccc1"),
        "mercapto" | "sulfanyl" => Some("S"),
        "acetyl" => Some("C(=O)C"),
        _ => None,
    }
}

fn parse_substituted_benzene(norm: &str) -> Result<Option<MolGraph>, MolError> {
    // Handle patterns like "4-chlorophenol", "p-nitrophenol", "2,4-dichlorotoluene", "1,4-dichlorobenzene", "3-aminobenzoic acid"
    let bases: &[(&str, &str)] = &[
        ("benzoic acid", "O=C(O)c1ccccc1"),
        ("benzaldehyde", "O=Cc1ccccc1"),
        ("benzonitrile", "N#Cc1ccccc1"),
        ("benzene", "c1ccccc1"),
        ("toluene", "Cc1ccccc1"),
        ("phenol", "Oc1ccccc1"),
        ("aniline", "Nc1ccccc1"),
        ("anisole", "COc1ccccc1"),
    ];
    for &(base_name, base_smi) in bases {
        if let Some(prefix_raw) = norm.strip_suffix(base_name) {
            let prefix_raw = prefix_raw.trim().trim_end_matches('-');
            if prefix_raw.is_empty() {
                continue;
            }
            let subs = parse_substituent_list(prefix_raw);
            if subs.is_empty() {
                continue;
            }
            let mut mol = parse_smiles(base_smi)?;
            // Find the 6-membered aromatic ring atom indices starting at C1 (substituted atom if any)
            let Some(ring) = mol.rings.iter().find(|r| r.len() == 6).cloned() else {
                continue;
            };
            // Rotate ring so index 0 is the substituted carbon (degree 3) if present
            let c1_pos = ring
                .iter()
                .position(|&a| mol.adj[a].len() >= 3)
                .unwrap_or(0);
            let mut ordered = Vec::with_capacity(6);
            for k in 0..6 {
                ordered.push(ring[(c1_pos + k) % 6]);
            }
            for (loc, smi) in subs {
                let idx = loc.saturating_sub(1);
                if idx >= 6 {
                    return Ok(None);
                }
                let attach = ordered[idx];
                attach_smiles_frag(&mut mol, attach, smi)?;
            }
            mol.perceive()?;
            return Ok(Some(mol));
        }
    }
    Ok(None)
}

fn parse_substituent_list(raw: &str) -> Vec<(usize, &'static str)> {
    let mut out = Vec::new();
    let s = raw
        .replace("ortho-", "2-")
        .replace("meta-", "3-")
        .replace("para-", "4-")
        .replace("o-", "2-")
        .replace("m-", "3-")
        .replace("p-", "4-");

    // Tokenize e.g. "2,4-dichloro-6-methyl" -> groups starting at locants
    let mut rem = s.as_str();
    while !rem.is_empty() {
        rem = rem.trim_start_matches('-');
        // Read optional locant list "2,4-"
        let mut locs: Vec<usize> = Vec::new();
        let mut idx = 0usize;
        let bytes = rem.as_bytes();
        while idx < bytes.len() && (bytes[idx].is_ascii_digit() || bytes[idx] == b',') {
            idx += 1;
        }
        if idx > 0 && idx < bytes.len() && bytes[idx] == b'-' {
            for tok in rem[..idx].split(',') {
                if let Ok(n) = tok.parse::<usize>() {
                    locs.push(n);
                }
            }
            rem = &rem[idx + 1..];
        } else if locs.is_empty() {
            locs.push(1);
        }

        // Strip Greek multiplier prefixes di/tri/tetra
        for mult in ["tetra", "tri", "di"] {
            if let Some(rest) = rem.strip_prefix(mult) {
                rem = rest;
                break;
            }
        }

        // Match longest known substituent prefix
        let known = [
            "trifluoromethyl",
            "dimethylamino",
            "tert-butyl",
            "isopropyl",
            "t-butyl",
            "sulfanyl",
            "mercapto",
            "methoxy",
            "hydroxy",
            "ethoxy",
            "phenyl",
            "propyl",
            "methyl",
            "fluoro",
            "chloro",
            "acetyl",
            "ethyl",
            "butyl",
            "bromo",
            "amino",
            "nitro",
            "cyano",
            "iodo",
        ];
        let mut matched = false;
        for &kw in &known {
            if let Some(rest) = rem.strip_prefix(kw)
                && let Some(smi) = subst_smiles(kw)
            {
                for &l in &locs {
                    out.push((l, smi));
                }
                rem = rest;
                matched = true;
                break;
            }
        }
        if !matched {
            return Vec::new();
        }
    }
    out
}

fn attach_smiles_frag(mol: &mut MolGraph, anchor: usize, frag_smi: &str) -> Result<(), MolError> {
    let frag = parse_smiles(frag_smi)?;
    let offset = mol.atoms.len();
    for fa in &frag.atoms {
        mol.add_atom(fa.clone());
    }
    for fb in &frag.bonds {
        mol.add_bond(offset + fb.u, offset + fb.v, fb.order, fb.stereo)?;
    }
    mol.add_bond(anchor, offset, BondOrder::Single, BondStereo::None)?;
    Ok(())
}

fn parse_iupac_aliphatic(norm: &str) -> Result<Option<MolGraph>, MolError> {
    // Match stem + unsaturation + functional suffix: e.g. "2-methylpropan-1-ol", "hexanoic acid", "cyclopentanone", "2-butene", "3-chloropentanenitrile"
    let stems: &[(&str, usize)] = &[
        ("dodeca", 12),
        ("dodec", 12),
        ("undeca", 11),
        ("undec", 11),
        ("deca", 10),
        ("dec", 10),
        ("nona", 9),
        ("non", 9),
        ("octa", 8),
        ("oct", 8),
        ("hepta", 7),
        ("hept", 7),
        ("hexa", 6),
        ("hex", 6),
        ("penta", 5),
        ("pent", 5),
        ("buta", 4),
        ("but", 4),
        ("propa", 3),
        ("prop", 3),
        ("etha", 2),
        ("eth", 2),
        ("metha", 1),
        ("meth", 1),
    ];

    // Check functional suffix
    let (body, suffix_kind, suffix_loc_opt) = extract_functional_suffix(norm);

    struct StemMatch<'a> {
        chain_len: usize,
        bond_mult: usize,
        is_cyclo: bool,
        unsat_loc: Option<usize>,
        prefix_str: &'a str,
    }
    let mut found: Option<StemMatch<'_>> = None;
    for &(stem_str, chain_len) in stems {
        if let Some(pos) = body.rfind(stem_str) {
            let after = &body[pos + stem_str.len()..];
            let (is_unsat_ord, unsat_loc) = parse_unsat_tail(after);
            if is_unsat_ord.is_none() {
                continue;
            }
            let before = &body[..pos];
            let (prefix_str, is_cyclo) = if let Some(b) = before.strip_suffix("cyclo") {
                (b, true)
            } else {
                (before, false)
            };
            found = Some(StemMatch {
                chain_len,
                bond_mult: is_unsat_ord.unwrap_or(1),
                is_cyclo,
                unsat_loc,
                prefix_str,
            });
            break;
        }
    }

    let Some(StemMatch {
        chain_len,
        bond_mult,
        is_cyclo,
        unsat_loc: unsat_loc_inline,
        prefix_str,
    }) = found
    else {
        return Ok(None);
    };
    if is_cyclo && chain_len < 3 {
        return Ok(None);
    }

    // Parse leading locant (e.g. "2-pentanol" -> leading loc = 2; "2-methylbutane" -> substituent)
    let mut leading_loc: Option<usize> = None;
    let mut subs: Vec<(usize, &'static str)> = Vec::new();
    let p_trimmed = prefix_str.trim().trim_end_matches('-');
    if !p_trimmed.is_empty() {
        if let Ok(n) = p_trimmed.parse::<usize>() {
            leading_loc = Some(n);
        } else {
            // Check if prefix ends with a bare locant for the main functional group/unsat, e.g. "2-methyl-1-"
            let (sub_part, tail_loc) = split_trailing_locant(p_trimmed);
            if tail_loc.is_some() {
                leading_loc = tail_loc;
            }
            if !sub_part.is_empty() {
                subs = parse_substituent_list(sub_part);
                if subs.is_empty() {
                    return Ok(None);
                }
            }
        }
    }

    // Build carbon backbone
    let mut mol = MolGraph::new();
    for _ in 0..chain_len {
        mol.add_atom(MolAtom::new(0, "C")?);
    }
    let unsat_pos = unsat_loc_inline
        .or(if bond_mult > 1 && suffix_kind == SuffixKind::None {
            leading_loc
        } else {
            None
        })
        .unwrap_or(1)
        .saturating_sub(1);
    for i in 0..chain_len.saturating_sub(1) {
        let ord = if bond_mult == 2 && i == unsat_pos {
            BondOrder::Double
        } else if bond_mult == 3 && i == unsat_pos {
            BondOrder::Triple
        } else {
            BondOrder::Single
        };
        mol.add_bond(i, i + 1, ord, BondStereo::None)?;
    }
    if is_cyclo {
        mol.add_bond(chain_len - 1, 0, BondOrder::Single, BondStereo::None)?;
    }

    // Attach functional suffix
    let func_pos = suffix_loc_opt
        .or(leading_loc)
        .unwrap_or(1)
        .saturating_sub(1);
    if func_pos >= chain_len {
        return Ok(None);
    }
    match suffix_kind {
        SuffixKind::None => {}
        SuffixKind::Ol => {
            let o = mol.add_atom(MolAtom::new(0, "O")?);
            mol.add_bond(func_pos, o, BondOrder::Single, BondStereo::None)?;
        }
        SuffixKind::Al => {
            let o = mol.add_atom(MolAtom::new(0, "O")?);
            mol.add_bond(0, o, BondOrder::Double, BondStereo::None)?;
        }
        SuffixKind::One => {
            let pos = if func_pos == 0 && chain_len >= 3 && !is_cyclo {
                1
            } else {
                func_pos
            };
            let o = mol.add_atom(MolAtom::new(0, "O")?);
            mol.add_bond(pos, o, BondOrder::Double, BondStereo::None)?;
        }
        SuffixKind::OicAcid => {
            let o1 = mol.add_atom(MolAtom::new(0, "O")?);
            let o2 = mol.add_atom(MolAtom::new(0, "O")?);
            mol.add_bond(0, o1, BondOrder::Double, BondStereo::None)?;
            mol.add_bond(0, o2, BondOrder::Single, BondStereo::None)?;
        }
        SuffixKind::Amine => {
            let n = mol.add_atom(MolAtom::new(0, "N")?);
            mol.add_bond(func_pos, n, BondOrder::Single, BondStereo::None)?;
        }
        SuffixKind::Amide => {
            let o = mol.add_atom(MolAtom::new(0, "O")?);
            let n = mol.add_atom(MolAtom::new(0, "N")?);
            mol.add_bond(0, o, BondOrder::Double, BondStereo::None)?;
            mol.add_bond(0, n, BondOrder::Single, BondStereo::None)?;
        }
        SuffixKind::Nitrile => {
            let n = mol.add_atom(MolAtom::new(0, "N")?);
            mol.add_bond(0, n, BondOrder::Triple, BondStereo::None)?;
        }
        SuffixKind::Thiol => {
            let s = mol.add_atom(MolAtom::new(0, "S")?);
            mol.add_bond(func_pos, s, BondOrder::Single, BondStereo::None)?;
        }
    }

    // Attach substituents
    for (loc, smi) in subs {
        let pos = loc.saturating_sub(1);
        if pos >= chain_len {
            return Ok(None);
        }
        attach_smiles_frag(&mut mol, pos, smi)?;
    }

    mol.perceive()?;
    Ok(Some(mol))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SuffixKind {
    None,
    Ol,
    Al,
    One,
    OicAcid,
    Amine,
    Amide,
    Nitrile,
    Thiol,
}

fn extract_functional_suffix(s: &str) -> (&str, SuffixKind, Option<usize>) {
    let pairs: &[(&str, SuffixKind)] = &[
        ("oic acid", SuffixKind::OicAcid),
        ("nitrile", SuffixKind::Nitrile),
        ("amine", SuffixKind::Amine),
        ("amide", SuffixKind::Amide),
        ("thiol", SuffixKind::Thiol),
        ("one", SuffixKind::One),
        ("ol", SuffixKind::Ol),
        ("al", SuffixKind::Al),
    ];
    for &(suf, kind) in pairs {
        if let Some(rest) = s.strip_suffix(suf) {
            let (before, loc) = split_trailing_locant(rest.trim_end_matches('-'));
            return (before, kind, loc);
        }
    }
    (s, SuffixKind::None, None)
}

fn split_trailing_locant(s: &str) -> (&str, Option<usize>) {
    if let Some(dash) = s.rfind('-') {
        let tail = &s[dash + 1..];
        if let Ok(n) = tail.parse::<usize>() {
            return (&s[..dash], Some(n));
        }
    }
    (s, None)
}

fn parse_unsat_tail(after: &str) -> (Option<usize>, Option<usize>) {
    let clean = after.trim_matches('-');
    // Handles "", "an", "ane", "en", "ene", "yn", "yne", "2-en", "2-ene", "1-yne"
    match clean {
        "" | "n" | "an" | "ane" | "e" => (Some(1), None),
        "en" | "ene" => (Some(2), None),
        "yn" | "yne" => (Some(3), None),
        other => {
            if let Some(dash) = other.find('-') {
                let loc = other[..dash].parse::<usize>().ok();
                let tail = &other[dash + 1..];
                let mult = match tail {
                    "en" | "ene" => Some(2),
                    "yn" | "yne" => Some(3),
                    "an" | "ane" => Some(1),
                    _ => None,
                };
                (mult, loc)
            } else {
                (None, None)
            }
        }
    }
}
