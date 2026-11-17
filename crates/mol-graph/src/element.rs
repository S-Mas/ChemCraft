/// Chemical element metadata used for valence checks, formula, and mass.
#[derive(Debug, Clone, Copy)]
pub struct ElementMeta {
    pub atomic_number: u8,
    pub symbol: &'static str,
    pub atomic_weight: f64,
    pub electronegativity: f64,
    pub default_valences: &'static [u8],
}

static ELEMENTS: &[ElementMeta] = &[
    ElementMeta {
        atomic_number: 0,
        symbol: "*",
        atomic_weight: 0.0,
        electronegativity: 0.0,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 1,
        symbol: "H",
        atomic_weight: 1.008,
        electronegativity: 2.20,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 3,
        symbol: "Li",
        atomic_weight: 6.94,
        electronegativity: 0.98,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 5,
        symbol: "B",
        atomic_weight: 10.81,
        electronegativity: 2.04,
        default_valences: &[3],
    },
    ElementMeta {
        atomic_number: 6,
        symbol: "C",
        atomic_weight: 12.011,
        electronegativity: 2.55,
        default_valences: &[4],
    },
    ElementMeta {
        atomic_number: 7,
        symbol: "N",
        atomic_weight: 14.007,
        electronegativity: 3.04,
        default_valences: &[3, 5],
    },
    ElementMeta {
        atomic_number: 8,
        symbol: "O",
        atomic_weight: 15.999,
        electronegativity: 3.44,
        default_valences: &[2],
    },
    ElementMeta {
        atomic_number: 9,
        symbol: "F",
        atomic_weight: 18.998,
        electronegativity: 3.98,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 11,
        symbol: "Na",
        atomic_weight: 22.990,
        electronegativity: 0.93,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 12,
        symbol: "Mg",
        atomic_weight: 24.305,
        electronegativity: 1.31,
        default_valences: &[2],
    },
    ElementMeta {
        atomic_number: 13,
        symbol: "Al",
        atomic_weight: 26.982,
        electronegativity: 1.61,
        default_valences: &[3],
    },
    ElementMeta {
        atomic_number: 14,
        symbol: "Si",
        atomic_weight: 28.085,
        electronegativity: 1.90,
        default_valences: &[4],
    },
    ElementMeta {
        atomic_number: 15,
        symbol: "P",
        atomic_weight: 30.974,
        electronegativity: 2.19,
        default_valences: &[3, 5],
    },
    ElementMeta {
        atomic_number: 16,
        symbol: "S",
        atomic_weight: 32.06,
        electronegativity: 2.58,
        default_valences: &[2, 4, 6],
    },
    ElementMeta {
        atomic_number: 17,
        symbol: "Cl",
        atomic_weight: 35.45,
        electronegativity: 3.16,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 19,
        symbol: "K",
        atomic_weight: 39.098,
        electronegativity: 0.82,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 20,
        symbol: "Ca",
        atomic_weight: 40.078,
        electronegativity: 1.00,
        default_valences: &[2],
    },
    ElementMeta {
        atomic_number: 25,
        symbol: "Mn",
        atomic_weight: 54.938,
        electronegativity: 1.55,
        default_valences: &[2, 4, 7],
    },
    ElementMeta {
        atomic_number: 26,
        symbol: "Fe",
        atomic_weight: 55.845,
        electronegativity: 1.83,
        default_valences: &[2, 3, 10],
    },
    ElementMeta {
        atomic_number: 29,
        symbol: "Cu",
        atomic_weight: 63.546,
        electronegativity: 1.90,
        default_valences: &[1, 2],
    },
    ElementMeta {
        atomic_number: 30,
        symbol: "Zn",
        atomic_weight: 65.38,
        electronegativity: 1.65,
        default_valences: &[2],
    },
    ElementMeta {
        atomic_number: 34,
        symbol: "Se",
        atomic_weight: 78.971,
        electronegativity: 2.55,
        default_valences: &[2, 4, 6],
    },
    ElementMeta {
        atomic_number: 35,
        symbol: "Br",
        atomic_weight: 79.904,
        electronegativity: 2.96,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 44,
        symbol: "Ru",
        atomic_weight: 101.07,
        electronegativity: 2.20,
        default_valences: &[2, 3, 4, 8],
    },
    ElementMeta {
        atomic_number: 45,
        symbol: "Rh",
        atomic_weight: 102.91,
        electronegativity: 2.28,
        default_valences: &[1, 3],
    },
    ElementMeta {
        atomic_number: 46,
        symbol: "Pd",
        atomic_weight: 106.42,
        electronegativity: 2.20,
        default_valences: &[2, 4],
    },
    ElementMeta {
        atomic_number: 47,
        symbol: "Ag",
        atomic_weight: 107.868,
        electronegativity: 1.93,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 53,
        symbol: "I",
        atomic_weight: 126.904,
        electronegativity: 2.66,
        default_valences: &[1],
    },
    ElementMeta {
        atomic_number: 56,
        symbol: "Ba",
        atomic_weight: 137.327,
        electronegativity: 0.89,
        default_valences: &[2],
    },
    ElementMeta {
        atomic_number: 78,
        symbol: "Pt",
        atomic_weight: 195.084,
        electronegativity: 2.28,
        default_valences: &[2, 4],
    },
];

pub fn element_by_symbol(sym: &str) -> Option<&'static ElementMeta> {
    ELEMENTS.iter().find(|e| e.symbol.eq_ignore_ascii_case(sym))
}

pub fn element_by_number(z: u8) -> Option<&'static ElementMeta> {
    ELEMENTS.iter().find(|e| e.atomic_number == z)
}

pub fn max_neutral_valence(symbol: &str, formal_charge: i8) -> u8 {
    let base = match symbol {
        "*" | "H" | "F" | "Cl" | "Br" | "I" | "Li" | "Na" | "K" => 1,
        "O" => 2,
        "N" => 3,
        "B" => 3,
        "C" | "Si" => 4,
        "P" => 5,
        "S" | "Se" => 6,
        _ => 12,
    };
    match (symbol, formal_charge) {
        ("N", c) if c > 0 => 4,
        ("N", c) if c < 0 => 2,
        ("O", c) if c > 0 => 3,
        ("O", c) if c < 0 => 1,
        ("S", c) if c < 0 => 1,
        ("B", c) if c < 0 => 4,
        ("C", c) if c < 0 => 3,
        ("F" | "Cl" | "Br" | "I", c) if c < 0 => 0,
        _ => base,
    }
}
