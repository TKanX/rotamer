use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5285,
        bond_angle: 109.4957,
        torsion: TorsionSrc::Fixed(-120.0138),
    },
    AtomSpec {
        name: "SG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.8141,
        bond_angle: 109.4981,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "HB2",
        refs: ["SG", "CA", "CB"],
        bond_length: 1.0901,
        bond_angle: 109.4724,
        torsion: TorsionSrc::Fixed(119.9825),
    },
    AtomSpec {
        name: "HB3",
        refs: ["SG", "CA", "CB"],
        bond_length: 1.0896,
        bond_angle: 109.4416,
        torsion: TorsionSrc::Fixed(-120.0181),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "CYM",
    type_name: "Cym",
    n_chi: 1,
    n_polar_h: 0,
    atoms: ATOMS,
};
