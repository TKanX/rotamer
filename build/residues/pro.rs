use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5434,
        bond_angle: 104.7219,
        torsion: TorsionSrc::Fixed(-118.8407),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5426,
        bond_angle: 105.0594,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CD",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.5437,
        bond_angle: 105.0631,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0903,
        bond_angle: 110.3550,
        torsion: TorsionSrc::Fixed(118.8174),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0895,
        bond_angle: 110.3543,
        torsion: TorsionSrc::Fixed(-118.8892),
    },
    AtomSpec {
        name: "HG2",
        refs: ["CD", "CB", "CG"],
        bond_length: 1.0900,
        bond_angle: 110.3647,
        torsion: TorsionSrc::Fixed(118.8354),
    },
    AtomSpec {
        name: "HG3",
        refs: ["CD", "CB", "CG"],
        bond_length: 1.0904,
        bond_angle: 110.3653,
        torsion: TorsionSrc::Fixed(-118.8865),
    },
    AtomSpec {
        name: "HD2",
        refs: ["N", "CG", "CD"],
        bond_length: 1.0903,
        bond_angle: 110.4602,
        torsion: TorsionSrc::Fixed(118.9023),
    },
    AtomSpec {
        name: "HD3",
        refs: ["N", "CG", "CD"],
        bond_length: 1.0899,
        bond_angle: 110.4192,
        torsion: TorsionSrc::Fixed(-118.6596),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "PRO",
    type_name: "Pro",
    n_chi: 2,
    n_polar_h: 0,
    atoms: ATOMS,
};
