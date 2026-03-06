use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5294,
        bond_angle: 109.4274,
        torsion: TorsionSrc::Fixed(120.0367),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5284,
        bond_angle: 109.5449,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "SD",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.8137,
        bond_angle: 109.5064,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "CE",
        refs: ["CB", "CG", "SD"],
        bond_length: 1.8135,
        bond_angle: 100.0339,
        torsion: TorsionSrc::Chi(2),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0901,
        bond_angle: 109.4389,
        torsion: TorsionSrc::Fixed(120.0398),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0905,
        bond_angle: 109.4620,
        torsion: TorsionSrc::Fixed(-120.0543),
    },
    AtomSpec {
        name: "HG2",
        refs: ["SD", "CB", "CG"],
        bond_length: 1.0903,
        bond_angle: 109.4601,
        torsion: TorsionSrc::Fixed(119.9798),
    },
    AtomSpec {
        name: "HG3",
        refs: ["SD", "CB", "CG"],
        bond_length: 1.0895,
        bond_angle: 109.4834,
        torsion: TorsionSrc::Fixed(-120.0199),
    },
    AtomSpec {
        name: "HE1",
        refs: ["CG", "SD", "CE"],
        bond_length: 1.0891,
        bond_angle: 109.5482,
        torsion: TorsionSrc::Fixed(-179.9723),
    },
    AtomSpec {
        name: "HE2",
        refs: ["CG", "SD", "CE"],
        bond_length: 1.0893,
        bond_angle: 109.4602,
        torsion: TorsionSrc::Fixed(60.0028),
    },
    AtomSpec {
        name: "HE3",
        refs: ["CG", "SD", "CE"],
        bond_length: 1.0899,
        bond_angle: 109.4479,
        torsion: TorsionSrc::Fixed(-60.0303),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "MET",
    type_name: "Met",
    n_chi: 3,
    n_polar_h: 0,
    atoms: ATOMS,
};
