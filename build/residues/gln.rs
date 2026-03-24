use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5288,
        bond_angle: 109.4592,
        torsion: TorsionSrc::Fixed(-120.0663),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5284,
        bond_angle: 109.5344,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CD",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.5066,
        bond_angle: 109.5426,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "OE1",
        refs: ["CB", "CG", "CD"],
        bond_length: 1.2122,
        bond_angle: 119.9367,
        torsion: TorsionSrc::Chi(2),
    },
    AtomSpec {
        name: "NE2",
        refs: ["OE1", "CG", "CD"],
        bond_length: 1.3471,
        bond_angle: 120.0933,
        torsion: TorsionSrc::Fixed(-179.9569),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0906,
        bond_angle: 109.4211,
        torsion: TorsionSrc::Fixed(119.9936),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0896,
        bond_angle: 109.4409,
        torsion: TorsionSrc::Fixed(-120.1141),
    },
    AtomSpec {
        name: "HG2",
        refs: ["CD", "CB", "CG"],
        bond_length: 1.0906,
        bond_angle: 109.4595,
        torsion: TorsionSrc::Fixed(119.9623),
    },
    AtomSpec {
        name: "HG3",
        refs: ["CD", "CB", "CG"],
        bond_length: 1.0896,
        bond_angle: 109.5022,
        torsion: TorsionSrc::Fixed(-120.0840),
    },
    AtomSpec {
        name: "HE21",
        refs: ["OE1", "CD", "NE2"],
        bond_length: 0.9694,
        bond_angle: 120.1198,
        torsion: TorsionSrc::Fixed(0.0314),
    },
    AtomSpec {
        name: "HE22",
        refs: ["OE1", "CD", "NE2"],
        bond_length: 0.9704,
        bond_angle: 119.9555,
        torsion: TorsionSrc::Fixed(-179.9603),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "GLN",
    type_name: "Gln",
    n_chi: 3,
    n_polar_h: 0,
    atoms: ATOMS,
};
