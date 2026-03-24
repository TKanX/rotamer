use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5302,
        bond_angle: 109.4824,
        torsion: TorsionSrc::Fixed(-119.9555),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5306,
        bond_angle: 109.4016,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CD",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.5076,
        bond_angle: 109.4303,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "OE1",
        refs: ["CB", "CG", "CD"],
        bond_length: 1.2084,
        bond_angle: 120.0030,
        torsion: TorsionSrc::Chi(2),
    },
    AtomSpec {
        name: "OE2",
        refs: ["OE1", "CG", "CD"],
        bond_length: 1.3425,
        bond_angle: 119.9977,
        torsion: TorsionSrc::Fixed(-179.9373),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0901,
        bond_angle: 109.5095,
        torsion: TorsionSrc::Fixed(119.9635),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0891,
        bond_angle: 109.5036,
        torsion: TorsionSrc::Fixed(-119.9773),
    },
    AtomSpec {
        name: "HG2",
        refs: ["CD", "CB", "CG"],
        bond_length: 1.0901,
        bond_angle: 109.5026,
        torsion: TorsionSrc::Fixed(119.9653),
    },
    AtomSpec {
        name: "HG3",
        refs: ["CD", "CB", "CG"],
        bond_length: 1.0894,
        bond_angle: 109.4463,
        torsion: TorsionSrc::Fixed(-119.9539),
    },
    AtomSpec {
        name: "HE2",
        refs: ["CG", "CD", "OE2"],
        bond_length: 0.9663,
        bond_angle: 116.9918,
        torsion: TorsionSrc::PolarH(0, 0.0),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "GLH",
    type_name: "Glh",
    n_chi: 3,
    n_polar_h: 1,
    atoms: ATOMS,
};
