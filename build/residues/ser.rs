use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5287,
        bond_angle: 109.4711,
        torsion: TorsionSrc::Fixed(120.0204),
    },
    AtomSpec {
        name: "OG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.4283,
        bond_angle: 109.5115,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "HB2",
        refs: ["OG", "CA", "CB"],
        bond_length: 1.0902,
        bond_angle: 109.4341,
        torsion: TorsionSrc::Fixed(-120.0069),
    },
    AtomSpec {
        name: "HB3",
        refs: ["OG", "CA", "CB"],
        bond_length: 1.0896,
        bond_angle: 109.4805,
        torsion: TorsionSrc::Fixed(120.0365),
    },
    AtomSpec {
        name: "HG",
        refs: ["CA", "CB", "OG"],
        bond_length: 0.9670,
        bond_angle: 106.8145,
        torsion: TorsionSrc::PolarH(0, 0.0),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "SER",
    type_name: "Ser",
    n_chi: 1,
    n_polar_h: 1,
    atoms: ATOMS,
};
