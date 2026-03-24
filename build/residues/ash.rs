use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5301,
        bond_angle: 109.4797,
        torsion: TorsionSrc::Fixed(-120.0113),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5075,
        bond_angle: 109.4630,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "OD1",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.2080,
        bond_angle: 119.9590,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "OD2",
        refs: ["OD1", "CB", "CG"],
        bond_length: 1.3415,
        bond_angle: 119.9993,
        torsion: TorsionSrc::Fixed(-179.9359),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0898,
        bond_angle: 109.4853,
        torsion: TorsionSrc::Fixed(119.9593),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0901,
        bond_angle: 109.4853,
        torsion: TorsionSrc::Fixed(-120.0111),
    },
    AtomSpec {
        name: "HD2",
        refs: ["OD1", "CG", "OD2"],
        bond_length: 0.9664,
        bond_angle: 117.0276,
        torsion: TorsionSrc::PolarH(0, 0.0),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "ASH",
    type_name: "Ash",
    n_chi: 2,
    n_polar_h: 1,
    atoms: ATOMS,
};
