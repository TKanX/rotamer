use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5286,
        bond_angle: 109.4163,
        torsion: TorsionSrc::Fixed(119.9728),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5303,
        bond_angle: 109.4950,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CD1",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.5300,
        bond_angle: 109.5000,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "CD2",
        refs: ["CD1", "CB", "CG"],
        bond_length: 1.5285,
        bond_angle: 109.5008,
        torsion: TorsionSrc::Fixed(-120.0918),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0898,
        bond_angle: 109.4621,
        torsion: TorsionSrc::Fixed(119.9648),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0902,
        bond_angle: 109.5405,
        torsion: TorsionSrc::Fixed(-120.0807),
    },
    AtomSpec {
        name: "HG",
        refs: ["CD1", "CB", "CG"],
        bond_length: 1.0896,
        bond_angle: 109.4031,
        torsion: TorsionSrc::Fixed(119.8613),
    },
    AtomSpec {
        name: "HD11",
        refs: ["CD2", "CG", "CD1"],
        bond_length: 1.0892,
        bond_angle: 109.4929,
        torsion: TorsionSrc::Fixed(-59.8967),
    },
    AtomSpec {
        name: "HD12",
        refs: ["CD2", "CG", "CD1"],
        bond_length: 1.0903,
        bond_angle: 109.5083,
        torsion: TorsionSrc::Fixed(179.9599),
    },
    AtomSpec {
        name: "HD13",
        refs: ["CD2", "CG", "CD1"],
        bond_length: 1.0898,
        bond_angle: 109.4476,
        torsion: TorsionSrc::Fixed(60.0409),
    },
    AtomSpec {
        name: "HD21",
        refs: ["CD1", "CG", "CD2"],
        bond_length: 1.0898,
        bond_angle: 109.4268,
        torsion: TorsionSrc::Fixed(-60.0282),
    },
    AtomSpec {
        name: "HD22",
        refs: ["CD1", "CG", "CD2"],
        bond_length: 1.0899,
        bond_angle: 109.4971,
        torsion: TorsionSrc::Fixed(179.9761),
    },
    AtomSpec {
        name: "HD23",
        refs: ["CD1", "CG", "CD2"],
        bond_length: 1.0903,
        bond_angle: 109.4847,
        torsion: TorsionSrc::Fixed(59.9357),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "LEU",
    type_name: "Leu",
    n_chi: 2,
    n_polar_h: 0,
    atoms: ATOMS,
};
