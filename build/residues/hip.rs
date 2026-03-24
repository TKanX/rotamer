use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5337,
        bond_angle: 111.1252,
        torsion: TorsionSrc::Fixed(-122.7794),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5100,
        bond_angle: 112.9791,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "ND1",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.3513,
        bond_angle: 120.3285,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "CD2",
        refs: ["ND1", "CB", "CG"],
        bond_length: 1.3376,
        bond_angle: 129.9283,
        torsion: TorsionSrc::Fixed(179.8457),
    },
    AtomSpec {
        name: "CE1",
        refs: ["CB", "CG", "ND1"],
        bond_length: 1.3369,
        bond_angle: 107.8621,
        torsion: TorsionSrc::Fixed(179.9049),
    },
    AtomSpec {
        name: "NE2",
        refs: ["CB", "CG", "CD2"],
        bond_length: 1.3739,
        bond_angle: 105.3317,
        torsion: TorsionSrc::Fixed(-179.8641),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0987,
        bond_angle: 110.3817,
        torsion: TorsionSrc::Fixed(121.1196),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0982,
        bond_angle: 110.1994,
        torsion: TorsionSrc::Fixed(-123.0905),
    },
    AtomSpec {
        name: "HD1",
        refs: ["CE1", "CG", "ND1"],
        bond_length: 1.0163,
        bond_angle: 127.0804,
        torsion: TorsionSrc::Fixed(-179.9889),
    },
    AtomSpec {
        name: "HD2",
        refs: ["NE2", "CG", "CD2"],
        bond_length: 1.0723,
        bond_angle: 137.1496,
        torsion: TorsionSrc::Fixed(179.8416),
    },
    AtomSpec {
        name: "HE1",
        refs: ["NE2", "ND1", "CE1"],
        bond_length: 1.0777,
        bond_angle: 126.1756,
        torsion: TorsionSrc::Fixed(-179.9653),
    },
    AtomSpec {
        name: "HE2",
        refs: ["CE1", "CD2", "NE2"],
        bond_length: 1.0156,
        bond_angle: 125.4757,
        torsion: TorsionSrc::Fixed(179.9562),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "HIP",
    type_name: "Hip",
    n_chi: 2,
    n_polar_h: 0,
    atoms: ATOMS,
};
