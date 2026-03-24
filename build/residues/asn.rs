use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5309,
        bond_angle: 109.4540,
        torsion: TorsionSrc::Fixed(-119.9959),
    },
    AtomSpec {
        name: "CG",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5066,
        bond_angle: 109.4843,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "OD1",
        refs: ["CA", "CB", "CG"],
        bond_length: 1.2133,
        bond_angle: 119.9743,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "ND2",
        refs: ["OD1", "CB", "CG"],
        bond_length: 1.3476,
        bond_angle: 120.0120,
        torsion: TorsionSrc::Fixed(-179.9262),
    },
    AtomSpec {
        name: "HB2",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0903,
        bond_angle: 109.4142,
        torsion: TorsionSrc::Fixed(119.9676),
    },
    AtomSpec {
        name: "HB3",
        refs: ["CG", "CA", "CB"],
        bond_length: 1.0895,
        bond_angle: 109.4697,
        torsion: TorsionSrc::Fixed(-120.0326),
    },
    AtomSpec {
        name: "HD21",
        refs: ["OD1", "CG", "ND2"],
        bond_length: 0.9701,
        bond_angle: 119.9863,
        torsion: TorsionSrc::Fixed(-179.9777),
    },
    AtomSpec {
        name: "HD22",
        refs: ["OD1", "CG", "ND2"],
        bond_length: 0.9696,
        bond_angle: 120.0576,
        torsion: TorsionSrc::Fixed(0.0640),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "ASN",
    type_name: "Asn",
    n_chi: 2,
    n_polar_h: 0,
    atoms: ATOMS,
};
