use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5294,
        bond_angle: 109.4645,
        torsion: TorsionSrc::Fixed(-119.9991),
    },
    AtomSpec {
        name: "HB1",
        refs: ["N", "CA", "CB"],
        bond_length: 1.0907,
        bond_angle: 109.4905,
        torsion: TorsionSrc::Fixed(-59.9713),
    },
    AtomSpec {
        name: "HB2",
        refs: ["N", "CA", "CB"],
        bond_length: 1.0899,
        bond_angle: 109.4323,
        torsion: TorsionSrc::Fixed(59.9976),
    },
    AtomSpec {
        name: "HB3",
        refs: ["N", "CA", "CB"],
        bond_length: 1.0901,
        bond_angle: 109.5239,
        torsion: TorsionSrc::Fixed(-179.9612),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "ALA",
    type_name: "Ala",
    n_chi: 0,
    n_polar_h: 0,
    atoms: ATOMS,
};
