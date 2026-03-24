use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5288,
        bond_angle: 109.4301,
        torsion: TorsionSrc::Fixed(-120.0733),
    },
    AtomSpec {
        name: "CG1",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5294,
        bond_angle: 109.5474,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CG2",
        refs: ["CG1", "CA", "CB"],
        bond_length: 1.5303,
        bond_angle: 109.4577,
        torsion: TorsionSrc::Fixed(-119.9719),
    },
    AtomSpec {
        name: "CD1",
        refs: ["CA", "CB", "CG1"],
        bond_length: 1.5288,
        bond_angle: 109.5474,
        torsion: TorsionSrc::Chi(1),
    },
    AtomSpec {
        name: "HB",
        refs: ["CG1", "CA", "CB"],
        bond_length: 1.0893,
        bond_angle: 109.4848,
        torsion: TorsionSrc::Fixed(120.0904),
    },
    AtomSpec {
        name: "HG12",
        refs: ["CD1", "CB", "CG1"],
        bond_length: 1.0898,
        bond_angle: 109.4338,
        torsion: TorsionSrc::Fixed(-119.9715),
    },
    AtomSpec {
        name: "HG13",
        refs: ["CD1", "CB", "CG1"],
        bond_length: 1.0895,
        bond_angle: 109.4756,
        torsion: TorsionSrc::Fixed(120.0331),
    },
    AtomSpec {
        name: "HG21",
        refs: ["CA", "CB", "CG2"],
        bond_length: 1.0885,
        bond_angle: 109.4850,
        torsion: TorsionSrc::Fixed(59.9726),
    },
    AtomSpec {
        name: "HG22",
        refs: ["CA", "CB", "CG2"],
        bond_length: 1.0900,
        bond_angle: 109.4497,
        torsion: TorsionSrc::Fixed(-179.9483),
    },
    AtomSpec {
        name: "HG23",
        refs: ["CA", "CB", "CG2"],
        bond_length: 1.0898,
        bond_angle: 109.4798,
        torsion: TorsionSrc::Fixed(-60.0297),
    },
    AtomSpec {
        name: "HD11",
        refs: ["CB", "CG1", "CD1"],
        bond_length: 1.0891,
        bond_angle: 109.6037,
        torsion: TorsionSrc::Fixed(-179.9919),
    },
    AtomSpec {
        name: "HD12",
        refs: ["CB", "CG1", "CD1"],
        bond_length: 1.0903,
        bond_angle: 109.4936,
        torsion: TorsionSrc::Fixed(-59.9014),
    },
    AtomSpec {
        name: "HD13",
        refs: ["CB", "CG1", "CD1"],
        bond_length: 1.0897,
        bond_angle: 109.4579,
        torsion: TorsionSrc::Fixed(59.9796),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "ILE",
    type_name: "Ile",
    n_chi: 2,
    n_polar_h: 0,
    atoms: ATOMS,
};
