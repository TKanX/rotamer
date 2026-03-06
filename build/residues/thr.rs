use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5290,
        bond_angle: 109.4115,
        torsion: TorsionSrc::Fixed(120.0002),
    },
    AtomSpec {
        name: "OG1",
        refs: ["N", "CA", "CB"],
        bond_length: 1.4280,
        bond_angle: 109.5053,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CG2",
        refs: ["OG1", "CA", "CB"],
        bond_length: 1.5301,
        bond_angle: 109.5255,
        torsion: TorsionSrc::Fixed(120.0308),
    },
    AtomSpec {
        name: "HB",
        refs: ["OG1", "CA", "CB"],
        bond_length: 1.0898,
        bond_angle: 109.4276,
        torsion: TorsionSrc::Fixed(-119.9922),
    },
    AtomSpec {
        name: "HG1",
        refs: ["CA", "CB", "OG1"],
        bond_length: 0.9666,
        bond_angle: 106.8126,
        torsion: TorsionSrc::PolarH(0, 0.0),
    },
    AtomSpec {
        name: "HG21",
        refs: ["OG1", "CB", "CG2"],
        bond_length: 1.0893,
        bond_angle: 109.4779,
        torsion: TorsionSrc::Fixed(-59.9996),
    },
    AtomSpec {
        name: "HG22",
        refs: ["OG1", "CB", "CG2"],
        bond_length: 1.0901,
        bond_angle: 109.4679,
        torsion: TorsionSrc::Fixed(-179.9932),
    },
    AtomSpec {
        name: "HG23",
        refs: ["OG1", "CB", "CG2"],
        bond_length: 1.0888,
        bond_angle: 109.4591,
        torsion: TorsionSrc::Fixed(60.0221),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "THR",
    type_name: "Thr",
    n_chi: 1,
    n_polar_h: 1,
    atoms: ATOMS,
};
