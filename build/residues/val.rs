use super::{AtomSpec, ResidueSpec, TorsionSrc};

const ATOMS: &[AtomSpec] = &[
    AtomSpec {
        name: "CB",
        refs: ["C", "N", "CA"],
        bond_length: 1.5287,
        bond_angle: 109.4452,
        torsion: TorsionSrc::Fixed(120.0037),
    },
    AtomSpec {
        name: "CG1",
        refs: ["N", "CA", "CB"],
        bond_length: 1.5299,
        bond_angle: 109.5086,
        torsion: TorsionSrc::Chi(0),
    },
    AtomSpec {
        name: "CG2",
        refs: ["CG1", "CA", "CB"],
        bond_length: 1.5292,
        bond_angle: 109.4895,
        torsion: TorsionSrc::Fixed(-120.0257),
    },
    AtomSpec {
        name: "HB",
        refs: ["CG1", "CA", "CB"],
        bond_length: 1.0909,
        bond_angle: 109.4560,
        torsion: TorsionSrc::Fixed(119.9712),
    },
    AtomSpec {
        name: "HG11",
        refs: ["CA", "CB", "CG1"],
        bond_length: 1.0905,
        bond_angle: 109.5198,
        torsion: TorsionSrc::Fixed(179.9893),
    },
    AtomSpec {
        name: "HG12",
        refs: ["CA", "CB", "CG1"],
        bond_length: 1.0896,
        bond_angle: 109.4569,
        torsion: TorsionSrc::Fixed(59.9808),
    },
    AtomSpec {
        name: "HG13",
        refs: ["CA", "CB", "CG1"],
        bond_length: 1.0889,
        bond_angle: 109.4590,
        torsion: TorsionSrc::Fixed(-60.0332),
    },
    AtomSpec {
        name: "HG21",
        refs: ["CA", "CB", "CG2"],
        bond_length: 1.0899,
        bond_angle: 109.4732,
        torsion: TorsionSrc::Fixed(60.0212),
    },
    AtomSpec {
        name: "HG22",
        refs: ["CA", "CB", "CG2"],
        bond_length: 1.0897,
        bond_angle: 109.4964,
        torsion: TorsionSrc::Fixed(-59.9747),
    },
    AtomSpec {
        name: "HG23",
        refs: ["CA", "CB", "CG2"],
        bond_length: 1.0903,
        bond_angle: 109.5304,
        torsion: TorsionSrc::Fixed(-179.9528),
    },
];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "VAL",
    type_name: "Val",
    n_chi: 1,
    n_polar_h: 0,
    atoms: ATOMS,
};
