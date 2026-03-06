use super::{AtomSpec, ResidueSpec};

const ATOMS: &[AtomSpec] = &[];

pub const SPEC: ResidueSpec = ResidueSpec {
    name: "GLY",
    type_name: "Gly",
    n_chi: 0,
    n_polar_h: 0,
    atoms: ATOMS,
};
