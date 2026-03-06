/// Torsion source for an atom's dihedral angle.
#[derive(Debug, Clone, Copy)]
pub enum TorsionSrc {
    /// Fixed dihedral (degrees).
    Fixed(f64),
    /// Equals chi[i] (runtime input).
    Chi(usize),
    /// Equals polar_h[i] + offset (degrees).
    PolarH(usize, f64),
}

/// Internal-coordinate specification for one atom.
#[derive(Debug, Clone, Copy)]
pub struct AtomSpec {
    pub name: &'static str,
    pub refs: [&'static str; 3],
    pub bond_length: f64,
    pub bond_angle: f64,
    pub torsion: TorsionSrc,
}

/// Complete specification of one residue type.
#[derive(Debug, Clone)]
pub struct ResidueSpec {
    pub name: &'static str,
    pub type_name: &'static str,
    pub n_chi: usize,
    pub n_polar_h: usize,
    pub atoms: &'static [AtomSpec],
}

mod ala;
mod arg;
mod arn;
mod ash;
mod asn;
mod asp;
mod cym;
mod cys;
mod cyx;
mod glh;
mod gln;
mod glu;
mod gly;
mod hid;
mod hie;
mod hip;
mod ile;
mod leu;
mod lyn;
mod lys;
mod met;
mod phe;
mod pro;
mod ser;
mod thr;
mod trp;
mod tym;
mod tyr;
mod val;

pub const ALL_RESIDUES: &[ResidueSpec] = &[
    ala::SPEC,
    arg::SPEC,
    arn::SPEC,
    ash::SPEC,
    asn::SPEC,
    asp::SPEC,
    cym::SPEC,
    cys::SPEC,
    cyx::SPEC,
    glh::SPEC,
    gln::SPEC,
    glu::SPEC,
    gly::SPEC,
    hid::SPEC,
    hie::SPEC,
    hip::SPEC,
    ile::SPEC,
    leu::SPEC,
    lyn::SPEC,
    lys::SPEC,
    met::SPEC,
    phe::SPEC,
    pro::SPEC,
    ser::SPEC,
    thr::SPEC,
    trp::SPEC,
    tym::SPEC,
    tyr::SPEC,
    val::SPEC,
];
