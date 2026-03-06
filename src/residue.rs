use crate::math::Vec3;
use crate::sealed::Sealed;

/// Provides a flat slice view of sidechain atom coordinates.
///
/// Each concrete `Coords` struct is `#[repr(C)]` with `N` fields of type
/// [`Vec3`], so it is memory-equivalent to `[Vec3; N]`. The [`as_slice`](Self::as_slice)
/// method performs a zero-copy reinterpretation.
pub trait SidechainCoords: Copy + 'static {
    /// Number of non-backbone atoms (heavy + hydrogen) in this sidechain.
    const N: usize;

    /// Zero-copy view of the coordinates as a contiguous slice.
    ///
    /// The returned slice has exactly [`N`](Self::N) elements, in the same
    /// order as the named fields of the concrete `Coords` struct.
    fn as_slice(&self) -> &[Vec3];
}

/// Sidechain geometry builder.
///
/// Given the three backbone anchor coordinates (N, Cα, C) and a set of χ
/// dihedral angles plus optional polar-hydrogen rotation angles, the
/// inherent `build` method on each concrete type returns the full set of
/// non-backbone atom coordinates using the NERF algorithm.
///
/// The trait is `sealed` — only the 29 residue types
/// defined in this crate implement it.
pub trait Sidechain: Sealed + Copy + 'static {
    /// Number of χ dihedral angles from the rotamer library.
    const N_CHI: usize;

    /// Number of independently rotatable polar-hydrogen torsions.
    ///
    /// Only hydrogens on heteroatoms (O, N, S) that freely rotate around
    /// a single bond are counted. Hydrogens whose positions are fixed by
    /// sp² or ring geometry (e.g. backbone N–H, imidazole N–H, amide NH₂)
    /// are placed with `Fixed` torsions and are **not** included here.
    const N_POLAR_H: usize;

    /// Three-letter residue name (uppercase ASCII, e.g. `"SER"`, `"ARG"`).
    const NAME: &'static str;

    /// Concrete coordinate type for this residue's sidechain.
    type Coords: SidechainCoords;
}

macro_rules! residue_type {
    ($(#[$meta:meta])* $t:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $t;
        impl Sealed for $t {}
    };
}

residue_type!(
    /// Alanine — 0 χ, 0 polar H.
    Ala
);
residue_type!(
    /// Arginine (protonated, +1) — 4 χ, 0 polar H.
    ///
    /// Guanidinium NH groups are coplanar (sp²); no rotatable polar H.
    Arg
);
residue_type!(
    /// Arginine (neutral) — 4 χ, 0 polar H.
    ///
    /// Guanidine NH groups remain coplanar (sp²); no rotatable polar H.
    Arn
);
residue_type!(
    /// Aspartate (protonated, neutral) — 2 χ, 1 polar H.
    ///
    /// OD2–HD2 hydroxyl rotation.
    Ash
);
residue_type!(
    /// Asparagine — 2 χ, 0 polar H.
    ///
    /// Amide NH₂ is coplanar with the CG=OD1 plane.
    Asn
);
residue_type!(
    /// Aspartate (deprotonated, −1) — 2 χ, 0 polar H.
    ///
    /// Both carboxylate oxygens (OD1, OD2) are deprotonated; no O–H hydrogen.
    Asp
);
residue_type!(
    /// Cysteine (deprotonated thiolate, −1) — 1 χ, 0 polar H.
    ///
    /// SG bears no hydrogen (thiolate S⁻).
    Cym
);
residue_type!(
    /// Cysteine (free thiol) — 1 χ, 1 polar H.
    ///
    /// SG–HG thiol rotation.
    Cys
);
residue_type!(
    /// Cystine (disulfide-bonded) — 1 χ, 0 polar H.
    ///
    /// SG is bonded to the partner residue's SG; bears no hydrogen.
    Cyx
);
residue_type!(
    /// Glutamate (protonated, neutral) — 3 χ, 1 polar H.
    ///
    /// OE2–HE2 hydroxyl rotation.
    Glh
);
residue_type!(
    /// Glutamine — 3 χ, 0 polar H.
    ///
    /// Amide NH₂ is coplanar with the CD=OE1 plane.
    Gln
);
residue_type!(
    /// Glutamate (deprotonated, −1) — 3 χ, 0 polar H.
    ///
    /// Both γ-carboxylate oxygens (OE1, OE2) are deprotonated; no O–H hydrogen.
    Glu
);
residue_type!(
    /// Glycine — 0 χ, 0 polar H.
    Gly
);
residue_type!(
    /// Histidine (Nδ-protonated) — 2 χ, 0 polar H.
    ///
    /// HD1 is fixed in the imidazole ring plane.
    Hid
);
residue_type!(
    /// Histidine (Nε-protonated) — 2 χ, 0 polar H.
    ///
    /// HE2 is fixed in the imidazole ring plane.
    Hie
);
residue_type!(
    /// Histidine (doubly protonated, +1) — 2 χ, 0 polar H.
    ///
    /// Both HD1 and HE2 are fixed in the imidazole ring plane.
    Hip
);
residue_type!(
    /// Isoleucine — 2 χ, 0 polar H.
    Ile
);
residue_type!(
    /// Leucine — 2 χ, 0 polar H.
    Leu
);
residue_type!(
    /// Lysine (neutral, −NH₂) — 4 χ, 1 polar H.
    ///
    /// NZ amine rotation (single torsion for the NH₂ group).
    Lyn
);
residue_type!(
    /// Lysine (protonated, +1, −NH₃⁺) — 4 χ, 1 polar H.
    ///
    /// NZ ammonium rotation (single torsion for the NH₃⁺ group).
    Lys
);
residue_type!(
    /// Methionine — 3 χ, 0 polar H.
    Met
);
residue_type!(
    /// Phenylalanine — 2 χ, 0 polar H.
    Phe
);
residue_type!(
    /// Proline — 2 χ, 0 polar H.
    ///
    /// No backbone N–H; pyrrolidine ring constrains geometry.
    Pro
);
residue_type!(
    /// Serine — 1 χ, 1 polar H.
    ///
    /// OG–HG hydroxyl rotation.
    Ser
);
residue_type!(
    /// Threonine — 1 χ, 1 polar H.
    ///
    /// OG1–HG1 hydroxyl rotation.
    Thr
);
residue_type!(
    /// Tryptophan — 2 χ, 0 polar H.
    ///
    /// NE1–HE1 is fixed in the indole ring plane.
    Trp
);
residue_type!(
    /// Tyrosine (deprotonated phenolate, −1) — 2 χ, 0 polar H.
    ///
    /// OH bears no hydrogen (phenolate O⁻).
    Tym
);
residue_type!(
    /// Tyrosine (neutral phenol) — 2 χ, 1 polar H.
    ///
    /// OH–HH phenol hydroxyl rotation.
    Tyr
);
residue_type!(
    /// Valine — 1 χ, 0 polar H.
    Val
);
