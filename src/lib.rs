//! # Rotamer
//!
//! **Runtime sidechain coordinate placement via compile-time-generated NERF chains.**
//!
//! Generates all sidechain atom coordinates (heavy atoms + hydrogens) for 29 amino acid sidechain types from three backbone anchor points (N, Cα, C) and a set of χ / polar-H torsion angles. All bond lengths, bond angles, and fixed torsion sin/cos values are hardcoded as `f32` literals at compile time by `build.rs`; only the runtime-variable χ and polar-H angles incur `sincosf` calls. Zero heap allocation, zero runtime dependencies, `#![no_std]`.
//!
//! [Features](#features) • [Installation](#installation) • [Usage](#usage) • [Sidechain Types](#sidechain-types) • [Performance](#performance) • [Verification](#verification)
//!
//! ---
//!
//! ## Features
//!
//! - **Compile-time code generation.** `build.rs` reads 29 residue specifications from `build/residues/` and emits a dedicated `fn build()` per type — a straight-line sequence of `place()` calls with all bond geometry baked in as `f32` literals. No loops, no branches, no table lookups at runtime.
//! - **Precomputed torsion sin/cos.** Fixed dihedral angles (e.g. CB placement at 120°, methyl H at ±120°) are converted to `TrigPair { cos, sin }` constants at build time. Only χ and polar-H torsions call `sincosf` at runtime, eliminating 78–100% of trig work per build depending on the residue.
//! - **Zero heap allocation.** Each `build()` returns a `#[repr(C)]` struct of named `Vec3` fields (e.g. `SerCoords { cb, og, hb2, hb3, hg }`), stack-allocated. The struct can be reinterpreted as `&[Vec3]` via [`SidechainCoords::as_slice()`] at zero cost.
//! - **`#![no_std]` compatible.** No standard library, no libm linkage. Custom `sincosf` (minimax polynomial, max error < 5×10⁻⁷) and `rsqrtf` (Quake III + 3 Newton–Raphson steps). Usable in embedded, kernel, and WASM environments.
//! - **Build-time validation.** `build.rs` asserts bond lengths (0.8–2.5 Å), bond angles (90°–180°), torsion source consistency (Chi/PolarH indices within bounds), and atom placement ordering (all reference atoms must precede the placed atom) for every atom in every residue. Compilation fails on invalid geometry.
//! - **Sealed type-safe API.** Each residue is a zero-sized type implementing the sealed [`Sidechain`] trait with compile-time constants `N_CHI`, `N_POLAR_H`, `NAME`, and associated `Coords` type. The `build()` signature varies per type — no unused parameters, no Option wrapping.
//! - **`for_all_sidechains!` macro.** A generated declarative macro invoking `$callback!(Type, N_CHI, N_POLAR_H, N)` for all 29 types. Drives tests, benchmarks, and generic dispatch with zero boilerplate.
//!
//! ---
//!
//! ## Installation
//!
//! ```toml
//! [dependencies]
//! rotamer = "0.1.0"
//! ```
//!
//! **Note:** `build.rs` generates all 29 `build()` functions from the residue specifications in `build/residues/`. Initial compilation takes a few seconds.
//!
//! ---
//!
//! ## Usage
//!
//! ### Basic Build
//!
//! ```
//! use rotamer::{Ser, Vec3};
//!
//! let n  = Vec3::new(0.0, 1.458, 0.0);
//! let ca = Vec3::new(0.0, 0.0, 0.0);
//! let c  = Vec3::new(1.525, 0.0, 0.0);
//!
//! // Serine: 1 χ angle, 1 polar-H angle
//! let chi = [1.0_f32];       // χ₁ in radians
//! let polar_h = [0.0_f32];   // OG–HG rotation in radians
//! let coords = Ser::build(n, ca, c, chi, polar_h);
//!
//! // coords.cb, coords.og, coords.hb2, coords.hb3, coords.hg — all Vec3
//! println!("OG = ({:.3}, {:.3}, {:.3})", coords.og.x, coords.og.y, coords.og.z);
//! ```
//!
//! ### Slice Access
//!
//! All coordinate structs implement [`SidechainCoords`], providing a flat `&[Vec3]` view:
//!
//! ```
//! use rotamer::{Ala, SidechainCoords, Vec3};
//!
//! let coords = Ala::build(
//!     Vec3::new(0.0, 1.458, 0.0),
//!     Vec3::new(0.0, 0.0, 0.0),
//!     Vec3::new(1.525, 0.0, 0.0),
//! );
//! let atoms: &[Vec3] = coords.as_slice();
//! assert_eq!(atoms.len(), 4); // CB, HB1, HB2, HB3
//! ```
//!
//! ### Build Signature Variants
//!
//! The `build()` signature adapts to each residue's degrees of freedom:
//!
//! ```
//! use rotamer::*;
//!
//! let n  = Vec3::new(0.0, 1.458, 0.0);
//! let ca = Vec3::new(0.0, 0.0, 0.0);
//! let c  = Vec3::new(1.525, 0.0, 0.0);
//!
//! // 0 χ, 0 polar H → build(N, CA, C)
//! let _ = Gly::build(n, ca, c);
//! let _ = Ala::build(n, ca, c);
//!
//! // N_CHI χ, 0 polar H → build(N, CA, C, chi)
//! let _ = Val::build(n, ca, c, [1.0]);
//! let _ = Phe::build(n, ca, c, [1.0, 0.5]);
//!
//! // N_CHI χ, N_PH polar H → build(N, CA, C, chi, polar_h)
//! let _ = Ser::build(n, ca, c, [1.0], [0.0]);
//! let _ = Lys::build(n, ca, c, [1.0, 0.5, 0.0, -0.5], [0.0]);
//! ```
//!
//! ### `for_all_sidechains!` Macro
//!
//! Invokes `$callback!(Type, N_CHI, N_POLAR_H, N)` for all 29 sidechain types:
//!
//! ```
//! use rotamer::*;
//!
//! macro_rules! print_info {
//!     ($T:ident, $nc:literal, $np:literal, $n:literal) => {
//!         println!("{}: {} χ, {} polar H, {} atoms",
//!             <$T as Sidechain>::NAME, $nc, $np, $n);
//!     };
//! }
//!
//! for_all_sidechains!(print_info);
//! ```
//!
//! ### Integration with `dunbrack`
//!
//! ```
//! use dunbrack::Residue as _;
//! use rotamer::*;
//!
//! let n  = Vec3::new(0.0, 1.458, 0.0);
//! let ca = Vec3::new(0.0, 0.0, 0.0);
//! let c  = Vec3::new(1.525, 0.0, 0.0);
//!
//! // Query Dunbrack rotamer library → build NERF coordinates
//! for rot in dunbrack::Val::rotamers(-60.0, -40.0) {
//!     let chi = [rot.chi_mean[0].to_radians()];
//!     let coords = Val::build(n, ca, c, chi);
//!     println!("CB = ({:.3}, {:.3}, {:.3})", coords.cb.x, coords.cb.y, coords.cb.z);
//! }
//! ```
//!
//! ---
//!
//! ## Sidechain Types
//!
//! All 29 sidechain types, including protonation state variants:
//!
//! | Type    | N_CHI | N_PH |   N | Notes                             |
//! | :------ | :---: | :--: | --: | :-------------------------------- |
//! | [`Gly`] |   0   |  0   |   0 | No sidechain atoms                |
//! | [`Ala`] |   0   |  0   |   4 | CB + 3 HB (all Fixed torsions)    |
//! | [`Val`] |   1   |  0   |  10 |                                   |
//! | [`Cym`] |   1   |  0   |   4 | Cysteine thiolate (S⁻)            |
//! | [`Cyx`] |   1   |  0   |   4 | Cystine (disulfide)               |
//! | [`Cys`] |   1   |  1   |   5 | Free cysteine (SG–HG)             |
//! | [`Ser`] |   1   |  1   |   5 | Serine (OG–HG)                    |
//! | [`Thr`] |   1   |  1   |   8 |                                   |
//! | [`Pro`] |   2   |  0   |   9 | Pyrrolidine ring, no backbone N–H |
//! | [`Asp`] |   2   |  0   |   6 | Aspartate (COO⁻)                  |
//! | [`Asn`] |   2   |  0   |   8 |                                   |
//! | [`Ile`] |   2   |  0   |  13 |                                   |
//! | [`Leu`] |   2   |  0   |  13 |                                   |
//! | [`Phe`] |   2   |  0   |  14 | Aromatic ring                     |
//! | [`Tym`] |   2   |  0   |  14 | Tyrosine phenolate (O⁻)           |
//! | [`Hid`] |   2   |  0   |  11 | Histidine (Nδ-protonated)         |
//! | [`Hie`] |   2   |  0   |  11 | Histidine (Nε-protonated)         |
//! | [`Hip`] |   2   |  0   |  12 | Histidine (doubly protonated, +1) |
//! | [`Trp`] |   2   |  0   |  18 | Indole ring                       |
//! | [`Ash`] |   2   |  1   |   7 | Aspartate protonated (OD2–HD2)    |
//! | [`Tyr`] |   2   |  1   |  15 | Tyrosine phenol (OH–HH)           |
//! | [`Met`] |   3   |  0   |  11 |                                   |
//! | [`Glu`] |   3   |  0   |   9 | Glutamate (COO⁻)                  |
//! | [`Gln`] |   3   |  0   |  11 |                                   |
//! | [`Glh`] |   3   |  1   |  10 | Glutamate protonated (OE2–HE2)    |
//! | [`Arg`] |   4   |  0   |  18 | Arginine (protonated, +1)         |
//! | [`Arn`] |   4   |  0   |  17 | Arginine (neutral)                |
//! | [`Lyn`] |   4   |  1   |  15 | Lysine neutral (NH₂)              |
//! | [`Lys`] |   4   |  1   |  16 | Lysine protonated (NH₃⁺)          |
//!
//! - **N_CHI** — number of χ dihedral angles from the rotamer library (runtime `sincosf`)
//! - **N_PH** — number of polar-hydrogen torsions (runtime `sincosf`)
//! - **N** — total sidechain atoms placed (heavy atoms + all sidechain hydrogens; backbone HN and HA excluded)
//!
//! Each type implements [`Sidechain`] `+ Copy + PartialEq + Eq + Hash + Debug`.
//!
//! ---
//!
//! ## Performance
//!
//! Benchmarked with `Criterion.rs` on an Intel® Core™ i7-13620H (Raptor Lake, 4.90 GHz turbo, AVX2), Linux, `opt-level=3, lto=true, codegen-units=1`.
//!
//! **Single build** — time to call `Type::build(N, CA, C [, chi] [, polar_h])` once:
//!
//! | Type    |   N | Time (ns) | ns/atom |
//! | :------ | --: | --------: | ------: |
//! | [`Gly`] |   0 |      0.78 |       — |
//! | [`Ala`] |   4 |      41.5 |    10.4 |
//! | [`Val`] |  10 |     118.9 |    11.9 |
//! | [`Asp`] |   6 |     130.8 |    21.8 |
//! | [`Pro`] |   9 |     150.5 |    16.7 |
//! | [`Leu`] |  13 |     168.9 |    13.0 |
//! | [`Met`] |  11 |     175.4 |    16.0 |
//! | [`Hid`] |  11 |     230.2 |    20.9 |
//! | [`Lys`] |  16 |     242.5 |    15.2 |
//! | [`Phe`] |  14 |     264.5 |    18.9 |
//! | [`Arg`] |  18 |     291.3 |    16.2 |
//! | [`Tyr`] |  15 |     304.2 |    20.3 |
//! | [`Trp`] |  18 |     329.3 |    18.3 |
//!
//! The cost has two components: **N × ~10 ns** (NERF `place()` per atom) + **(N_CHI + N_PH) × ~20 ns** (runtime `sincosf` per torsion angle). Fixed-torsion atoms (Ala, methyl H's, ring atoms) require no `sincosf` — their trig values are compile-time `f32` literals.
//!
//! **Full pipeline sweep** (Dunbrack query + NERF build, 37×37 = 1,369 grid points):
//!
//! | Type    | N_ROT | Time (ms) | Per-Point (µs) |
//! | :------ | :---: | --------: | -------------: |
//! | [`Val`] |   3   |     0.310 |           0.23 |
//! | [`Phe`] |  18   |     5.338 |           3.90 |
//! | [`Trp`] |  36   |    13.740 |          10.04 |
//! | [`Gln`] |  108  |    20.255 |          14.80 |
//! | [`Arg`] |  75   |    23.538 |          17.19 |
//!
//! For full data including all 29 types and detailed analysis, see [BENCHMARKS.md](https://github.com/TKanX/rotamer/blob/main/BENCHMARKS.md).
//!
//! ### Why it's fast
//!
//! | Optimization                            | Savings                                                               |
//! | :-------------------------------------- | :-------------------------------------------------------------------- |
//! | Precomputed torsion `TrigPair` literals | Eliminates 78–100% of `sincosf` calls (e.g. Trp: 16/18 are Fixed)     |
//! | Dedicated `fn build()` per type         | No loop, no branch, no array indexing — straight-line `place()` chain |
//! | Custom branchless `sincosf` + `rsqrtf`  | No libm; polynomial + bit-trick; fully inlined                        |
//! | `#[inline(always)]` everywhere          | Entire build compiles to a single basic block in the caller           |
//! | `#[repr(C)]` coordinate struct          | Zero-copy `as_slice()`; no indirection, no allocator                  |
//!
//! ---
//!
//! ## Verification
//!
//! The library is verified at three levels:
//!
//! **Compile time (`build.rs` assertions)** — compilation aborts if any of the following fail:
//!
//! - Atom reference ordering: all three reference atoms placed before the current atom
//! - Bond length ∈ [0.8, 2.5] Å for every atom in every residue
//! - Bond angle ∈ (90°, 180°) for every atom in every residue
//! - Chi/PolarH torsion indices within declared `n_chi` / `n_polar_h`
//! - `#[repr(C)]` layout: `size_of::<Coords>() == N * size_of::<Vec3>()`
//!
//! **Unit tests** (37 tests in `src/`):
//!
//! - `sincosf` accuracy: max error < 5×10⁻⁷ over [−2π, 2π], Pythagorean identity, quadrant signs
//! - `rsqrtf` accuracy across float range
//! - `Vec3` arithmetic: add, sub, cross, dot, normalize, `#[repr(C)]` layout
//! - NERF `place()`: bond length preservation, bond angle preservation, dihedral round-trip
//!
//! **Integration tests** (227 tests in `tests/`):
//!
//! | File           | Tests | What is verified                                                                                                     |
//! | :------------- | :---: | :------------------------------------------------------------------------------------------------------------------- |
//! | `smoke.rs`     |  29   | All 29 types build without panic; atom count matches; all coordinates finite and non-NaN                             |
//! | `constants.rs` |  29   | `N_CHI`, `N_POLAR_H`, `N`, `NAME` match declared values for all 29 types                                             |
//! | `bond.rs`      |  30   | All sidechain atoms have at least one neighbor within 2.1 Å across the full 37×37 (φ, ψ) grid with Dunbrack χ values |
//! | `chirality.rs` |  60   | L-amino acid CB chirality sign correct; spot checks + 37×37 grid with Dunbrack χ values for all 28 non-Gly types     |
//! | `geometry.rs`  |  62   | χ dihedral round-trip (< 5×10⁻⁴ rad) across uniform χ sampling + 37×37 grid; polar-H dihedral round-trip             |
//! | `ring.rs`      |  17   | Ring closure gap (placed vs expected distance) < 0.05 Å aromatic, < 0.50 Å proline; all 8 ring-bearing residues      |
//!
//! Total: **264 tests** + **58 criterion benchmarks**.
//!
//! Run the full suite:
//!
//! ```bash
//! cargo test
//! cargo bench --bench sidechain
//! ```

#![no_std]

mod math;
mod nerf;
mod residue;
mod sealed;

include!(concat!(env!("OUT_DIR"), "/generated.rs"));

pub use math::Vec3;
pub use residue::{
    Ala, Arg, Arn, Ash, Asn, Asp, Cym, Cys, Cyx, Glh, Gln, Glu, Gly, Hid, Hie, Hip, Ile, Leu, Lyn,
    Lys, Met, Phe, Pro, Ser, Thr, Trp, Tym, Tyr, Val,
};
pub use residue::{Sidechain, SidechainCoords};
