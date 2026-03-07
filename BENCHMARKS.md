# Performance Benchmarks

Sidechain coordinate placement via compile-time-generated NERF chains. All benchmarks with Criterion.rs, 100 samples.

## Executive Summary

- **Trivial**: Gly — **0.78 ns** (0 atoms, empty struct, optimized to no-op)
- **Smallest**: Ala — **41.5 ns** (4 atoms, 0 `sincosf`, pure `place` chain)
- **Largest**: Trp — **329 ns** (18 atoms, 2 `sincosf`)
- **Full pipeline** (Dunbrack query + NERF build, 37×37 grid): Val **310 µs** → Arg **23.5 ms**
- **Zero allocation**: coordinates returned as `#[repr(C)]` stack struct (max 216 bytes for Trp)
- **Zero runtime constants**: bond lengths, bond angles, and fixed torsion sin/cos are `f32` literals hardcoded at compile time

## Detailed Results

### Single Build (`build/*`)

Time to call `Type::build(N, CA, C [, chi] [, polar_h])` once. Each call runs a fully-inlined, straight-line NERF chain that places all sidechain atoms from three backbone coordinates and up to N_CHI + N_PH runtime torsion angles.

| Type |  N  | sincosf | Time (ns) |  ns/atom |
| :--- | :-: | :-----: | --------: | -------: |
| Gly  |  0  |    0    |  **0.78** |        — |
| Ala  |  4  |    0    |  **41.5** | **10.4** |
| Cym  |  4  |    1    |      61.0 |     15.2 |
| Cyx  |  4  |    1    |      61.0 |     15.2 |
| Cys  |  5  |    2    |      94.7 |     18.9 |
| Ser  |  5  |    2    |      95.1 |     19.0 |
| Val  | 10  |    1    |     118.9 |     11.9 |
| Asp  |  6  |    2    |     130.8 |     21.8 |
| Thr  |  8  |    2    |     132.4 |     16.6 |
| Pro  |  9  |    2    |     150.5 |     16.7 |
| Asn  |  8  |    2    |     156.3 |     19.5 |
| Ile  | 13  |    2    |     165.3 |     12.7 |
| Ash  |  7  |    3    |     168.3 |     24.0 |
| Leu  | 13  |    2    |     168.9 |     13.0 |
| Glu  |  9  |    3    |     174.2 |     19.4 |
| Met  | 11  |    3    |     175.4 |     16.0 |
| Gln  | 11  |    3    |     200.9 |     18.3 |
| Glh  | 10  |    4    |     201.8 |     20.2 |
| Hid  | 11  |    2    |     230.2 |     20.9 |
| Hie  | 11  |    2    |     231.0 |     21.0 |
| Hip  | 12  |    2    |     235.9 |     19.7 |
| Lyn  | 15  |    5    |     240.6 |     16.0 |
| Lys  | 16  |    5    |     242.5 |     15.2 |
| Phe  | 14  |    2    |     264.5 |     18.9 |
| Tym  | 14  |    2    |     272.4 |     19.5 |
| Arg  | 18  |    4    |     291.3 |     16.2 |
| Arn  | 17  |    4    |     295.7 |     17.4 |
| Tyr  | 15  |    3    |     304.2 |     20.3 |
| Trp  | 18  |    2    | **329.3** | **18.3** |

- **N** — total sidechain atoms placed per build
- **sincosf** — runtime `sincosf` calls = N_CHI + N_PH; all other torsions use compile-time `TrigPair` literals

### Full Pipeline Sweep (`sweep/*`)

Time to query all 37 × 37 = 1,369 (φ, ψ) grid points via `dunbrack`, then call `build()` for every rotamer at every grid point. Measures the end-to-end Dunbrack query + NERF build pipeline.

| Type      | Source | N_ROT | Time (ms) | Per-Point (µs) |
| :-------- | :----: | :---: | --------: | -------------: |
| Cym       |  Cyh   |   3   |     0.122 |           0.09 |
| Cyx       |  Cyd   |   3   |     0.121 |           0.09 |
| Cys       |  Cyh   |   3   |     0.196 |           0.14 |
| Ser       |  Ser   |   3   |     0.196 |           0.14 |
| Thr       |  Thr   |   3   |     0.292 |           0.21 |
| Val       |  Val   |   3   |     0.310 |           0.23 |
| Pro(tpr)  |  Tpr   |   2   |     0.353 |           0.26 |
| Pro(cpr)  |  Cpr   |   2   |     0.353 |           0.26 |
| Pro(pool) |  Pro   |   2   |     0.353 |           0.26 |
| Ile       |  Ile   |   9   |     1.466 |           1.07 |
| Leu       |  Leu   |   9   |     1.522 |           1.11 |
| Asp       |  Asp   |  18   |     1.881 |           1.37 |
| Ash       |  Asp   |  18   |     2.536 |           1.85 |
| Met       |  Met   |  27   |     4.632 |           3.38 |
| Phe       |  Phe   |  18   |     5.338 |           3.90 |
| Asn       |  Asn   |  36   |     5.401 |           3.94 |
| Tym       |  Tyr   |  18   |     5.724 |           4.18 |
| Tyr       |  Tyr   |  18   |     6.305 |           4.61 |
| Hid       |  His   |  36   |     8.721 |           6.37 |
| Hie       |  His   |  36   |     8.741 |           6.39 |
| Glu       |  Glu   |  54   |     9.123 |           6.66 |
| Hip       |  His   |  36   |     9.313 |           6.80 |
| Glh       |  Glu   |  54   |    10.475 |           7.65 |
| Trp       |  Trp   |  36   |    13.740 |          10.04 |
| Lyn       |  Lys   |  73   |    20.279 |          14.81 |
| Lys       |  Lys   |  73   |    20.398 |          14.90 |
| Gln       |  Gln   |  108  |    20.255 |          14.80 |
| Arg       |  Arg   |  75   |    23.538 |          17.19 |
| Arn       |  Arg   |  75   |    23.525 |          17.18 |

- **Source** — Dunbrack residue type queried for rotamer library _lookup_
- **N_ROT** — rotamers per grid point from `-180°` to `+180°` in 10° steps
- **Per-Point** — Time / 1,369 = amortized cost to build all rotamers at one (φ, ψ)

## Test Environment

- **CPU**: Intel® Core™ i7-13620H (Raptor Lake), 6P+4E, 4.90 GHz turbo, AVX2
- **OS**: Linux (Arch Linux, Kernel 6.12.63-1)
- **Profile**: `opt-level=3, lto=true, codegen-units=1`
- **Date**: March 2026

## Performance Analysis

### 1. Cost Model: `place()` + `sincosf()`

Each `build()` executes two kinds of work:

- **N × `place()`** — one NERF placement per atom: two `normalize()` (each containing `rsqrtf` via Quake III bit-trick + 3 Newton–Raphson steps), one `cross()`, three `mul_scalar()`, and one final `add()`. Roughly 30 FP ops per atom.
- **(N_CHI + N_PH) × `sincosf()`** — custom branchless minimax polynomial (degree 9 sin, degree 8 cos, quadrant folding via bit manipulation). Called once per runtime torsion angle.

Approximate model: `time ≈ N × 10 ns + (N_CHI + N_PH) × 20 ns`

```text
Ala  (N=4,  S=0):   4×10 +  0×20 =  40 ns → actual  41.5 ns ✓
Val  (N=10, S=1):  10×10 +  1×20 = 120 ns → actual 118.9 ns ✓
Met  (N=11, S=3):  11×10 +  3×20 = 170 ns → actual 175.4 ns ✓
Lys  (N=16, S=5):  16×10 +  5×20 = 260 ns → actual 242.5 ns ≈
```

Ring residues (Phe, Tyr, Trp, His, Pro) run ~10–20% above the model due to longer dependency chains in ring closure atoms and higher register pressure.

### 2. Compile-Time Precomputation

The build script (`build/codegen.rs`) computes at **build time**:

- **Bond lengths** → `f32` literals (e.g., `1.5287_f32`), zero runtime cost
- **Bond angles** → `TrigPair { cos: ..., sin: ... }` literals, zero runtime `sincosf`
- **Fixed torsion angles** → `TrigPair` literals (sin/cos computed at `f64` precision, truncated to `f32`)
- **Range validation** — bond length (0.8–2.5 Å) and bond angle (90°–180°) assertions at compile time

Only χ and polar-H torsions require runtime `sincosf` calls. The ratio of eliminated trig operations:

```text
Ala  (N=4):   4 Fixed,  0 runtime → 4 sincosf eliminated (100%)
Val  (N=10):  9 Fixed,  1 runtime → 9 sincosf eliminated (90%)
Trp  (N=18): 16 Fixed,  2 runtime → 16 sincosf eliminated (89%)
Arg  (N=18): 14 Fixed,  4 runtime → 14 sincosf eliminated (78%)
```

Estimated speedup from precomputation (Trp): 16 eliminated × ~20 ns = **320 ns saved**, nearly equal to the entire measured build time (329 ns). Without codegen, Trp would cost ~650 ns — **2× slower**.

### 3. Monomorphized Build Functions

Each of the 29 types gets a **dedicated** `fn build()` generated by the build script. Properties:

- **No loop** — N sequential `place()` calls, fully unrolled at source level
- **No branch** — each atom's torsion source (Fixed/Chi/PolarH) is resolved at compile time
- **No array indexing** — bond parameters are inline literals, not table lookups
- **`#[inline(always)]`** on `place()`, `sincosf()`, `rsqrtf()` → entire build compiles to a single basic block
- **LLVM constant propagation** — backbone coords + precomputed TrigPairs enable further arithmetic simplification

The generated code for Ala (4 atoms, 0 runtime torsions) is 4 straight `place()` calls with all 6 parameters as compile-time constants except the three backbone `Vec3`s.

### 4. Sweep Pipeline Decomposition

Each sweep iteration: **Dunbrack query** (rotamer lookup + bilinear interpolation) → **chi extraction** (`to_radians()`) → **NERF build** (backbone + chi → 3D coordinates).

The total cost per grid point scales with `N_ROT × N_ATOMS`:

```text
Val  (3 rot × 10 atoms = 30):      0.23 µs/point
Phe  (18 rot × 14 atoms = 252):    3.90 µs/point
Gln  (108 rot × 11 atoms = 1188):  14.80 µs/point
Arg  (75 rot × 18 atoms = 1350):   17.19 µs/point
```

Arg is slower than Gln despite fewer rotamers because 75 × 18 = 1,350 > 108 × 11 = 1,188 total atom placements per grid point. The dominant cost is NERF `place()`, not the Dunbrack query.

Types sharing the same Dunbrack source show the NERF cost difference directly: Asp (6 atoms, 1.37 µs/point) vs Ash (7 atoms + 1 polar-H sincosf, 1.85 µs/point) — the extra atom + extra sincosf adds 35%.

### 5. Memory and Allocation

All output is a `#[repr(C)]` struct of `N` × `Vec3` (12 bytes each), returned on the stack:

```text
Gly:   0 bytes    Ala:  48 bytes    Val: 120 bytes
Trp: 216 bytes    Arg: 216 bytes    (largest)
```

No heap allocation, no `Vec`, no `Box`. The `#![no_std]` crate has zero runtime dependencies.

## Reproducibility

```bash
cargo bench --bench sidechain
```

Full HTML reports are generated in `target/criterion/`.
