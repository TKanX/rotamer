use criterion::{Criterion, black_box, criterion_group, criterion_main};
use dunbrack::Residue as _;
use rotamer::*;

const N: Vec3 = Vec3::new(0.0, 1.458, 0.0);
const CA: Vec3 = Vec3::new(0.0, 0.0, 0.0);
const C: Vec3 = Vec3::new(1.525, 0.0, 0.0);

const GRID_COUNT: usize = 37;
const GRID_MIN: f32 = -180.0;
const GRID_STEP: f32 = 10.0;

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("build");

    macro_rules! add_bench {
        ($T:ident, 0, 0, $n:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$T:lower>]), |b| {
                    b.iter(|| {
                        black_box($T::build(
                            black_box(N),
                            black_box(CA),
                            black_box(C),
                        ))
                    });
                });
            }
        };
        ($T:ident, $nc:literal, 0, $n:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$T:lower>]), |b| {
                    b.iter(|| {
                        black_box($T::build(
                            black_box(N),
                            black_box(CA),
                            black_box(C),
                            black_box([0.0_f32; $nc]),
                        ))
                    });
                });
            }
        };
        ($T:ident, $nc:literal, $np:literal, $n:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$T:lower>]), |b| {
                    b.iter(|| {
                        black_box($T::build(
                            black_box(N),
                            black_box(CA),
                            black_box(C),
                            black_box([0.0_f32; $nc]),
                            black_box([0.0_f32; $np]),
                        ))
                    });
                });
            }
        };
    }
    for_all_sidechains!(add_bench);

    group.finish();
}

fn bench_sweep(c: &mut Criterion) {
    let mut group = c.benchmark_group("sweep");

    macro_rules! add_bench {
        ($rot:ident, $dun:ident, $nc:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$rot:lower>]), |b| {
                    b.iter(|| {
                        for phi_idx in 0..GRID_COUNT {
                            for psi_idx in 0..GRID_COUNT {
                                let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                                let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;
                                for rot in dunbrack::$dun::rotamers(phi, psi) {
                                    let chi: [f32; $nc] =
                                        core::array::from_fn(|i| rot.chi_mean[i].to_radians());
                                    black_box($rot::build(N, CA, C, chi));
                                }
                            }
                        }
                    });
                });
            }
        };
    }

    macro_rules! add_bench_ph {
        ($rot:ident, $dun:ident, $nc:literal) => {
            paste::paste! {
                group.bench_function(stringify!([<$rot:lower>]), |b| {
                    b.iter(|| {
                        for phi_idx in 0..GRID_COUNT {
                            for psi_idx in 0..GRID_COUNT {
                                let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                                let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;
                                for rot in dunbrack::$dun::rotamers(phi, psi) {
                                    let chi: [f32; $nc] =
                                        core::array::from_fn(|i| rot.chi_mean[i].to_radians());
                                    black_box($rot::build(N, CA, C, chi, [0.0_f32]));
                                }
                            }
                        }
                    });
                });
            }
        };
    }

    add_bench!(Arg, Arg, 4);
    add_bench!(Arn, Arg, 4);
    add_bench!(Asn, Asn, 2);
    add_bench!(Asp, Asp, 2);
    add_bench!(Cym, Cyh, 1);
    add_bench!(Cyx, Cyd, 1);
    add_bench!(Gln, Gln, 3);
    add_bench!(Glu, Glu, 3);
    add_bench!(Hid, His, 2);
    add_bench!(Hie, His, 2);
    add_bench!(Hip, His, 2);
    add_bench!(Ile, Ile, 2);
    add_bench!(Leu, Leu, 2);
    add_bench!(Met, Met, 3);
    add_bench!(Phe, Phe, 2);
    add_bench!(Trp, Trp, 2);
    add_bench!(Tym, Tyr, 2);
    add_bench!(Val, Val, 1);

    add_bench_ph!(Ash, Asp, 2);
    add_bench_ph!(Cys, Cyh, 1);
    add_bench_ph!(Glh, Glu, 3);
    add_bench_ph!(Lyn, Lys, 4);
    add_bench_ph!(Lys, Lys, 4);
    add_bench_ph!(Ser, Ser, 1);
    add_bench_ph!(Thr, Thr, 1);
    add_bench_ph!(Tyr, Tyr, 2);

    group.bench_function("pro_tpr", |b| {
        b.iter(|| {
            for phi_idx in 0..GRID_COUNT {
                for psi_idx in 0..GRID_COUNT {
                    let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                    let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;
                    for rot in dunbrack::Tpr::rotamers(phi, psi) {
                        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                        black_box(Pro::build(N, CA, C, chi));
                    }
                }
            }
        });
    });

    group.bench_function("pro_cpr", |b| {
        b.iter(|| {
            for phi_idx in 0..GRID_COUNT {
                for psi_idx in 0..GRID_COUNT {
                    let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                    let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;
                    for rot in dunbrack::Cpr::rotamers(phi, psi) {
                        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                        black_box(Pro::build(N, CA, C, chi));
                    }
                }
            }
        });
    });

    group.bench_function("pro_pool", |b| {
        b.iter(|| {
            for phi_idx in 0..GRID_COUNT {
                for psi_idx in 0..GRID_COUNT {
                    let phi = GRID_MIN + phi_idx as f32 * GRID_STEP;
                    let psi = GRID_MIN + psi_idx as f32 * GRID_STEP;
                    for rot in dunbrack::Pro::rotamers(phi, psi) {
                        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                        black_box(Pro::build(N, CA, C, chi));
                    }
                }
            }
        });
    });

    group.finish();
}

criterion_group!(benches, bench_build, bench_sweep);
criterion_main!(benches);
