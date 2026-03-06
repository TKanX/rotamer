use dunbrack::Residue as _;
use rotamer::*;

mod common;
use common::*;

fn assert_l_chirality(cb: Vec3, name: &str) {
    let sign = chirality_sign(cb);
    assert!(
        sign > 0.0,
        "{name}: wrong chirality at Cα (sign = {sign:.6}, expected positive for L-amino acid)"
    );
}

macro_rules! chirality {
    ($T:ident, $nc:literal, $np:literal, 0) => {};
    ($T:ident, 0, 0, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<chirality_ $T:lower>]() {
                let coords = $T::build(N, CA, C);
                assert_l_chirality(coords.as_slice()[0], stringify!($T));
            }
        }
    };
    ($T:ident, $nc:literal, 0, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<chirality_ $T:lower>]() {
                let coords = $T::build(N, CA, C, [0.0_f32; $nc]);
                assert_l_chirality(coords.as_slice()[0], stringify!($T));
            }
        }
    };
    ($T:ident, $nc:literal, $np:literal, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<chirality_ $T:lower>]() {
                let coords = $T::build(N, CA, C, [0.0_f32; $nc], [0.0_f32; $np]);
                assert_l_chirality(coords.as_slice()[0], stringify!($T));
            }
        }
    };
}

for_all_sidechains!(chirality);

macro_rules! chirality_grid {
    ($T:ident, $DunT:path, $nc:literal) => {
        paste::paste! {
            #[test]
            fn [<chirality_ $T:lower _grid>]() {
                for phi_i in 0..37_u32 {
                    let phi = -180.0 + phi_i as f32 * 10.0;
                    for psi_i in 0..37_u32 {
                        let psi = -180.0 + psi_i as f32 * 10.0;
                        for rot in <$DunT>::rotamers(phi, psi) {
                            let mut chi = [0.0_f32; $nc];
                            for (i, c) in chi.iter_mut().enumerate() {
                                *c = rot.chi_mean[i].to_radians();
                            }
                            let coords = $T::build(N, CA, C, chi);
                            assert_l_chirality(coords.as_slice()[0], stringify!($T));
                        }
                    }
                }
            }
        }
    };
    ($T:ident, $DunT:path, $nc:literal, ph) => {
        paste::paste! {
            #[test]
            fn [<chirality_ $T:lower _grid>]() {
                for phi_i in 0..37_u32 {
                    let phi = -180.0 + phi_i as f32 * 10.0;
                    for psi_i in 0..37_u32 {
                        let psi = -180.0 + psi_i as f32 * 10.0;
                        for rot in <$DunT>::rotamers(phi, psi) {
                            let mut chi = [0.0_f32; $nc];
                            for (i, c) in chi.iter_mut().enumerate() {
                                *c = rot.chi_mean[i].to_radians();
                            }
                            let coords = $T::build(N, CA, C, chi, [0.0_f32]);
                            assert_l_chirality(coords.as_slice()[0], stringify!($T));
                        }
                    }
                }
            }
        }
    };
}

chirality_grid!(Val, dunbrack::Val, 1);
chirality_grid!(Ile, dunbrack::Ile, 2);
chirality_grid!(Leu, dunbrack::Leu, 2);
chirality_grid!(Phe, dunbrack::Phe, 2);
chirality_grid!(Tym, dunbrack::Tyr, 2);
chirality_grid!(Trp, dunbrack::Trp, 2);
chirality_grid!(Cym, dunbrack::Cyh, 1);
chirality_grid!(Cyx, dunbrack::Cyd, 1);

#[test]
fn chirality_cys_pool_grid() {
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Cys::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians()];
                assert_l_chirality(Cys::build(N, CA, C, chi, [0.0]).as_slice()[0], "CYS(pool)");
                assert_l_chirality(Cym::build(N, CA, C, chi).as_slice()[0], "CYM(pool)");
                assert_l_chirality(Cyx::build(N, CA, C, chi).as_slice()[0], "CYX(pool)");
            }
        }
    }
}

chirality_grid!(Met, dunbrack::Met, 3);
chirality_grid!(Asn, dunbrack::Asn, 2);
chirality_grid!(Gln, dunbrack::Gln, 3);
chirality_grid!(Glu, dunbrack::Glu, 3);
chirality_grid!(Glh, dunbrack::Glu, 3, ph);
chirality_grid!(Asp, dunbrack::Asp, 2);
chirality_grid!(Ash, dunbrack::Asp, 2, ph);
chirality_grid!(Arg, dunbrack::Arg, 4);
chirality_grid!(Arn, dunbrack::Arg, 4);
chirality_grid!(Lyn, dunbrack::Lys, 4, ph);
chirality_grid!(Hid, dunbrack::His, 2);
chirality_grid!(Hie, dunbrack::His, 2);
chirality_grid!(Hip, dunbrack::His, 2);

#[test]
fn chirality_pro_tpr_grid() {
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Tpr::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                assert_l_chirality(Pro::build(N, CA, C, chi).as_slice()[0], "PRO(Tpr)");
            }
        }
    }
}

#[test]
fn chirality_pro_cpr_grid() {
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Cpr::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                assert_l_chirality(Pro::build(N, CA, C, chi).as_slice()[0], "PRO(Cpr)");
            }
        }
    }
}

#[test]
fn chirality_pro_pool_grid() {
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Pro::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                assert_l_chirality(Pro::build(N, CA, C, chi).as_slice()[0], "PRO(pool)");
            }
        }
    }
}

chirality_grid!(Ser, dunbrack::Ser, 1, ph);
chirality_grid!(Thr, dunbrack::Thr, 1, ph);
chirality_grid!(Cys, dunbrack::Cyh, 1, ph);
chirality_grid!(Tyr, dunbrack::Tyr, 2, ph);
chirality_grid!(Lys, dunbrack::Lys, 4, ph);

#[test]
fn chirality_invariant_under_chi() {
    let angles: [f32; 7] = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0];
    for &d in &angles {
        let chi = [d.to_radians()];
        let ph = [0.0];
        let c = Ser::build(N, CA, C, chi, ph);
        let sign = chirality_sign(c.cb);
        assert!(
            sign > 0.0,
            "SER chirality flipped at chi1={d}° (sign={sign:.4})"
        );
    }
}

#[test]
fn chirality_rotated_backbone() {
    use rotamer::Vec3;

    let n = Vec3::new(1.458, 0.0, 0.0);
    let ca = Vec3::new(0.0, 0.0, 0.0);
    let c = Vec3::new(-0.553, 1.424, 0.0);

    let coords = Ala::build(n, ca, c);
    let cb = coords.as_slice()[0];
    let sign = (n - ca).cross(c - ca).dot(cb - ca);

    let ser = Ser::build(n, ca, c, [1.0], [0.0]);
    let cb_ser = ser.as_slice()[0];
    let sign_ser = (n - ca).cross(c - ca).dot(cb_ser - ca);

    assert!(
        sign * sign_ser > 0.0,
        "Ala and Ser disagree on chirality with rotated backbone"
    );
}
