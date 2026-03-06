use dunbrack::Residue as _;
use rotamer::*;

mod common;
use common::*;

const AROMATIC_TOL: f32 = 0.05;
const PROLINE_TOL: f32 = 0.50;

const C_R_C_R: f32 = 1.39;
const C_R_N_R: f32 = 1.35;
const C_3_N_3: f32 = 1.47;

#[test]
fn ring_phe() {
    for rot in dunbrack::Phe::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Phe::build(N, CA, C, chi);
        let gap = dist(coords.cz, coords.ce2);
        let err = (gap - C_R_C_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "PHE ring closure CZ–CE2: |{gap:.4} − {C_R_C_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_phe_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Phe::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Phe::build(N, CA, C, chi);
                let err = (dist(c.cz, c.ce2) - C_R_C_R).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < AROMATIC_TOL,
        "PHE max ring closure error across grid: {max_err:.4} Å"
    );
}

#[test]
fn ring_tyr() {
    for rot in dunbrack::Tyr::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Tyr::build(N, CA, C, chi, [0.0]);
        let gap = dist(coords.cz, coords.ce2);
        let err = (gap - C_R_C_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "TYR ring closure CZ–CE2: |{gap:.4} − {C_R_C_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_tym() {
    for rot in dunbrack::Tyr::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Tym::build(N, CA, C, chi);
        let gap = dist(coords.cz, coords.ce2);
        let err = (gap - C_R_C_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "TYM ring closure CZ–CE2: |{gap:.4} − {C_R_C_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_tyr_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Tyr::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Tyr::build(N, CA, C, chi, [0.0]);
                let err = (dist(c.cz, c.ce2) - C_R_C_R).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < AROMATIC_TOL,
        "TYR max ring closure error across grid: {max_err:.4} Å"
    );
}

#[test]
fn ring_tym_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Tyr::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Tym::build(N, CA, C, chi);
                let err = (dist(c.cz, c.ce2) - C_R_C_R).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < AROMATIC_TOL,
        "TYM max ring closure error across grid: {max_err:.4} Å"
    );
}

#[test]
fn ring_hid() {
    for rot in dunbrack::His::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Hid::build(N, CA, C, chi);
        let gap = dist(coords.ce1, coords.ne2);
        let err = (gap - C_R_N_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "HID ring closure CE1–NE2: |{gap:.4} − {C_R_N_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_hie() {
    for rot in dunbrack::His::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Hie::build(N, CA, C, chi);
        let gap = dist(coords.ce1, coords.ne2);
        let err = (gap - C_R_N_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "HIE ring closure CE1–NE2: |{gap:.4} − {C_R_N_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_hip() {
    for rot in dunbrack::His::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Hip::build(N, CA, C, chi);
        let gap = dist(coords.ce1, coords.ne2);
        let err = (gap - C_R_N_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "HIP ring closure CE1–NE2: |{gap:.4} − {C_R_N_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_hid_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::His::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Hid::build(N, CA, C, chi);
                let err = (dist(c.ce1, c.ne2) - C_R_N_R).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < AROMATIC_TOL,
        "HID max ring closure error across grid: {max_err:.4} Å"
    );
}

#[test]
fn ring_hie_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::His::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Hie::build(N, CA, C, chi);
                let err = (dist(c.ce1, c.ne2) - C_R_N_R).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < AROMATIC_TOL,
        "HIE max ring closure error across grid: {max_err:.4} Å"
    );
}

#[test]
fn ring_hip_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::His::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Hip::build(N, CA, C, chi);
                let err = (dist(c.ce1, c.ne2) - C_R_N_R).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < AROMATIC_TOL,
        "HIP max ring closure error across grid: {max_err:.4} Å"
    );
}

#[test]
fn ring_trp_pyrrole() {
    for rot in dunbrack::Trp::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Trp::build(N, CA, C, chi);
        let gap = dist(coords.ne1, coords.ce2);
        let err = (gap - C_R_N_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "TRP pyrrole ring closure NE1–CE2: |{gap:.4} − {C_R_N_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_trp_benzene() {
    for rot in dunbrack::Trp::rotamers(-60.0, -40.0) {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Trp::build(N, CA, C, chi);
        let gap = dist(coords.ch2, coords.cz3);
        let err = (gap - C_R_C_R).abs();
        assert!(
            err < AROMATIC_TOL,
            "TRP benzene ring closure CH2–CZ3: |{gap:.4} − {C_R_C_R}| = {err:.4} Å"
        );
    }
}

#[test]
fn ring_trp_grid() {
    let mut max_pyrrole: f32 = 0.0;
    let mut max_benzene: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Trp::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Trp::build(N, CA, C, chi);
                let e_p = (dist(c.ne1, c.ce2) - C_R_N_R).abs();
                let e_b = (dist(c.ch2, c.cz3) - C_R_C_R).abs();
                if e_p > max_pyrrole {
                    max_pyrrole = e_p;
                }
                if e_b > max_benzene {
                    max_benzene = e_b;
                }
            }
        }
    }
    assert!(
        max_pyrrole < AROMATIC_TOL,
        "TRP max pyrrole closure error: {max_pyrrole:.4} Å"
    );
    assert!(
        max_benzene < AROMATIC_TOL,
        "TRP max benzene closure error: {max_benzene:.4} Å"
    );
}

#[test]
fn ring_pro() {
    for rot in dunbrack::Tpr::rotamers(-60.0, -40.0) {
        let chi: [f32; 2] = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let coords = Pro::build(N, CA, C, chi);
        let gap = dist(coords.cd, N);
        let err = (gap - C_3_N_3).abs();
        assert!(
            err < PROLINE_TOL,
            "PRO ring closure CD–N: |{gap:.4} − {C_3_N_3}| = {err:.4} Å (>{PROLINE_TOL})"
        );
    }
}

#[test]
fn ring_pro_grid() {
    let mut max_err: f32 = 0.0;
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Tpr::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Pro::build(N, CA, C, chi);
                let err = (dist(c.cd, N) - C_3_N_3).abs();
                if err > max_err {
                    max_err = err;
                }
            }
            for rot in dunbrack::Cpr::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Pro::build(N, CA, C, chi);
                let err = (dist(c.cd, N) - C_3_N_3).abs();
                if err > max_err {
                    max_err = err;
                }
            }
            for rot in dunbrack::Pro::rotamers(phi, psi) {
                let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let c = Pro::build(N, CA, C, chi);
                let err = (dist(c.cd, N) - C_3_N_3).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
    }
    assert!(
        max_err < PROLINE_TOL,
        "PRO max ring closure error across grid: {max_err:.4} Å"
    );
}
