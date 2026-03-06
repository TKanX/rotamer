use dunbrack::Residue as _;
use rotamer::*;

mod common;
use common::*;

/// Maximum acceptable angular error (radians) for chi round-trip.
const CHI_TOL: f32 = 5e-4;

/// Chi-1 angles to test (degrees).
const CHI1_DEGS: [f32; 7] = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0];

/// Chi angles for scanning deeper torsions (degrees).
const CHI_DEEP: [f32; 5] = [-180.0, -60.0, 0.0, 60.0, 180.0];

/// Polar-H angles to test (degrees).
const PH_DEGS: [f32; 5] = [-120.0, -60.0, 0.0, 60.0, 120.0];

// -----------------------------------------------------------------------
// Serine: 1 χ, 1 polar H
//   chi1 = dihedral(N, CA, CB, OG)
//   ph0  = dihedral(CA, CB, OG, HG)
// -----------------------------------------------------------------------

#[test]
fn chi1_ser() {
    for &d in &CHI1_DEGS {
        let chi = [d.to_radians()];
        let ph = [0.0_f32];
        let c = Ser::build(N, CA, C, chi, ph);
        let got = dihedral(N, CA, c.cb, c.og);
        let err = angle_diff(got, chi[0]);
        assert!(err.abs() < CHI_TOL, "SER chi1={d}°: err = {:.2e} rad", err);
    }
}

#[test]
fn polar_h_ser() {
    for &ph_deg in &PH_DEGS {
        let chi = [60.0_f32.to_radians()]; // gauche+
        let ph = [ph_deg.to_radians()];
        let c = Ser::build(N, CA, C, chi, ph);
        let got = dihedral(CA, c.cb, c.og, c.hg);
        let err = angle_diff(got, ph[0]);
        assert!(
            err.abs() < CHI_TOL,
            "SER ph0={ph_deg}° at chi1=60°: err = {:.2e} rad",
            err
        );
    }
}

// -----------------------------------------------------------------------
// Valine: 1 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG1)
// -----------------------------------------------------------------------

#[test]
fn chi1_val() {
    for &d in &CHI1_DEGS {
        let chi = [d.to_radians()];
        let c = Val::build(N, CA, C, chi);
        let got = dihedral(N, CA, c.cb, c.cg1);
        let err = angle_diff(got, chi[0]);
        assert!(err.abs() < CHI_TOL, "VAL chi1={d}°: err = {:.2e} rad", err);
    }
}

// -----------------------------------------------------------------------
// Phenylalanine: 2 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD1)
// -----------------------------------------------------------------------

#[test]
fn chi12_phe() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Phe::build(N, CA, C, chi);

            let got1 = dihedral(N, CA, c.cb, c.cg);
            let err1 = angle_diff(got1, chi[0]);
            assert!(
                err1.abs() < CHI_TOL,
                "PHE chi1={d1}°: err = {:.2e} rad",
                err1
            );

            let got2 = dihedral(CA, c.cb, c.cg, c.cd1);
            let err2 = angle_diff(got2, chi[1]);
            assert!(
                err2.abs() < CHI_TOL,
                "PHE chi2={d2}° (chi1={d1}°): err = {:.2e} rad",
                err2
            );
        }
    }
}

// -----------------------------------------------------------------------
// Threonine: 1 χ, 1 polar H
//   chi1 = dihedral(N, CA, CB, OG1)
//   ph0  = dihedral(CA, CB, OG1, HG1)
// -----------------------------------------------------------------------

#[test]
fn chi1_thr() {
    for &d in &CHI1_DEGS {
        let chi = [d.to_radians()];
        let ph = [0.0];
        let c = Thr::build(N, CA, C, chi, ph);
        let got = dihedral(N, CA, c.cb, c.og1);
        let err = angle_diff(got, chi[0]);
        assert!(err.abs() < CHI_TOL, "THR chi1={d}°: err = {:.2e} rad", err);
    }
}

#[test]
fn polar_h_thr() {
    for &ph_deg in &PH_DEGS {
        let chi = [60.0_f32.to_radians()];
        let ph = [ph_deg.to_radians()];
        let c = Thr::build(N, CA, C, chi, ph);
        let got = dihedral(CA, c.cb, c.og1, c.hg1);
        let err = angle_diff(got, ph[0]);
        assert!(
            err.abs() < CHI_TOL,
            "THR ph0={ph_deg}° at chi1=60°: err = {:.2e} rad",
            err
        );
    }
}

// -----------------------------------------------------------------------
// Methionine: 3 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, SD)
//   chi3 = dihedral(CB, CG, SD, CE)
// -----------------------------------------------------------------------

#[test]
fn chi123_met() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            for &d3 in &CHI_DEEP {
                let chi = [d1.to_radians(), d2.to_radians(), d3.to_radians()];
                let c = Met::build(N, CA, C, chi);

                let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                assert!(err1.abs() < CHI_TOL, "MET chi1={d1}°: err={err1:.2e}");

                let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.sd), chi[1]);
                assert!(err2.abs() < CHI_TOL, "MET chi2={d2}°: err={err2:.2e}");

                let err3 = angle_diff(dihedral(c.cb, c.cg, c.sd, c.ce), chi[2]);
                assert!(err3.abs() < CHI_TOL, "MET chi3={d3}°: err={err3:.2e}");
            }
        }
    }
}

// -----------------------------------------------------------------------
// Arginine: 4 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD)
//   chi3 = dihedral(CB, CG, CD, NE)
//   chi4 = dihedral(CG, CD, NE, CZ)
// -----------------------------------------------------------------------

#[test]
fn chi1234_arg() {
    let angles: [f32; 3] = [-60.0, 60.0, 180.0];
    for &d1 in &angles {
        for &d2 in &angles {
            for &d3 in &angles {
                for &d4 in &angles {
                    let chi = [
                        d1.to_radians(),
                        d2.to_radians(),
                        d3.to_radians(),
                        d4.to_radians(),
                    ];
                    let c = Arg::build(N, CA, C, chi);

                    let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                    assert!(err1.abs() < CHI_TOL, "ARG chi1={d1}°: err={err1:.2e}");

                    let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                    assert!(err2.abs() < CHI_TOL, "ARG chi2={d2}°: err={err2:.2e}");

                    let err3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ne), chi[2]);
                    assert!(err3.abs() < CHI_TOL, "ARG chi3={d3}°: err={err3:.2e}");

                    let err4 = angle_diff(dihedral(c.cg, c.cd, c.ne, c.cz), chi[3]);
                    assert!(err4.abs() < CHI_TOL, "ARG chi4={d4}°: err={err4:.2e}");
                }
            }
        }
    }
}

// -----------------------------------------------------------------------
// Tryptophan: 2 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD1)
// -----------------------------------------------------------------------

#[test]
fn chi12_trp() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Trp::build(N, CA, C, chi);

            let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            assert!(err1.abs() < CHI_TOL, "TRP chi1={d1}°: err={err1:.2e}");

            let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
            assert!(err2.abs() < CHI_TOL, "TRP chi2={d2}°: err={err2:.2e}");
        }
    }
}

// -----------------------------------------------------------------------
// Tyrosine: 2 χ, 1 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD1)
//   ph0  = dihedral(CE1, CZ, OH, HH)
// -----------------------------------------------------------------------

#[test]
fn chi12_tyr() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let ph = [0.0_f32];
            let c = Tyr::build(N, CA, C, chi, ph);

            let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            assert!(err1.abs() < CHI_TOL, "TYR chi1={d1}°: err={err1:.2e}");

            let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
            assert!(err2.abs() < CHI_TOL, "TYR chi2={d2}°: err={err2:.2e}");
        }
    }
}

#[test]
fn polar_h_tyr() {
    for &ph_deg in &PH_DEGS {
        let chi = [60.0_f32.to_radians(), 90.0_f32.to_radians()];
        let ph = [ph_deg.to_radians()];
        let c = Tyr::build(N, CA, C, chi, ph);
        let got = dihedral(c.ce1, c.cz, c.oh, c.hh);
        let err = angle_diff(got, ph[0]);
        assert!(
            err.abs() < CHI_TOL,
            "TYR ph0={ph_deg}°: err = {:.2e} rad",
            err
        );
    }
}

// -----------------------------------------------------------------------
// Lysine: 4 χ, 1 polar H
//   chi[0..4] = N-CA-CB-CG, CA-CB-CG-CD, CB-CG-CD-CE, CG-CD-CE-NZ
//   ph0       = dihedral(CD, CE, NZ, HZ1)
// -----------------------------------------------------------------------

#[test]
fn chi1234_lys() {
    let angles: [f32; 3] = [-60.0, 60.0, 180.0];
    for &d1 in &angles {
        for &d2 in &angles {
            for &d3 in &angles {
                for &d4 in &angles {
                    let chi = [
                        d1.to_radians(),
                        d2.to_radians(),
                        d3.to_radians(),
                        d4.to_radians(),
                    ];
                    let ph = [0.0];
                    let c = Lys::build(N, CA, C, chi, ph);

                    let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                    assert!(err1.abs() < CHI_TOL, "LYS chi1={d1}°: err={err1:.2e}");

                    let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                    assert!(err2.abs() < CHI_TOL, "LYS chi2={d2}°: err={err2:.2e}");

                    let err3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ce), chi[2]);
                    assert!(err3.abs() < CHI_TOL, "LYS chi3={d3}°: err={err3:.2e}");

                    let err4 = angle_diff(dihedral(c.cg, c.cd, c.ce, c.nz), chi[3]);
                    assert!(err4.abs() < CHI_TOL, "LYS chi4={d4}°: err={err4:.2e}");
                }
            }
        }
    }
}

#[test]
fn polar_h_lys() {
    let chi = [60.0_f32.to_radians(); 4];
    for &ph_deg in &PH_DEGS {
        let ph = [ph_deg.to_radians()];
        let c = Lys::build(N, CA, C, chi, ph);
        let got = dihedral(c.cd, c.ce, c.nz, c.hz1);
        let err = angle_diff(got, ph[0]);
        assert!(
            err.abs() < CHI_TOL,
            "LYS ph0={ph_deg}°: err = {:.2e} rad",
            err
        );
    }
}

// -----------------------------------------------------------------------
// Isoleucine: 2 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG1)
//   chi2 = dihedral(CA, CB, CG1, CD1)
// -----------------------------------------------------------------------

#[test]
fn chi12_ile() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Ile::build(N, CA, C, chi);

            let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg1), chi[0]);
            assert!(err1.abs() < CHI_TOL, "ILE chi1={d1}°: err={err1:.2e}");

            let err2 = angle_diff(dihedral(CA, c.cb, c.cg1, c.cd1), chi[1]);
            assert!(err2.abs() < CHI_TOL, "ILE chi2={d2}°: err={err2:.2e}");
        }
    }
}

// -----------------------------------------------------------------------
// Leucine: 2 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD1)
// -----------------------------------------------------------------------

#[test]
fn chi12_leu() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Leu::build(N, CA, C, chi);

            let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            assert!(err1.abs() < CHI_TOL, "LEU chi1={d1}°: err={err1:.2e}");

            let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
            assert!(err2.abs() < CHI_TOL, "LEU chi2={d2}°: err={err2:.2e}");
        }
    }
}

// -----------------------------------------------------------------------
// Glutamine: 3 χ, 0 polar H
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD)
//   chi3 = dihedral(CB, CG, CD, OE1)
// -----------------------------------------------------------------------

#[test]
fn chi123_gln() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            for &d3 in &CHI_DEEP {
                let chi = [d1.to_radians(), d2.to_radians(), d3.to_radians()];
                let c = Gln::build(N, CA, C, chi);

                let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                assert!(err1.abs() < CHI_TOL, "GLN chi1={d1}°: err={err1:.2e}");

                let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                assert!(err2.abs() < CHI_TOL, "GLN chi2={d2}°: err={err2:.2e}");

                let err3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.oe1), chi[2]);
                assert!(err3.abs() < CHI_TOL, "GLN chi3={d3}°: err={err3:.2e}");
            }
        }
    }
}

// -----------------------------------------------------------------------
// Cysteine (free thiol): 1 χ, 1 polar H
//   chi1 = dihedral(N, CA, CB, SG)
//   ph0  = dihedral(CA, CB, SG, HG)
// -----------------------------------------------------------------------

#[test]
fn chi1_cys() {
    for &d in &CHI1_DEGS {
        let chi = [d.to_radians()];
        let ph = [0.0];
        let c = Cys::build(N, CA, C, chi, ph);
        let got = dihedral(N, CA, c.cb, c.sg);
        let err = angle_diff(got, chi[0]);
        assert!(err.abs() < CHI_TOL, "CYS chi1={d}°: err={err:.2e}");
    }
}

#[test]
fn polar_h_cys() {
    for &ph_deg in &PH_DEGS {
        let chi = [60.0_f32.to_radians()];
        let ph = [ph_deg.to_radians()];
        let c = Cys::build(N, CA, C, chi, ph);
        let got = dihedral(CA, c.cb, c.sg, c.hg);
        let err = angle_diff(got, ph[0]);
        assert!(err.abs() < CHI_TOL, "CYS ph0={ph_deg}°: err={err:.2e}");
    }
}

// -----------------------------------------------------------------------
// Proline: 2 χ → round-trip (fixed angles)
//   chi1 = dihedral(N, CA, CB, CG)
//   chi2 = dihedral(CA, CB, CG, CD)
// -----------------------------------------------------------------------

#[test]
fn chi12_pro() {
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Pro::build(N, CA, C, chi);

            let err1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            assert!(err1.abs() < CHI_TOL, "PRO chi1={d1}°: err={err1:.2e}");

            let err2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
            assert!(err2.abs() < CHI_TOL, "PRO chi2={d2}°: err={err2:.2e}");
        }
    }
}

// -----------------------------------------------------------------------
// Residues covered only by full-grid (protonation variants, His family,
// Asn, Asp, Glu). Add uniform-sample round-trips to exercise 0° and
// other angles not present in Dunbrack rotamer means.
// -----------------------------------------------------------------------

#[test]
fn chi12_asn() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, OD1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Asn::build(N, CA, C, chi);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.od1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "ASN chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "ASN chi2={d2}°: err={e2:.2e}");
        }
    }
}

#[test]
fn chi12_asp() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, OD1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Asp::build(N, CA, C, chi);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.od1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "ASP chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "ASP chi2={d2}°: err={e2:.2e}");
        }
    }
}

#[test]
fn chi12_ash() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, OD1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Ash::build(N, CA, C, chi, [0.0]);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.od1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "ASH chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "ASH chi2={d2}°: err={e2:.2e}");
        }
    }
}

#[test]
fn chi123_glu() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, CD)
    //   chi3 = dihedral(CB, CG, CD, OE1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            for &d3 in &CHI_DEEP {
                let chi = [d1.to_radians(), d2.to_radians(), d3.to_radians()];
                let c = Glu::build(N, CA, C, chi);
                let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.oe1), chi[2]);
                assert!(e1.abs() < CHI_TOL, "GLU chi1={d1}°: err={e1:.2e}");
                assert!(e2.abs() < CHI_TOL, "GLU chi2={d2}°: err={e2:.2e}");
                assert!(e3.abs() < CHI_TOL, "GLU chi3={d3}°: err={e3:.2e}");
            }
        }
    }
}

#[test]
fn chi123_glh() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, CD)
    //   chi3 = dihedral(CB, CG, CD, OE1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            for &d3 in &CHI_DEEP {
                let chi = [d1.to_radians(), d2.to_radians(), d3.to_radians()];
                let c = Glh::build(N, CA, C, chi, [0.0]);
                let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.oe1), chi[2]);
                assert!(e1.abs() < CHI_TOL, "GLH chi1={d1}°: err={e1:.2e}");
                assert!(e2.abs() < CHI_TOL, "GLH chi2={d2}°: err={e2:.2e}");
                assert!(e3.abs() < CHI_TOL, "GLH chi3={d3}°: err={e3:.2e}");
            }
        }
    }
}

#[test]
fn chi1234_arn() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, CD)
    //   chi3 = dihedral(CB, CG, CD, NE)
    //   chi4 = dihedral(CG, CD, NE, CZ)
    let angles: [f32; 3] = [-60.0, 60.0, 180.0];
    for &d1 in &angles {
        for &d2 in &angles {
            for &d3 in &angles {
                for &d4 in &angles {
                    let chi = [
                        d1.to_radians(),
                        d2.to_radians(),
                        d3.to_radians(),
                        d4.to_radians(),
                    ];
                    let c = Arn::build(N, CA, C, chi);
                    let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                    let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                    let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ne), chi[2]);
                    let e4 = angle_diff(dihedral(c.cg, c.cd, c.ne, c.cz), chi[3]);
                    assert!(e1.abs() < CHI_TOL, "ARN chi1={d1}°: err={e1:.2e}");
                    assert!(e2.abs() < CHI_TOL, "ARN chi2={d2}°: err={e2:.2e}");
                    assert!(e3.abs() < CHI_TOL, "ARN chi3={d3}°: err={e3:.2e}");
                    assert!(e4.abs() < CHI_TOL, "ARN chi4={d4}°: err={e4:.2e}");
                }
            }
        }
    }
}

#[test]
fn chi1234_lyn() {
    //   chi1-4: same atoms as LYS but NZ is neutral
    let angles: [f32; 3] = [-60.0, 60.0, 180.0];
    for &d1 in &angles {
        for &d2 in &angles {
            for &d3 in &angles {
                for &d4 in &angles {
                    let chi = [
                        d1.to_radians(),
                        d2.to_radians(),
                        d3.to_radians(),
                        d4.to_radians(),
                    ];
                    let c = Lyn::build(N, CA, C, chi, [0.0]);
                    let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
                    let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
                    let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ce), chi[2]);
                    let e4 = angle_diff(dihedral(c.cg, c.cd, c.ce, c.nz), chi[3]);
                    assert!(e1.abs() < CHI_TOL, "LYN chi1={d1}°: err={e1:.2e}");
                    assert!(e2.abs() < CHI_TOL, "LYN chi2={d2}°: err={e2:.2e}");
                    assert!(e3.abs() < CHI_TOL, "LYN chi3={d3}°: err={e3:.2e}");
                    assert!(e4.abs() < CHI_TOL, "LYN chi4={d4}°: err={e4:.2e}");
                }
            }
        }
    }
}

#[test]
fn chi1_cym() {
    //   chi1 = dihedral(N, CA, CB, SG)
    for &d in &CHI1_DEGS {
        let chi = [d.to_radians()];
        let c = Cym::build(N, CA, C, chi);
        let err = angle_diff(dihedral(N, CA, c.cb, c.sg), chi[0]);
        assert!(err.abs() < CHI_TOL, "CYM chi1={d}°: err={err:.2e}");
    }
}

#[test]
fn chi1_cyx() {
    //   chi1 = dihedral(N, CA, CB, SG)
    for &d in &CHI1_DEGS {
        let chi = [d.to_radians()];
        let c = Cyx::build(N, CA, C, chi);
        let err = angle_diff(dihedral(N, CA, c.cb, c.sg), chi[0]);
        assert!(err.abs() < CHI_TOL, "CYX chi1={d}°: err={err:.2e}");
    }
}

#[test]
fn chi12_tym() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, CD1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Tym::build(N, CA, C, chi);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "TYM chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "TYM chi2={d2}°: err={e2:.2e}");
        }
    }
}

#[test]
fn chi12_hid() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, ND1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Hid::build(N, CA, C, chi);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.nd1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "HID chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "HID chi2={d2}°: err={e2:.2e}");
        }
    }
}

#[test]
fn chi12_hie() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, ND1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Hie::build(N, CA, C, chi);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.nd1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "HIE chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "HIE chi2={d2}°: err={e2:.2e}");
        }
    }
}

#[test]
fn chi12_hip() {
    //   chi1 = dihedral(N, CA, CB, CG)
    //   chi2 = dihedral(CA, CB, CG, ND1)
    for &d1 in &CHI_DEEP {
        for &d2 in &CHI_DEEP {
            let chi = [d1.to_radians(), d2.to_radians()];
            let c = Hip::build(N, CA, C, chi);
            let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
            let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.nd1), chi[1]);
            assert!(e1.abs() < CHI_TOL, "HIP chi1={d1}°: err={e1:.2e}");
            assert!(e2.abs() < CHI_TOL, "HIP chi2={d2}°: err={e2:.2e}");
        }
    }
}

// -----------------------------------------------------------------------
// Full Dunbrack grid scans
// Scan the complete 37×37 (φ,ψ) grid × all Dunbrack rotamers and verify
// that each chi angle round-trips exactly through the NERF geometry.
// -----------------------------------------------------------------------

/// Iterate the full 37×37 (φ,ψ) grid and all Dunbrack rotamers.
macro_rules! grid_scan {
    ($DunT:path, |$rot:ident| $body:block) => {{
        for phi_i in 0..37_u32 {
            let phi = -180.0 + phi_i as f32 * 10.0;
            for psi_i in 0..37_u32 {
                let psi = -180.0 + psi_i as f32 * 10.0;
                #[allow(unused_variables)]
                for $rot in <$DunT>::rotamers(phi, psi) {
                    $body
                }
            }
        }
    }};
}

#[test]
fn chi1_ser_grid() {
    grid_scan!(dunbrack::Ser, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let c = Ser::build(N, CA, C, chi, [0.0]);
        let err = angle_diff(dihedral(N, CA, c.cb, c.og), chi[0]);
        assert!(err.abs() < CHI_TOL, "SER chi1 full-grid: err={err:.2e}");
    });
}

#[test]
fn chi1_val_grid() {
    grid_scan!(dunbrack::Val, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let c = Val::build(N, CA, C, chi);
        let err = angle_diff(dihedral(N, CA, c.cb, c.cg1), chi[0]);
        assert!(err.abs() < CHI_TOL, "VAL chi1 full-grid: err={err:.2e}");
    });
}

#[test]
fn chi12_phe_grid() {
    grid_scan!(dunbrack::Phe, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Phe::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "PHE chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "PHE chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi1_thr_grid() {
    grid_scan!(dunbrack::Thr, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let c = Thr::build(N, CA, C, chi, [0.0]);
        let err = angle_diff(dihedral(N, CA, c.cb, c.og1), chi[0]);
        assert!(err.abs() < CHI_TOL, "THR chi1 full-grid: err={err:.2e}");
    });
}

#[test]
fn chi123_met_grid() {
    grid_scan!(dunbrack::Met, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
        ];
        let c = Met::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.sd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.sd, c.ce), chi[2]);
        assert!(e1.abs() < CHI_TOL, "MET chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "MET chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "MET chi3 full-grid: err={e3:.2e}");
    });
}

#[test]
fn chi1234_arg_grid() {
    grid_scan!(dunbrack::Arg, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
            rot.chi_mean[3].to_radians(),
        ];
        let c = Arg::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ne), chi[2]);
        let e4 = angle_diff(dihedral(c.cg, c.cd, c.ne, c.cz), chi[3]);
        assert!(e1.abs() < CHI_TOL, "ARG chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "ARG chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "ARG chi3 full-grid: err={e3:.2e}");
        assert!(e4.abs() < CHI_TOL, "ARG chi4 full-grid: err={e4:.2e}");
    });
}

#[test]
fn chi12_trp_grid() {
    grid_scan!(dunbrack::Trp, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Trp::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "TRP chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "TRP chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi12_tyr_grid() {
    grid_scan!(dunbrack::Tyr, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Tyr::build(N, CA, C, chi, [0.0]);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "TYR chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "TYR chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi1234_lys_grid() {
    grid_scan!(dunbrack::Lys, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
            rot.chi_mean[3].to_radians(),
        ];
        let c = Lys::build(N, CA, C, chi, [0.0]);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ce), chi[2]);
        let e4 = angle_diff(dihedral(c.cg, c.cd, c.ce, c.nz), chi[3]);
        assert!(e1.abs() < CHI_TOL, "LYS chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "LYS chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "LYS chi3 full-grid: err={e3:.2e}");
        assert!(e4.abs() < CHI_TOL, "LYS chi4 full-grid: err={e4:.2e}");
    });
}

#[test]
fn chi12_ile_grid() {
    grid_scan!(dunbrack::Ile, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Ile::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg1), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg1, c.cd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "ILE chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "ILE chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi12_leu_grid() {
    grid_scan!(dunbrack::Leu, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Leu::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "LEU chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "LEU chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi123_gln_grid() {
    grid_scan!(dunbrack::Gln, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
        ];
        let c = Gln::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.oe1), chi[2]);
        assert!(e1.abs() < CHI_TOL, "GLN chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "GLN chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "GLN chi3 full-grid: err={e3:.2e}");
    });
}

#[test]
fn chi1_cys_grid() {
    grid_scan!(dunbrack::Cyh, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let c = Cys::build(N, CA, C, chi, [0.0]);
        let err = angle_diff(dihedral(N, CA, c.cb, c.sg), chi[0]);
        assert!(err.abs() < CHI_TOL, "CYS chi1 full-grid: err={err:.2e}");
    });
}

#[test]
fn chi12_pro_tpr_grid() {
    grid_scan!(dunbrack::Tpr, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Pro::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        assert!(e1.abs() < CHI_TOL, "PRO chi1 full-grid (Tpr): err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "PRO chi2 full-grid (Tpr): err={e2:.2e}");
    });
}

#[test]
fn chi12_asp_grid() {
    grid_scan!(dunbrack::Asp, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Asp::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.od1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "ASP chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "ASP chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi12_ash_grid() {
    grid_scan!(dunbrack::Asp, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Ash::build(N, CA, C, chi, [0.0]);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.od1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "ASH chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "ASH chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi123_glu_grid() {
    grid_scan!(dunbrack::Glu, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
        ];
        let c = Glu::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.oe1), chi[2]);
        assert!(e1.abs() < CHI_TOL, "GLU chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "GLU chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "GLU chi3 full-grid: err={e3:.2e}");
    });
}

#[test]
fn chi123_glh_grid() {
    grid_scan!(dunbrack::Glu, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
        ];
        let c = Glh::build(N, CA, C, chi, [0.0]);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.oe1), chi[2]);
        assert!(e1.abs() < CHI_TOL, "GLH chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "GLH chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "GLH chi3 full-grid: err={e3:.2e}");
    });
}

#[test]
fn chi12_hid_grid() {
    grid_scan!(dunbrack::His, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Hid::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.nd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "HID chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "HID chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi12_hie_grid() {
    grid_scan!(dunbrack::His, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Hie::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.nd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "HIE chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "HIE chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi12_hip_grid() {
    grid_scan!(dunbrack::His, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Hip::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.nd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "HIP chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "HIP chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi1234_arn_grid() {
    grid_scan!(dunbrack::Arg, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
            rot.chi_mean[3].to_radians(),
        ];
        let c = Arn::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ne), chi[2]);
        let e4 = angle_diff(dihedral(c.cg, c.cd, c.ne, c.cz), chi[3]);
        assert!(e1.abs() < CHI_TOL, "ARN chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "ARN chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "ARN chi3 full-grid: err={e3:.2e}");
        assert!(e4.abs() < CHI_TOL, "ARN chi4 full-grid: err={e4:.2e}");
    });
}

#[test]
fn chi1234_lyn_grid() {
    grid_scan!(dunbrack::Lys, |rot| {
        let chi = [
            rot.chi_mean[0].to_radians(),
            rot.chi_mean[1].to_radians(),
            rot.chi_mean[2].to_radians(),
            rot.chi_mean[3].to_radians(),
        ];
        let c = Lyn::build(N, CA, C, chi, [0.0]);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        let e3 = angle_diff(dihedral(c.cb, c.cg, c.cd, c.ce), chi[2]);
        let e4 = angle_diff(dihedral(c.cg, c.cd, c.ce, c.nz), chi[3]);
        assert!(e1.abs() < CHI_TOL, "LYN chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "LYN chi2 full-grid: err={e2:.2e}");
        assert!(e3.abs() < CHI_TOL, "LYN chi3 full-grid: err={e3:.2e}");
        assert!(e4.abs() < CHI_TOL, "LYN chi4 full-grid: err={e4:.2e}");
    });
}

#[test]
fn chi1_cym_grid() {
    grid_scan!(dunbrack::Cyh, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let c = Cym::build(N, CA, C, chi);
        let err = angle_diff(dihedral(N, CA, c.cb, c.sg), chi[0]);
        assert!(err.abs() < CHI_TOL, "CYM chi1 full-grid: err={err:.2e}");
    });
}

#[test]
fn chi1_cyx_grid() {
    grid_scan!(dunbrack::Cyd, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let c = Cyx::build(N, CA, C, chi);
        let err = angle_diff(dihedral(N, CA, c.cb, c.sg), chi[0]);
        assert!(
            err.abs() < CHI_TOL,
            "CYX chi1 full-grid (Cyd): err={err:.2e}"
        );
    });
}

#[test]
fn chi1_cys_pool_grid() {
    grid_scan!(dunbrack::Cys, |rot| {
        let chi = [rot.chi_mean[0].to_radians()];
        let cs = Cys::build(N, CA, C, chi, [0.0]);
        let ec = angle_diff(dihedral(N, CA, cs.cb, cs.sg), chi[0]);
        let cm = Cym::build(N, CA, C, chi);
        let em = angle_diff(dihedral(N, CA, cm.cb, cm.sg), chi[0]);
        let cx = Cyx::build(N, CA, C, chi);
        let ex = angle_diff(dihedral(N, CA, cx.cb, cx.sg), chi[0]);
        assert!(ec.abs() < CHI_TOL, "CYS chi1 (pool): err={ec:.2e}");
        assert!(em.abs() < CHI_TOL, "CYM chi1 (pool): err={em:.2e}");
        assert!(ex.abs() < CHI_TOL, "CYX chi1 (pool): err={ex:.2e}");
    });
}

#[test]
fn chi12_tym_grid() {
    grid_scan!(dunbrack::Tyr, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Tym::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "TYM chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "TYM chi2 full-grid: err={e2:.2e}");
    });
}

#[test]
fn chi12_pro_cpr_grid() {
    grid_scan!(dunbrack::Cpr, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Pro::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        assert!(e1.abs() < CHI_TOL, "PRO chi1 full-grid (Cpr): err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "PRO chi2 full-grid (Cpr): err={e2:.2e}");
    });
}

#[test]
fn chi12_pro_pool_grid() {
    grid_scan!(dunbrack::Pro, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Pro::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.cd), chi[1]);
        assert!(
            e1.abs() < CHI_TOL,
            "PRO chi1 full-grid (Pro pool): err={e1:.2e}"
        );
        assert!(
            e2.abs() < CHI_TOL,
            "PRO chi2 full-grid (Pro pool): err={e2:.2e}"
        );
    });
}

#[test]
fn chi12_asn_grid() {
    grid_scan!(dunbrack::Asn, |rot| {
        let chi = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
        let c = Asn::build(N, CA, C, chi);
        let e1 = angle_diff(dihedral(N, CA, c.cb, c.cg), chi[0]);
        let e2 = angle_diff(dihedral(CA, c.cb, c.cg, c.od1), chi[1]);
        assert!(e1.abs() < CHI_TOL, "ASN chi1 full-grid: err={e1:.2e}");
        assert!(e2.abs() < CHI_TOL, "ASN chi2 full-grid: err={e2:.2e}");
    });
}
