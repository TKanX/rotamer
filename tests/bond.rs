use dunbrack::Residue as _;
use rotamer::*;

mod common;
use common::*;

fn verify_chi_build(atoms: &[Vec3], name: &str) {
    assert_coords_valid(atoms, name);
    assert_neighbors(atoms, name);
}

macro_rules! bond_chi {
    ($rot:ident, $dun:ident, $nc:literal) => {
        paste::paste! {
            #[test]
            fn [<bond_ $rot:lower>]() {
                for phi_i in 0..37_u32 {
                    let phi = -180.0 + phi_i as f32 * 10.0;
                    for psi_i in 0..37_u32 {
                        let psi = -180.0 + psi_i as f32 * 10.0;
                        for rot in dunbrack::$dun::rotamers(phi, psi) {
                            let chi: [f32; $nc] =
                                core::array::from_fn(|i| rot.chi_mean[i].to_radians());
                            let coords = $rot::build(N, CA, C, chi);
                            verify_chi_build(coords.as_slice(), stringify!($rot));
                        }
                    }
                }
            }
        }
    };
}

macro_rules! bond_chi_ph {
    ($rot:ident, $dun:ident, $nc:literal, $np:literal) => {
        paste::paste! {
            #[test]
            fn [<bond_ $rot:lower>]() {
                const PH_SAMPLES: [f32; 3] = [0.0, 2.094_395_1, -2.094_395_1];
                for phi_i in 0..37_u32 {
                    let phi = -180.0 + phi_i as f32 * 10.0;
                    for psi_i in 0..37_u32 {
                        let psi = -180.0 + psi_i as f32 * 10.0;
                        for rot in dunbrack::$dun::rotamers(phi, psi) {
                            let chi: [f32; $nc] =
                                core::array::from_fn(|i| rot.chi_mean[i].to_radians());
                            for &ph in &PH_SAMPLES {
                                let polar_h = [ph; $np];
                                let coords = $rot::build(N, CA, C, chi, polar_h);
                                verify_chi_build(coords.as_slice(), stringify!($rot));
                            }
                        }
                    }
                }
            }
        }
    };
}

bond_chi!(Arg, Arg, 4);
bond_chi!(Asn, Asn, 2);
bond_chi!(Asp, Asp, 2);
bond_chi!(Gln, Gln, 3);
bond_chi!(Glu, Glu, 3);
bond_chi!(Ile, Ile, 2);
bond_chi!(Leu, Leu, 2);
bond_chi!(Met, Met, 3);
bond_chi!(Phe, Phe, 2);
bond_chi!(Trp, Trp, 2);
bond_chi!(Val, Val, 1);

bond_chi!(Arn, Arg, 4);
bond_chi!(Cym, Cyh, 1);
bond_chi!(Cyx, Cyd, 1);
bond_chi!(Hid, His, 2);
bond_chi!(Hie, His, 2);
bond_chi!(Hip, His, 2);
bond_chi!(Tym, Tyr, 2);

bond_chi_ph!(Ser, Ser, 1, 1);
bond_chi_ph!(Thr, Thr, 1, 1);
bond_chi_ph!(Tyr, Tyr, 2, 1);
bond_chi_ph!(Lys, Lys, 4, 1);

bond_chi_ph!(Ash, Asp, 2, 1);
bond_chi_ph!(Cys, Cyh, 1, 1);
bond_chi_ph!(Glh, Glu, 3, 1);
bond_chi_ph!(Lyn, Lys, 4, 1);

#[test]
fn bond_pro() {
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Tpr::rotamers(phi, psi) {
                let chi: [f32; 2] = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let coords = Pro::build(N, CA, C, chi);
                verify_chi_build(coords.as_slice(), "Pro/Tpr");
            }
            for rot in dunbrack::Cpr::rotamers(phi, psi) {
                let chi: [f32; 2] = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let coords = Pro::build(N, CA, C, chi);
                verify_chi_build(coords.as_slice(), "Pro/Cpr");
            }
            for rot in dunbrack::Pro::rotamers(phi, psi) {
                let chi: [f32; 2] = [rot.chi_mean[0].to_radians(), rot.chi_mean[1].to_radians()];
                let coords = Pro::build(N, CA, C, chi);
                verify_chi_build(coords.as_slice(), "Pro/pool");
            }
        }
    }
}

#[test]
fn bond_cys_pool() {
    for phi_i in 0..37_u32 {
        let phi = -180.0 + phi_i as f32 * 10.0;
        for psi_i in 0..37_u32 {
            let psi = -180.0 + psi_i as f32 * 10.0;
            for rot in dunbrack::Cys::rotamers(phi, psi) {
                let chi: [f32; 1] = [rot.chi_mean[0].to_radians()];
                verify_chi_build(Cys::build(N, CA, C, chi, [0.0]).as_slice(), "Cys/pool");
                verify_chi_build(Cym::build(N, CA, C, chi).as_slice(), "Cym/pool");
                verify_chi_build(Cyx::build(N, CA, C, chi).as_slice(), "Cyx/pool");
            }
        }
    }
}

#[test]
fn bond_ala() {
    let coords = Ala::build(N, CA, C);
    verify_chi_build(coords.as_slice(), "Ala");
}

#[test]
fn bond_gly() {
    let coords = Gly::build(N, CA, C);
    verify_chi_build(coords.as_slice(), "Gly");
}
