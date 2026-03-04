use crate::math::{TrigPair, Vec3};

/// Places atom D given reference frame A–B–C.
///
/// # Geometry
///
/// Returns the unique point D such that:
///
/// - `|C − D| = d`
/// - `∠BCD = arccos(theta.cos)` where `theta.sin = sin(∠BCD) > 0`
/// - `dihedral(A, B, C, D) = arctan2(phi.sin, phi.cos)`
///
/// # Preconditions
///
/// - A, B, C must not be collinear (zero cross product → undefined frame).
/// - `theta.sin > 0` (the bond angle must not be 0° or 180°).
/// - `d > 0`.
///
/// Violating any precondition yields an unspecified but defined (non-UB) result.
#[inline(always)]
pub fn place(a: Vec3, b: Vec3, c: Vec3, d: f32, theta: TrigPair, phi: TrigPair) -> Vec3 {
    let bc = (c - b).normalize();
    let n = (b - a).cross(bc).normalize();
    let m = n.cross(bc);
    c + (bc * -theta.cos + m * (theta.sin * phi.cos) - n * (theta.sin * phi.sin)) * d
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use core::f32::consts::PI;

    fn dihedral(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f32 {
        let b2 = c - b;
        let n1 = (b - a).cross(b2);
        let n2 = b2.cross(d - c);
        let m1 = n1.cross(b2.normalize());
        f32::atan2(m1.dot(n2), n1.dot(n2))
    }

    const A: Vec3 = Vec3::new(0.0, 1.0, 0.0);
    const B: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    const C: Vec3 = Vec3::new(1.0, 0.0, 0.0);
    const D_BOND: f32 = 1.522; // C_3–C_3 bond length (Å)
    const COS_TET: f32 = -0.333_329_7; // cos(109.471°)
    const SIN_TET: f32 = 0.942_809_5; // sin(109.471°)

    #[test]
    fn place_bond_length() {
        for phi_deg in [0.0_f32, 60.0, 120.0, 180.0, -120.0, -60.0] {
            let phi = TrigPair {
                cos: phi_deg.to_radians().cos(),
                sin: phi_deg.to_radians().sin(),
            };
            let d = place(
                A,
                B,
                C,
                D_BOND,
                TrigPair {
                    cos: COS_TET,
                    sin: SIN_TET,
                },
                phi,
            );
            let got = (d - C).len();
            assert_relative_eq!(got, D_BOND, max_relative = 1e-5, epsilon = 1e-5);
        }
    }

    #[test]
    fn place_bond_angle() {
        for phi_deg in [0.0_f32, 60.0, 120.0, 180.0, -120.0, -60.0] {
            let phi = TrigPair {
                cos: phi_deg.to_radians().cos(),
                sin: phi_deg.to_radians().sin(),
            };
            let d = place(
                A,
                B,
                C,
                D_BOND,
                TrigPair {
                    cos: COS_TET,
                    sin: SIN_TET,
                },
                phi,
            );
            let cb_hat = (B - C).normalize();
            let cd_hat = (d - C).normalize();
            let cos_angle = cb_hat.dot(cd_hat);
            assert_relative_eq!(cos_angle, COS_TET, max_relative = 1e-5, epsilon = 1e-5);
        }
    }

    #[test]
    fn place_dihedral_roundtrip() {
        for phi_deg in [0.0_f32, 60.0, 120.0, 180.0, -120.0, -60.0, -1.0, 37.3] {
            let phi = phi_deg.to_radians();
            let d = place(
                A,
                B,
                C,
                D_BOND,
                TrigPair {
                    cos: COS_TET,
                    sin: SIN_TET,
                },
                TrigPair {
                    cos: phi.cos(),
                    sin: phi.sin(),
                },
            );
            let phi_out = dihedral(A, B, C, d);
            let diff = (phi_out - phi + 3.0 * PI).rem_euclid(2.0 * PI) - PI;
            assert!(
                diff.abs() < 2e-4,
                "dihedral round-trip failed at φ_in = {phi_deg:.1}°: \
                 got {:.4}°, diff = {:.2e} rad",
                phi_out.to_degrees(),
                diff
            );
        }
    }

    #[test]
    fn place_cis_trans() {
        let d_cis = place(
            A,
            B,
            C,
            D_BOND,
            TrigPair {
                cos: COS_TET,
                sin: SIN_TET,
            },
            TrigPair { cos: 1.0, sin: 0.0 },
        );
        let d_trans = place(
            A,
            B,
            C,
            D_BOND,
            TrigPair {
                cos: COS_TET,
                sin: SIN_TET,
            },
            TrigPair {
                cos: -1.0,
                sin: 0.0,
            },
        );

        let a_perp_y = A.y - C.y;
        let cis_perp_y = d_cis.y - C.y;
        let trans_perp_y = d_trans.y - C.y;

        assert!(
            a_perp_y * cis_perp_y > 0.0,
            "φ=0 (cis): D should be on the same y-side as A, \
             but a_perp_y={a_perp_y:.4}, cis_perp_y={cis_perp_y:.4}"
        );
        assert!(
            a_perp_y * trans_perp_y < 0.0,
            "φ=π (trans): D should be on the opposite y-side from A, \
             but a_perp_y={a_perp_y:.4}, trans_perp_y={trans_perp_y:.4}"
        );
    }

    #[test]
    fn place_ser_carbonyl_o_roundtrip() {
        let n_atom = Vec3::new(1.525, 0.493, -0.608);
        let ca = Vec3::new(0.100, 0.469, -0.252);
        let c_atom = Vec3::new(-0.053, 0.004, 1.173);
        let o_ref = Vec3::new(0.751, -0.760, 1.649);

        let d = (o_ref - c_atom).len();
        let cb_hat = (ca - c_atom).normalize();
        let cd_hat = (o_ref - c_atom).normalize();
        let cos_theta = cb_hat.dot(cd_hat);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let phi = dihedral(n_atom, ca, c_atom, o_ref);

        let o_calc = place(
            n_atom,
            ca,
            c_atom,
            d,
            TrigPair {
                cos: cos_theta,
                sin: sin_theta,
            },
            TrigPair {
                cos: phi.cos(),
                sin: phi.sin(),
            },
        );
        let err = (o_calc - o_ref).len();
        assert!(
            err < 5e-5,
            "SER carbonyl O round-trip error {err:.2e} Å exceeds 5e-5"
        );
    }

    #[test]
    fn place_ser_cb_roundtrip() {
        let c_atom = Vec3::new(-0.053, 0.004, 1.173);
        let n_atom = Vec3::new(1.525, 0.493, -0.608);
        let ca = Vec3::new(0.100, 0.469, -0.252);
        let cb_ref = Vec3::new(-0.642, -0.489, -1.184);

        let d = (cb_ref - ca).len();
        let nb_hat = (n_atom - ca).normalize();
        let cd_hat = (cb_ref - ca).normalize();
        let cos_theta = nb_hat.dot(cd_hat);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let phi = dihedral(c_atom, n_atom, ca, cb_ref);

        let cb_calc = place(
            c_atom,
            n_atom,
            ca,
            d,
            TrigPair {
                cos: cos_theta,
                sin: sin_theta,
            },
            TrigPair {
                cos: phi.cos(),
                sin: phi.sin(),
            },
        );
        let err = (cb_calc - cb_ref).len();
        assert!(
            err < 5e-5,
            "SER Cβ round-trip error {err:.2e} Å exceeds 5e-5"
        );
    }
}
