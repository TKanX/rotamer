use crate::math::Vec3;

/// Places atom D given reference frame A–B–C.
///
/// # Geometry
///
/// Returns the unique point D such that:
///
/// - `|C − D| = d`
/// - `∠BCD = arccos(cos_theta)` where `sin_theta = sin(∠BCD) > 0`
/// - `dihedral(A, B, C, D) = arctan2(sin_phi, cos_phi)`
///
/// # Preconditions
///
/// - A, B, C must not be collinear (zero cross product → undefined frame).
/// - `sin_theta > 0` (the bond angle must not be 0° or 180°).
/// - `d > 0`.
///
/// Violating any precondition yields an unspecified but defined (non-UB) result.
#[inline(always)]
pub fn place(
    a: Vec3,
    b: Vec3,
    c: Vec3,
    d: f32,
    cos_theta: f32,
    sin_theta: f32,
    cos_phi: f32,
    sin_phi: f32,
) -> Vec3 {
    let bc = (c - b).normalize();
    let n = (b - a).cross(bc).normalize();
    let m = n.cross(bc);
    c + (bc * -cos_theta + m * (sin_theta * cos_phi) - n * (sin_theta * sin_phi)) * d
}
