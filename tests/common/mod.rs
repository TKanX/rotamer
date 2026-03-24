#![allow(dead_code)]

use rotamer::Vec3;

pub const N: Vec3 = Vec3::new(0.0, 1.458, 0.0);
pub const CA: Vec3 = Vec3::new(0.0, 0.0, 0.0);
pub const C: Vec3 = Vec3::new(1.525, 0.0, 0.0);

#[inline]
pub fn dist(a: Vec3, b: Vec3) -> f32 {
    (b - a).len()
}

#[inline]
pub fn dihedral(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f32 {
    let b2 = c - b;
    let n1 = (b - a).cross(b2);
    let n2 = b2.cross(d - c);
    f32::atan2(b2.len() * (b - a).dot(n2), n1.dot(n2))
}

#[inline]
pub fn angle_diff(a: f32, b: f32) -> f32 {
    let d = a - b;
    (d + 3.0 * core::f32::consts::PI).rem_euclid(2.0 * core::f32::consts::PI)
        - core::f32::consts::PI
}

pub fn assert_coords_valid(atoms: &[Vec3], name: &str) {
    for (i, v) in atoms.iter().enumerate() {
        assert!(
            v.x.is_finite() && v.y.is_finite() && v.z.is_finite(),
            "{name}: atom {i} has non-finite coordinate: ({}, {}, {})",
            v.x,
            v.y,
            v.z
        );
    }
    for (i, v) in atoms.iter().enumerate() {
        let d = dist(*v, CA);
        assert!(
            d < 12.0,
            "{name}: atom {i} is {d:.2} Å from CA (exceeds 12 Å limit)"
        );
    }
}

pub fn assert_neighbors(atoms: &[Vec3], name: &str) {
    let anchors = [N, CA, C];
    for (i, v) in atoms.iter().enumerate() {
        let mut min_d = f32::INFINITY;
        for a in &anchors {
            let d = dist(*v, *a);
            if d < min_d {
                min_d = d;
            }
        }
        for (j, u) in atoms.iter().enumerate() {
            if i == j {
                continue;
            }
            let d = dist(*v, *u);
            if d < min_d {
                min_d = d;
            }
        }
        assert!(
            min_d < 2.1,
            "{name}: atom {i} has no neighbor within 2.1 Å (nearest = {min_d:.4} Å)"
        );
    }
}

#[inline]
pub fn chirality_sign(cb: Vec3) -> f32 {
    (N - CA).cross(C - CA).dot(cb - CA)
}
