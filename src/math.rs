use core::f32::consts::{FRAC_2_PI, FRAC_PI_2};

/// Computes `1 / √x` via the Quake III bit-cast trick followed by three
/// Newton–Raphson refinement steps.
#[inline(always)]
fn rsqrtf(x: f32) -> f32 {
    const MAGIC: u32 = 0x5f375a86;
    let x2 = x * 0.5;
    let y = f32::from_bits(MAGIC - (x.to_bits() >> 1));
    let y = y * (1.5 - x2 * y * y);
    let y = y * (1.5 - x2 * y * y);
    y * (1.5 - x2 * y * y)
}

/// Computes `(sin(x), cos(x))` simultaneously for any finite `x` (radians).
///
/// # Algorithm
///
/// 1. Round to nearest quadrant index `j`, reduce `r = x − j·π/2`.
/// 2. Evaluate degree-9 (sin) and degree-8 (cos) minimax polynomials in `r²`.
/// 3. Swap and negate via branchless bit manipulation on the quadrant index.
///
/// Maximum absolute error over `[−2π, 2π]` is below `5 × 10⁻⁷`.
#[inline(always)]
pub fn sincosf(x: f32) -> (f32, f32) {
    let t = x * FRAC_2_PI;
    let j = (t + f32::copysign(0.5, t)) as i32;
    let r = x - j as f32 * FRAC_PI_2;

    let u = r * r;
    #[rustfmt::skip]
    let sin_r = r * (1.0 + u * (-1.666_666_7e-1
                         + u * ( 8.333_334e-3
                         + u * (-1.984_127e-4
                         + u *   2.755_731_7e-6))));
    #[rustfmt::skip]
    let cos_r =     1.0 + u * (-5.0e-1
                         + u * ( 4.166_666_7e-2
                         + u * (-1.388_888_9e-3
                         + u *   2.480_158_7e-5)));

    let j_u = j as u32;
    let swap = j_u & 1;
    let s_neg = (j_u >> 1) & 1;
    let c_neg = ((j_u.wrapping_add(1)) >> 1) & 1;

    let (a, b) = (sin_r.to_bits(), cos_r.to_bits());
    let xab = a ^ b;
    let sel = 0u32.wrapping_sub(swap);
    let s_bits = a ^ (sel & xab);
    let c_bits = b ^ (sel & xab);

    const SIGN: u32 = 0x8000_0000;
    let sin_out = f32::from_bits(s_bits ^ (0u32.wrapping_sub(s_neg) & SIGN));
    let cos_out = f32::from_bits(c_bits ^ (0u32.wrapping_sub(c_neg) & SIGN));

    (sin_out, cos_out)
}

/// A three-component single-precision floating-point vector.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    /// Construct a vector from its three components.
    #[inline(always)]
    pub const fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    /// The zero vector `(0, 0, 0)`.
    #[inline(always)]
    pub const fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Euclidean inner product `self · other`.
    #[inline(always)]
    pub fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Euclidean cross product `self × other` (right-hand rule).
    #[inline(always)]
    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Squared Euclidean length `|self|²`.
    #[inline(always)]
    pub fn len_sq(self) -> f32 {
        self.dot(self)
    }

    /// Euclidean length `|self|`.
    ///
    /// # Precondition
    ///
    /// `self` must not be the zero vector. Passing the zero vector yields an
    /// unspecified (but not UB) result.
    #[inline(always)]
    pub fn len(self) -> f32 {
        let q = self.len_sq();
        q * rsqrtf(q)
    }

    /// Return the unit vector in the direction of `self`.
    ///
    /// # Precondition
    ///
    /// `self` must not be the zero vector. Passing the zero vector yields an
    /// unspecified (but not UB) result.
    #[inline(always)]
    pub fn normalize(self) -> Self {
        self * rsqrtf(self.len_sq())
    }
}

impl core::ops::Add for Vec3 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl core::ops::Sub for Vec3 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl core::ops::Neg for Vec3 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl core::ops::Mul<f32> for Vec3 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, s: f32) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn rsqrtf_one() {
        assert_relative_eq!(rsqrtf(1.0), 1.0, max_relative = 1.5e-7);
    }

    #[test]
    fn rsqrtf_four() {
        assert_relative_eq!(rsqrtf(4.0), 0.5, max_relative = 1.5e-7);
    }

    #[test]
    fn rsqrtf_nine() {
        assert_relative_eq!(rsqrtf(9.0), 1.0 / 3.0, max_relative = 1.5e-7);
    }

    #[test]
    fn rsqrtf_max_error() {
        let mut max_rel_err = 0.0_f64;
        for exp in -60_i32..=60 {
            for k in 0_u32..1000 {
                let x = (2.0_f64).powi(exp) * (1.0 + k as f64 / 1000.0);
                let x_f = x as f32;
                if x_f <= 0.0 || !x_f.is_finite() {
                    continue;
                }
                let got = rsqrtf(x_f) as f64;
                let expected = 1.0 / x.sqrt();
                let rel = (got - expected).abs() / expected;
                if rel > max_rel_err {
                    max_rel_err = rel;
                }
            }
        }
        assert!(
            max_rel_err < 1.5e-7,
            "max rsqrtf relative error {max_rel_err:.2e} exceeds 1.5e-7"
        );
    }

    #[test]
    fn sincosf_zero() {
        let (s, c) = sincosf(0.0);
        assert_relative_eq!(s, 0.0, max_relative = 5e-7, epsilon = 1e-7);
        assert_relative_eq!(c, 1.0, max_relative = 5e-7);
    }

    #[test]
    fn sincosf_pi_over_6() {
        let (s, c) = sincosf(core::f32::consts::FRAC_PI_6);
        assert_relative_eq!(s, 0.5, max_relative = 5e-7);
        assert_relative_eq!(c, 3.0_f32.sqrt() / 2.0, max_relative = 5e-7);
    }

    #[test]
    fn sincosf_pi_over_4() {
        let (s, c) = sincosf(core::f32::consts::FRAC_PI_4);
        let root2_over2 = core::f32::consts::FRAC_1_SQRT_2;
        assert_relative_eq!(s, root2_over2, max_relative = 5e-7);
        assert_relative_eq!(c, root2_over2, max_relative = 5e-7);
    }

    #[test]
    fn sincosf_pi_over_2() {
        let (s, c) = sincosf(FRAC_PI_2);
        assert_relative_eq!(s, 1.0, max_relative = 5e-7);
        assert_relative_eq!(c, 0.0, max_relative = 5e-7, epsilon = 1e-7);
    }

    #[test]
    fn sincosf_pi() {
        let (s, c) = sincosf(core::f32::consts::PI);
        assert_relative_eq!(s, 0.0, max_relative = 5e-7, epsilon = 1e-7);
        assert_relative_eq!(c, -1.0, max_relative = 5e-7);
    }

    #[test]
    fn sincosf_3pi_over_2() {
        let (s, c) = sincosf(3.0 * FRAC_PI_2);
        assert_relative_eq!(s, -1.0, max_relative = 5e-7);
        assert_relative_eq!(c, 0.0, max_relative = 5e-7, epsilon = 1e-7);
    }

    #[test]
    fn sincosf_sin_odd() {
        for &x in &[0.3_f32, 1.0, 2.5, core::f32::consts::PI] {
            let (s_pos, _) = sincosf(x);
            let (s_neg, _) = sincosf(-x);
            assert_relative_eq!(s_neg, -s_pos, max_relative = 1e-6, epsilon = 1e-7);
        }
    }

    #[test]
    fn sincosf_cos_even() {
        for &x in &[0.3_f32, 1.0, 2.5, core::f32::consts::PI] {
            let (_, c_pos) = sincosf(x);
            let (_, c_neg) = sincosf(-x);
            assert_relative_eq!(c_neg, c_pos, max_relative = 1e-6, epsilon = 1e-7);
        }
    }

    #[test]
    fn sincosf_pythagorean_identity() {
        let test_angles: &[f32] = &[
            0.0,
            0.1,
            0.5,
            1.0,
            1.5,
            2.0,
            2.5,
            3.0,
            -0.1,
            -1.0,
            -2.0,
            -3.0,
            FRAC_PI_2,
            core::f32::consts::PI,
            3.0 * FRAC_PI_2,
        ];
        for &x in test_angles {
            let (s, c) = sincosf(x);
            let norm_sq = s * s + c * c;
            assert!(
                (norm_sq - 1.0).abs() < 1e-6,
                "sin²+cos² = {norm_sq:.8} at x = {x:.4}, deviation {:.2e}",
                (norm_sq - 1.0).abs()
            );
        }
    }

    #[test]
    fn sincosf_quadrant_signs() {
        let (s, c) = sincosf(core::f32::consts::FRAC_PI_4);
        assert!(s > 0.0 && c > 0.0);
        let (s, c) = sincosf(3.0 * core::f32::consts::FRAC_PI_4);
        assert!(s > 0.0 && c < 0.0);
        let (s, c) = sincosf(-3.0 * core::f32::consts::FRAC_PI_4);
        assert!(s < 0.0 && c < 0.0);
        let (s, c) = sincosf(-core::f32::consts::FRAC_PI_4);
        assert!(s < 0.0 && c > 0.0);
    }

    #[test]
    fn sincosf_max_error() {
        let mut max_err_sin = 0.0_f32;
        let mut max_err_cos = 0.0_f32;
        let n = 100_000_u32;
        let range = 2.0 * core::f32::consts::PI;
        for i in 0..=n {
            let x = -range + 2.0 * range * i as f32 / n as f32;
            let (s, c) = sincosf(x);
            let err_s = (s - x.sin()).abs();
            let err_c = (c - x.cos()).abs();
            if err_s > max_err_sin {
                max_err_sin = err_s;
            }
            if err_c > max_err_cos {
                max_err_cos = err_c;
            }
        }
        assert!(
            max_err_sin < 5e-7,
            "max sin error {max_err_sin:.2e} over [-2π,2π] exceeds 5e-7"
        );
        assert!(
            max_err_cos < 5e-7,
            "max cos error {max_err_cos:.2e} over [-2π,2π] exceeds 5e-7"
        );
    }

    #[test]
    fn vec3_add() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let c = a + b;
        assert_eq!(c, Vec3::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn vec3_sub() {
        let a = Vec3::new(5.0, 7.0, 9.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_eq!(a - b, Vec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn vec3_neg() {
        let v = Vec3::new(1.0, -2.0, 3.0);
        assert_eq!(-v, Vec3::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn vec3_mul_scalar() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v * 2.0, Vec3::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn vec3_zero() {
        const Z: Vec3 = Vec3::zero();
        assert_eq!(Z, Vec3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn vec3_dot_perpendicular() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.dot(y), 0.0);
    }

    #[test]
    fn vec3_dot_self() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert_relative_eq!(v.dot(v), v.len_sq(), max_relative = 1e-6);
    }

    #[test]
    fn vec3_cross_basis() {
        let ex = Vec3::new(1.0, 0.0, 0.0);
        let ey = Vec3::new(0.0, 1.0, 0.0);
        let ez = Vec3::new(0.0, 0.0, 1.0);
        assert_eq!(ex.cross(ey), ez);
        assert_eq!(ey.cross(ez), ex);
        assert_eq!(ez.cross(ex), ey);
    }

    #[test]
    fn vec3_cross_anticommutative() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let ab = a.cross(b);
        let ba = b.cross(a);
        assert_relative_eq!(ab.x, -ba.x, max_relative = 1e-6, epsilon = 1e-7);
        assert_relative_eq!(ab.y, -ba.y, max_relative = 1e-6, epsilon = 1e-7);
        assert_relative_eq!(ab.z, -ba.z, max_relative = 1e-6, epsilon = 1e-7);
    }

    #[test]
    fn vec3_cross_parallel_is_zero() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let c = v.cross(v * 3.7);
        assert!(c.len_sq() < 1e-10, "cross(v, kv) = {c:?} should be zero");
    }

    #[test]
    fn vec3_cross_perpendicular_to_inputs() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-1.0, 4.0, 0.5);
        let c = a.cross(b);
        assert!(
            a.dot(c).abs() < 1e-5,
            "a · (a×b) = {} should be 0",
            a.dot(c)
        );
        assert!(
            b.dot(c).abs() < 1e-5,
            "b · (a×b) = {} should be 0",
            b.dot(c)
        );
    }

    #[test]
    fn vec3_len_345() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert_relative_eq!(v.len(), 5.0, max_relative = 1e-6);
    }

    #[test]
    fn vec3_normalize_unit_length() {
        let v = Vec3::new(3.0, -7.5, 2.1);
        let n = v.normalize();
        assert_relative_eq!(n.len(), 1.0, max_relative = 1.5e-7);
    }

    #[test]
    fn vec3_normalize_direction_preserved() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let n = v.normalize();
        assert!(v.dot(n) > 0.0);
    }

    #[test]
    fn vec3_repr_c_layout() {
        assert_eq!(
            core::mem::size_of::<Vec3>(),
            3 * core::mem::size_of::<f32>()
        );
        assert_eq!(core::mem::align_of::<Vec3>(), core::mem::align_of::<f32>());
    }

    #[test]
    fn vec3_f32_slice_reinterpret() {
        let vs = [Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 5.0, 6.0)];
        let flat: &[f32] = unsafe { core::slice::from_raw_parts(vs.as_ptr() as *const f32, 6) };
        assert_eq!(flat, &[1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }
}
