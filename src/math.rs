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
