// Copyright (c) 2024-2025 Ken Barker

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

//! The geodesic module contains functions for calculating the geodesic segment
//! between two points on the surface of an ellipsoid.
//!
//! See CFF Karney: [Algorithms for geodesics](https://arxiv.org/pdf/1109.4448.pdf).

#![allow(clippy::float_cmp)]
#![allow(clippy::many_single_char_names)]

use crate::{ellipsoid, Ellipsoid, Metres};
use angle_sc::{trig, trig::UnitNegRange, Angle, Radians};
use unit_sphere::{great_circle, LatLong};

/// Estimate omega12 by solving the astroid problem.
/// Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
/// See CFF Karney section 5.
/// * `x`, `y` - astroid parameter
///
/// returns the solution to the astroid problem.
#[must_use]
fn calculate_astroid(x: f64, y: f64) -> f64 {
    let p = x * x;
    let q = y * y;
    let r = (p + q - 1.0) / 6.0;

    // y = 0 with |x| <= 1
    // for y small, positive root is k = abs(y)/sqrt(1-x^2)
    if (q <= 0.0) && (r <= 0.0) {
        0.0
    } else {
        let s = p * q / 4.0;
        let r2 = r * r;
        let r3 = r * r2;
        let mut u = r;

        // The discriminant of the quadratic equation for T3.
        // This is zero on the evolute curve p^(1/3)+q^(1/3) = 1
        let discriminant = s * (s + 2.0 * r3);
        if 0.0 <= discriminant {
            let mut t3 = s + r3;
            // Pick the sign on the sqrt to maximize abs(T3), to minimise loss
            // of precision due to cancellation.
            t3 += if t3 < 0.0 {
                -libm::sqrt(discriminant)
            } else {
                libm::sqrt(discriminant)
            };
            let t = libm::cbrt(t3);
            u += if t == 0.0 { 0.0 } else { t + r2 / t };
        } else {
            // T is complex, but the way u is defined the result is real.
            let angle = libm::atan2(libm::sqrt(-discriminant), -(s + r3));
            // There are three possible cube roots.  We choose the root which
            // avoids cancellation.  Note: discriminant < 0 implies that r < 0.
            u += 2.0 * r * libm::cos(angle / 3.0);
        }

        let v = libm::sqrt(u * u + q); // guaranteed positive
        let uv = if u < 0.0 { q / (v - u) } else { u + v }; // u+v, guaranteed positive
        let w = (uv - q) / (2.0 * v); // positive?

        // Rearrange expression for k to avoid loss of accuracy due to subtraction.
        // Division by 0 not possible because uv > 0, w >= 0.
        uv / (libm::sqrt(uv + w * w) + w) // guaranteed positive
    }
}

/// Calculate: m12b = (reduced length)/_b
#[allow(clippy::similar_names)]
#[must_use]
fn calculate_reduced_length(
    eps: f64,
    sigma12: Radians,
    sigma1: Angle,
    dn1: f64,
    sigma2: Angle,
    dn2: f64,
) -> f64 {
    let a1 = ellipsoid::coefficients::evaluate_a1(eps);
    let a2 = ellipsoid::coefficients::evaluate_a2(eps);
    let m0x = a1 - a2;

    let a1p1 = 1.0 + a1;
    let a2p1 = 1.0 + a2;

    let ca = ellipsoid::coefficients::evaluate_coeffs_c1(eps);
    let mut cb = ellipsoid::coefficients::evaluate_coeffs_c2(eps);

    // Assume here that ca.len() >= cb.len()
    for i in 1..cb.len() {
        cb[i] = a1p1 * ca[i] - a2p1 * cb[i];
    }

    let j12 = m0x
        * (sigma12
            + (ellipsoid::coefficients::sin_cos_series(&cb, sigma2)
                - ellipsoid::coefficients::sin_cos_series(&cb, sigma1)))
        .0;
    dn2 * (sigma1.cos().0 * sigma2.sin().0)
        - dn1 * (sigma1.sin().0 * sigma2.cos().0)
        - sigma1.cos().0 * sigma2.cos().0 * j12
}

/// Estimate the initial azimuth on the auxiliary sphere for a nearly antipodal arc.
/// It calculates and solves the astroid problem.
/// * `beta1`, `beta2` - the parametric latitudes of the start and finish points
///   on the auxiliary sphere.
/// * `lambda12` - Longitude difference between start and finish points.
///
/// returns the estimate of the initial azimuth on the auxiliary sphere.
#[must_use]
fn estimate_antipodal_initial_azimuth(
    beta1: Angle,
    beta2: Angle,
    lambda12: Angle,
    ellipsoid: &Ellipsoid,
) -> Angle {
    const Y_TOLERANCE: f64 = 200.0 * f64::EPSILON;
    let x_threshold: f64 = 1000.0 * libm::sqrt(f64::EPSILON);

    // Calculate the integration parameter for geodesic
    let clairaut = beta1.cos(); // Note: assumes sin_alpha_1 = 1
    let eps = ellipsoid.calculate_epsilon(clairaut);
    let a3f = ellipsoid.calculate_a3f(eps);

    let lamscale = ellipsoid.f() * beta1.cos().0 * a3f * core::f64::consts::PI;
    let betscale = lamscale * beta1.cos().0;

    // Solve astroid problem
    let x = Radians::from(lambda12.opposite()).0 / lamscale;
    let y = trig::sine_sum(beta1.sin(), beta1.cos(), beta2.sin(), beta2.cos()).0 / betscale;

    // Test x and y params
    if (y > -Y_TOLERANCE) && (x > -1.0 - x_threshold) {
        let sin_alpha = UnitNegRange::clamp(-x);
        Angle::new(sin_alpha, trig::swap_sin_cos(sin_alpha)).negate_cos()
    } else {
        let k = calculate_astroid(x, y);
        let omg12a = lamscale * (-x * k / (1.0 + k));

        let omega12 = Radians(core::f64::consts::PI - omg12a);
        great_circle::calculate_gc_azimuth(beta1, beta2, Angle::from(omega12))
    }
}

/// Calculate the cosine of the longitude difference from the equator crossing.
/// * `beta` the `parametric` latitude
/// * `cos_azimuth` the cosine of the azimuth at the `parametric` latitude
///
/// returns the cosine of the longitude difference
#[must_use]
pub fn calculate_cos_omega(beta: Angle, cos_azimuth: UnitNegRange) -> UnitNegRange {
    UnitNegRange(cos_azimuth.0 * beta.cos().0)
}

/// Calculate the azimuth on the auxiliary sphere at `parametric` latitude
/// beta2 given the `parametric` latitude beta1 and azimuth, `alpha1`.
/// * `beta1`, `beta2` - the parametric latitudes of the start and finish points
///   on the auxiliary sphere.
/// * `alpha1` - start point azimuth.
///
/// returns the finish point azimuth.
#[must_use]
fn calculate_end_azimuth(beta1: Angle, beta2: Angle, alpha1: Angle) -> Angle {
    let clairaut = UnitNegRange(alpha1.sin().0 * beta1.cos().0);

    let sin_alpha2 = if beta2.cos() == beta1.cos() {
        alpha1.sin()
    } else {
        UnitNegRange::clamp(clairaut.0 / beta2.cos().0)
    };

    // Karney's method to calculate the cosine of the end azimuth
    let cos_alpha2 = if (beta2.cos() != beta1.cos()) || (beta2.sin().abs().0 != -beta1.sin().0) {
        let temp1 = alpha1.cos().0 * beta1.cos().0;
        let temp2 = if beta1.cos().0 < beta1.abs().sin().0 {
            (beta2.cos().0 - beta1.cos().0) * (beta1.cos().0 + beta2.cos().0)
        } else {
            (beta1.sin().0 - beta2.sin().0) * (beta1.sin().0 + beta2.sin().0)
        };
        let temp3 = temp1 * temp1 + temp2;
        let temp4 = if 0.0 < temp3 {
            libm::sqrt(temp3) / beta2.cos().0
        } else {
            0.0
        };
        UnitNegRange::clamp(temp4)
    } else {
        alpha1.cos().abs()
    };

    Angle::new(sin_alpha2, cos_alpha2)
}

/// Calculate the longitude difference between the auxiliary sphere and
/// ellipsoid.
#[allow(clippy::similar_names)]
#[must_use]
fn delta_omega12(
    clairaut: UnitNegRange,
    eps: f64,
    sigma12: Radians,
    sigma1: Angle,
    sigma2: Angle,
    ellipsoid: &Ellipsoid,
) -> Radians {
    let c3 = ellipsoid.calculate_c3y(eps);
    let b31 = ellipsoid::coefficients::sin_cos_series(&c3, sigma1);
    let b32 = ellipsoid::coefficients::sin_cos_series(&c3, sigma2);

    let a3c = ellipsoid.calculate_a3c(clairaut, eps);
    Radians(a3c * (sigma12 + (b32 - b31)).0)
}

/// Estimate the initial azimuth on the auxiliary sphere for normal arcs,
/// i.e. NOT nearly antipodal points.
///
/// * `beta1`, `beta2` the parametric latitudes on the auxiliary sphere.
/// * `lambda12` longitude difference between points.
/// * `ellipsoid` - the `Ellipsoid`.
///
/// @return an estimate of the initial azimuth on the auxiliary sphere.
#[allow(clippy::similar_names)]
#[must_use]
fn estimate_initial_azimuth(
    beta1: Angle,
    beta2: Angle,
    lambda12: Angle,
    ellipsoid: &Ellipsoid,
) -> Angle {
    // Calculate the great circle azimuth between parametric latitudes
    let alpha1 = great_circle::calculate_gc_azimuth(beta1, beta2, lambda12);

    // Calculate Clairaut's constant
    let clairaut = UnitNegRange(alpha1.sin().0 * beta1.cos().0);
    let eps = ellipsoid.calculate_epsilon(clairaut);

    // Calculate sigma1
    let cos_omega1 = calculate_cos_omega(beta1, alpha1.cos());
    let sigma1 = Angle::from_y_x(beta1.sin().0, cos_omega1.0);

    // Calculate sigma2
    let alpha2 = calculate_end_azimuth(beta1, beta2, alpha1);
    let cos_omega2 = calculate_cos_omega(beta2, alpha2.cos());
    let sigma2 = Angle::from_y_x(beta2.sin().0, cos_omega2.0);

    // Calculate omega12
    let sigma12 = Radians::from(sigma2) - Radians::from(sigma1);
    let domg12 = delta_omega12(clairaut, eps, sigma12, sigma1, sigma2, ellipsoid);
    let omega12 = lambda12 + Angle::from(domg12);

    // Recalculate the great circle azimuth using the longitude difference
    // estimate: omega12
    great_circle::calculate_gc_azimuth(beta1, beta2, omega12)
}

/// Find the azimuth and great circle length on the auxiliary sphere.
///
/// It uses Newton's method to solve:
///   f(alp1) = lambda12(alp1) - lam12 = 0
///
/// * `beta1`, `beta2` - the start and end `parametric` latitudes.
/// * `abs_lambda12` - Longitude difference between start and finish points.
/// * `alpha` - the initial azimuth.
/// * `gc_length` - the auxiliary sphere great circle length.
/// * `ellipsoid` - the `Ellipsoid`.
/// * `tolerance` - the tolerance to perform the calculation to in Radians.
///
/// returns the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic segment and the number of iterations required to
/// calculate them.
#[allow(clippy::similar_names)]
#[must_use]
fn find_azimuth_length_newtons_method(
    beta1: Angle,
    beta2: Angle,
    abs_lambda12: Angle,
    alpha: Angle,
    gc_length: Radians,
    ellipsoid: &Ellipsoid,
    tolerance: Radians,
) -> (Angle, Radians, u32) {
    // The first iteration threshold
    const MAX_ITER1: u32 = 20;
    // The maximum number of iterations to attempt.
    const MAX_ITERS: u32 = MAX_ITER1 + f64::DIGITS + 10;

    let dn1 = libm::sqrt(1.0 + ellipsoid.ep_2() * beta1.sin().0 * beta1.sin().0);
    let dn2 = libm::sqrt(1.0 + ellipsoid.ep_2() * beta2.sin().0 * beta2.sin().0);

    let mut alpha1 = alpha;
    let mut sigma12_rad = gc_length;
    let mut iterations: u32 = 1;
    for i in 1..=MAX_ITERS {
        iterations = i;
        // Calculate Clairaut's constant
        let clairaut = UnitNegRange(alpha1.sin().0 * beta1.cos().0);
        let eps = ellipsoid.calculate_epsilon(clairaut);

        // Calculate first longitude (omega1) and distance (sigma1) from the
        // Northbound equator crossing
        let sin_omega1 = UnitNegRange(clairaut.0 * beta1.sin().0);
        let cos_omega1 = calculate_cos_omega(beta1, alpha1.cos());
        let sigma1 = Angle::from_y_x(beta1.sin().0, cos_omega1.0);
        let omega1 = Angle::new(sin_omega1, cos_omega1);

        // Calculate azimuth at the end point
        let alpha2 = calculate_end_azimuth(beta1, beta2, alpha1);

        // Calculate second longitude (omega2) and distance (sigma2) from the
        // Northbound equator crossing
        let sin_omega2 = UnitNegRange(clairaut.0 * beta2.sin().0);
        let cos_omega2 = calculate_cos_omega(beta2, alpha2.cos());
        let sigma2 = Angle::from_y_x(beta2.sin().0, cos_omega2.0);
        let omega2 = Angle::new(sin_omega2, cos_omega2);

        // Calculate great circle length on the auxiliary sphere
        let sigma12 = sigma2 - sigma1;
        // clamp to range 0 to Pi
        sigma12_rad = if sigma12.sin().0 > 0.0 {
            Radians::from(sigma12)
        } else {
            Radians(libm::atan2(0.0, sigma12.cos().0))
        };
        let domg12 = delta_omega12(clairaut, eps, sigma12_rad, sigma1, sigma2, ellipsoid);

        // Calculate Longitude difference on the auxiliary sphere
        let mut omega12 = omega2 - omega1;
        // clamp to range 0 to Pi
        if omega12.sin().0.is_sign_negative() {
            omega12 = Angle::from_y_x(0.0, omega12.cos().0);
        }
        let eta = Radians::from(omega12 - abs_lambda12);

        // Difference between differences
        let v = eta.0 - domg12.0;
        if libm::fabs(v) < tolerance.0 {
            break;
        }

        // Calculate the denominator for Newton's method
        let dv = if alpha2.cos().0 == 0.0 {
            -2.0 * ellipsoid.one_minus_f() * dn1 / beta1.sin().0
        } else {
            let m12 = calculate_reduced_length(eps, sigma12_rad, sigma1, dn1, sigma2, dn2);
            ellipsoid.one_minus_f() * m12 / (alpha2.cos().0 * beta2.cos().0)
        };
        if libm::fabs(dv) < tolerance.0 {
            break;
        }

        // Calculate the change in initial azimuth
        let dalpha1 = UnitNegRange::clamp(-v / dv);
        if dalpha1.abs().0 <= great_circle::MIN_VALUE {
            break;
        }

        // Adjust the azimuth by dalpha1
        alpha1 += Angle::from(Radians(dalpha1.0));
    }

    (alpha1, sigma12_rad, iterations)
}

/// Find the azimuths and great circle length on the auxiliary sphere.
///
/// It adjusts the latitudes and longitude difference so that the azimuth of the
/// geodesic segment lies between 0° and 90°.
/// It calls `find_azimuth_length_newtons_method` and then changes the resulting
/// azimuth to match the orientation of the geodesic segment.
///
/// * `beta_a`, `beta_b` - the `parametric` latitudes of the start and finish points.
/// * `lambda12` - Longitude difference between start and finish points.
/// * `gc_length` - the auxiliary sphere great circle length.
/// * `tolerance` - the tolerance to perform the calculation to in Radians.
/// * `ellipsoid` - the `Ellipsoid`.
///
/// returns the azimuth at the start of the geodesic segment, the great circle
/// arc length on the auxiliary sphere, the azimuth at the end of
/// geodesic segment and the number of iterations required to calculate them.
#[allow(clippy::similar_names)]
#[must_use]
fn find_azimuths_and_arc_length(
    beta_a: Angle,
    beta_b: Angle,
    lambda12: Angle,
    gc_length: Radians,
    tolerance: Radians,
    ellipsoid: &Ellipsoid,
) -> (Angle, Radians, Angle, u32) {
    let antipodal_arc_threshold: f64 = core::f64::consts::PI * ellipsoid.one_minus_f();

    // Start at the latitude furthest from the Equator,
    let abs_beta_a = beta_a.abs();
    let abs_beta_b = beta_b.abs();
    let swap_latitudes = (abs_beta_a.sin() < abs_beta_b.sin())
        || ((abs_beta_a.sin() == abs_beta_b.sin()) && (abs_beta_a.cos() > abs_beta_b.cos()));
    let mut beta1 = if swap_latitudes { beta_b } else { beta_a };
    let mut beta2 = if swap_latitudes { beta_a } else { beta_b };

    // Start South of the Equator
    let negate_latitude = beta1.sin().0.is_sign_positive();
    if negate_latitude {
        beta1 = -beta1;
        beta2 = -beta2;
    }

    // Use positive lambda12, so all azimuths are positive
    let abs_lambda12 = lambda12.abs();

    // Estimate the azimuth at the start of the geodesic segment
    let alpha0 = if gc_length.0 >= antipodal_arc_threshold {
        estimate_antipodal_initial_azimuth(beta1, beta2, abs_lambda12, ellipsoid)
    } else {
        estimate_initial_azimuth(beta1, beta2, abs_lambda12, ellipsoid)
    };

    // Use Newton's method to calculate the initial azimuth and aux length
    let (alpha, sigma12, iterations) = find_azimuth_length_newtons_method(
        beta1,
        beta2,
        abs_lambda12,
        alpha0,
        gc_length,
        ellipsoid,
        tolerance,
    );

    // Calculate the correct azimuths for the start and finish points
    let mut alpha1 = alpha;
    let mut alpha2 = calculate_end_azimuth(beta1, beta2, alpha1);
    if swap_latitudes {
        alpha1 = alpha2;
        alpha2 = alpha;
    }

    if swap_latitudes != negate_latitude {
        alpha1 = alpha1.negate_cos();
        alpha2 = alpha2.negate_cos();
    }

    if lambda12.sin().0.is_sign_negative() {
        alpha1 = -alpha1;
        alpha2 = -alpha2;
    }

    (alpha1, sigma12, alpha2, iterations)
}

/// Calculate the initial azimuth and great circle length between a pair
/// of points on the auxiliary sphere.
/// * `beta1`, `beta2` - the `parametric` latitudes of the start and finish
///   points on the auxiliary sphere.
/// * `delta_long` - the longitude difference on the auxiliary sphere.
/// * `tolerance` - the tolerance to perform the calculation to.
/// * `ellipsoid` - the `Ellipsoid`.
///
/// returns the azimuth at the start of the geodesic segment, the great circle
/// arc length on the auxiliary sphere, the azimuth at the end of
/// geodesic segment and the number of iterations required to calculate them.
#[must_use]
pub fn aux_sphere_azimuths_length(
    beta1: Angle,
    beta2: Angle,
    delta_long: Angle,
    tolerance: Radians,
    ellipsoid: &Ellipsoid,
) -> (Angle, Radians, Angle, u32) {
    let max_equatorial_length = Radians(core::f64::consts::PI * ellipsoid.one_minus_f());

    let gc_azimuth = great_circle::calculate_gc_azimuth(beta1, beta2, delta_long);
    let gc_length = great_circle::calculate_gc_distance(beta1, beta2, delta_long);

    // Determine whether on a meridian, i.e. a great circle which passes through the North and South poles
    if gc_azimuth.abs().sin().0 < great_circle::MIN_VALUE {
        // gc_azimuth is 0° or 180°

        // Use opposite azimuth if points on opposite meridians
        let end_azimuth = if delta_long.cos().0 < 0.0 {
            gc_azimuth.opposite()
        } else {
            gc_azimuth
        };
        (gc_azimuth, gc_length, end_azimuth, 0)
    } else {
        // Determine whether on an equatorial path, i.e. the circle around the equator.
        // gc_azimuth is +/-90°, the ends are NOT nearly antipodal
        // and both latitudes are very close to the equator
        if (gc_azimuth.cos().0 < f64::EPSILON)
            && (gc_length < max_equatorial_length)
            && (beta1.abs().sin().0 < f64::EPSILON)
            && (beta2.abs().sin().0 < f64::EPSILON)
        {
            // Calculate the distance around the equator on the axillary sphere
            let equatorial_length = Radians(gc_length.0 * ellipsoid.recip_one_minus_f());
            (gc_azimuth, equatorial_length, gc_azimuth, 0)
        } else {
            // Iterate to find the azimuth and length on the axillary sphere
            find_azimuths_and_arc_length(beta1, beta2, delta_long, gc_length, tolerance, ellipsoid)
        }
    }
}

/// Calculate the `geodesic` azimuth and great circle length on the auxiliary sphere
/// between a pair of positions.
/// * `a`, `b` - the start and finish positions in geodetic coordinates.
/// * `tolerance` - the tolerance to perform the calculation to.
/// * `ellipsoid` - the `Ellipsoid`.
///
/// returns the azimuth at the start of the geodesic segment, the great circle
/// arc length on the auxiliary sphere, the azimuth at the end of
/// geodesic segment and the number of iterations required to calculate them.
#[must_use]
pub fn calculate_azimuths_aux_length(
    a: &LatLong,
    b: &LatLong,
    tolerance: Radians,
    ellipsoid: &Ellipsoid,
) -> (Angle, Radians, Angle, u32) {
    // calculate the parametric latitudes on the auxiliary sphere
    let beta_a = ellipsoid.calculate_parametric_latitude(Angle::from(a.lat()));
    let beta_b = ellipsoid.calculate_parametric_latitude(Angle::from(b.lat()));

    // calculate the longitude difference
    let delta_long = Angle::from((b.lon(), a.lon()));
    aux_sphere_azimuths_length(beta_a, beta_b, delta_long, tolerance, ellipsoid)
}

/// Convert a great circle distance on the auxiliary sphere in radians to
/// metres on the ellipsoid.
/// * `beta1`, the start parametric Latitude on the auxiliary sphere.
/// * `alpha1`, the azimuth at the start point.
/// * `gc_distance`, the great circle distance on the auxiliary sphere in radians.
/// * `ellipsoid` - the `Ellipsoid`.
///
/// returns the geodesic distance in metres.
#[must_use]
pub fn convert_radians_to_metres(
    beta1: Angle,
    alpha1: Angle,
    gc_distance: Radians,
    ellipsoid: &Ellipsoid,
) -> Metres {
    // Calculate the distance from the first equator crossing
    let cos_omega1 = calculate_cos_omega(beta1, alpha1.cos());
    let sigma1 = Angle::from_y_x(beta1.sin().0, cos_omega1.0);
    let sigma_sum = sigma1 + Angle::from(gc_distance);

    // Calculate the ellipsoid coefficients
    let clairaut = UnitNegRange(alpha1.sin().0 * beta1.cos().0);
    let eps = ellipsoid.calculate_epsilon(clairaut);
    let a1 = ellipsoid::coefficients::evaluate_a1(eps) + 1.0;
    let c1 = ellipsoid::coefficients::evaluate_coeffs_c1(eps);
    let b11 = ellipsoid::coefficients::sin_cos_series(&c1, sigma1);
    let b12 = ellipsoid::coefficients::sin_cos_series(&c1, sigma_sum);

    Metres(ellipsoid.b().0 * a1 * (gc_distance + b12 - b11).0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeodesicSegment, WGS84_ELLIPSOID};
    use angle_sc::{is_within_tolerance, Degrees};

    #[test]
    fn test_calculate_astroid() {
        const Y_TOLERANCE: f64 = 200.0 * f64::EPSILON;
        let x_threshold: f64 = 1000.0 * libm::sqrt(f64::EPSILON);

        assert_eq!(0.0, calculate_astroid(0.0, 0.0));
        assert_eq!(0.0, calculate_astroid(1.0, 0.0));

        // 0.0, 0.0 to 0.5, 179.5
        assert_eq!(
            0.91583665308532092,
            calculate_astroid(-0.82852367684428574, -0.82576675584253256)
        );
        // 0.0, 0.0 to 1.0, 179.0
        assert_eq!(
            1.9858096632693705,
            calculate_astroid(-1.6572357126833825, -1.6518470456464789)
        );
        // -30.0, 0.0 to 30.0, 179.0
        assert_eq!(
            0.9121190093974804,
            calculate_astroid(-1.9121190093974805, 0.0)
        );
        // -30.0, 0.0 to 30.5, 179.5
        assert_eq!(
            1.2324261949931818,
            calculate_astroid(-0.96091919533424308, -1.1124132048023443)
        );

        assert_eq!(
            4.463096559488633e-14,
            calculate_astroid(x_threshold, -Y_TOLERANCE - f64::EPSILON)
        );
    }

    #[test]
    fn test_calculate_end_azimuth() {
        let angle_50 = Angle::from(Degrees(50.0));
        let angle_20 = Angle::from(Degrees(20.0));

        let result: Angle = calculate_end_azimuth(angle_20, angle_50, angle_20);
        assert!(is_within_tolerance(
            30.0,
            Degrees::from(result).0,
            32.0 * f64::EPSILON
        ));

        let result: Angle = calculate_end_azimuth(-angle_50, angle_50, angle_20);
        assert_eq!(20.0, Degrees::from(result).0);
    }

    #[test]
    fn test_delta_omega12() {
        // 0.0, 0.0 to 30.0, 90.0
        let clairaut_30_90 = Angle::from(Degrees(60.0)).sin();
        let eps_30_90 = WGS84_ELLIPSOID.calculate_epsilon(clairaut_30_90);
        let lam12_30_90 = delta_omega12(
            clairaut_30_90,
            eps_30_90,
            Radians(std::f64::consts::FRAC_PI_2),
            Angle::from_y_x(0.0, 1.0),
            Angle::from_y_x(1.0, 0.0),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(0.0045600360192803542, lam12_30_90.0);

        // 0.0, 0.0 to 45.0, 90.0
        let clairaut_45_90 = Angle::from(Degrees(45.0)).sin();
        let eps_45_90 = WGS84_ELLIPSOID.calculate_epsilon(clairaut_45_90);
        let lam12_45_90 = delta_omega12(
            clairaut_45_90,
            eps_45_90,
            Radians(std::f64::consts::FRAC_PI_2),
            Angle::from_y_x(0.0, 1.0),
            Angle::from_y_x(1.0, 0.0),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(0.0037224722989948442, lam12_45_90.0);

        // 0.0, 0.0 to 60.0, 90.0
        let clairaut_60_90 = Angle::from(Degrees(30.0)).sin();
        let eps_60_90 = WGS84_ELLIPSOID.calculate_epsilon(clairaut_60_90);
        let lam12_60_90 = delta_omega12(
            clairaut_60_90,
            eps_60_90,
            Radians(std::f64::consts::FRAC_PI_2),
            Angle::from_y_x(0.0, 1.0),
            Angle::from_y_x(1.0, 0.0),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(0.0026316334829412581, lam12_60_90.0);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_meridian() {
        let latlon1 = LatLong::new(Degrees(-70.0), Degrees(40.0));
        let latlon2 = LatLong::new(Degrees(80.0), Degrees(40.0));

        let tolerance = Radians(great_circle::MIN_VALUE);

        // Northbound geodesic segment along a meridian
        let result = calculate_azimuths_aux_length(&latlon1, &latlon2, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(0.0, Degrees::from(result.0).0);
        assert_eq!(2.6163378712682306, (result.1).0);
        assert_eq!(0.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);

        // Southbound geodesic segment along a meridian
        let result = calculate_azimuths_aux_length(&latlon2, &latlon1, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(180.0, Degrees::from(result.0).0);
        assert_eq!(2.6163378712682306, (result.1).0);
        assert_eq!(180.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);

        // Northbound geodesic segment past the North pole
        let latlon3: LatLong = LatLong::new(Degrees(80.0), Degrees(-140.0));
        let result = calculate_azimuths_aux_length(&latlon2, &latlon3, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(0.0, Degrees::from(result.0).0);
        assert_eq!(0.3502163200513691, (result.1).0);
        assert_eq!(180.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_equator() {
        let latlon1 = LatLong::new(Degrees(0.0), Degrees(-40.0));
        let latlon2 = LatLong::new(Degrees(0.0), Degrees(50.0));

        let tolerance = Radians(great_circle::MIN_VALUE);

        // Eastbound geodesic segment along the equator
        let result = calculate_azimuths_aux_length(&latlon1, &latlon2, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(90.0, Degrees::from(result.0).0);
        assert_eq!(1.5760806267286946, (result.1).0);
        assert_eq!(90.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);

        // Westbound geodesic segment along the equator
        let result = calculate_azimuths_aux_length(&latlon2, &latlon1, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(-90.0, Degrees::from(result.0).0);
        assert_eq!(1.5760806267286946, (result.1).0);
        assert_eq!(-90.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);

        // Long Eastbound geodesic segment along the equator
        let latlon3 = LatLong::new(Degrees(0.0), Degrees(135.0));
        let result = calculate_azimuths_aux_length(&latlon1, &latlon3, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(90.0, Degrees::from(result.0).0);
        assert_eq!(3.0646012186391296, (result.1).0);
        assert_eq!(90.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_equator_antipodal() {
        let latlon1 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(0.0), Degrees(180.0));

        let tolerance = Radians(great_circle::MIN_VALUE);

        // Northbound geodesic segment along the equator
        let result = calculate_azimuths_aux_length(&latlon1, &latlon2, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(0.0, Degrees::from(result.0).0);
        assert_eq!(core::f64::consts::PI, (result.1).0);
        assert_eq!(180.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_equator_nearly_antipodal() {
        let latlon1 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(0.0), Degrees(179.5));

        let tolerance = Radians(great_circle::MIN_VALUE);

        // Northbound geodesic segment along the equator
        let result = calculate_azimuths_aux_length(&latlon1, &latlon2, tolerance, &WGS84_ELLIPSOID);
        assert_eq!(55.96649514015865, Degrees::from(result.0).0);
        assert_eq!(core::f64::consts::PI, (result.1).0);
        assert_eq!(124.03350485984134, Degrees::from(result.2).0);
        assert_eq!(3, result.3);
    }
    #[test]
    fn test_calculate_azimuths_aux_length_normal_01() {
        // North West bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(-40.0), Degrees(70.0));
        let latlon2 = LatLong::new(Degrees(30.0), Degrees(0.0));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(-55.00473169905792, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428877, (result.1).0);
        assert_eq!(-46.47061016713593, Degrees::from(result.2).0);
        assert_eq!(5, result.3);

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(-40.0)));
        let g = GeodesicSegment::new(
            beta1,
            Angle::from(Degrees(70.0)),
            result.0,
            result.1,
            &WGS84_ELLIPSOID,
        );
        let latlon3 = g.aux_lat_long(result.1);
        assert_eq!(
            Degrees::from(latlon3.lat()).0,
            Degrees::from(latlon2.lat()).0
        );
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lon()).0,
            Degrees::from(latlon3.lon()).0,
            64.0 * f64::EPSILON
        ));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_02() {
        // South West bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(30.0), Degrees(70.0));
        let latlon2 = LatLong::new(Degrees(-40.0), Degrees(0.0));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(-133.52938983286407, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428877, (result.1).0);
        assert_eq!(-124.99526830094209, Degrees::from(result.2).0);
        assert_eq!(5, result.3);

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(30.0)));
        let g = GeodesicSegment::new(
            beta1,
            Angle::from(Degrees(70.0)),
            result.0,
            result.1,
            &WGS84_ELLIPSOID,
        );
        let latlon3 = g.aux_lat_long(result.1);
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lat()).0,
            Degrees::from(latlon3.lat()).0,
            64.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lon()).0,
            Degrees::from(latlon3.lon()).0,
            64.0 * f64::EPSILON
        ));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_03() {
        // South East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(30.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(-40.0), Degrees(70.0));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(133.52938983286407, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428877, (result.1).0);
        assert_eq!(124.99526830094209, Degrees::from(result.2).0);
        assert_eq!(5, result.3);

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(30.0)));
        let g = GeodesicSegment::new(
            beta1,
            Angle::default(),
            result.0,
            result.1,
            &WGS84_ELLIPSOID,
        );
        let latlon3 = g.aux_lat_long(result.1);
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lat()).0,
            Degrees::from(latlon3.lat()).0,
            64.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lon()).0,
            Degrees::from(latlon3.lon()).0,
            64.0 * f64::EPSILON
        ));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_04() {
        // North East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(-40.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(30.0), Degrees(70.0));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(55.00473169905792, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428877, (result.1).0);
        assert_eq!(46.47061016713593, Degrees::from(result.2).0);
        assert_eq!(5, result.3);

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(-40.0)));
        let g = GeodesicSegment::new(
            beta1,
            Angle::default(),
            result.0,
            result.1,
            &WGS84_ELLIPSOID,
        );
        let latlon3 = g.aux_lat_long(result.1);
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lat()).0,
            Degrees::from(latlon3.lat()).0,
            64.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            Degrees::from(latlon2.lon()).0,
            Degrees::from(latlon3.lon()).0,
            64.0 * f64::EPSILON
        ));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_05() {
        // North East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(0.5), Degrees(179.98));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(1.0420381519981552, Degrees::from(result.0).0);
        assert_eq!(3.132893826005981, (result.1).0);
        assert_eq!(178.9579224301469, Degrees::from(result.2).0);
        assert_eq!(5, result.3);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_06() {
        // GeodTest.dat line 460107
        // 89.985810803742 0 90.033923043742 -89.985810803761488692 179.999716989078075251 89.966210133068275597 20003931.4528694 179.999999966908046132 .0036837809003 -47969483155.576793
        let latlon1 = LatLong::new(Degrees(89.985810803742), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-89.985810803761488692),
            Degrees(179.999716989078075251),
        );

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(90.03393799266541, Degrees::from(result.0).0); // 90.033923043742
        assert!(is_within_tolerance(
            179.999999966908046132_f64.to_radians(),
            (result.1).0,
            2.0 * f64::EPSILON
        ));
        assert_eq!(89.96619518414488, Degrees::from(result.2).0); // 89.966210133068275597
        assert_eq!(2, result.3);

        let beta_1 =
            WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(89.985810803742)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(20003931.452869404, distance.0); // 20003931.4528694
        assert!(is_within_tolerance(20003931.4528694, distance.0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_07() {
        // GeodTest.dat line 4701132
        // 85.89252711453 0 90.028477847874 -85.892527114530046541 179.956663887832388079 89.971522153429881464 20003758.1089151 179.999999999906050673 .0000961077487 -40347877014.28062
        let latlon1 = LatLong::new(Degrees(85.89252711453), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-85.892527114530046541),
            Degrees(179.956663887832388079),
        );

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(90.02817870320018, Degrees::from(result.0).0); // 90.028477847874
        assert_eq!(89.97182129799214, Degrees::from(result.2).0); // 89.971522153429881464
        assert_eq!(2, result.3);

        let beta_1 =
            WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(85.89252711453)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(20003758.108915098, distance.0); // 20003758.1089151
        assert!(is_within_tolerance(20003758.1089151, distance.0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_08() {
        // GeodTest.dat line 470905
        // 56.706063494255 0 89.959815697468 -56.706063494254993895 179.668131859151492609 90.0401843025452288 19993766.626695 179.999999999991283078 .0100047302652 56858075837.296755
        let latlon1 = LatLong::new(Degrees(56.706063494255), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-56.706063494254993895),
            Degrees(179.668131859151492609),
        );

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(89.9598195917545, Degrees::from(result.0).0); // 89.959815697468
        assert_eq!(90.0401804082455, Degrees::from(result.2).0); // 90.0401843025452288
        assert_eq!(2, result.3);

        let beta_1 =
            WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(56.706063494255)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(19993766.626695, distance.0); // 19993766.626695
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_09() {
        // GeodTest.dat line 497725
        // 85.224117973184 0 89.957303913327 -85.2241179731839905 179.949627233487121769 90.042696086825586442 20003697.2437342 179.999999999987208748 .000261547236 60492019939.566295
        let latlon1 = LatLong::new(Degrees(85.224117973184), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-85.2241179731839905),
            Degrees(179.949627233487121769),
        );

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        // GeodTest.dat azimuths are swapped around
        assert_eq!(89.9575382708637, Degrees::from(result.0).0); // 89.957303913327
        assert!(is_within_tolerance(
            179.999999999987208748_f64.to_radians(),
            (result.1).0,
            1024.0 * f64::EPSILON
        ));
        assert_eq!(90.0424617291363, Degrees::from(result.2).0); // 90.042696086825586442
        assert_eq!(2, result.3);

        let beta_1 =
            WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(85.224117973184)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(20003697.2437342, distance.0); // 20003697.2437342
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_10() {
        // GeodTest.dat line 451464. but antipodal
        let lat1d = 30.815985336295;
        let lat2d = -30.8159853362949972;
        let lon2d = 180.0;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(0.0, Degrees::from(result.0).0);
        assert_eq!(3.1415926237874707, (result.1).0);
        assert_eq!(180.0, Degrees::from(result.2).0);
        assert_eq!(0, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(20003931.26901283, distance.0);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_11() {
        // GeodTest.dat line 451464. but nearly antipodal
        let lat1d = 30.815985336295;
        let lat2d = -30.8159853362949972;
        let lon2d = 179.99;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(1.1054776533957336, Degrees::from(result.0).0);
        assert_eq!(3.141592653589793, (result.1).0);
        assert_eq!(178.89452234660428, Degrees::from(result.2).0);
        assert_eq!(3, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(20003922.22814904, distance.0);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_normal_12() {
        // GeodTest.dat line 451464. but slightly closer to antiopodal
        let lat1d = 30.815985336295;
        let lat2d = -30.8159853362949972;
        let lon2d = 179.5;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(74.60015893697746, Degrees::from(result.0).0);
        assert_eq!(3.1415926535897927, (result.1).0);
        assert_eq!(105.39984106302255, Degrees::from(result.2).0);
        assert_eq!(3, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(19980861.90889096, distance.0);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_451464() {
        // GeodTest.dat line 451464
        // 30.815985336295 0 89.999989151475 -30.8159853362949972 179.481356807121660669 90.000010857299389702 19979110.018652 179.999999985240671654 .0016389343224 15326161.345917
        let lat1d = 30.815985336295;
        let lat2d = -30.8159853362949972;
        let lon2d = 179.481356807121660669;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(89.9999712763982, Degrees::from(result.0).0); // 89.999989151475
        assert_eq!(90.00002873237621, Degrees::from(result.2).0); // 90.00000000877435
        assert_eq!(11, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(19979110.018652, distance.0);
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_457539() {
        // GeodTest.dat line 457539
        // 50.853729012223 0 89.99999358451 -50.853729012222999629 179.61842826953446462 90.00000641955036513 19990494.4151505 179.999999996683663418 .0003686784385 9076399.666959
        let lat1d = 50.853729012223;
        let lat2d = -50.853729012222999629;
        let lon2d = 179.61842826953446462;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(90.00000000203018, Degrees::from(result.0).0); // 89.99999358451
        assert!(is_within_tolerance(
            179.999999996683663418_f64.to_radians(),
            (result.1).0,
            f64::EPSILON
        ));
        assert_eq!(90.00000000203018, Degrees::from(result.2).0); // 90.00000641955036513
        assert_eq!(2, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        // assert_eq!(19990494.415150497, distance.0); // 19990494.4151505
        assert!(is_within_tolerance(19990494.4151505, distance.0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_459042() {
        // GeodTest.dat line 459042
        // 80.614649787777 0 89.999981005382 -80.614649787776995822 179.901342365578970087 90.000019070703465263 20003033.0359723 179.999999987381812901 .0014045245045 26963180.980693
        let lat1d = 80.614649787777;
        let lat2d = -80.614649787776995822;
        let lon2d = 179.901342365578970087;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(90.00000003804273, Degrees::from(result.0).0); // 89.999981005382
        assert!(is_within_tolerance(
            179.999999987381812901_f64.to_radians(),
            (result.1).0,
            f64::EPSILON
        ));
        assert_eq!(90.00000003804273, Degrees::from(result.2).0); // 90.000019070703465263
        assert_eq!(2, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        // assert_eq!(20003033.035972305, distance.0); // 20003033.0359723
        assert!(is_within_tolerance(20003033.0359723, distance.0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_461124() {
        // GeodTest.dat line 461124
        // 83.829512252973 0 90.006690097427 -83.829512252973003785 179.934969072220767866 89.993309902872831362 20003541.1016439 179.999999999967475093 .0000142648344 -9478364844.203696
        let lat1d = 83.829512252973;
        let lat2d = -83.829512252973003785;
        let lon2d = 179.934969072220767866;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(89.99621690421323, Degrees::from(result.0).0); // 90.006690097427
        assert_eq!(90.00378309578677, Degrees::from(result.2).0); // 89.993309902872831362
        assert_eq!(1, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        // assert_eq!(20003541.101643898, distance.0); // 20003541.1016439
        assert!(is_within_tolerance(20003541.1016439, distance.0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_472656() {
        // GeodTest.dat line 472656
        // 89.725731151764 0 89.99999412823 -89.725731151763992167 179.997094777788067591 90.000014899078757451 20003930.6795515 179.999999956641523093 .0048266431982 14714576.102069
        let lat1d = 89.725731151764;
        let lat2d = -89.725731151763992167;
        let lon2d = 179.997094777788067591;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(90.00000451365437, Degrees::from(result.0).0); // 89.99999412823
        assert!(is_within_tolerance(
            179.999999956641523093_f64.to_radians(),
            (result.1).0,
            f64::EPSILON
        ));
        assert_eq!(90.00000451365437, Degrees::from(result.2).0); // 90.000014899078757451
        assert_eq!(2, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        assert_eq!(20003930.6795515, distance.0); // 20003930.6795515
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_475646() {
        // GeodTest.dat line 475646
        // 38.464554001342 0 90.002808565642 -38.464554001342002684 179.526999021189467903 89.997191434401322223 19983285.4392258 179.999999999945283834 .0001053731912 -3968382822.969946
        let lat1d = 38.464554001342;
        let lat2d = -38.464554001342002684;
        let lon2d = 179.526999021189467903;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(89.99736901205219, Degrees::from(result.0).0); // 90.002808565642
        assert_eq!(90.00263098794781, Degrees::from(result.2).0); // 89.997191434401322223
        assert_eq!(2, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        // assert_eq!(19983285.439225797, distance.0); // 19983285.4392258
        assert!(is_within_tolerance(19983285.4392258, distance.0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuths_aux_length_geodtest_470131() {
        // GeodTest.dat line 470131
        // 85.89252711453 0 90.028477847874 -85.892527114530046541 179.956663887832388079 89.971522153429881464 20003758.1089151 179.999999999906050673 .0000961077487 -40347877014.28062
        let lat1d = 85.89252711453;
        let lat2d = -85.892527114530046541;
        let lon2d = 179.956663887832388079;

        let latlon1 = LatLong::new(Degrees(lat1d), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(lat2d), Degrees(lon2d));

        let result = calculate_azimuths_aux_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        assert_eq!(90.02817870320018, Degrees::from(result.0).0); // 90.028477847874
        assert_eq!(89.97182129799214, Degrees::from(result.2).0); // 89.97152215342988146
        assert_eq!(2, result.3);

        let beta_1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(Degrees(lat1d)));
        let distance = convert_radians_to_metres(beta_1, result.0, result.1, &WGS84_ELLIPSOID);
        // assert_eq!(20003758.108915098, distance.0); // 20003758.1089151
        assert!(is_within_tolerance(20003758.1089151, distance.0, 1e-8));
    }
}
