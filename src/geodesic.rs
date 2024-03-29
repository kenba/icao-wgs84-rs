// Copyright (c) 2024 Ken Barker

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

//! The geodesic module contains functions for calculating the geodesic path
//! between two points on the surface of an ellipsoid.

#![allow(clippy::float_cmp)]
#![allow(clippy::many_single_char_names)]

use crate::ellipsoid::coefficients::{
    evaluate_a1, evaluate_a2, evaluate_coeffs_c1, evaluate_coeffs_c2, evaluate_coeffs_c3y,
    evaluate_polynomial, sin_cos_series,
};
use crate::ellipsoid::{calculate_epsilon, calculate_parametric_latitude, Metres};
use crate::Ellipsoid;
use angle_sc::trig::{cosine_from_sine, UnitNegRange};
use angle_sc::{is_small, Angle, Radians};
use unit_sphere::great_circle;
use unit_sphere::LatLong;

/// The maximum precision, in Radians.
pub const MAX_PRECISION: Radians = Radians(2.0 * core::f64::EPSILON);

/// Estimate omega12 by solving the astroid problem.
/// Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
/// * `x`, `y` - astroid parameters, see Karney section 7.
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
    let a1 = evaluate_a1(eps);
    let a2 = evaluate_a2(eps);
    let m0x = a1 - a2;

    let a1p1 = 1.0 + a1;
    let a2p1 = 1.0 + a2;

    let ca = evaluate_coeffs_c1(eps);
    let mut cb = evaluate_coeffs_c2(eps);

    // Assume here that ca.len() >= cb.len()
    for i in 1..cb.len() {
        cb[i] = a1p1 * ca[i] - a2p1 * cb[i];
    }

    let j12 = m0x * (sigma12 + (sin_cos_series(&cb, sigma2) - sin_cos_series(&cb, sigma1))).0;
    dn2 * (sigma1.cos().0 * sigma2.sin().0)
        - dn1 * (sigma1.sin().0 * sigma2.cos().0)
        - sigma1.cos().0 * sigma2.cos().0 * j12
}

/// Estimate the initial azimuth on the auxiliary sphere for a nearly antipodal arc.  
/// It calculates and solves the astroid problem.
/// * `beta1`, `beta2` - the parametric latitudes of the start and finish points
/// on the auxiliary sphere.
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
    const Y_TOLERANCE: f64 = 200.0 * core::f64::EPSILON;
    const X_TOLERANCE: f64 = 2000.0 / core::f64::consts::FRAC_2_SQRT_PI;

    // Calculate the integration parameter for geodesic
    let clairaut = beta1.cos(); // Note: assumes sin_alpha_1 = 1
    let eps = calculate_epsilon(clairaut, ellipsoid.ep_2());
    let a3f = evaluate_polynomial(&ellipsoid.a3(), eps);

    let lamscale = ellipsoid.f() * beta1.cos().0 * a3f * core::f64::consts::PI;
    let betscale = lamscale * beta1.cos().0;

    // Solve astroid problem
    let x = Radians::from(lambda12.opposite()).0 / lamscale;
    let y = (beta1 + beta2).sin().0 / betscale;

    // Test x and y params
    if (x <= -(1.0 + X_TOLERANCE)) || (y < -Y_TOLERANCE) {
        let k = calculate_astroid(x, y);
        let omg12a = lamscale * (-x * k / (1.0 + k));

        let omega12 = Radians(core::f64::consts::PI - omg12a);
        great_circle::calculate_gc_azimuth(beta1, beta2, Angle::from(omega12))
    } else {
        let sin_alpha = UnitNegRange(if -x < 1.0 { -x } else { 1.0 });
        Angle::new(sin_alpha, cosine_from_sine(sin_alpha, -1.0))
    }
}

/// Calculate the cosine of the longitude difference from the equator crossing.
/// * `beta` the parametric latitude
/// * `cos_azimuth` the cosine of the azimuth at the parametric latitude
///
/// returns the cosine of the longitude difference, zero if the parametric
/// latitude is close to the equator.
#[must_use]
pub fn calculate_cos_omega(beta: Angle, cos_azimuth: UnitNegRange) -> UnitNegRange {
    if is_small(libm::fabs(beta.sin().0), core::f64::EPSILON) {
        UnitNegRange(1.0)
    } else {
        UnitNegRange(cos_azimuth.0 * beta.cos().0)
    }
}

/// Calculate the azimuth on the auxiliary sphere at latitude beta2 given the
/// latitude beta1 and the azimuth at that latitude, alpha1.
/// * `beta1`, `beta2` - the parametric latitudes of the start and finish points
/// on the auxiliary sphere.
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
    let cos_alpha2 =
        if (beta2.cos() != beta1.cos()) || (libm::fabs(beta2.sin().0) != -beta1.sin().0) {
            let temp1 = alpha1.cos().0 * beta1.cos().0;
            let temp2 = if beta1.cos().0 < libm::fabs(beta1.sin().0) {
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
            UnitNegRange(libm::fabs(alpha1.cos().0))
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
) -> f64 {
    let a3f = evaluate_polynomial(&ellipsoid.a3(), eps);
    let a3c = ellipsoid.f() * clairaut.0 * a3f;

    let c3 = evaluate_coeffs_c3y(&ellipsoid.c3x(), eps);
    let b31 = sin_cos_series(&c3, sigma1);
    let b32 = sin_cos_series(&c3, sigma2);

    a3c * (sigma12 + (b32 - b31)).0
}

/// Find the azimuth and great circle length on the auxiliary sphere.
/// It uses Newton's method to solve:
///   f(alp1) = lambda12(alp1) - lam12 = 0
/// * `lat_a`, `lat_b` - the geodetic latitudes of the start and finish points.
/// * `lambda12` - Longitude difference between start and finish points.
///
/// returns the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic.
#[allow(clippy::similar_names)]
#[must_use]
fn find_azimuth_and_aux_length(
    lat_a: Angle,
    lat_b: Angle,
    lambda12: Angle,
    gc_length: Radians,
    ellipsoid: &Ellipsoid,
) -> (Angle, Radians) {
    // The maximum number of iterations to attempt.
    const MAX_ITERS: u32 = 20;
    let antipodal_arc_threshold: f64 = core::f64::consts::PI * ellipsoid.one_minus_f();

    // Start at the latitude furthest from the Equator
    let swap_latitudes = libm::fabs(lat_a.sin().0) < libm::fabs(lat_b.sin().0);
    let mut lat1 = if swap_latitudes { lat_b } else { lat_a };
    let mut lat2 = if swap_latitudes { lat_a } else { lat_b };

    // Start South of the Equator
    let negate_latitude = 0.0 < lat1.sin().0;
    if negate_latitude {
        lat1 = -lat1;
        lat2 = -lat2;
    }

    // project latitudes onto the auxiliary sphere
    let beta1 = calculate_parametric_latitude(lat1, ellipsoid.one_minus_f());
    let beta2 = calculate_parametric_latitude(lat2, ellipsoid.one_minus_f());

    let dn1 = libm::sqrt(1.0 + ellipsoid.ep_2() * beta1.sin().0 * beta1.sin().0);
    let dn2 = libm::sqrt(1.0 + ellipsoid.ep_2() * beta2.sin().0 * beta2.sin().0);

    // Use positive lambda12, so all azimuths are positive
    let abs_lambda12 = lambda12.abs();

    // Estimate the azimuth at the start of the geodesic
    let mut alpha1 = if antipodal_arc_threshold < gc_length.0 {
        estimate_antipodal_initial_azimuth(beta1, beta2, abs_lambda12, ellipsoid)
    } else {
        // Use great circle azimuth at the start
        great_circle::calculate_gc_azimuth(lat1, lat2, abs_lambda12)
    };
    let mut alpha2 = alpha1;

    let mut sigma12_rad = gc_length;

    for _i in 0..MAX_ITERS {
        // Calculate Clairaut's constant
        let clairaut = UnitNegRange(alpha1.sin().0 * beta1.cos().0);
        let eps = calculate_epsilon(clairaut, ellipsoid.ep_2());

        // Calculate first longitude (omega1) and distance (sigma1) from the
        // Northbound equator crossing
        let sin_omega1 = UnitNegRange(clairaut.0 * beta1.sin().0);
        let cos_omega1 = calculate_cos_omega(beta1, alpha1.cos());
        let omega1 = Angle::from_y_x(sin_omega1.0, cos_omega1.0);
        let sigma1 = Angle::from_y_x(beta1.sin().0, cos_omega1.0);

        // Calculate azimuth at the end point
        alpha2 = calculate_end_azimuth(beta1, beta2, alpha1);

        // Calculate second longitude (omega2) and distance (sigma2) from the
        // Northbound equator crossing
        let sin_omega2 = UnitNegRange(clairaut.0 * beta2.sin().0);
        let cos_omega2 = calculate_cos_omega(beta2, alpha2.cos());
        let omega2 = Angle::from_y_x(sin_omega2.0, cos_omega2.0);
        let sigma2 = Angle::from_y_x(beta2.sin().0, cos_omega2.0);

        // Calculate Longitude difference on the auxiliary sphere
        let mut omega12 = omega2 - omega1;
        // clamp to range 0 to Pi
        if omega12.sin() < UnitNegRange(0.0) {
            omega12 = Angle::from_y_x(0.0, omega12.cos().0);
        }

        // Calculate great circle length on the auxiliary sphere
        let mut sigma12 = sigma2 - sigma1;
        // clamp to range 0 to Pi
        if sigma12.sin() < UnitNegRange(0.0) {
            sigma12 = Angle::from_y_x(0.0, sigma12.cos().0);
        }

        // Calculate difference between geodesic and great circle longitudes
        let eta = Radians::from(omega12 - abs_lambda12);

        sigma12_rad = Radians::from(sigma12);
        let domg12 = delta_omega12(clairaut, eps, sigma12_rad, sigma1, sigma2, ellipsoid);

        // Difference between differences
        let v = eta.0 - domg12;
        if is_small(libm::fabs(v), MAX_PRECISION.0) {
            break;
        }

        // Calculate the denominator for Newton's method
        let dv = if is_small(libm::fabs(alpha2.cos().0), core::f64::EPSILON) {
            -2.0 * ellipsoid.one_minus_f() * dn1 / beta1.sin().0
        } else {
            let m12 = calculate_reduced_length(eps, sigma12_rad, sigma1, dn1, sigma2, dn2);
            ellipsoid.one_minus_f() * m12 / (alpha2.cos().0 * beta2.cos().0)
        };
        if is_small(libm::fabs(dv), MAX_PRECISION.0) {
            break;
        }

        // Calculate the change in initial azimuth
        let dalpha1 = UnitNegRange::clamp(-v / dv);
        if is_small(libm::fabs(dalpha1.0), MAX_PRECISION.0) {
            break;
        }

        // Adjust the azimuth by dalpha1
        alpha1 = alpha1 + Angle::from(Radians(dalpha1.0));
    }

    if swap_latitudes {
        alpha1 = alpha2;
    }

    if swap_latitudes != negate_latitude {
        alpha1 = alpha1.negate_cos();
    }

    let lambda12_negative = lambda12.sin().0 < 0.0;
    if lambda12_negative {
        alpha1 = -alpha1;
    }

    (alpha1, sigma12_rad)
}

/// Calculate the initial azimuth and great circle length between a pair
/// of points on the auxiliary sphere.
/// * `lat1`, `lat2` - the geodesic latitudes of the start and finish points
/// on the auxiliary sphere.
/// * `delta_long` - the geodesic longitude difference.
///
/// returns the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic.
#[must_use]
pub fn aux_sphere_azimuth_length(
    lat1: Angle,
    lat2: Angle,
    delta_long: Angle,
    ellipsoid: &Ellipsoid,
) -> (Angle, Radians) {
    const MIN_VALUE: UnitNegRange = UnitNegRange(2.0 * core::f64::EPSILON);
    const MAX_LENGTH: Radians = Radians(core::f64::consts::PI - 2.0 * MIN_VALUE.0);

    let gc_azimuth = great_circle::calculate_gc_azimuth(lat1, lat2, delta_long);
    let gc_length = great_circle::calculate_gc_distance(lat1, lat2, delta_long);
    // Ensure that the points are far enough apart
    if gc_length.0 <= MIN_VALUE.0 {
        return (gc_azimuth, Radians(0.0));
    }

    // Determine whether a meridional path
    let abs_delta_long = Radians::from(delta_long.abs());
    if (abs_delta_long.0 <= MIN_VALUE.0)
        || (MAX_LENGTH <= abs_delta_long)
        || (lat1.cos() <= MIN_VALUE)
        || (lat2.cos() <= MIN_VALUE)
    {
        let meridian_length = if MAX_LENGTH <= gc_length {
            Radians(core::f64::consts::PI)
        } else {
            let beta1 = calculate_parametric_latitude(lat1, ellipsoid.one_minus_f());
            let beta2 = calculate_parametric_latitude(lat2, ellipsoid.one_minus_f());
            great_circle::calculate_gc_distance(beta1, beta2, delta_long)
        };
        return (gc_azimuth, meridian_length);
    }

    // Determine whether an equitorial path, i.e. both latitudes on the equator.
    if (lat1.abs().sin() <= MIN_VALUE) && (lat2.abs().sin() <= MIN_VALUE) {
        let equator_length = Radians(gc_length.0 * ellipsoid.recip_one_minus_f());
        return (gc_azimuth, equator_length);
    }

    // Iterate using Newton's method to find the azimuth and length
    find_azimuth_and_aux_length(lat1, lat2, delta_long, gc_length, ellipsoid)
}

/// Calculate the geodesic azimuth and great circle length on the auxiliary sphere
/// between a pair of positions .
/// * `a`, `b` - the start and finish positions in geodetic coordinates.
///
/// returns the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic.
#[must_use]
pub fn calculate_azimuth_aux_length(
    a: &LatLong,
    b: &LatLong,
    ellipsoid: &Ellipsoid,
) -> (Angle, Radians) {
    let a_lat = Angle::from(a.lat());
    let b_lat = Angle::from(b.lat());
    let delta_long = Angle::from(b.lon() - a.lon());
    aux_sphere_azimuth_length(a_lat, b_lat, delta_long, ellipsoid)
}

/// Convert a great circle distance on the auxiliary sphere in radians to
/// metres on the ellipsoid.
/// * `beta1`, the start parametric Latitude on the auxiliary sphere.
/// * `alpha1`, the azimuth at the start point.
/// * `gc_distance`, the great circle distance on the auxiliary sphere in radians.
/// * `ellipsoid`, the Ellipsoid
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
    let eps = calculate_epsilon(clairaut, ellipsoid.ep_2());
    let a1 = evaluate_a1(eps) + 1.0;
    let c1 = evaluate_coeffs_c1(eps);
    let b11 = sin_cos_series(&c1, sigma1);
    let b12 = sin_cos_series(&c1, sigma_sum);

    Metres(ellipsoid.b().0 * a1 * (gc_distance + b12 - b11).0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Ellipsoid;
    use angle_sc::{is_within_tolerance, Degrees};

    #[test]
    fn test_calculate_astroid() {
        const Y_TOLERANCE: f64 = 200.0 * core::f64::EPSILON;
        const X_TOLERANCE: f64 = 2000.0 / core::f64::consts::FRAC_2_SQRT_PI;

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
            1771.453850905516,
            calculate_astroid(X_TOLERANCE, -Y_TOLERANCE - core::f64::EPSILON)
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
            32.0 * core::f64::EPSILON
        ));

        let result: Angle = calculate_end_azimuth(-angle_50, angle_50, angle_20);
        assert_eq!(20.0, Degrees::from(result).0);
    }

    #[test]
    fn test_delta_omega12() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // 0.0, 0.0 to 30.0, 90.0
        let clairaut_30_90 = Angle::from(Degrees(60.0)).sin();
        let eps_30_90 = calculate_epsilon(clairaut_30_90, wgs84_ellipsoid.ep_2());
        let lam12_30_90 = delta_omega12(
            clairaut_30_90,
            eps_30_90,
            Radians(std::f64::consts::FRAC_PI_2),
            Angle::from_y_x(0.0, 1.0),
            Angle::from_y_x(1.0, 0.0),
            &wgs84_ellipsoid,
        );
        assert_eq!(0.0045600360192803542, lam12_30_90);

        // 0.0, 0.0 to 45.0, 90.0
        let clairaut_45_90 = Angle::from(Degrees(45.0)).sin();
        let eps_45_90 = calculate_epsilon(clairaut_45_90, wgs84_ellipsoid.ep_2());
        let lam12_45_90 = delta_omega12(
            clairaut_45_90,
            eps_45_90,
            Radians(std::f64::consts::FRAC_PI_2),
            Angle::from_y_x(0.0, 1.0),
            Angle::from_y_x(1.0, 0.0),
            &wgs84_ellipsoid,
        );
        assert_eq!(0.0037224722989948442, lam12_45_90);

        // 0.0, 0.0 to 60.0, 90.0
        let clairaut_60_90 = Angle::from(Degrees(30.0)).sin();
        let eps_60_90 = calculate_epsilon(clairaut_60_90, wgs84_ellipsoid.ep_2());
        let lam12_60_90 = delta_omega12(
            clairaut_60_90,
            eps_60_90,
            Radians(std::f64::consts::FRAC_PI_2),
            Angle::from_y_x(0.0, 1.0),
            Angle::from_y_x(1.0, 0.0),
            &wgs84_ellipsoid,
        );
        assert_eq!(0.0026316334829412581, lam12_60_90);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_meridian() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        let latlon1 = LatLong::new(Degrees(-70.0), Degrees(40.0));
        let latlon2 = LatLong::new(Degrees(80.0), Degrees(40.0));

        // Northbound geodesic along a meridian
        let result = calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(0.0, Degrees::from(result.0).0);
        assert_eq!(2.6163378712682306, (result.1).0);

        // Southbound geodesic along a meridian
        let result = calculate_azimuth_aux_length(&latlon2, &latlon1, &wgs84_ellipsoid);
        assert_eq!(180.0, Degrees::from(result.0).0);
        assert_eq!(2.6163378712682306, (result.1).0);

        // Northbound geodesic past the North pole
        let latlon3 = LatLong::new(Degrees(80.0), Degrees(-140.0));
        let result = calculate_azimuth_aux_length(&latlon2, &latlon3, &wgs84_ellipsoid);
        assert_eq!(0.0, Degrees::from(result.0).0);
        assert_eq!(0.3502163200513691, (result.1).0);
    }
    #[test]
    fn test_calculate_azimuth_aux_length_equator() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        let latlon1 = LatLong::new(Degrees(0.0), Degrees(-40.0));
        let latlon2 = LatLong::new(Degrees(0.0), Degrees(50.0));

        // Eastbound geodesic along the equator
        let result = calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(90.0, Degrees::from(result.0).0);
        assert_eq!(1.5760806267286946, (result.1).0);

        // Westbound geodesic along the equator
        let result = calculate_azimuth_aux_length(&latlon2, &latlon1, &wgs84_ellipsoid);
        assert_eq!(-90.0, Degrees::from(result.0).0);
        assert_eq!(1.5760806267286946, (result.1).0);

        // Long Eastbound geodesic along the equator
        let latlon3 = LatLong::new(Degrees(0.0), Degrees(135.0));
        let result = calculate_azimuth_aux_length(&latlon1, &latlon3, &wgs84_ellipsoid);
        assert_eq!(90.0, Degrees::from(result.0).0);
        assert_eq!(3.0646012186391296, (result.1).0);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_normal_01() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // North West bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(-40.0), Degrees(70.0));
        let latlon2 = LatLong::new(Degrees(30.0), Degrees(0.0));

        let result = calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(-55.00473169905792, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428875, (result.1).0);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_normal_02() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // South West bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(30.0), Degrees(70.0));
        let latlon2 = LatLong::new(Degrees(-40.0), Degrees(0.0));

        let result = calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(-133.52938983286407, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428875, (result.1).0);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_normal_03() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // South East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(30.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(-40.0), Degrees(70.0));

        let result = calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(133.52938983286407, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428875, (result.1).0);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_normal_04() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // North East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(-40.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(30.0), Degrees(70.0));

        let result = calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(55.00473169905792, Degrees::from(result.0).0);
        assert_eq!(1.6656790467428875, (result.1).0);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_normal_05() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // North East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let latlon2 = LatLong::new(Degrees(0.5), Degrees(179.98));

        let result: (Angle, Radians) =
            calculate_azimuth_aux_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(1.042038151998155, Degrees::from(result.0).0);
        assert_eq!(3.132893826005981, (result.1).0);
    }
}
