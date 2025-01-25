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

//! The ellipsoid module contains types and functions for defining an ellipsoid
//! given its Semimajor axis (the equivalent of its radius) and flattening ratio.

#![allow(clippy::suboptimal_flops)]

pub mod coefficients;
pub mod wgs84;

use crate::Metres;
use angle_sc::{trig, Angle};

/// Calculate the Semiminor axis of an ellipsoid.
/// * `a` - the Semimajor axis of an ellipsoid.
/// * `f` - the flattening ratio.
/// # Examples
/// ```
/// use icao_wgs84::Metres;
/// use icao_wgs84::ellipsoid::{calculate_minor_axis, wgs84};
///
/// // The WGS 84 Semiminor axis measured in metres.
/// let b : Metres = Metres(6_356_752.314_245_179);
/// assert_eq!(b, calculate_minor_axis(wgs84::A, wgs84::F));
/// ```
#[must_use]
pub fn calculate_minor_axis(a: Metres, f: f64) -> Metres {
    Metres(a.0 * (1.0 - f))
}

/// Calculate the square of the Eccentricity of an ellipsoid.
/// * `f` - the flattening ratio.
/// # Examples
/// ```
/// use icao_wgs84::ellipsoid::{calculate_sq_eccentricity, wgs84};
///
/// // The WGS 84 sq_eccentricity.
/// assert_eq!(0.0066943799901413165, calculate_sq_eccentricity(wgs84::F));
/// ```
#[must_use]
pub fn calculate_sq_eccentricity(f: f64) -> f64 {
    f * (2.0 - f)
}

/// Calculate the square of the second Eccentricity of an ellipsoid.
/// * `f` - the flattening ratio.
/// # Examples
/// ```
/// use icao_wgs84::ellipsoid::{calculate_sq_2nd_eccentricity, wgs84};
///
/// // The WGS 84 sq 2nd eccentricity.
/// assert_eq!(0.006739496742276434, calculate_sq_2nd_eccentricity(wgs84::F));
/// ```
#[must_use]
pub fn calculate_sq_2nd_eccentricity(f: f64) -> f64 {
    let one_minus_f = 1.0 - f;
    calculate_sq_eccentricity(f) / (one_minus_f * one_minus_f)
}

/// Calculate the third flattening of an ellipsoid.
/// * `f` - the flattening ratio.
/// # Examples
/// ```
/// use icao_wgs84::ellipsoid::{calculate_3rd_flattening, wgs84};
///
/// // The WGS 84 3rd flattening.
/// assert_eq!(0.0016792203863837047, calculate_3rd_flattening(wgs84::F));
/// ```
#[must_use]
pub fn calculate_3rd_flattening(f: f64) -> f64 {
    f / (2.0 - f)
}

/// Function to calculate `epsilon`, the variable used in series expansions,
/// derived from Clairaut's constant.
///
/// Note: `epsilon` is positive and small.
/// CFF Karney, [Algorithms for geodesics](https://arxiv.org/pdf/1109.4448.pdf)
/// Eqs 9 & 16.
/// * `clairaut` - Clairaut's constant.
/// * `ep_2` - the square of the second Eccentricity of the ellipsoid.
#[must_use]
pub fn calculate_epsilon(clairaut: trig::UnitNegRange, ep_2: f64) -> f64 {
    // Clairaut's constant is sin alpha0; sq_cos_alpha0 is 1 - clairaut^2
    let sq_cos_alpha0 = (1.0 - clairaut.0) * (1.0 + clairaut.0);
    let k2 = ep_2 * sq_cos_alpha0; // square of Karney equation 9
    let sqrt_k2_1 = libm::sqrt(1.0 + k2) + 1.0;
    k2 / (sqrt_k2_1 * sqrt_k2_1) // Karney equation 16
}

/// Function to convert a `geodetic` Latitude to a `parametric` Latitude on the
/// auxiliary sphere.
/// * `lat` - the `geodetic` Latitude
/// * `one_minus_f` - one minus the flattening ratio.
#[must_use]
pub fn calculate_parametric_latitude(lat: Angle, one_minus_f: f64) -> Angle {
    Angle::from_y_x(one_minus_f * lat.sin().0, lat.cos().0)
}

/// Function to convert a `parametric` Latitude on the auxiliary sphere to a
/// `geodetic` Latitude.
/// * `lat` - the `parametric` Latitude
/// * `one_minus_f` - one minus the flattening ratio.
#[must_use]
pub fn calculate_geodetic_latitude(lat: Angle, one_minus_f: f64) -> Angle {
    Angle::from_y_x(lat.sin().0 / one_minus_f, lat.cos().0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use angle_sc::{is_within_tolerance, Degrees, Radians};

    #[test]
    fn test_calculate_epsilon() {
        let wgs84_ep2 = calculate_sq_2nd_eccentricity(wgs84::F);
        assert_eq!(
            0.0016792203863837047,
            calculate_epsilon(trig::UnitNegRange(0.0), wgs84_ep2)
        );
        assert_eq!(
            0.0015745990877544997,
            calculate_epsilon(trig::UnitNegRange(0.25), wgs84_ep2)
        );
        assert_eq!(
            0.0012604720416530619,
            calculate_epsilon(trig::UnitNegRange(0.5), wgs84_ep2)
        );
        assert_eq!(
            0.0007360477262034019,
            calculate_epsilon(trig::UnitNegRange(0.75), wgs84_ep2)
        );
        assert_eq!(0.0, calculate_epsilon(trig::UnitNegRange(1.0), wgs84_ep2));
    }

    #[test]
    fn test_calculate_parametric_and_geodetic_latitude() {
        let one_minus_f = 1.0 - wgs84::F;

        for i in -90..91 {
            let latitude = i as f64;
            let lat = Angle::from(Degrees(latitude));
            let parametric_lat = calculate_parametric_latitude(lat, one_minus_f);
            let result = calculate_geodetic_latitude(parametric_lat, one_minus_f);

            assert!(is_within_tolerance(
                Radians::from(lat).0,
                Radians::from(result).0,
                f64::EPSILON
            ));
        }
    }
}
