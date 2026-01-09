// Copyright (c) 2024-2026 Ken Barker

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

//! This module contains coefficients and functions for performing calculations
//! on the surface of an ellipsoid.
//!
//! It uses the equations given by CFF Karney in
//! [Algorithms for geodesics](https://arxiv.org/pdf/1109.4448.pdf) and
//! the equation for A2 in
//! [Geodesics on an arbitrary ellipsoid of revolution](https://arxiv.org/pdf/2208.00492.pdf).

use angle_sc::{Angle, Radians};

/// The scale factor `A1`.
/// CFF Karney, Eq. 17.
/// * `eps` - epsilon the integration variable derived from Clairaut's constant.
/// # Examples
/// ```
/// use icao_wgs84::ellipsoid::{calculate_sq_2nd_eccentricity, wgs84};
/// use icao_wgs84::ellipsoid::coefficients::evaluate_a1;
///
/// // evaluate_a1 for WGS 84 latitude 45.0
/// let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
/// let a1 = evaluate_a1(eps45);
///
/// assert_eq!(0.0033839903702120875, a1);
/// ```
#[must_use]
pub fn evaluate_a1(eps: f64) -> f64 {
    let eps2 = eps * eps;
    let t = eps2 * (eps2 * (eps2 + 4.0) + 64.0) / 256.0;
    (t + eps) / (1.0 - eps)
}

/// The scale factor `A2`.
///
/// CFF Karney [Geodesics on an arbitrary ellipsoid of revolution](https://arxiv.org/pdf/2208.00492.pdf),
/// Eq. A1.
/// * `eps` - epsilon the integration variable derived from Clairaut's constant.
/// # Examples
/// ```
/// use icao_wgs84::ellipsoid::{calculate_sq_2nd_eccentricity, wgs84};
/// use icao_wgs84::ellipsoid::coefficients::evaluate_a2;
///
/// // evaluate_a2 for WGS 84 latitude 45.0
/// let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
/// let a2 = evaluate_a2(eps45);
///
/// assert_eq!(-0.0033669191180908161, a2);
/// ```
#[must_use]
pub fn evaluate_a2(eps: f64) -> f64 {
    let eps2 = eps * eps;
    let t = eps2 * ((-11. * eps2 - 28.) * eps2 - 192.) / 256.;
    (t - eps) / (1. + eps)
}

/// The coefficients `A3`.
/// CFF Karney, Eq. 23.
/// * `n` - the third flattening of the ellipsoid.
#[must_use]
pub fn evaluate_coeffs_a3(n: f64) -> [f64; 6] {
    [
        1.,
        (n - 1.) / 2.,
        (n * (3. * n - 1.) - 2.) / 8.,
        ((-n - 3.) * n - 1.) / 16.,
        (-2. * n - 3.) / 64.,
        -3. / 128.,
    ]
}

/// The coefficients `C1[l]` in the Fourier expansion of `B1`.
/// CFF Karney, Eq. 18.
/// * `eps` - epsilon the integration variable derived from Clairaut's constant.
#[must_use]
pub fn evaluate_coeffs_c1(eps: f64) -> [f64; 7] {
    let eps2 = eps * eps;
    let eps4 = (eps2 * eps) * eps; // Note: not the same as eps2 * eps2!
    let eps6 = (eps4 * eps) * eps;

    [
        0.,
        eps * ((6. - eps2) * eps2 - 16.) / 32.,
        eps2 * ((64. - 9. * eps2) * eps2 - 128.) / 2048.,
        eps * eps2 * (9. * eps2 - 16.) / 768.,
        eps4 * (3. * eps2 - 5.) / 512.,
        eps * eps4 * (-7. / 1280.),
        eps6 * (-7. / 2048.),
    ]
}

/// The coefficients `C1p[l]` in the Fourier expansion of `B1p`.
/// CFF Karney, Eq. 21.
/// * `eps` - epsilon the integration variable derived from Clairaut's constant.
#[must_use]
pub fn evaluate_coeffs_c1p(eps: f64) -> [f64; 6] {
    let eps2 = eps * eps;
    let eps4 = (eps2 * eps) * eps; // Note: not the same as eps2 * eps2!

    [
        0.,
        eps * (eps2 * (205. * eps2 - 432.) + 768.) / 1536.,
        eps2 * (30. - 37. * eps2) / 96.,
        eps * eps2 * (116. - 225. * eps2) / 384.,
        eps4 * 539. / 1536.,
        (eps * eps4) * 3467. / 7680.,
    ]
}

/// The coefficients `C2[l]` in the Fourier expansion of `B2`.
/// CFF Karney, Eq. 42.
/// * `eps` - epsilon the integration variable derived from Clairaut's constant.
#[must_use]
pub fn evaluate_coeffs_c2(eps: f64) -> [f64; 7] {
    let eps2 = eps * eps;
    let eps4 = (eps2 * eps) * eps; // Note: not the same as eps2 * eps2!
    let eps6 = (eps4 * eps) * eps;

    [
        0.,
        eps * (eps2 * (eps2 + 2.) + 16.) / 32.,
        eps2 * (eps2 * (35. * eps2 + 64.) + 384.) / 2048.,
        eps * eps2 * (15. * eps2 + 80.) / 768.,
        eps4 * (7. * eps2 + 35.) / 512.,
        eps * eps4 * 63. / 1280.,
        eps6 * 77. / 2048.,
    ]
}

/// The coefficients `C3x[l]` in the Fourier expansion of `C3`.
/// CFF Karney, Eq. 25.
/// * `n` - the third flattening of the ellipsoid.
#[must_use]
pub fn evaluate_coeffs_c3x(n: f64) -> [f64; 15] {
    [
        (1. - n) / 4.,
        (1. - n * n) / 8.,
        (n * ((-5. * n - 1.) * n + 3.) + 3.) / 64.,
        (n * ((2. - 2. * n) * n + 2.) + 5.) / 128.,
        (n * (3. * n + 11.) + 12.) / 512.,
        ((n - 3.) * n + 2.) / 32.,
        (n * (n * (2. * n - 3.) - 2.) + 3.) / 64.,
        (n * ((-6. * n - 9.) * n + 2.) + 6.) / 256.,
        ((1. - 2. * n) * n + 5.) / 256.,
        (n * ((5. - n) * n - 9.) + 5.) / 192.,
        (n * (n * (10. * n - 6.) - 10.) + 9.) / 384.,
        ((-77. * n - 8.) * n + 42.) / 3072.,
        (n * ((20. - 7. * n) * n - 28.) + 14.) / 1024.,
        ((-7. * n - 40.) * n + 28.) / 2048.,
        (n * (75. * n - 90.) + 42.) / 5120.,
    ]
}

/// Evaluate a first degree polynomial in x using
/// [Estrin's scheme](https://en.wikipedia.org/wiki/Estrin%27s_scheme).
/// * `coeffs` - the polynomial coefficients.
/// * `x` - the variable.
#[must_use]
fn evaluate_2_coeffs(coeffs: &[f64], x: f64) -> f64 {
    x.mul_add(coeffs[1], coeffs[0])
}

/// Evaluate the polynomial in x using
/// [Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method).
/// * `coeffs` - the polynomial coefficients.
/// * `x` - the variable.
#[must_use]
pub fn evaluate_polynomial(coeffs: &[f64], x: f64) -> f64 {
    let mut result: f64 = 0.;

    match coeffs.len() {
        // Use Estrin's scheme for 2 coefficients, since same result as Horner's method
        2 => result = evaluate_2_coeffs(coeffs, x),
        _ => {
            if let Some((last, elements)) = coeffs.split_last() {
                result = *last;
                for element in elements.iter().rev() {
                    result = result.mul_add(x, *element);
                }
            }
        }
    }

    result
}

/// The coefficients `C3[l]` in the Fourier expansion of `C3`.
/// CFF Karney, Eq. 26.
/// * `coeffs` - the polynomial coefficients from `evaluate_coeffs_C3x`.
/// * `eps` - epsilon the integration variable derived from Clairaut's constant.
#[must_use]
pub fn evaluate_coeffs_c3y(coeffs: &[f64], eps: f64) -> [f64; 6] {
    let c1 = eps * evaluate_polynomial(&coeffs[0..5], eps);
    let eps_2 = eps * eps;
    let c2 = eps_2 * evaluate_polynomial(&coeffs[5..9], eps);
    let eps_3 = eps * eps_2;
    let c3 = eps_3 * evaluate_polynomial(&coeffs[9..12], eps);
    let eps_4 = eps * eps_3;
    let c4 = eps_4 * evaluate_polynomial(&coeffs[12..14], eps);
    let eps_5: f64 = eps * eps_4;
    let c5 = eps_5 * evaluate_polynomial(&coeffs[14..15], eps);
    [0.0, c1, c2, c3, c4, c5]
}

/// Evaluate the following:
///   `y = sum(c[i] * sin(2*i * angle), i, 1, n)`
/// using [Clenshaw summation](https://en.wikipedia.org/wiki/Clenshaw_algorithm).
/// * `coeffs` - the polynomial coefficients.
/// * `angle` - the Angle.
#[must_use]
pub fn sin_cos_series(coeffs: &[f64], angle: Angle) -> Radians {
    let angle2x = angle.double();

    if angle2x.sin().abs().0 < f64::EPSILON {
        Radians(0.0)
    } else {
        // the Clenshaw ak(theta) parameter, beta(k) = -1
        let ar = 2.0 * angle2x.cos().0;

        let mut index = coeffs.len() - 1;
        let coeffs_length_is_odd = 0 != (index & 1);
        let mut k1 = if coeffs_length_is_odd {
            0.0
        } else {
            coeffs[index]
        };
        if !coeffs_length_is_odd {
            index -= 1;
        }
        let mut k0 = ar.mul_add(k1, coeffs[index]);
        index -= 1;

        // Unroll loop x 2, so accumulators return to their original role.
        while 0 < index {
            k1 = coeffs[index] + ar.mul_add(k0, -k1);
            index -= 1;
            k0 = coeffs[index] + ar.mul_add(k1, -k0);
            index -= 1;
        }
        Radians(angle2x.sin().0 * k0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ellipsoid::{calculate_3rd_flattening, calculate_sq_2nd_eccentricity, wgs84};
    use angle_sc::{Angle, Radians};

    #[test]
    fn test_evaluate_coeffs_a3() {
        // evaluate_coeffs_a3 for WGS 84 flattening
        let n = calculate_3rd_flattening(wgs84::F);
        let a3 = evaluate_coeffs_a3(n);

        assert_eq!(1.0, a3[0]);
        assert_eq!(-0.49916038980680816, a3[1]);
        assert_eq!(-0.2502088451303832, a3[2]);
        assert_eq!(-0.06281503005876607, a3[3]);
        assert_eq!(-0.046927475637074494, a3[4]);
        assert_eq!(-0.0234375, a3[5]);
    }

    #[test]
    fn test_evaluate_coeffs_c1() {
        // evaluate_coeffs_c1 for WGS 84 latitude 45.0
        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let c1 = evaluate_coeffs_c1(eps45);

        assert_eq!(0.0, c1[0]);
        assert_eq!(-0.0016848670110488485, c1[1]);
        assert_eq!(-7.09696225910107e-07, c1[2]);
        assert_eq!(-7.971653346618919e-10, c1[3]);
        assert_eq!(-1.259177551940401e-12, c1[4]);
        assert_eq!(-2.3761586316497056e-15, c1[5]);
        assert_eq!(-5.004410424104756e-18, c1[6]);
    }
    #[test]
    fn test_evaluate_coeffs_c1p() {
        // evaluate_coeffs_c1 for WGS 84 latitude 45.0
        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let c1p = evaluate_coeffs_c1p(eps45);

        assert_eq!(0.0, c1p[0]);
        assert_eq!(0.0016848634238263412, c1p[1]);
        assert_eq!(3.548451580617158e-06, c1p[2]);
        assert_eq!(1.1558716594815811e-08, c1p[3]);
        assert_eq!(4.52467549714072e-11, c1p[4]);
        assert_eq!(1.9614623752213165e-13, c1p[5]);
    }
    #[test]
    fn test_evaluate_coeffs_c2() {
        // evaluate_coeffs_c1 for WGS 84 latitude 45.0
        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let c2 = evaluate_coeffs_c2(eps45);

        assert_eq!(0.0, c2[0]);
        assert_eq!(0.0016848765770939658, c2[1]);
        assert_eq!(2.129104795318516e-06, c2[2]);
        assert_eq!(3.985860618432769e-09, c2[3]);
        assert_eq!(8.814322934149593e-12, c2[4]);
        assert_eq!(2.138542768484735e-14, c2[5]);
        assert_eq!(5.5048514665152315e-17, c2[6]);
    }

    #[test]
    fn test_evaluate_coeffs_c3x() {
        // evaluate_coeffs_c3x for WGS 84 flattening
        let n = calculate_3rd_flattening(wgs84::F);
        let c3x = evaluate_coeffs_c3x(n);

        assert_eq!(0.24958019490340408, c3x[0]);
        assert_eq!(0.12499964752736174, c3x[1]);
        assert_eq!(0.04695366902660743, c3x[2]);
        assert_eq!(0.03908878180363212, c3x[3]);
        assert_eq!(0.02347359352264363, c3x[4]);
        assert_eq!(0.062342661206936094, c3x[5]);
        assert_eq!(0.046822392333655975, c3x[6]);
        assert_eq!(0.023450519665361755, c3x[7]);
        assert_eq!(0.01953778742509442, c3x[8]);
        assert_eq!(0.02596302661819293, c3x[9]);
        assert_eq!(0.023393726366666315, c3x[10]);
        assert_eq!(0.013667431352136642, c3x[11]);
        assert_eq!(0.013626013859041152, c3x[12]);
        assert_eq!(0.01363906808890474, c3x[13]);
        assert_eq!(0.008173648759532895, c3x[14]);
    }

    #[test]
    fn test_evaluate_coeffs_c3y() {
        // evaluate_coeffs_c3x for WGS 84 flattening
        let n = calculate_3rd_flattening(wgs84::F);
        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let c3x = evaluate_coeffs_c3x(n);
        let c3y = evaluate_coeffs_c3y(&c3x, eps45);

        assert_eq!(0.0, c3y[0]);
        assert_eq!(0.0008424436534463023, c3y[1]);
        assert_eq!(7.09708293878426e-07, c3y[2]);
        assert_eq!(9.964762770100143e-10, c3y[3]);
        assert_eq!(1.7628733595825342e-12, c3y[4]);
        assert_eq!(3.5514305923724795e-15, c3y[5]);
    }

    #[test]
    fn test_sin_cos_series_c1() {
        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let c1 = evaluate_coeffs_c1(eps45);

        let angle = Angle::from(Radians(0.1 * std::f64::consts::PI));
        let sin_cos_c1 = sin_cos_series(&c1, angle);

        assert_eq!(-0.0009910157012782634, sin_cos_c1.0);
    }

    #[test]
    fn test_sin_cos_series_c3() {
        // evaluate_coeffs_c3x for WGS 84 flattening
        let n = calculate_3rd_flattening(wgs84::F);
        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let c3x = evaluate_coeffs_c3x(n);
        let c3y = evaluate_coeffs_c3y(&c3x, eps45);

        let angle = Angle::from(Radians(0.1 * std::f64::consts::PI));
        let sin_cos_c3 = sin_cos_series(&c3y, angle);

        assert_eq!(0.00049585187682213688, sin_cos_c3.0);
    }

    #[test]
    fn test_evaluate_poynomial_a3() {
        // evaluate_coeffs_a3 for WGS 84 flattening
        let n = calculate_3rd_flattening(wgs84::F);
        let a3 = evaluate_coeffs_a3(n);

        let eps45 = calculate_sq_2nd_eccentricity(wgs84::F) / 2.0;
        let a3_eps = evaluate_polynomial(&a3, eps45);
        assert_eq!(0.9983151115073848, a3_eps);

        let small = evaluate_polynomial(&a3, 0.0);
        assert_eq!(1.0, small);

        let result = evaluate_polynomial(&a3, 0.25);
        assert_eq!(0.85838416624767966, result);

        let result = evaluate_polynomial(&a3, 0.5);
        assert_eq!(0.67635032595433719, result);

        let result = evaluate_polynomial(&a3, 1.0);
        assert_eq!(0.11745075936696803, result);

        let empty: &[f64] = &[];
        let zero = evaluate_polynomial(&empty, eps45);
        assert_eq!(0.0, zero);
    }
}
