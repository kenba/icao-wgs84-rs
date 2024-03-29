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

//! icao-wgs84
//!
//! [![crates.io](https://img.shields.io/crates/v/icao-wgs84.svg)](https://crates.io/crates/icao-wgs84)
//! [![docs.io](https://docs.rs/icao-wgs84/badge.svg)](https://docs.rs/icao-wgs84/)
//! [![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
//! [![Rust](https://github.com/kenba/icao-wgs84-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/kenba/icao-wgs84-rs/actions)
//! [![codecov](https://codecov.io/gh/kenba/icao-wgs84-rs/graph/badge.svg?token=85TJX5VAHF)](https://codecov.io/gh/kenba/icao-wgs84-rs)
//!
//! A library for performing geometric calculations on the
//! [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) ellipsoid,
//! see *Figure 1*.
//!
//! <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/3e/WGS84_mean_Earth_radius.svg/800px-WGS84_mean_Earth_radius.svg.png" width="400">
//!
//! *Figure 1 The WGS84 Ellipsoid (not to scale)*
//!
//! WGS84 has become the de facto standard for satellite navigation since its adoption
//! by the Navstar [Global Positioning System](https://en.wikipedia.org/wiki/Global_Positioning_System)
//! (GPS) and US president Ronald Reagan's 1983 decision to make GPS available for civilian use
//! after airliner [KAL 007](https://en.wikipedia.org/wiki/Korean_Air_Lines_Flight_007)
//! was shot down by Soviet interceptor aircraft when it strayed into
//! prohibited airspace due to navigational errors.
//!
//! This library uses the WGS84 primary parameters defined in Tab. 3-1 of the
//! [ICAO WGS 84 Implementation Manual](https://www.icao.int/safety/pbn/Documentation/EUROCONTROL/Eurocontrol%20WGS%2084%20Implementation%20Manual.pdf).
//!
//! ## Geodesic navigation
//!
//! The shortest path between two points on the surface of an ellipsoid is a
//! [geodesic](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid) -
//! the equivalent of straight line segments in planar geometry or
//! [great circles](https://en.wikipedia.org/wiki/Great_circle) on the surface of a
//! sphere, see *Figure 2*.
//!
//! <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Geodesic_problem_on_an_ellipsoid.svg/1024px-Geodesic_problem_on_an_ellipsoid.svg.png" width="400">
//!
//! *Figure 2 A geodesic between points A and B*
//!
//! [Karney(2013)](https://link.springer.com/article/10.1007/s00190-012-0578-z)
//! solves geodesic problems by mapping a geodesic onto the auxiliary sphere
//! and then solving the corresponding problem in
//! [great-circle navigation](https://en.wikipedia.org/wiki/Great-circle_navigation).
//!
//! [Baselga and Martinez-Llario(2017)](https://www.researchgate.net/publication/321358300_Intersection_and_point-to-line_solutions_for_geodesics_on_the_ellipsoid)
//! solve geodesic intersection and point-to-line problems by using the
//! correspondence between geodesics on an ellipsoid and great-circles on the
//! auxiliary sphere.
//!
//! Nick Korbey and I [Barker and Korbey(2019)](https://www.researchgate.net/publication/335749834_Geodesic_Geometry) developed Baselga and Martinez-Llario's algorithms
//! by using vectors to solve geodesic intersection and point-to-line problems on
//! the auxiliary sphere.
//!
//! This library also uses the correspondence between geodesics on an ellipsoid
//! and great-circles on the auxiliary sphere together with 3D vectors to calculate:
//!
//! - the initial azimuth and length of a geodesic between two positions;
//! - the along track distance and across track distance of a position relative to a geodesic;
//! - and the intersection of a pair of geodesics.
//!
//! ## Design
//!
//! The library is based on Charles Karney's [GeographicLib](https://geographiclib.sourceforge.io/) library.
//! Like `GeographicLib`, it models geodesic paths as great circles on
//! the surface of an auxiliary sphere. However, it also uses vectors to
//! calculate along track distances, across track distances and
//! intersections between geodesics.
//!
//! The library depends upon the following crates:
//!
//! - [angle-sc](https://crates.io/crates/angle-sc) - to define `Angle`,
//! `Degrees` and `Radians` and perform trigonometric calculations;
//! - [unit-sphere](https://crates.io/crates/unit-sphere) - to define `LatLong`
//! and perform great-circle and vector calculations.
//! - [icao_units](https://crates.io/crates/icao-units) - to define `Metres` and
//! `NauticalMiles` and perform conversions between them.
//!
//! The library is declared [no_std](https://docs.rust-embedded.org/book/intro/no-std.html)
//! so it can be used in embedded applications.

#![cfg_attr(not(test), no_std)]

extern crate angle_sc;
extern crate icao_units;
extern crate unit_sphere;

pub mod ellipsoid;
pub mod geodesic;

pub use angle_sc::{Angle, Degrees, Radians, Validate};
pub use icao_units::non_si::NauticalMiles;
pub use icao_units::si::Metres;
pub use unit_sphere::LatLong;

use angle_sc::{is_small, trig};
use unit_sphere::Vector3d;

/// The parameters of an Ellipsoid.  
/// The default value is the WGS84 Ellipsoid.
#[derive(Clone, Debug, PartialEq)]
pub struct Ellipsoid {
    /// The Semimajor axis of the ellipsoid.
    a: Metres,
    /// The flattening of the ellipsoid, a ratio.
    f: f64,

    /// The Semiminor axis of the ellipsoid.
    b: Metres,
    /// One minus the flattening ratio.
    one_minus_f: f64,
    /// The reciprocal of one minus the flattening ratio.
    recip_one_minus_f: f64,
    /// The square of the Eccentricity of the ellipsoid.
    e_2: f64,
    /// The square of the second Eccentricity of the ellipsoid.
    ep_2: f64,
    /// The third flattening of the ellipsoid.
    n: f64,

    /// The A3 series ellipsoid::coefficients of the ellipsoid.
    a3: [f64; 6],
    /// The C3x series ellipsoid::coefficients of the ellipsoid.
    c3x: [f64; 15],
}

impl Ellipsoid {
    /// Constructor.
    /// * `a` - the Semimajor axis of the ellipsoid.
    /// * `f` - the flattening of the ellipsoid, a ratio.
    #[must_use]
    pub fn new(a: Metres, f: f64) -> Self {
        let one_minus_f = 1.0 - f;
        let n = ellipsoid::calculate_3rd_flattening(f);
        Self {
            a,
            f,
            b: ellipsoid::calculate_minor_axis(a, f),
            one_minus_f,
            recip_one_minus_f: 1.0 / one_minus_f,
            e_2: ellipsoid::calculate_sq_eccentricity(f),
            ep_2: ellipsoid::calculate_sq_2nd_eccentricity(f),
            n,
            a3: ellipsoid::coefficients::evaluate_coeffs_a3(n),
            c3x: ellipsoid::coefficients::evaluate_coeffs_c3x(n),
        }
    }

    /// Construct an Ellipsoid with the WGS 84 parameters.
    #[must_use]
    pub fn wgs84() -> Self {
        Self::new(ellipsoid::wgs84::A, ellipsoid::wgs84::F)
    }

    /// The Semimajor axis of the ellipsoid.
    #[must_use]
    pub const fn a(&self) -> Metres {
        self.a
    }

    /// The flattening of the ellipsoid, a ratio.
    #[must_use]
    pub const fn f(&self) -> f64 {
        self.f
    }

    /// The Semiminor axis of the ellipsoid.
    #[must_use]
    pub const fn b(&self) -> Metres {
        self.b
    }

    /// One minus the flattening ratio.
    #[must_use]
    pub const fn one_minus_f(&self) -> f64 {
        self.one_minus_f
    }

    /// The reciprocal of one minus the flattening ratio.
    #[must_use]
    pub const fn recip_one_minus_f(&self) -> f64 {
        self.recip_one_minus_f
    }

    /// The square of the Eccentricity of the ellipsoid.
    #[must_use]
    pub const fn e_2(&self) -> f64 {
        self.e_2
    }

    /// The square of the second Eccentricity of the ellipsoid.
    #[must_use]
    pub const fn ep_2(&self) -> f64 {
        self.ep_2
    }

    /// The third flattening of the ellipsoid.
    #[must_use]
    pub const fn n(&self) -> f64 {
        self.n
    }

    /// The A3 series `ellipsoid::coefficients` of the ellipsoid.
    #[must_use]
    pub const fn a3(&self) -> [f64; 6] {
        self.a3
    }

    /// The C3x series `ellipsoid::coefficients` of the ellipsoid.
    #[must_use]
    pub const fn c3x(&self) -> [f64; 15] {
        self.c3x
    }

    /// Calculate epsilon, the variable used in series expansions.  
    /// Note: epsilon is positive and small.
    /// * `clairaut` - Clairaut's constant.
    #[must_use]
    pub fn calculate_epsilon(&self, clairaut: trig::UnitNegRange) -> f64 {
        ellipsoid::calculate_epsilon(clairaut, self.ep_2)
    }

    /// Calculate a3c from the A3 series `ellipsoid::coefficients` of the ellipsoid.
    /// * `clairaut` - Clairaut's constant.
    /// * `eps` - epsilon
    #[must_use]
    pub fn calculate_a3c(&self, clairaut: trig::UnitNegRange, eps: f64) -> f64 {
        self.f * clairaut.0 * ellipsoid::coefficients::evaluate_polynomial(&self.a3, eps)
    }

    /// Calculate the coefficients `C3[l]` in the Fourier expansion of `C3`.
    /// * `eps` - epsilon
    #[must_use]
    pub fn calculate_c3y(&self, eps: f64) -> [f64; 6] {
        ellipsoid::coefficients::evaluate_coeffs_c3y(&self.c3x, eps)
    }

    /// Convert a geodetic Latitude to a parametric Latitude on the
    /// auxiliary sphere.
    /// * `lat` - the geodetic Latitude
    #[must_use]
    pub fn calculate_parametric_latitude(&self, lat: Angle) -> Angle {
        ellipsoid::calculate_parametric_latitude(lat, self.one_minus_f)
    }

    /// Convert a parametric Latitude on the auxiliary sphere to a
    /// geodetic Latitude.
    /// * `beta` - the parametric Latitude
    #[must_use]
    pub fn calculate_geodetic_latitude(&self, beta: Angle) -> Angle {
        ellipsoid::calculate_geodetic_latitude(beta, self.one_minus_f)
    }
}

/// A default Ellipsoid: WGS 84.
impl Default for Ellipsoid {
    /// Construct an Ellipsoid with the WGS 84 parameters.
    #[must_use]
    fn default() -> Self {
        Self::wgs84()
    }
}

/// Calculate the azimuth and geodesic length (in metres) between a pair
/// of positions on the ellipsoid.
///  * `a`, `b` - the start and finish positions in geodetic coordinates.
///
/// returns the azimuth at the start position and the length of the geodesic
/// on the ellipsoid in metres.
///
/// # Examples
/// ```
/// use icao_wgs84::*;
///
/// let wgs84_ellipsoid = Ellipsoid::wgs84();
///
/// let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
/// let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
/// let (azimuth, length) = calculate_azimuth_and_geodesic_length(&istanbul, &washington, &wgs84_ellipsoid);
///
/// let azimuth_degrees = Degrees::from(azimuth);
/// println!("Istanbul-Washington initial azimuth: {:?}", azimuth_degrees.0);
///
/// let distance_nm = NauticalMiles::from(length);
/// println!("Istanbul-Washington distance: {:?}", distance_nm);
/// ```
#[must_use]
pub fn calculate_azimuth_and_geodesic_length(
    a: &LatLong,
    b: &LatLong,
    ellipsoid: &Ellipsoid,
) -> (Angle, Metres) {
    let (alpha1, gc_distance) = geodesic::calculate_azimuth_aux_length(a, b, ellipsoid);
    let beta1 =
        ellipsoid::calculate_parametric_latitude(Angle::from(a.lat()), ellipsoid.one_minus_f());
    (
        alpha1,
        geodesic::convert_radians_to_metres(beta1, alpha1, gc_distance, ellipsoid),
    )
}

/// A geodesic path on the surface of an ellipsoid.  
#[derive(Clone, Debug, PartialEq)]
pub struct Geodesic<'a> {
    /// The parametric start latitude on the auxiliary sphere.
    beta: Angle,
    /// The start longitude.
    lon: Angle,
    /// The start azimuth.
    azi: Angle,
    /// Azimuth at the Equator.
    azi0: Angle,
    /// Great circle distance to the first Equator crossing.
    sigma1: Angle,
    /// The Great Circle length on the auxiliary sphere in radians.
    aux_length: Radians,
    /// integration constant epsilon, derived from Clairaut's constant.
    eps: f64,
    /// constant used to convert geodesic/great circle distances.
    a1: f64,
    /// constant used to convert geodesic/great circle longitudes.
    a3c: f64,
    /// start point geodesic/great circle distance difference.
    b11: Radians,
    /// A reference to the underlying Ellipsoid.
    ellipsoid: &'a Ellipsoid,
}

impl<'a> Validate for Geodesic<'a> {
    /// Test whether a `Geodesic` is valid.
    /// Whether latitude <= 90.0 and `aux_length` is positive and less than PI.
    fn is_valid(&self) -> bool {
        self.beta.cos().0 >= 0.0 && (0.0..core::f64::consts::PI).contains(&self.aux_length.0)
    }
}

impl<'a> Geodesic<'a> {
    /// Construct a `Geodesic`
    /// * `beta` - the start point parametric latitude on the auxiliary sphere.
    /// * `lon` - the start point longitude.
    /// * `azi` - the start azimuth.
    /// * `aux_length` - the Great Circle length on the auxiliary sphere in radians.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    pub fn new(
        beta: Angle,
        lon: Angle,
        azi: Angle,
        aux_length: Radians,
        ellipsoid: &'a Ellipsoid,
    ) -> Geodesic<'a> {
        // Calculate the azimuth at the first Equator crossing
        let clairaut = trig::UnitNegRange(azi.sin().0 * beta.cos().0);
        let cos_azi0 = trig::UnitNegRange(libm::hypot(azi.cos().0, azi.sin().0 * beta.sin().0));
        let azi0 = Angle::new(clairaut, cos_azi0);

        // Calculate the distance to the first Equator crossing
        let cos_omega1 = geodesic::calculate_cos_omega(beta, azi.cos());
        let sigma1 = Angle::from_y_x(beta.sin().0, cos_omega1.0);

        // Calculate eps and c1 for calculating coefficients
        let eps = ellipsoid.calculate_epsilon(clairaut);
        let c1 = ellipsoid::coefficients::evaluate_coeffs_c1(eps);
        Geodesic {
            beta,
            lon,
            azi,
            azi0,
            sigma1,
            aux_length,
            eps,
            a1: ellipsoid::coefficients::evaluate_a1(eps) + 1.0,
            a3c: ellipsoid.calculate_a3c(clairaut, eps),
            b11: ellipsoid::coefficients::sin_cos_series(&c1, sigma1),
            ellipsoid,
        }
    }

    /// Construct a `Geodesic` using the "direct" method.  
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `aux_length` - the Great Circle length on the auxiliary sphere in radians.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    pub fn from_lat_lon_azi_aux_length(
        a: &LatLong,
        azimuth: Angle,
        aux_length: Radians,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        let a_lat = Angle::from(a.lat());
        let a_lon = Angle::from(a.lon());
        Geodesic::new(
            ellipsoid.calculate_parametric_latitude(a_lat),
            a_lon,
            azimuth,
            aux_length,
            ellipsoid,
        )
    }

    /// Construct a `Geodesic` using the "direct" method with the length in metres.  
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `length` - the length on the Ellipsoid in metres.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    pub fn from_lat_lon_azi_length(
        a: &LatLong,
        azimuth: Angle,
        length: Metres,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        let mut arc = Geodesic::from_lat_lon_azi_aux_length(a, azimuth, Radians(0.0), ellipsoid);
        arc.set_aux_length(arc.metres_to_radians(length));
        arc
    }

    /// Construct a `Geodesic` between a pair of positions, the "indirect" method.  
    /// @pre |lat| <= 90.0 degrees.
    /// * `a`, `b` - the start and finish positions in geodetic coordinates.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    pub fn between_positions(a: &LatLong, b: &LatLong, ellipsoid: &'a Ellipsoid) -> Self {
        Self::from((a, b, ellipsoid))
    }

    /// Accessor for the start latitude on the auxiliary sphere.
    #[must_use]
    pub const fn beta(&self) -> Angle {
        self.beta
    }

    /// Accessor for the start longitude.
    #[must_use]
    pub const fn lon(&self) -> Angle {
        self.lon
    }

    /// Accessor for the start azimuth.
    #[must_use]
    pub const fn azi(&self) -> Angle {
        self.azi
    }

    /// Set the `aux_length` of a `Geodesic`
    /// * `aux_length` - the aux length of the `Geodesic`.
    pub fn set_aux_length(&mut self, aux_length: Radians) -> &mut Self {
        self.aux_length = aux_length;
        self
    }

    /// Accessor for the arc length on the auxiliary sphere in radians.
    #[must_use]
    pub const fn aux_length(&self) -> Radians {
        self.aux_length
    }

    /// Accessor for the reference to the underlying Ellipsoid.
    #[must_use]
    pub const fn ellipsoid(&self) -> &Ellipsoid {
        self.ellipsoid
    }

    /// Convert a distance in metres on the ellipsoid to radians on the
    /// auxiliary sphere.
    /// * `distance_m` - the distance in metres along the Geodesic.
    ///
    /// returns the distance in radians on the auxiliary sphere.
    #[must_use]
    pub fn metres_to_radians(&self, distance_m: Metres) -> Radians {
        if is_small(libm::fabs(distance_m.0), core::f64::EPSILON) {
            Radians(0.0)
        } else {
            let tau12 = Radians(distance_m.0 / (self.ellipsoid.b().0 * self.a1));
            let tau_sum = Angle::from(self.b11 + tau12);
            let c1p = ellipsoid::coefficients::evaluate_coeffs_c1p(self.eps);
            let b12 = ellipsoid::coefficients::sin_cos_series(&c1p, self.sigma1 + tau_sum);

            tau12 + b12 + self.b11
        }
    }

    /// Convert a distance in radians on the auxiliary sphere to metres
    /// on the ellipsoid.
    /// * `gc_distance` - the distance in radians on the auxiliary sphere.
    ///
    /// returns the distance in metres on the ellipsoid.
    #[must_use]
    pub fn radians_to_metres(&self, gc_distance: Radians) -> Metres {
        if is_small(libm::fabs(gc_distance.0), core::f64::EPSILON) {
            Metres(0.0)
        } else {
            let sigma_sum = self.sigma1 + Angle::from(gc_distance);
            let c1 = ellipsoid::coefficients::evaluate_coeffs_c1(self.eps);
            let b12 = ellipsoid::coefficients::sin_cos_series(&c1, sigma_sum);
            Metres(self.ellipsoid.b().0 * self.a1 * (gc_distance + b12 - self.b11).0)
        }
    }

    /// Accessor for the length of the Geodesic in metres.
    #[must_use]
    pub fn length(&self) -> Metres {
        self.radians_to_metres(self.aux_length)
    }

    /// Calculate the parametric latitude at the great circle length.
    /// * `length` - the length on the auxiliary sphere as an Angle.
    ///
    /// return the parametric latitude of the position at length.
    #[must_use]
    fn aux_beta(&self, length: Angle) -> Angle {
        let sin_beta = trig::UnitNegRange(
            length.cos().0 * self.beta.sin().0
                + length.sin().0 * self.beta.cos().0 * self.azi.cos().0,
        );
        Angle::new(sin_beta, trig::swap_sin_cos(sin_beta))
    }

    /// Calculate the geodetic latitude at the great circle length.
    /// * `gc_length` - the length in radians on the auxiliary sphere.
    ///
    /// return the geodetic latitude of the position at `gc_length`.
    #[must_use]
    pub fn aux_latitude(&self, gc_length: Radians) -> Angle {
        let length = Angle::from(gc_length);
        self.ellipsoid
            .calculate_geodetic_latitude(self.aux_beta(length))
    }

    /// Calculate the geodetic latitude at the length along the geodesic.
    /// * `length` - the length in metres along the Geodesic.
    ///
    /// return the geodetic latitude of the position at length.
    #[must_use]
    pub fn latitude(&self, length: Metres) -> Angle {
        let gc_length = self.metres_to_radians(length);
        self.aux_latitude(gc_length)
    }

    /// Calculate the azimuth at the great circle length.
    /// * `gc_length` - the length in radians on the auxiliary sphere.
    ///
    /// return the azimuth of the geodesic/great circle at `gc_length`.
    #[must_use]
    pub fn aux_azimuth(&self, gc_length: Radians) -> Angle {
        const MAX_LAT: f64 = 1.0 - 2.0 * core::f64::EPSILON;

        if is_small(libm::fabs(gc_length.0), 2.0 * core::f64::EPSILON) {
            self.azi
        } else {
            let length = Angle::from(gc_length);
            let sigma_sum = self.sigma1 + length;
            let sin_beta = self.azi0.cos().0 * sigma_sum.sin().0;

            // if at North pole, only valid azimuth is due South
            if MAX_LAT < sin_beta {
                Angle::new(trig::UnitNegRange(0.0), trig::UnitNegRange(-1.0))
            } else {
                Angle::from_y_x(self.azi0.sin().0, self.azi0.cos().0 * sigma_sum.cos().0)
            }
        }
    }

    /// Calculate the azimuth at the length along the geodesic.
    /// * `length` - the length in metres along the Geodesic.
    ///
    /// return the azimuth of the geodesic/great circle at length.
    #[must_use]
    pub fn azimuth(&self, length: Metres) -> Angle {
        let gc_length = self.metres_to_radians(length);
        self.aux_azimuth(gc_length)
    }

    /// Calculate the geodesic longitude difference at a great circle length
    /// along the auxiliary sphere.
    /// * `gc_length` - the length in radians on the auxiliary sphere.
    ///
    /// return the longitude difference of the geodesic at `gc_length`.
    #[allow(clippy::similar_names)]
    #[must_use]
    pub fn delta_longitude(&self, gc_length: Radians) -> Angle {
        if libm::fabs(gc_length.0) <= core::f64::EPSILON {
            return Angle::default();
        }
        // The great circle distance from Northward Equator crossing.
        let sigma_sum = self.sigma1 + Angle::from(gc_length);

        // The longitude difference on the auxiliary sphere, omega12.
        let omega12 = Angle::from_y_x(self.azi0.sin().0 * sigma_sum.sin().0, sigma_sum.cos().0)
            - Angle::from_y_x(
                self.azi0.sin().0 * self.beta.sin().0,
                geodesic::calculate_cos_omega(self.beta, self.azi.cos()).0,
            );

        let c3 = self.ellipsoid.calculate_c3y(self.eps);
        let b31 = ellipsoid::coefficients::sin_cos_series(&c3, self.sigma1);
        let b32 = ellipsoid::coefficients::sin_cos_series(&c3, sigma_sum);

        omega12 - Angle::from(Radians(self.a3c * (gc_length.0 + (b32.0 - b31.0))))
    }

    /// Calculate the geodesic longitude at the great circle length along
    /// the auxiliary sphere.
    /// * `gc_length` - the length in radians on the auxiliary sphere.
    ///
    /// return the longitude of the geodesic at `gc_length`.
    #[must_use]
    pub fn aux_longitude(&self, gc_length: Radians) -> Angle {
        self.lon + self.delta_longitude(gc_length)
    }

    /// Calculate the geodesic longitude at the length along the geodesic.
    /// * `length` - the length in metres along the Geodesic.
    ///
    /// return the longitude of the geodesic at length.
    #[must_use]
    pub fn longitude(&self, length: Metres) -> Angle {
        let gc_length = self.metres_to_radians(length);
        self.aux_longitude(gc_length)
    }

    /// Calculate the geodesic `LatLong` at the great circle length along
    /// the auxiliary sphere.
    /// * `gc_length` - the length on the auxiliary sphere, in `Radians`.
    ///
    /// return the `LatLong` of the geodesic position at `gc_length`.
    #[must_use]
    pub fn aux_lat_long(&self, gc_length: Radians) -> LatLong {
        LatLong::new(
            angle_sc::Degrees::from(self.aux_latitude(gc_length)),
            angle_sc::Degrees::from(self.aux_longitude(gc_length)),
        )
    }

    /// Calculate the geodesic `LatLong` at the length along the `Geogesic`.
    /// * `length` - the length in `Metres`.
    ///
    /// return the `LatLong` of the geodesic position at `length`.
    #[must_use]
    pub fn lat_long(&self, length: Metres) -> LatLong {
        let gc_length = self.metres_to_radians(length);
        self.aux_lat_long(gc_length)
    }

    /// Calculate the geodesic point and pole at the length along the geodesic.
    /// * `gc_length` the great circle length on the auxiliary sphere
    ///
    /// returns the point and pole on the auxiliary sphere at `gc_length`.
    #[must_use]
    pub fn aux_point_and_pole(&self, gc_length: Radians) -> (Vector3d, Vector3d) {
        let point = unit_sphere::vector::to_point(self.beta, self.lon);
        let pole = unit_sphere::vector::calculate_pole(self.beta, self.lon, self.azi);
        // If at the start of the geodesic
        if is_small(libm::fabs(gc_length.0), 2.0 * core::f64::EPSILON) {
            (point, pole)
        } else {
            let length = Angle::from(gc_length);
            let beta: Angle = self.aux_beta(length);
            let lon = self.aux_longitude(gc_length);
            let point = unit_sphere::vector::to_point(beta, lon);

            // if point is on a meridional Geodesic use auxiliary sphere point and pole
            if is_small(libm::fabs(self.azi0.sin().0), 2.0 * core::f64::EPSILON) {
                (point, pole)
            } else {
                // Note: point cannot be at North pole, since it is not on a meridional Geodesic
                // Use Karney's method to calculate azimuth.
                let sigma_sum = self.sigma1 + length;
                let azimuth =
                    Angle::from_y_x(self.azi0.sin().0, self.azi0.cos().0 * sigma_sum.cos().0);
                (
                    point,
                    unit_sphere::vector::calculate_pole(beta, lon, azimuth),
                )
            }
        }
    }

    /// Calculate along and across track distances to a position from a geodesic.
    /// * `position` the position as a `LatLong`
    /// * `precision` the required precision
    ///
    /// returns the along and across track distances to the position in `Radians`.
    #[allow(clippy::similar_names)]
    #[must_use]
    pub fn calculate_aux_atd_and_xtd(
        &self,
        position: &LatLong,
        precision: Radians,
    ) -> (Radians, Radians, u32) {
        const MAX_ITERATIONS: u32 = 10;

        // calculate the position as a point on the auxiliary sphere
        let lat = Angle::from(position.lat());
        let beta = self.ellipsoid.calculate_parametric_latitude(lat);
        let lon = Angle::from(position.lon());
        let point = unit_sphere::vector::to_point(beta, lon);

        // calculate the start point and pole of the geodesic on the auxiliary sphere
        let (a, pole) = self.aux_point_and_pole(Radians(0.0));
        let gc_d =
            unit_sphere::great_circle::e2gc_distance(unit_sphere::vector::distance(&a, &point));

        // if the point is close to the start point of the Geodesic
        if is_small(gc_d, precision) {
            (Radians(0.0), Radians(0.0), 0)
        } else {
            // estimate initial along track distance on the auxiliary sphere
            let (mut atd, mut xtd) = unit_sphere::vector::calculate_atd_and_xtd(&a, &pole, &point);

            let mut iterations = 1;
            while iterations < MAX_ITERATIONS {
                // calculate the position and azimuth at atd along the Geodesic
                let lat_x = self.aux_latitude(atd);
                let lon_x = self.aux_longitude(atd);
                let azi_x = self.aux_azimuth(atd);

                // calculate the geodesic azimuth and length to the point from the Geodesic position at atd
                let (azi_p, length) =
                    geodesic::aux_sphere_azimuth_length(lat_x, lat, lon - lon_x, self.ellipsoid);
                let delta_azi = azi_x - azi_p;
                let delta_atd = trig::spherical_cosine_rule(delta_azi.cos(), length);
                atd = atd + delta_atd;
                xtd = length;

                if is_small(libm::fabs(delta_atd.0), precision.0) {
                    break;
                }

                iterations += 1;
            }
            // get the cross track distance (and sign) at the along track distance
            xtd = if is_small(xtd, precision) {
                Radians(0.0)
            } else {
                let (_a, pole) = self.aux_point_and_pole(atd);
                let sign = pole.dot(&point);
                Radians(libm::copysign(xtd.0, sign))
            };
            (atd, xtd, iterations)
        }
    }

    /// Calculate along and across track distances to a position from a geodesic.
    /// * `position` the position as a `LatLong`
    /// * `precision` the required precision
    ///
    /// returns the along and across track distances to the position in `Metres`.
    /// # Examples
    /// ```
    /// use icao_wgs84::*;
    /// use angle_sc::is_within_tolerance;
    ///
    /// let wgs84_ellipsoid = Ellipsoid::wgs84();
    ///
    /// let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
    /// let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
    /// let g1 = Geodesic::between_positions(&istanbul, &washington, &wgs84_ellipsoid);
    ///
    /// let azimuth_degrees = Degrees::from(g1.azimuth(Metres(0.0)));
    /// println!("Istanbul-Washington initial azimuth: {:?}", azimuth_degrees.0);
    ///
    /// let distance_nm = NauticalMiles::from(g1.length());
    /// println!("Istanbul-Washington distance: {:?}", distance_nm);
    ///
    /// let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
    ///
    /// // Calculate geodesic along track and across track distances to 1mm precision.
    /// let (atd, xtd, iterations) = g1.calculate_atd_and_xtd(&reyjavik, Metres(1e-3));
    /// assert!(is_within_tolerance(3928788.572, atd.0, 1e-3));
    /// assert!(is_within_tolerance(-1010585.9988368, xtd.0, 1e-3));
    /// println!("calculate_atd_and_xtd iterations: {:?}", iterations);
    ///
    /// // The expected latitude and longitude are from:
    /// // <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#8a93>
    /// let position = g1.lat_long(atd);
    /// assert!(is_within_tolerance(
    ///     54.92853149711691,
    ///     Degrees::from(position.lat()).0,
    ///     128.0 * core::f64::EPSILON
    /// ));
    /// assert!(is_within_tolerance(
    ///     -21.93729106604878,
    ///     Degrees::from(position.lon()).0,
    ///     2048.0 * core::f64::EPSILON
    /// ));
    /// ```
    #[allow(clippy::similar_names)]
    #[must_use]
    pub fn calculate_atd_and_xtd(
        &self,
        position: &LatLong,
        precision: Metres,
    ) -> (Metres, Metres, u32) {
        // convert precision to Radians
        let precision = Radians(precision.0 / self.ellipsoid.a().0);

        let (atd, xtd, iterations) = self.calculate_aux_atd_and_xtd(position, precision);

        let beta = self.aux_beta(Angle::from(atd));
        let alpha = self.aux_azimuth(atd).quarter_turn_ccw();
        (
            self.radians_to_metres(atd),
            geodesic::convert_radians_to_metres(beta, alpha, xtd, self.ellipsoid),
            iterations,
        )
    }
}

impl<'a> From<(&LatLong, Angle, Radians, &'a Ellipsoid)> for Geodesic<'a> {
    /// Construct a `Geodesic` using the "direct" method.  
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `aux_length` - the Great Circle length on the auxiliary sphere in radians.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    fn from(params: (&LatLong, Angle, Radians, &'a Ellipsoid)) -> Geodesic<'a> {
        Geodesic::from_lat_lon_azi_aux_length(params.0, params.1, params.2, params.3)
    }
}

impl<'a> From<(&LatLong, Angle, Metres, &'a Ellipsoid)> for Geodesic<'a> {
    /// Construct a `Geodesic` using the "direct" method with the length in metres.  
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `length` - the length on the Ellipsoid in metres.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    fn from(params: (&LatLong, Angle, Metres, &'a Ellipsoid)) -> Geodesic<'a> {
        Geodesic::from_lat_lon_azi_length(params.0, params.1, params.2, params.3)
    }
}

impl<'a> From<(&LatLong, &LatLong, &'a Ellipsoid)> for Geodesic<'a> {
    /// Construct a `Geodesic` between a pair of positions, the "indirect" method.  
    /// @pre |lat| <= 90.0 degrees.
    /// * `a`, `b` - the start and finish positions in geodetic coordinates.
    /// * `ellipsoid` - a reference to the Ellipsoid.
    #[must_use]
    fn from(params: (&LatLong, &LatLong, &'a Ellipsoid)) -> Geodesic<'a> {
        let (azimuth, aux_length) =
            geodesic::calculate_azimuth_aux_length(params.0, params.1, params.2);
        Geodesic::from((params.0, azimuth, aux_length, params.2))
    }
}

/// Calculate the distances along a pair of Geodesics (in Radians) to their
/// closest intersection or reference points.
/// * `g1`, `g2` the Geodesics.
///
/// returns the distances along the Geodesics to the intersection point or to
/// their closest (reference) points if the Geodesics do not intersect.
///
/// # Panics
///
/// The function will panic if the Geodesics are **not** on the same `Ellipsoid`.
///
/// # Examples
/// ```
/// use icao_wgs84::*;
/// use angle_sc::is_within_tolerance;
///
/// let wgs84_ellipsoid = Ellipsoid::wgs84();
///
/// let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
/// let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
/// let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
/// let accra = LatLong::new(Degrees(6.0), Degrees(0.0));
///
/// let g1 = Geodesic::from((&istanbul, &washington, &wgs84_ellipsoid));
/// let g2 = Geodesic::from((&reyjavik, &accra, &wgs84_ellipsoid));
///
/// // Calculate distances from the geodesic start points to the intersection point to 1mm precision.
/// let (distance1, _distance2, iterations) =
///     calculate_intersection_distances(&g1, &g2, Metres(1e-3));
/// println!(
///     "calculate_intersection_distances iterations: {:?}",
///     iterations
/// );
///
/// // Get the intersection point position
/// // The expected latitude and longitude are from:
/// // <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a>
/// let lat_lon = g1.aux_lat_long(distance1);
/// assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));
/// assert!(is_within_tolerance(-14.56385574430775, lat_lon.lon().0, 1e-6));
/// ```
#[must_use]
pub fn calculate_intersection_distances(
    g1: &Geodesic,
    g2: &Geodesic,
    precision: Metres,
) -> (Radians, Radians, u32) {
    const MAX_ITERATIONS: u32 = 10;

    // Geodesic MUST be on the same Ellipsoid
    assert!(g1.ellipsoid() == g2.ellipsoid());

    // Convert the precision to Radians
    let precision = Radians(precision.0 / g1.ellipsoid().a().0);

    let a1_lat = g1.beta();
    let a1 = unit_sphere::vector::to_point(a1_lat, g1.lon());

    let a2_lat = g2.beta();
    let a2 = unit_sphere::vector::to_point(a2_lat, g2.lon());
    let gc_d = unit_sphere::great_circle::e2gc_distance(unit_sphere::vector::distance(&a1, &a2));

    // if the start points are within precision of each other
    if is_small(gc_d, precision) {
        (Radians(0.0), Radians(0.0), 0)
    } else {
        let pole1 = unit_sphere::vector::calculate_pole(a1_lat, g1.lon(), g1.azi());
        let pole2 = unit_sphere::vector::calculate_pole(a2_lat, g2.lon(), g2.azi());

        let length1 = g1.aux_length();
        let length2 = g2.aux_length();

        let c = unit_sphere::vector::intersection::calculate_intersection_point(&pole1, &pole2);
        c.map_or_else(
            || {
                let distances =
                    unit_sphere::vector::intersection::calculate_same_gc_reference_distances(
                        &a1, &pole1, length1, &a2, &pole2, length2, gc_d,
                    );
                (distances.0, distances.1, 0)
            },
            |c| {
                let (mut distance1, mut distance2, use_other_intersection) =
                    unit_sphere::vector::intersection::calculate_closest_intersection_distances(
                        &a1, &pole1, length1, &a2, &pole2, length2, &c,
                    );

                if use_other_intersection {
                    distance1 = distance1 + Radians(core::f64::consts::PI);
                    distance2 = distance2 + Radians(core::f64::consts::PI);
                }

                // Convert precision in radians to the square of Euclidean precision.
                let e_precision = unit_sphere::great_circle::gc2e_distance(precision);
                let sq_precision = e_precision * e_precision;

                let mut iterations = 1;

                while iterations < MAX_ITERATIONS {
                    let (pos1, pole1) = g1.aux_point_and_pole(distance1);
                    let (pos2, pole2) = g2.aux_point_and_pole(distance2);

                    let sq_d = unit_sphere::vector::sq_distance(&pos1, &pos2);
                    if is_small(sq_d, sq_precision) {
                        break;
                    }

                    iterations += 1;

                    let x = if use_other_intersection {
                        unit_sphere::vector::intersection::calculate_intersection_point(
                            &pole2, &pole1,
                        )
                    } else {
                        unit_sphere::vector::intersection::calculate_intersection_point(
                            &pole1, &pole2,
                        )
                    };
                    match x {
                        Some(x) => {
                            let (delta1, delta2) =
                                unit_sphere::vector::intersection::calculate_intersection_distances(
                                    &pos1, &pole1, &pos2, &pole2, &x,
                                );
                            distance1 = distance1 + delta1;
                            distance2 = distance2 + delta2;
                        }
                        None => break,
                    }
                }

                (distance1, distance2, iterations)
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use angle_sc::is_within_tolerance;
    use core::mem::size_of;
    use unit_sphere::LatLong;

    #[test]
    fn test_ellipsoid_wgs84() {
        let geoid = Ellipsoid::wgs84();
        assert_eq!(ellipsoid::wgs84::A, geoid.a());
        assert_eq!(ellipsoid::wgs84::F, geoid.f());
        assert_eq!(
            ellipsoid::calculate_minor_axis(ellipsoid::wgs84::A, ellipsoid::wgs84::F),
            geoid.b()
        );
        assert_eq!(1.0 - ellipsoid::wgs84::F, geoid.one_minus_f());
        assert_eq!(1.0 / (1.0 - ellipsoid::wgs84::F), geoid.recip_one_minus_f());
        assert_eq!(
            ellipsoid::calculate_sq_eccentricity(ellipsoid::wgs84::F),
            geoid.e_2()
        );
        assert_eq!(
            ellipsoid::calculate_sq_2nd_eccentricity(ellipsoid::wgs84::F),
            geoid.ep_2()
        );
        assert_eq!(
            ellipsoid::calculate_3rd_flattening(ellipsoid::wgs84::F),
            geoid.n()
        );
    }

    #[test]
    fn test_ellipsoid_traits() {
        let geoid = Ellipsoid::default();

        let geoid_clone = geoid.clone();
        assert!(geoid_clone == geoid);

        println!("Ellipsoid: {:?}", geoid);
    }

    #[test]
    fn test_calculate_azimuth_aux_length_normal_05() {
        // GeodTest.dat line 2874
        // 5.421025561218 0 84.846843174846
        // 3.027329237478900117 109.666857465735641205 96.826992198613537236
        // 12161089.9991805 109.607910081857488806 5988906.6319258056178 8449589948776.249238
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // North East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(5.421025561218), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(3.027329237478900117),
            Degrees(109.666857465735641205),
        );

        let result = calculate_azimuth_and_geodesic_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert_eq!(84.846843174846, Degrees::from(result.0).0);
        assert!(is_within_tolerance(12161089.9991805, (result.1).0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuth_aux_length_nearly_antipodal_1() {
        // GeodTest.dat line 100001
        // 8.226828747671 0 111.1269645725
        // -8.516119211674268968 178.688979582629224039 68.982798544955243193
        // 19886305.6710041 179.197987814300505446 97496.4436255989712 -29736790544759.340534
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        let latlon1 = LatLong::new(Degrees(8.226828747671), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-8.516119211674268968),
            Degrees(178.688979582629224039),
        );

        let result = calculate_azimuth_and_geodesic_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert!(is_within_tolerance(
            111.1269645725,
            Degrees::from(result.0).0,
            1e-9
        ));
        assert!(is_within_tolerance(19886305.6710041, (result.1).0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuth_aux_length_nearly_antipodal_2() {
        // GeodTest.dat line 100017
        // .322440123063 0 100.319048368176
        // -.367465171996537868 179.160624688175359763 79.682430612745621077
        // 19943611.6727803 179.749470297545372441 29954.0028615773743 -14555544282075.683105
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        let latlon1 = LatLong::new(Degrees(0.322440123063), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-0.367465171996537868),
            Degrees(179.160624688175359763),
        );

        let result = calculate_azimuth_and_geodesic_length(&latlon1, &latlon2, &wgs84_ellipsoid);
        assert!(is_within_tolerance(
            100.319048368176,
            Degrees::from(result.0).0,
            1e-9
        ));
        assert!(is_within_tolerance(19943611.6727803, (result.1).0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuth_and_geodesic_length_karney() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let (azimuth, length) =
            calculate_azimuth_and_geodesic_length(&istanbul, &washington, &wgs84_ellipsoid);

        let azimuth_degrees = Degrees::from(azimuth);
        assert_eq!(-50.69375304113997, azimuth_degrees.0);
        assert_eq!(8339863.136005359, length.0);

        println!("Istanbul-Washington azimuth: {:?}", azimuth_degrees.0);

        let distance_nm = NauticalMiles::from(length);
        println!("Istanbul-Washington distance: {:?}", distance_nm);
    }

    #[test]
    fn test_geodesicarc_direct_constructors() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();
        let length = Metres(9_000_000.0);

        // Ensure that two Geodesics can fit on a cache line.
        assert_eq!(128, size_of::<Geodesic>());

        let a = LatLong::new(Degrees(45.0), Degrees(45.0));

        // Increase azimuth around compass from due South to due North
        for i in -180..180 {
            let azi = i as f64;
            let azimuth = Angle::from(Degrees(azi));

            let geodesic1 = Geodesic::from((&a, azimuth, length, &wgs84_ellipsoid));
            assert!(geodesic1.is_valid());
            let azi0 = geodesic1.azimuth(Metres(0.0));
            assert!(is_within_tolerance(
                Radians::from(azimuth).0,
                Radians::from(azi0).0,
                2.0 * std::f64::EPSILON
            ));

            let len0 = geodesic1.length();
            assert!(is_within_tolerance(length.0, len0.0, 1.0e-8));
        }
    }

    #[test]
    fn test_geodesicarc_between_positions() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));

        let g1 = Geodesic::between_positions(&istanbul, &washington, &wgs84_ellipsoid);
        assert!(g1.is_valid());

        let end_azimuth = Degrees::from(g1.azimuth(g1.length()));
        assert_eq!(-132.2646607116376, end_azimuth.0);

        let mut g1_clone = g1.clone();
        assert_eq!(g1_clone, g1);

        let g1_clone = g1_clone.set_aux_length(Radians(1.0));
        println!("Geodesic: {:?}", &g1_clone);

        // test start position
        assert!(is_within_tolerance(
            42.0,
            Degrees::from(wgs84_ellipsoid.calculate_geodetic_latitude(g1.beta())).0,
            32.0 * core::f64::EPSILON
        ));
        assert!(is_within_tolerance(
            29.0,
            Degrees::from(g1.lon()).0,
            16.0 * core::f64::EPSILON
        ));

        let lat_long: LatLong = g1.lat_long(Metres(0.0));
        // assert_eq!(istanbul, lat_long);
        assert!(is_within_tolerance(
            istanbul.lat().0,
            lat_long.lat().0,
            32.0 * core::f64::EPSILON
        ));
        assert!(is_within_tolerance(
            istanbul.lon().0,
            lat_long.lon().0,
            32.0 * core::f64::EPSILON
        ));

        let (atd, xtd, iterations) = g1.calculate_atd_and_xtd(&istanbul, Metres(1e-3));
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);
        assert_eq!(0.0, atd.0);
        assert_eq!(0.0, xtd.0);

        // test end position
        let aux_length = g1.aux_length();
        assert_eq!(1.309412846249522, aux_length.0);

        let length = g1.length();
        assert_eq!(8339863.136005359, length.0);

        // let end_position = g1.aux_lat_long(aux_length);
        assert!(is_within_tolerance(
            39.0,
            Degrees::from(g1.latitude(length)).0,
            32.0 * core::f64::EPSILON
        ));
        assert!(is_within_tolerance(
            -77.0,
            Degrees::from(g1.longitude(length)).0,
            64.0 * core::f64::EPSILON
        ));

        // test mid position
        let half_length = Metres(0.5 * length.0);
        let mid_position = g1.lat_long(half_length);

        assert!(is_within_tolerance(
            54.86379153725445,
            Degrees::from(mid_position.lat()).0,
            64.0 * core::f64::EPSILON
        ));

        assert!(is_within_tolerance(
            -25.694568908316413,
            Degrees::from(mid_position.lon()).0,
            32.0 * core::f64::EPSILON
        ));

        let mid_length = g1.metres_to_radians(half_length);
        assert_eq!(0.654673165141749, mid_length.0);

        let precision = Radians(1e-3 / wgs84_ellipsoid.a().0);
        let (atd, xtd, iterations) = g1.calculate_aux_atd_and_xtd(&mid_position, precision);
        assert!(is_within_tolerance(mid_length.0, atd.0, precision.0));
        assert_eq!(0.0, xtd.0);
        println!("calculate_aux_atd_and_xtd iterations: {:?}", iterations);

        let (atd, xtd, iterations) = g1.calculate_atd_and_xtd(&mid_position, Metres(1e-3));
        assert!(is_within_tolerance(half_length.0, atd.0, 1e-3));
        assert!(is_small(libm::fabs(xtd.0), 1e-3));
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);
    }

    #[test]
    fn test_meridonal_geodesicarc() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // A Geodesic along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(45.0), Degrees(180.0));
        let g1 = Geodesic::from((&a, &b, &wgs84_ellipsoid));
        let (_point, pole0) = g1.aux_point_and_pole(Radians(0.0));

        // Calculate the azimuth at the North pole
        let mid_length = Radians(0.5 * g1.aux_length().0);
        let azimuth = Degrees::from(g1.aux_azimuth(mid_length));
        assert_eq!(180.0, azimuth.0);

        // Calculate the point and great circle pole at the North pole
        let (point1, pole1) = g1.aux_point_and_pole(mid_length);
        assert_eq!(Vector3d::new(0.0, 0.0, 1.0), point1);
        assert_eq!(pole0, pole1);
    }

    #[test]
    fn test_calculate_atd_and_xtd() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // Karney's example
        // Istanbul, Washington and Reyjavik
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let g1 = Geodesic::from((&istanbul, &washington, &wgs84_ellipsoid));

        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));

        // Calculate geodesic along track and across track distances to 1mm precision.
        let (atd, xtd, iterations) = g1.calculate_atd_and_xtd(&reyjavik, Metres(1e-3));
        assert!(is_within_tolerance(3928788.572, atd.0, 1e-3));
        assert!(is_within_tolerance(-1010585.9988368, xtd.0, 1e-3));
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);

        // Karney's latitude and longitude from Final result at:
        // https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#8a93
        let position = g1.lat_long(atd);
        assert!(is_within_tolerance(
            54.92853149711691,
            Degrees::from(position.lat()).0,
            128.0 * core::f64::EPSILON
        ));
        assert!(is_within_tolerance(
            -21.93729106604878,
            Degrees::from(position.lon()).0,
            2048.0 * core::f64::EPSILON
        ));

        // Test delta_azimuth at interception, should be PI/2
        let azimuth_1 = g1.azimuth(atd);
        let g2 = Geodesic::from((&position, &reyjavik, &wgs84_ellipsoid));
        let azimuth_2 = g2.azi();
        let delta_azimuth = azimuth_2 - azimuth_1;
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            Radians::from(delta_azimuth).0,
            128.0 * core::f64::EPSILON
        ));

        // opposite geodesic
        let g1 = Geodesic::from((&washington, &istanbul, &wgs84_ellipsoid));
        let (atd, xtd, iterations) = g1.calculate_atd_and_xtd(&reyjavik, Metres(1e-3));
        assert!(is_within_tolerance(
            g1.length().0 - 3928788.572,
            atd.0,
            1e-3
        ));
        assert!(is_within_tolerance(1010585.9988368, xtd.0, 1e-3));
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);
    }

    #[test]
    fn test_intersection_point_distance() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // Karney's example
        // Istanbul, Washington, Reyjavik and Accra
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
        let accra = LatLong::new(Degrees(6.0), Degrees(0.0));

        let g1 = Geodesic::from((&istanbul, &washington, &wgs84_ellipsoid));
        let g2 = Geodesic::from((&reyjavik, &accra, &wgs84_ellipsoid));

        let (distance1, distance2, iterations) =
            calculate_intersection_distances(&g1, &g2, Metres(1e-3));
        println!(
            "calculate_intersection_distances iterations: {:?}",
            iterations
        );
        assert_eq!(0.5423772210660033, distance1.0);
        assert_eq!(0.17504915893240863, distance2.0);

        let lat_lon = g1.aux_lat_long(distance1);
        // Geodesic intersection latitude is 54.7170296089477
        assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));

        // Geodesic intersection longitude is -14.56385574430775
        assert!(is_within_tolerance(
            -14.56385574430775,
            lat_lon.lon().0,
            1e-6
        ));

        let lat_lon = g2.aux_lat_long(distance2);
        // Geodesic intersection latitude is 54.7170296089477
        assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));

        // Geodesic intersection longitude is -14.56385574430775
        assert!(is_within_tolerance(
            -14.56385574430775,
            lat_lon.lon().0,
            1e-6
        ));

        // Swap geodesics
        let (distance1, distance2, iterations) =
            calculate_intersection_distances(&g2, &g1, Metres(1e-3));
        println!(
            "calculate_intersection_distances iterations: {:?}",
            iterations
        );
        assert_eq!(0.1750491589330442, distance1.0);
        assert_eq!(0.5423772210601135, distance2.0);

        let lat_lon = g2.aux_lat_long(distance1);
        // Geodesic intersection latitude is 54.7170296089477
        assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));

        // Geodesic intersection longitude is -14.56385574430775
        assert!(is_within_tolerance(
            -14.56385574430775,
            lat_lon.lon().0,
            1e-6
        ));

        let lat_lon = g1.aux_lat_long(distance2);
        // Geodesic intersection latitude is 54.7170296089477
        assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));

        // Geodesic intersection longitude is -14.56385574430775
        assert!(is_within_tolerance(
            -14.56385574430775,
            lat_lon.lon().0,
            1e-6
        ));
    }

    #[test]
    fn test_intersection_same_start_point() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // A Geodesic along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(45.0), Degrees(180.0));
        let g1 = Geodesic::from((&a, &b, &wgs84_ellipsoid));

        let c = LatLong::new(Degrees(45.0), Degrees(45.0));
        let g2 = Geodesic::from((&a, &c, &wgs84_ellipsoid));

        let (distance1, distance2, iterations) =
            calculate_intersection_distances(&g1, &g2, Metres(1e-3));
        println!(
            "calculate_intersection_distances iterations: {:?}",
            iterations
        );
        assert_eq!(0.0, distance1.0);
        assert_eq!(0.0, distance2.0);
    }

    #[test]
    fn test_non_intersection_same_geodesic() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // A Geodesic along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(50.0), Degrees(0.0));
        let g1 = Geodesic::from((&a, &b, &wgs84_ellipsoid));

        let c = LatLong::new(Degrees(55.0), Degrees(0.0));
        let d = LatLong::new(Degrees(60.0), Degrees(0.0));
        let g2 = Geodesic::from((&c, &d, &wgs84_ellipsoid));

        let (distance1, distance2, iterations) =
            calculate_intersection_distances(&g1, &g2, Metres(1e-3));
        assert_eq!(0.17463328753266863, distance1.0);
        assert_eq!(0.0, distance2.0);
        assert_eq!(0, iterations);
    }
}
