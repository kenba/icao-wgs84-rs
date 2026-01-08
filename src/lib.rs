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

//! icao-wgs84
//!
//! [![crates.io](https://img.shields.io/crates/v/icao-wgs84.svg)](https://crates.io/crates/icao-wgs84)
//! [![docs.io](https://docs.rs/icao-wgs84/badge.svg)](https://docs.rs/icao-wgs84/)
//! [![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
//! [![Rust](https://github.com/kenba/icao-wgs84-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/kenba/icao-wgs84-rs/actions)
//! [![codecov](https://codecov.io/gh/kenba/icao-wgs84-rs/graph/badge.svg?token=85TJX5VAHF)](https://codecov.io/gh/kenba/icao-wgs84-rs)
//!
//! A library for performing geometric calculations on the
//! [WGS-84](https://www.icao.int/NACC/Documents/Meetings/2014/ECARAIM/REF08-Doc9674.pdf)
//! ellipsoid, see *Figure 1*.
//!
//! <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/3e/WGS84_mean_Earth_radius.svg/800px-WGS84_mean_Earth_radius.svg.png" width="400">
//!
//! *Figure 1 The WGS-84 Ellipsoid (not to scale)
//! [Cmglee](https://commons.wikimedia.org/wiki/User:Cmglee), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), via Wikimedia Commons*
//!
//! [WGS-84](https://www.icao.int/NACC/Documents/Meetings/2014/ECARAIM/REF08-Doc9674.pdf)
//! has become the de facto standard for satellite navigation since its adoption
//! by the Navstar Global Positioning System
//! ([GPS](https://www.gps.gov/systems/gps/performance/accuracy/))
//! and the USA making GPS available for civilian use in 1983.
//!
//! This library uses the WGS-84 primary parameters defined in Tab. 3-1 of the
//! [ICAO WGS-84 Implementation Manual](https://www.icao.int/NACC/Documents/Meetings/2014/ECARAIM/REF08-Doc9674.pdf).
//!
//! ## Geodesic navigation
//!
//! The shortest path between two points on the surface of an ellipsoid is a
//! [geodesic segment](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid).
//! It is the equivalent of a straight line segment in planar geometry or a
//! [great circle arc](https://en.wikipedia.org/wiki/Great_circle) on the
//! surface of a sphere, see *Figure 2*.
//!
//! <img src="https://via-technology.aero/img/navigation/ellipsoid/sphere_mercator_long_geodesic.png" width="400">
//!
//! *Figure 2 A geodesic segment (orange) and a great circle  arc (blue)*
//!
//! This library uses the correspondence between geodesic segments on an ellipsoid
//! and great-circle arcs on a unit sphere, together with 3D vectors to calculate:
//!
//! - the length and azimuths of a geodesic segment between two positions;
//! - the along track and across track distances of a point relative to a geodesic segment;
//! - and the intersection of two geodesic segments.
//!
//! See: [geodesic algorithms](https://via-technology.aero/navigation/geodesic-algorithms/).
//!
//! ## Design
//!
//! The library is based on Charles Karney's [GeographicLib](https://geographiclib.sourceforge.io/) library.
//!
//! Like `GeographicLib`, it models geodesic segments as great circle arcs on
//! the surface of a unit sphere. However, it also uses vectors to perform
//! calculations between geodesic segments.
//!
//! The `Ellipsoid` class represents an ellipsoid of revolution.
//! The static `WGS84_ELLIPSOID` represents the WGS-84 `Ellipsoid` which is used
//! by the `GeodesicSegment` `From` traits to create `GeodesicSegment`s on the WGS-84 `Ellipsoid`.
//!
//! The library depends upon the following crates:
//!
//! - [angle-sc](https://crates.io/crates/angle-sc) - to define `Angle`,
//!   `Degrees` and `Radians` and perform trigonometric calculations;
//! - [unit-sphere](https://crates.io/crates/unit-sphere) - to define `LatLong`
//!   and perform great-circle and vector calculations.
//! - [icao_units](https://crates.io/crates/icao-units) - to define `Metres` and
//!   `NauticalMiles` and perform conversions between them.
//!
//! <img src="https://via-technology.aero/img/software/ellipsoid_class_diagram.svg" width="400">
//!
//! *Figure 3 Class Diagram*
//!
//! The library is declared [no_std](https://docs.rust-embedded.org/book/intro/no-std.html)
//! so it can be used in embedded applications.

#![cfg_attr(not(test), no_std)]

extern crate angle_sc;
extern crate icao_units;
extern crate unit_sphere;

pub mod ellipsoid;
pub mod geodesic;
pub mod intersection;

pub use angle_sc::{Angle, Degrees, Radians, Validate};
pub use icao_units::non_si::NauticalMiles;
pub use icao_units::si::Metres;
pub use unit_sphere::LatLong;

use angle_sc::trig;
use once_cell::sync::Lazy;
use unit_sphere::{Vector3d, great_circle};

/// The parameters of an `Ellipsoid`.
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

    /// The A3 series `coefficients` of the ellipsoid.
    a3: [f64; 6],
    /// The C3x series `coefficients` of the ellipsoid.
    c3x: [f64; 15],
}

impl Ellipsoid {
    /// Constructor.
    /// * `a` - the Semimajor axis of the `Ellipsoid`.
    /// * `f` - the flattening of the `Ellipsoid`, a ratio.
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

    /// Construct an `Ellipsoid` with the WGS-84 parameters.
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

    /// Calculate epsilon, the variable used in series expansions.
    /// Note: epsilon is positive and small.
    /// * `clairaut` - Clairaut's constant.
    #[must_use]
    pub fn calculate_epsilon(&self, clairaut: trig::UnitNegRange) -> f64 {
        ellipsoid::calculate_epsilon(clairaut, self.ep_2)
    }

    /// Calculate a3f from the A3 series `coefficients` of the ellipsoid.
    /// * `eps` - epsilon
    #[must_use]
    pub fn calculate_a3f(&self, eps: f64) -> f64 {
        ellipsoid::coefficients::evaluate_polynomial(&self.a3, eps)
    }

    /// Calculate a3c from the A3 series `coefficients` of the ellipsoid.
    /// * `clairaut` - Clairaut's constant.
    /// * `eps` - epsilon
    #[must_use]
    pub fn calculate_a3c(&self, clairaut: trig::UnitNegRange, eps: f64) -> f64 {
        self.f * clairaut.0 * self.calculate_a3f(eps)
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

    /// Convert a geodetic latitude and longitude to a point on the
    /// auxiliary sphere.
    /// @pre |lat| <= 90.0 degrees.
    /// * `lat` - the latitude.
    /// * `lon` - the longitude.
    ///
    /// returns a `Vector3d` of the point on the auxiliary sphere.
    #[must_use]
    pub fn to_arc_point(&self, lat: Angle, lon: Angle) -> Vector3d {
        let beta = self.calculate_parametric_latitude(lat);
        unit_sphere::vector::to_point(beta, lon)
    }
}

/// A static instance of the WGS-84 `Ellipsoid`.
pub static WGS84_ELLIPSOID: Lazy<Ellipsoid> = Lazy::new(Ellipsoid::wgs84);

/// Calculate the azimuths and geodesic length (in metres) between a pair
/// of positions on the ellipsoid.
/// * `a`, `b` - the start and finish positions in geodetic coordinates.
/// * `tolerance` - the tolerance to perform the calculation to.
/// * `ellipsoid` - the `Ellipsoid`.
///
/// returns the azimuth at the start and end positions and the length of
/// the geodesic segment on the ellipsoid in metres.
///
/// # Examples
/// ```
/// use icao_wgs84::*;
/// use unit_sphere::great_circle;
///
/// let tolerance = Radians(great_circle::MIN_VALUE);
///
/// let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
/// let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
/// let (azimuth, length, end_azimuth) = calculate_azimuths_and_geodesic_length(&istanbul, &washington, tolerance, &WGS84_ELLIPSOID);
///
/// let azimuth_degrees = Degrees::from(azimuth);
/// println!("Istanbul-Washington initial azimuth: {:?}", azimuth_degrees.0);
///
/// let distance_nm = NauticalMiles::from(length);
/// println!("Istanbul-Washington distance: {:?}", distance_nm);
///
/// let azimuth_degrees = Degrees::from(end_azimuth.opposite());
/// println!("Washington-Istanbul initial azimuth: {:?}", azimuth_degrees.0);
#[must_use]
pub fn calculate_azimuths_and_geodesic_length(
    a: &LatLong,
    b: &LatLong,
    tolerance: Radians,
    ellipsoid: &Ellipsoid,
) -> (Angle, Metres, Angle) {
    let (alpha1, arc_length, alpha2, _) =
        geodesic::calculate_azimuths_arc_length(a, b, tolerance, ellipsoid);
    let beta1 =
        ellipsoid::calculate_parametric_latitude(Angle::from(a.lat()), ellipsoid.one_minus_f());
    (
        alpha1,
        geodesic::convert_radians_to_metres(beta1, alpha1, arc_length, ellipsoid),
        alpha2,
    )
}

/// A geodesic segment on the surface of an ellipsoid.
///
/// A geodesic segment on an ellipsoid is the shortest path between two points.
/// It is represented by a great circle arc on the auxiliary sphere.
#[derive(Clone, Debug, PartialEq)]
pub struct GeodesicSegment<'a> {
    /// The parametric start latitude on the auxiliary sphere.
    beta: Angle,
    /// The start longitude.
    lon: Angle,
    /// The start azimuth.
    azi: Angle,
    /// Azimuth at the Equator.
    azi0: Angle,
    /// Great circle arc distance to the first Equator crossing.
    sigma1: Angle,
    /// Great circle arc length on the auxiliary sphere in radians.
    arc_length: Radians,
    /// The half width of a Geodesic Rectangle in metres.
    half_width: Metres,
    /// Integration constant: epsilon, derived from Clairaut's constant.
    eps: f64,
    /// constant used to convert geodesic/great circle longitudes.
    a3c: f64,
    /// Start parameter for geodesic/great circle distance differences.
    b11: Radians,
    /// A reference to the underlying `Ellipsoid`.
    ellipsoid: &'a Ellipsoid,
}

impl Validate for GeodesicSegment<'_> {
    /// Test whether a `GeodesicSegment` is valid.
    /// Whether 0° <= `latitude` <= 90° and 0 <= `arc_length` <= π.
    fn is_valid(&self) -> bool {
        self.beta.cos().0 >= 0.0
            && (0.0..core::f64::consts::PI).contains(&self.arc_length.0)
            && self.half_width.0 >= 0.0
    }
}

impl<'a> GeodesicSegment<'a> {
    /// Construct a `GeodesicSegment`
    /// * `beta` - the start point parametric latitude on the auxiliary sphere.
    /// * `lon` - the start point longitude.
    /// * `azi` - the start azimuth.
    /// * `arc_length` - the great circle arc length on the auxiliary sphere in radians.
    /// * `half_width` - the `GeodesicSegment` half width in Metres.
    /// * `ellipsoid` - a reference to the `Ellipsoid`.
    #[must_use]
    pub fn new(
        beta: Angle,
        lon: Angle,
        azi: Angle,
        arc_length: Radians,
        half_width: Metres,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        // Calculate the azimuth at the first Equator crossing
        let clairaut = trig::UnitNegRange(azi.sin().0 * beta.cos().0);
        let azi0 = Angle::new(clairaut, trig::swap_sin_cos(clairaut));

        // Calculate the distance to the first Equator crossing
        let sigma1 = Angle::from_y_x(beta.sin().0, beta.cos().0 * azi.cos().0);

        // Calculate eps and c1 for calculating coefficients
        let eps = ellipsoid.calculate_epsilon(azi0.sin());
        let c1 = ellipsoid::coefficients::evaluate_coeffs_c1(eps);
        Self {
            beta,
            lon,
            azi,
            azi0,
            sigma1,
            arc_length,
            half_width,
            eps,
            a3c: ellipsoid.calculate_a3c(azi0.sin(), eps),
            b11: ellipsoid::coefficients::sin_cos_series(&c1, sigma1),
            ellipsoid,
        }
    }

    /// Construct a `GeodesicSegment` using the "direct" method.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `arc_length` - the Great Circle arc length on the auxiliary sphere in radians.
    /// * `half_width` - the `GeodesicSegment` half width in Metres.
    /// * `ellipsoid` - a reference to the `Ellipsoid`.
    #[must_use]
    pub fn from_lat_lon_azi_arc_length_half_width(
        a: &LatLong,
        azimuth: Angle,
        arc_length: Radians,
        half_width: Metres,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        let a_lat = Angle::from(a.lat());
        let a_lon = Angle::from(a.lon());
        GeodesicSegment::new(
            ellipsoid.calculate_parametric_latitude(a_lat),
            a_lon,
            azimuth,
            arc_length,
            half_width,
            ellipsoid,
        )
    }

    /// Construct a `GeodesicSegment` using the "direct" method.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `arc_length` - the Great Circle arc length on the auxiliary sphere in radians.
    /// * `ellipsoid` - a reference to the `Ellipsoid`.
    #[must_use]
    pub fn from_lat_lon_azi_arc_length(
        a: &LatLong,
        azimuth: Angle,
        arc_length: Radians,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        GeodesicSegment::from_lat_lon_azi_arc_length_half_width(
            a,
            azimuth,
            arc_length,
            Metres(0.0),
            ellipsoid,
        )
    }

    /// Construct a `GeodesicSegment` using the "direct" method with the length in metres.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `length` - the length on the `Ellipsoid` in metres.
    /// * `ellipsoid` - a reference to the `Ellipsoid`.
    #[must_use]
    pub fn from_lat_lon_azi_length(
        a: &LatLong,
        azimuth: Angle,
        length: Metres,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        let mut arc =
            GeodesicSegment::from_lat_lon_azi_arc_length(a, azimuth, Radians(0.0), ellipsoid);
        arc.set_arc_length(arc.metres_to_radians(length));
        arc
    }

    /// Construct a `GeodesicSegment` between a pair of positions, the "indirect" method.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a`, `b` - the start and finish positions in geodetic coordinates.
    /// * `half_width` - the `GeodesicSegment` half width in Metres.
    /// * `tolerance` - the tolerance to perform the calculation to.
    /// * `ellipsoid` - a reference to the `Ellipsoid`.
    #[must_use]
    pub fn between_positions(
        a: &LatLong,
        b: &LatLong,
        half_width: Metres,
        tolerance: Radians,
        ellipsoid: &'a Ellipsoid,
    ) -> Self {
        let (azimuth, arc_length, _, _) =
            geodesic::calculate_azimuths_arc_length(a, b, tolerance, ellipsoid);
        let a_lat = Angle::from(a.lat());
        // if a is at the North or South pole
        if a_lat.cos().0 < great_circle::MIN_VALUE {
            // use b's longitude
            Self::from_lat_lon_azi_arc_length_half_width(
                &LatLong::new(a.lat(), b.lon()),
                azimuth,
                arc_length,
                half_width,
                ellipsoid,
            )
        } else {
            Self::from_lat_lon_azi_arc_length_half_width(
                a, azimuth, arc_length, half_width, ellipsoid,
            )
        }
    }

    /// Accessor for the start parametric latitude on the auxiliary sphere.
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

    /// Set the `arc_length` of a `GeodesicSegment`
    /// * `arc_length` - the great circle arc length of the `GeodesicSegment`.
    pub const fn set_arc_length(&mut self, arc_length: Radians) -> &mut Self {
        self.arc_length = arc_length;
        self
    }

    /// Accessor for the arc length on the auxiliary sphere in radians.
    #[must_use]
    pub const fn arc_length(&self) -> Radians {
        self.arc_length
    }

    /// Method to set the half width in metres.
    /// Set the `half_width` of a `GeodesicSegment`
    /// * `half_width` - the half width of the `GeodesicSegment`.
    pub const fn set_half_width(&mut self, half_width: Metres) -> &mut Self {
        self.half_width = half_width;
        self
    }

    /// Accessor for the half width in metres.
    #[must_use]
    pub const fn half_width(&self) -> Metres {
        self.half_width
    }

    /// Accessor for the reference to the underlying `Ellipsoid`.
    #[must_use]
    pub const fn ellipsoid(&self) -> &Ellipsoid {
        self.ellipsoid
    }

    /// Accessor for the start point on the unit sphere.
    #[must_use]
    pub fn a(&self) -> Vector3d {
        unit_sphere::vector::to_point(self.beta, self.lon)
    }

    /// Convert a distance in metres on the ellipsoid to radians on the
    /// auxiliary sphere.
    /// * `distance` - the distance along the `GeodesicSegment` in metres.
    ///
    /// returns the distance along the great circle arc in radians.
    #[must_use]
    pub fn metres_to_radians(&self, distance: Metres) -> Radians {
        if distance.0.abs() < great_circle::MIN_VALUE {
            Radians(0.0)
        } else {
            let a1 = ellipsoid::coefficients::evaluate_a1(self.eps) + 1.0;
            let tau12 = Radians(distance.0 / (self.ellipsoid.b().0 * a1));
            let tau_sum = Angle::from(self.b11 + tau12);
            let c1p = ellipsoid::coefficients::evaluate_coeffs_c1p(self.eps);
            let b12 = ellipsoid::coefficients::sin_cos_series(&c1p, self.sigma1 + tau_sum);

            tau12 + b12 + self.b11
        }
    }

    /// Convert a great circle distance in radians on the auxiliary sphere to metres
    /// on the ellipsoid.
    /// * `arc_distance` - the great circle distance in radians on the auxiliary sphere.
    /// * `sigma` the `arc_distance` as an `Angle`.
    ///
    /// returns the distance in metres on the ellipsoid.
    #[must_use]
    pub fn radians_to_metres(&self, arc_distance: Radians, sigma: Angle) -> Metres {
        let sigma_sum = self.sigma1 + sigma;
        let c1 = ellipsoid::coefficients::evaluate_coeffs_c1(self.eps);
        let b12 = ellipsoid::coefficients::sin_cos_series(&c1, sigma_sum);
        let a1 = ellipsoid::coefficients::evaluate_a1(self.eps) + 1.0;
        Metres(self.ellipsoid.b().0 * a1 * (arc_distance + b12 - self.b11).0)
    }

    /// Accessor for the length of the `GeodesicSegment` in metres.
    #[must_use]
    pub fn length(&self) -> Metres {
        self.radians_to_metres(self.arc_length, Angle::from(self.arc_length))
    }

    /// Calculate the parametric latitude at the great circle length.
    /// * `sigma` - the arc distance on the auxiliary sphere as an Angle.
    ///
    /// return the parametric latitude of the position at sigma.
    #[must_use]
    pub fn arc_beta(&self, sigma: Angle) -> Angle {
        great_circle::calculate_latitude(self.beta, self.azi, sigma)
    }

    /// Calculate the geodetic latitude at the great circle arc distance.
    /// * `arc_distance` - the great circle arc distance on the auxiliary sphere.
    ///
    /// return the geodetic latitude of the position at `arc_distance`.
    #[must_use]
    pub fn arc_latitude(&self, arc_distance: Radians) -> Angle {
        let sigma = Angle::from(arc_distance);
        self.ellipsoid
            .calculate_geodetic_latitude(self.arc_beta(sigma))
    }

    /// Calculate the geodetic latitude at the length along the geodesic.
    /// * `distance` - the distance along the `GeodesicSegment`, in metres.
    ///
    /// return the geodetic latitude of the position at distance.
    #[must_use]
    pub fn latitude(&self, distance: Metres) -> Angle {
        let arc_distance = self.metres_to_radians(distance);
        self.arc_latitude(arc_distance)
    }

    /// Calculate the azimuth at the great circle length.
    /// * `sigma` - the arc distance on the auxiliary sphere as an Angle.
    ///
    /// return the azimuth at `sigma`.
    #[must_use]
    pub fn arc_azimuth(&self, sigma: Angle) -> Angle {
        const MAX_LAT: f64 = 1.0 - great_circle::MIN_VALUE;

        let sigma_sum = self.sigma1 + sigma;
        let sin_beta = self.azi0.cos().0 * sigma_sum.sin().0;

        // if at North pole, only valid azimuth is due South
        if MAX_LAT < sin_beta {
            Angle::new(trig::UnitNegRange(0.0), trig::UnitNegRange(-1.0))
        } else {
            Angle::from_y_x(self.azi0.sin().0, self.azi0.cos().0 * sigma_sum.cos().0)
        }
    }

    /// Calculate the azimuth at the length along the geodesic.
    /// * `distance` - the distance along the `GeodesicSegment`, in metres.
    ///
    /// return the azimuth of the geodesic/great circle at length.
    #[must_use]
    pub fn azimuth(&self, distance: Metres) -> Angle {
        let sigma = Angle::from(self.metres_to_radians(distance));
        self.arc_azimuth(sigma)
    }

    /// Calculate the geodesic longitude difference at arc distance
    /// along the auxiliary sphere.
    /// * `arc_distance` - the great circle arc distance on the auxiliary sphere.
    /// * `sigma` - the arc distance as an Angle.
    ///
    /// return the longitude difference from the start point.
    #[must_use]
    pub fn delta_longitude(&self, arc_distance: Radians, sigma: Angle) -> Angle {
        if arc_distance.abs().0 < great_circle::MIN_VALUE {
            Angle::default()
        } else {
            // The great circle distance from Northward Equator crossing.
            let sigma_sum = self.sigma1 + sigma;

            // The longitude difference on the auxiliary sphere, omega12.
            let omega12 = Angle::from_y_x(self.azi0.sin().0 * sigma_sum.sin().0, sigma_sum.cos().0)
                - Angle::from_y_x(
                    self.azi0.sin().0 * self.beta.sin().0,
                    self.beta.cos().0 * self.azi.cos().0,
                );

            let c3 = self.ellipsoid.calculate_c3y(self.eps);
            let b31 = ellipsoid::coefficients::sin_cos_series(&c3, self.sigma1);
            let b32 = ellipsoid::coefficients::sin_cos_series(&c3, sigma_sum);

            omega12 - Angle::from(Radians(self.a3c * (arc_distance.0 + (b32.0 - b31.0))))
        }
    }

    /// Calculate the geodesic longitude at the great circle length along
    /// the auxiliary sphere.
    /// * `arc_distance` - the great circle arc distance on the auxiliary sphere.
    ///
    /// return the longitude of the geodesic at `arc_distance`.
    #[must_use]
    pub fn arc_longitude(&self, arc_distance: Radians) -> Angle {
        let sigma = Angle::from(arc_distance);
        self.lon + self.delta_longitude(arc_distance, sigma)
    }

    /// Calculate the geodesic longitude at the distance along the geodesic.
    /// * `distance` - the distance along the `GeodesicSegment`, in metres.
    ///
    /// return the longitude of the geodesic at distance.
    #[must_use]
    pub fn longitude(&self, distance: Metres) -> Angle {
        self.arc_longitude(self.metres_to_radians(distance))
    }

    /// Calculate the parametric latitude and longitude at the arc distance.
    ///
    /// * `arc_distance` - the arc distance on the auxiliary sphere in Radians.
    ///
    /// return the parametric latitude and longitude at `arc_distance`.
    #[must_use]
    pub fn arc_beta_long(&self, arc_distance: Radians) -> (Angle, Angle) {
        let sigma = Angle::from(arc_distance);
        (
            self.arc_beta(sigma),
            self.lon + self.delta_longitude(arc_distance, sigma),
        )
    }

    /// Calculate the geodesic `LatLong` at the arc distance along
    /// the auxiliary sphere.
    /// * `arc_distance` - the great circle arc distance on the auxiliary sphere.
    ///
    /// return the `LatLong` of the geodesic position at `arc_distance`.
    #[must_use]
    pub fn arc_lat_long(&self, arc_distance: Radians) -> LatLong {
        let (beta, lon) = self.arc_beta_long(arc_distance);
        LatLong::new(
            angle_sc::Degrees::from(self.ellipsoid.calculate_geodetic_latitude(beta)),
            angle_sc::Degrees::from(lon),
        )
    }

    /// Calculate the geodesic `LatLong` at the distance along the `GeodesicSegment`.
    /// * `distance` - the distance in `Metres`.
    ///
    /// return the `LatLong` of the geodesic position at `distance`.
    #[must_use]
    pub fn lat_long(&self, distance: Metres) -> LatLong {
        let arc_distance = self.metres_to_radians(distance);
        self.arc_lat_long(arc_distance)
    }

    /// Calculate the parametric latitude, longitude and azimuth at the arc distance.
    ///
    /// * `arc_distance` - the arc distance on the auxiliary sphere in Radians.
    ///
    /// return the parametric latitude, longitude and azimuth at `arc_distance`.
    #[must_use]
    pub fn arc_angles(&self, arc_distance: Radians) -> (Angle, Angle, Angle) {
        let sigma = Angle::from(arc_distance);
        let beta: Angle = self.arc_beta(sigma);
        let lon = self.lon + self.delta_longitude(arc_distance, sigma);
        let azimuth = self.arc_azimuth(sigma);

        (beta, lon, azimuth)
    }

    /// Calculate the vector on the auxiliary sphere at `arc_distance` in `Radians`.
    /// * `arc_distance` the great circle distance on the auxiliary sphere
    ///
    /// returns the point on the auxiliary sphere at `arc_distance`.
    #[must_use]
    pub fn arc_point(&self, arc_distance: Radians) -> Vector3d {
        if arc_distance.abs().0 < great_circle::MIN_VALUE {
            unit_sphere::vector::to_point(self.beta, self.lon)
        } else {
            let (beta, lon) = self.arc_beta_long(arc_distance);
            unit_sphere::vector::to_point(beta, lon)
        }
    }

    /// Calculate the vector on the auxiliary sphere at the mid point of the `GeodesicSegment`.
    ///
    /// returns the mid point vector of the `GeodesicSegment`.
    #[must_use]
    pub fn mid_point(&self) -> Vector3d {
        self.arc_point(self.metres_to_radians(self.length().half()))
    }

    /// Calculate the geodesic and pole at the great circle arc distance.
    /// * `arc_distance` the great circle arc distance on the auxiliary sphere
    ///
    /// returns the point and pole on the auxiliary sphere at `arc_distance`.
    #[must_use]
    pub fn arc_pole(&self, arc_distance: Radians) -> Vector3d {
        // if point is on a meridional GeodesicSegment use auxiliary sphere point and pole
        if self.azi0.sin().abs().0 < great_circle::MIN_VALUE {
            unit_sphere::vector::calculate_pole(self.beta, self.lon, self.azi)
        } else {
            let (beta, lon, azimuth) = self.arc_angles(arc_distance);
            unit_sphere::vector::calculate_pole(beta, lon, azimuth)
        }
    }

    /// Calculate the geodesic point and pole at the arc distance along the great circle.
    /// * `arc_distance` the great circle arc distance on the auxiliary sphere
    ///
    /// returns the point and pole on the auxiliary sphere at `arc_distance`.
    #[must_use]
    pub fn arc_point_and_pole(&self, arc_distance: Radians) -> (Vector3d, Vector3d) {
        let (beta, lon, azimuth) = self.arc_angles(arc_distance);

        // if point is on a meridional GeodesicSegment use auxiliary sphere point and pole
        let pole = if self.azi0.sin().abs().0 < great_circle::MIN_VALUE {
            unit_sphere::vector::calculate_pole(self.beta, self.lon, self.azi)
        } else {
            unit_sphere::vector::calculate_pole(beta, lon, azimuth)
        };

        (unit_sphere::vector::to_point(beta, lon), pole)
    }

    /// The reverse `GeodesicSegment` from end to start.
    ///
    /// returns the reverse `GeodesicSegment` from end to start.
    #[must_use]
    pub fn reverse(&self) -> GeodesicSegment<'_> {
        let sigma = Angle::from(self.arc_length);
        let mut segment = GeodesicSegment::new(
            self.arc_beta(sigma),
            self.lon + self.delta_longitude(self.arc_length, sigma),
            self.arc_azimuth(sigma).opposite(),
            self.arc_length,
            self.half_width,
            self.ellipsoid,
        );
        segment.set_half_width(self.half_width);
        segment.clone()
    }

    /// Calculate along and across track distances to a position from a geodesic segment.
    /// * `beta` the latitude of the position
    /// * `lon` the longitude of the position
    /// * `precision` the required precision
    ///
    /// returns the along and across track distances to the position in `Radians`.
    #[allow(clippy::similar_names)]
    #[must_use]
    pub fn calculate_sphere_atd_and_xtd(
        &self,
        beta: Angle,
        lon: Angle,
        precision: Radians,
    ) -> (Radians, Radians, u32) {
        const MAX_ITERATIONS: u32 = 10;

        // calculate the position as a point on the unit sphere
        let point = unit_sphere::vector::to_point(beta, lon);

        // calculate the start point and pole of the geodesic on the unit sphere
        let (a, pole) = self.arc_point_and_pole(Radians(0.0));
        let gc_d =
            unit_sphere::great_circle::e2gc_distance(unit_sphere::vector::distance(&a, &point));

        // if the point is close to the start point of the GeodesicSegment
        if gc_d < precision {
            (Radians(0.0), Radians(0.0), 0)
        } else {
            // estimate initial along track distance on the unit sphere
            let (mut atd, mut xtd) = unit_sphere::vector::calculate_atd_and_xtd(&a, &pole, &point);
            let mut iterations = 1;
            while iterations < MAX_ITERATIONS {
                // calculate the position and azimuth at atd along the GeodesicSegment
                let (beta_x, lon_x, azi_x) = self.arc_angles(atd);

                // calculate the geodesic azimuth and length to the point from the GeodesicSegment position at atd
                let (azi_p, length, _, _) = geodesic::aux_sphere_azimuths_length(
                    beta_x,
                    beta,
                    lon - lon_x,
                    Radians(great_circle::MIN_VALUE),
                    self.ellipsoid,
                );
                let delta_azi = azi_x - azi_p;
                let delta_atd = trig::spherical_cosine_rule(delta_azi.cos(), length);
                atd += delta_atd;
                xtd = length;

                if delta_atd.abs().0 < precision.0 {
                    break;
                }

                iterations += 1;
            }
            // get the cross track distance (and sign) at the along track distance
            xtd = if xtd < precision {
                Radians(0.0)
            } else {
                let pole = self.arc_pole(atd);
                let sign = pole.dot(&point);
                Radians(xtd.0.copysign(sign))
            };
            (atd, xtd, iterations)
        }
    }

    /// Calculate the shortest geodesic distance of point from the `GeodesicSegment`.
    ///
    /// * `beta` - the parametric latitude of the point.
    /// * `lon` - the longitude of the point.
    /// * `precision` the required precision in `Radians`
    ///
    /// returns the shortest distance of the point from the `GeodesicSegment` in Metres.
    #[allow(clippy::similar_names)]
    #[must_use]
    pub fn calculate_sphere_shortest_distance(
        &self,
        beta: Angle,
        lon: Angle,
        precision: Radians,
    ) -> Metres {
        let (atd, xtd, _) = self.calculate_sphere_atd_and_xtd(beta, lon, precision);

        // if the position is beside the geodesic segment
        if (-precision <= atd) && (atd <= self.arc_length + precision) {
            if xtd.abs() < precision {
                Metres(0.0)
            } else {
                // convert cross track distance to Metres
                let (beta_x, _lon, azi) = self.arc_angles(atd);
                let alpha = azi.quarter_turn_ccw();
                let distance =
                    geodesic::convert_radians_to_metres(beta_x, alpha, xtd, self.ellipsoid);
                // return the abs cross track distance in Metres
                Metres(distance.0.abs())
            }
        } else {
            // adjust atd to measure the distance from the centre of the Arc to the point
            let atd_centre = atd - self.arc_length.half();
            if atd_centre.0.is_sign_negative() {
                // calculate the geodesic distance from the start of the segment
                let delta_long = lon - self.lon;
                let (alpha, distance, _, _) = geodesic::aux_sphere_azimuths_length(
                    self.beta,
                    beta,
                    delta_long,
                    precision,
                    self.ellipsoid,
                );
                geodesic::convert_radians_to_metres(self.beta, alpha, distance, self.ellipsoid)
            } else {
                // calculate the geodesic distance from the end of the segment
                let (arc_beta, arc_lon, _azi) = self.arc_angles(self.arc_length());
                let delta_long = lon - arc_lon;
                let (alpha, distance, _, _) = geodesic::aux_sphere_azimuths_length(
                    arc_beta,
                    beta,
                    delta_long,
                    precision,
                    self.ellipsoid,
                );
                geodesic::convert_radians_to_metres(arc_beta, alpha, distance, self.ellipsoid)
            }
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
    /// use unit_sphere::great_circle;
    ///
    /// let tolerance = Radians(great_circle::MIN_VALUE);
    ///
    /// let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
    /// let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
    /// let g_0 = GeodesicSegment::between_positions(&istanbul, &washington, Metres(0.0), tolerance, &WGS84_ELLIPSOID);
    ///
    /// let azimuth_degrees = Degrees::from(g_0.azimuth(Metres(0.0)));
    /// println!("Istanbul-Washington initial azimuth: {:?}", azimuth_degrees.0);
    ///
    /// let distance_nm = NauticalMiles::from(g_0.length());
    /// println!("Istanbul-Washington distance: {:?}", distance_nm);
    ///
    /// let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
    ///
    /// // Calculate geodesic along track and across track distances to 1mm precision.
    /// let (atd, xtd, iterations) = g_0.calculate_atd_and_xtd(&reyjavik, Metres(1e-3));
    /// assert!(is_within_tolerance(3928788.572, atd.0, 1e-3));
    /// assert!(is_within_tolerance(-1010585.9988368, xtd.0, 1e-3));
    /// println!("calculate_atd_and_xtd iterations: {:?}", iterations);
    ///
    /// // The expected latitude and longitude are from:
    /// // <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#8a93>
    /// let position = g_0.lat_long(atd);
    /// assert!(is_within_tolerance(
    ///     54.92853149711691,
    ///     Degrees::from(position.lat()).0,
    ///     128.0 * f64::EPSILON
    /// ));
    /// assert!(is_within_tolerance(
    ///     -21.93729106604878,
    ///     Degrees::from(position.lon()).0,
    ///     2048.0 * f64::EPSILON
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

        // calculate the parametric latitude and longitude of the position
        let beta = self
            .ellipsoid
            .calculate_parametric_latitude(Angle::from(position.lat()));
        let lon = Angle::from(position.lon());

        let (atd, xtd, iterations) = self.calculate_sphere_atd_and_xtd(beta, lon, precision);

        // calculate the parametric latitude and azimuth at the abeam point
        let atd_angle = Angle::from(atd);
        let beta = self.arc_beta(atd_angle);
        let alpha = self.arc_azimuth(atd_angle).quarter_turn_ccw();
        (
            self.radians_to_metres(atd, atd_angle),
            geodesic::convert_radians_to_metres(beta, alpha, xtd, self.ellipsoid),
            iterations,
        )
    }

    /// Calculate the shortest geodesic distance of point from the `GeodesicSegment`.
    ///
    /// * `position` - the point.
    /// * `precision` the required precision in `Metres`
    ///
    /// returns the shortest distance of the point from the `GeodesicSegment` in Metres.
    #[allow(clippy::similar_names)]
    #[must_use]
    pub fn shortest_distance(&self, position: &LatLong, precision: Metres) -> Metres {
        // calculate the parametric latitude and longitude of the position
        let beta = self
            .ellipsoid
            .calculate_parametric_latitude(Angle::from(position.lat()));
        let lon = Angle::from(position.lon());
        // convert precision to Radians
        let precision = Radians(precision.0 / self.ellipsoid.a().0);
        self.calculate_sphere_shortest_distance(beta, lon, precision)
    }
}

impl From<(&LatLong, Angle, Radians)> for GeodesicSegment<'_> {
    /// Construct a `GeodesicSegment` on the WGS-84  `Ellipsoid` using the "direct"
    /// method with the length in `Radians`.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `arc_length` - the great circle arc length on the auxiliary sphere in radians.
    fn from(params: (&LatLong, Angle, Radians)) -> Self {
        GeodesicSegment::from_lat_lon_azi_arc_length(params.0, params.1, params.2, &WGS84_ELLIPSOID)
    }
}

impl From<(&LatLong, Angle, Metres)> for GeodesicSegment<'_> {
    /// Construct a `GeodesicSegment` on the WGS-84 `Ellipsoid` using the "direct"
    /// method with the length in metres.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a` - the start position in geodetic coordinates.
    /// * `azimuth` - the azimuth at the start position.
    /// * `length` - the length on the `Ellipsoid` in metres.
    fn from(params: (&LatLong, Angle, Metres)) -> Self {
        GeodesicSegment::from_lat_lon_azi_length(params.0, params.1, params.2, &WGS84_ELLIPSOID)
    }
}

impl From<(&LatLong, &LatLong)> for GeodesicSegment<'_> {
    /// Construct a `GeodesicSegment` between a pair of positions on the WGS-84
    /// `Ellipsoid`, the "indirect" method.
    /// @pre |lat| <= 90.0 degrees.
    /// * `a`, `b` - the start and finish positions in geodetic coordinates.
    fn from(params: (&LatLong, &LatLong)) -> Self {
        Self::between_positions(
            params.0,
            params.1,
            Metres(0.0),
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        )
    }
}

/// Calculate the distances along a pair of `GeodesicSegment`s (in Radians)
/// to their closest intersection or reference points.
/// * `g_0`, `g_1` the Geodesic`GeodesicSegment`s.
/// * `precision` the precision in `Metres`
///
/// returns the distances along the `GeodesicSegment`s to the intersection
/// point or to their closest (reference) points if the `GeodesicSegment`s
/// do not intersect and the number of iterations required.
///
/// # Panics
///
/// The function will panic if the Geodesics are **not** on the same `Ellipsoid`.
#[must_use]
pub fn calculate_intersection_distances(
    g_0: &GeodesicSegment,
    g_1: &GeodesicSegment,
    precision: Metres,
) -> (Radians, Radians) {
    let precision = Radians(precision.0 / g_0.ellipsoid().a().0);
    let (distance1, distance2, _, _) =
        intersection::calculate_arc_reference_distances_and_angle(g_0, g_1, precision);
    (
        distance1 + g_0.arc_length().half(),
        distance2 + g_1.arc_length().half(),
    )
}

/// Calculate the position (Latitude and Longitude) where a pair of `GeodesicSegment`s
/// intersect, or None if the `GeodesicSegment`s do not intersect.
/// * `g_0`, `g_1` the `GeodesicSegment`s.
/// * `precision` the precision in `Metres`
///
/// returns the distances along the `GeodesicSegment`s to the intersection point
/// or to their closest (reference) points if the `GeodesicSegment`s do not intersect.
///
/// # Panics
///
/// The function will panic if the `GeodesicSegment`s are **not** on the same `Ellipsoid`.
///
/// # Examples
/// ```
/// use icao_wgs84::*;
/// use angle_sc::is_within_tolerance;
///
/// let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
/// let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
/// let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
/// let accra = LatLong::new(Degrees(6.0), Degrees(0.0));
///
/// let g_0 = GeodesicSegment::from((&istanbul, &washington));
/// let g_1 = GeodesicSegment::from((&reyjavik, &accra));
///
/// // Calculate the intersection point position
/// // The expected latitude and longitude are from:
/// // <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a>
/// let result = calculate_intersection_point(&g_0, &g_1, Metres(1e-3));
/// let lat_lon = result.unwrap();
///
/// assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));
/// assert!(is_within_tolerance(-14.56385574430775, lat_lon.lon().0, 1e-6));
/// ```
#[must_use]
pub fn calculate_intersection_point(
    g_0: &GeodesicSegment,
    g_1: &GeodesicSegment,
    precision: Metres,
) -> Option<LatLong> {
    let precision = Radians(precision.0 / g_0.ellipsoid().a().0);
    let (distance1, distance2, angle, _) =
        intersection::calculate_arc_reference_distances_and_angle(g_0, g_1, precision);

    let segments_are_coincident = angle.sin().0 == 0.0;
    let segments_intersect_or_overlap = if segments_are_coincident {
        // do coincident segments overlap?
        distance1.abs() + distance2.abs()
            <= g_0.arc_length().half() + g_1.arc_length().half() + precision
    } else {
        // do geodesic paths intersect inside both segments
        (distance1.abs() <= g_0.arc_length().half() + precision)
            && distance2.abs() <= (g_1.arc_length().half() + precision)
    };

    if segments_intersect_or_overlap {
        let distance = (distance1 + g_0.arc_length().half()).clamp(g_0.arc_length());
        Some(g_0.arc_lat_long(distance))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use angle_sc::is_within_tolerance;
    use core::mem::size_of;
    use unit_sphere::{LatLong, great_circle};

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

        let point = geoid.to_arc_point(Angle::from(Degrees(45.0)), Angle::from(Degrees(45.0)));
        assert!(is_within_tolerance(
            44.903787849420226,
            Degrees::from(unit_sphere::vector::latitude(&point)).0,
            f64::EPSILON
        ));
        assert_eq!(
            45.0,
            Degrees::from(unit_sphere::vector::longitude(&point)).0
        );
    }

    #[test]
    fn test_ellipsoid_traits() {
        let geoid = Ellipsoid::wgs84();

        let geoid_clone = geoid.clone();
        assert!(geoid_clone == geoid);

        println!("Ellipsoid: {:?}", geoid);
    }

    #[test]
    fn test_calculate_azimuth_arc_length_normal_05() {
        // GeodTest.dat line 2874
        // 5.421025561218 0 84.846843174846
        // 3.027329237478900117 109.666857465735641205 96.826992198613537236
        // 12161089.9991805 109.607910081857488806 5988906.6319258056178 8449589948776.249238

        // North East bound, straddle Equator
        let latlon1 = LatLong::new(Degrees(5.421025561218), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(3.027329237478900117),
            Degrees(109.666857465735641205),
        );

        let result = calculate_azimuths_and_geodesic_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert_eq!(84.846843174846, Degrees::from(result.0).0);
        assert!(is_within_tolerance(12161089.9991805, (result.1).0, 1e-8));
        // assert_eq!(84.846843174846, Degrees::from(result.2).0);
    }

    #[test]
    fn test_calculate_azimuth_arc_length_nearly_antipodal_1() {
        // GeodTest.dat line 100001
        // 8.226828747671 0 111.1269645725
        // -8.516119211674268968 178.688979582629224039 68.982798544955243193
        // 19886305.6710041 179.197987814300505446 97496.4436255989712 -29736790544759.340534

        let latlon1 = LatLong::new(Degrees(8.226828747671), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-8.516119211674268968),
            Degrees(178.688979582629224039),
        );

        let result = calculate_azimuths_and_geodesic_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert!(is_within_tolerance(
            111.1269645725,
            Degrees::from(result.0).0,
            1e-9
        ));
        assert!(is_within_tolerance(19886305.6710041, (result.1).0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuth_arc_length_nearly_antipodal_2() {
        // GeodTest.dat line 100017
        // .322440123063 0 100.319048368176
        // -.367465171996537868 179.160624688175359763 79.682430612745621077
        // 19943611.6727803 179.749470297545372441 29954.0028615773743 -14555544282075.683105

        let latlon1 = LatLong::new(Degrees(0.322440123063), Degrees(0.0));
        let latlon2 = LatLong::new(
            Degrees(-0.367465171996537868),
            Degrees(179.160624688175359763),
        );

        let result = calculate_azimuths_and_geodesic_length(
            &latlon1,
            &latlon2,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        assert!(is_within_tolerance(
            100.319048368176,
            Degrees::from(result.0).0,
            1e-9
        ));
        assert!(is_within_tolerance(19943611.6727803, (result.1).0, 1e-8));
    }

    #[test]
    fn test_calculate_azimuth_and_geodesic_length_karney() {
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let (azimuth, length, end_azimuth) = calculate_azimuths_and_geodesic_length(
            &istanbul,
            &washington,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        let azimuth_degrees = Degrees::from(azimuth);
        assert_eq!(-50.69375304113997, azimuth_degrees.0);
        assert_eq!(8339863.136005359, length.0);

        println!("Istanbul-Washington azimuth: {:?}", azimuth_degrees.0);

        let distance_nm = NauticalMiles::from(length);
        println!("Istanbul-Washington distance: {:?}", distance_nm);

        let azimuth_degrees = Degrees::from(end_azimuth.opposite());
        assert_eq!(47.735339288362425, azimuth_degrees.0);
        println!("Washington-Istanbul azimuth: {:?}", azimuth_degrees.0);
    }

    #[test]
    fn test_geodesicarc_direct_constructors() {
        let length = Metres(9_000_000.0);
        let arc_length = Radians(core::f64::consts::FRAC_PI_2);

        // Ensure that two Geodesics can fit on a cache line.
        assert_eq!(128, size_of::<GeodesicSegment>());

        let a = LatLong::new(Degrees(45.0), Degrees(45.0));

        // Increase azimuth around compass from due South to due North
        for i in -180..180 {
            let azi = i as f64;
            let azimuth = Angle::from(Degrees(azi));

            let geodesic1 = GeodesicSegment::from((&a, azimuth, length));
            assert!(geodesic1.is_valid());
            assert_eq!(Metres(0.0), geodesic1.half_width());
            let azi0 = geodesic1.azimuth(Metres(0.0));
            assert!(is_within_tolerance(
                Radians::from(azimuth).0,
                Radians::from(azi0).0,
                2.0 * f64::EPSILON
            ));

            let len0 = geodesic1.length();
            assert!(is_within_tolerance(length.0, len0.0, 1.0e-8));

            let geodesic2 = GeodesicSegment::from((&a, azimuth, arc_length));
            assert!(geodesic2.is_valid());
            assert_eq!(Metres(0.0), geodesic2.half_width());
            let azi0 = geodesic2.azimuth(Metres(0.0));
            assert!(is_within_tolerance(
                Radians::from(azimuth).0,
                Radians::from(azi0).0,
                2.0 * f64::EPSILON
            ));

            let len0 = geodesic2.arc_length();
            assert!(is_within_tolerance(arc_length.0, len0.0, 1.0e-8));
        }
    }

    #[test]
    fn test_geodesicarc_between_positions() {
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));

        let tolerance = Radians(great_circle::MIN_VALUE);

        let g_0 = GeodesicSegment::between_positions(
            &istanbul,
            &washington,
            Metres(0.0),
            tolerance,
            &WGS84_ELLIPSOID,
        );
        assert!(g_0.is_valid());
        assert_eq!(Metres(0.0), g_0.half_width());

        let end_azimuth = Degrees::from(g_0.azimuth(g_0.length()));
        assert_eq!(-132.2646607116376, end_azimuth.0);

        let mut g1_clone = g_0.clone();
        assert_eq!(g1_clone, g_0);

        let g1_clone = g1_clone.set_arc_length(Radians(1.0));
        assert_eq!(Radians(1.0), g1_clone.arc_length());
        println!("GeodesicSegment: {:?}", &g1_clone);

        let g1_clone = g1_clone.set_half_width(Metres(2.0));
        assert_eq!(Metres(2.0), g1_clone.half_width());

        // test start position
        assert!(is_within_tolerance(
            42.0,
            Degrees::from(WGS84_ELLIPSOID.calculate_geodetic_latitude(g_0.beta())).0,
            32.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            29.0,
            Degrees::from(g_0.lon()).0,
            16.0 * f64::EPSILON
        ));

        let lat_long: LatLong = g_0.lat_long(Metres(0.0));
        // assert_eq!(istanbul, lat_long);
        assert!(is_within_tolerance(
            istanbul.lat().0,
            lat_long.lat().0,
            64.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            istanbul.lon().0,
            lat_long.lon().0,
            32.0 * f64::EPSILON
        ));

        // test start point
        let start_point = g_0.arc_point(Radians(0.0));
        assert!(is_within_tolerance(
            istanbul.lat().0,
            Degrees::from(
                g_0.ellipsoid()
                    .calculate_geodetic_latitude(unit_sphere::vector::latitude(&start_point))
            )
            .0,
            64.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            istanbul.lon().0,
            Degrees::from(unit_sphere::vector::longitude(&start_point)).0,
            16.0 * f64::EPSILON
        ));

        let precision = Metres(1e-3);
        let (atd, xtd, iterations) = g_0.calculate_atd_and_xtd(&istanbul, precision);
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);
        assert_eq!(0.0, atd.0);
        assert_eq!(0.0, xtd.0);

        let distance = g_0.shortest_distance(&istanbul, precision);
        assert_eq!(0.0, distance.0);

        // test end position
        let arc_length = g_0.arc_length();
        assert_eq!(1.309412846249522, arc_length.0);

        let length = g_0.length();
        assert_eq!(8339863.136005359, length.0);

        // let end_position = g_0.arc_lat_long(arc_length);
        assert!(is_within_tolerance(
            39.0,
            Degrees::from(g_0.latitude(length)).0,
            32.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            -77.0,
            Degrees::from(g_0.longitude(length)).0,
            256.0 * f64::EPSILON
        ));

        // test mid position
        let half_length = length.half();
        let mid_position = g_0.lat_long(half_length);

        assert!(is_within_tolerance(
            54.86379153725445,
            Degrees::from(mid_position.lat()).0,
            64.0 * f64::EPSILON
        ));

        assert!(is_within_tolerance(
            -25.694568908316413,
            Degrees::from(mid_position.lon()).0,
            128.0 * f64::EPSILON
        ));

        let mid_length = g_0.metres_to_radians(half_length);
        assert_eq!(0.654673165141749, mid_length.0);
        let mid_point = g_0.mid_point();
        let mid_beta = unit_sphere::vector::latitude(&mid_point);
        let mid_lat = g_0.ellipsoid().calculate_geodetic_latitude(mid_beta);
        assert!(is_within_tolerance(
            54.86379153725445,
            Degrees::from(mid_lat).0,
            64.0 * f64::EPSILON
        ));

        let mid_lon = unit_sphere::vector::longitude(&mid_point);
        assert!(is_within_tolerance(
            -25.694568908316413,
            Degrees::from(mid_lon).0,
            128.0 * f64::EPSILON
        ));

        let precision_r = Radians(1e-3 / WGS84_ELLIPSOID.a().0);
        let (atd, xtd, iterations) =
            g_0.calculate_sphere_atd_and_xtd(mid_beta, mid_lon, precision_r);
        assert!(is_within_tolerance(mid_length.0, atd.0, f64::EPSILON));
        assert_eq!(0.0, xtd.0);
        println!("calculate_sphere_atd_and_xtd iterations: {:?}", iterations);

        let (atd, xtd, iterations) = g_0.calculate_atd_and_xtd(&mid_position, precision);
        assert!(is_within_tolerance(half_length.0, atd.0, 1e-3));
        assert!(xtd.0.abs() < 1e-3);
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);

        let distance = g_0.shortest_distance(&mid_position, precision);
        assert_eq!(0.0, distance.0);

        let g_1 = g_0.reverse();
        assert_eq!(g_0.arc_length(), g_1.arc_length());
        assert_eq!(g_0.length(), g_1.length());
        assert_eq!(
            washington.lat().0,
            Degrees::from(g_1.arc_latitude(Radians::default())).0
        );
        assert_eq!(
            washington.lon().0,
            Degrees::from(g_1.arc_longitude(Radians(0.0))).0
        );
        assert!(is_within_tolerance(
            istanbul.lat().0,
            Degrees::from(g_1.arc_latitude(g_1.arc_length())).0,
            64.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            istanbul.lon().0,
            Degrees::from(g_1.arc_longitude(g_1.arc_length())).0,
            64.0 * f64::EPSILON
        ));
    }

    #[test]
    fn test_meridonal_geodesicarc() {
        // A GeodesicSegment along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(45.0), Degrees(180.0));
        let g_0 = GeodesicSegment::from((&a, &b));
        let pole0 = g_0.arc_pole(Radians(0.0));

        // Calculate the azimuth at the North pole
        let mid_length = g_0.arc_length().half();
        let azimuth = Degrees::from(g_0.arc_azimuth(Angle::from(mid_length)));
        assert_eq!(180.0, azimuth.0);

        // Calculate the point and great circle pole at the North pole
        let (point1, pole1) = g_0.arc_point_and_pole(mid_length);
        assert_eq!(Vector3d::new(0.5 * f64::EPSILON, 0.0, 1.0), point1);
        assert_eq!(pole0, pole1);
    }

    #[test]
    fn test_geodesicarc_90n_0n_0e() {
        let a = LatLong::new(Degrees(90.0), Degrees(0.0));
        let b = LatLong::new(Degrees(0.0), Degrees(0.0));
        let g_0 = GeodesicSegment::from((&a, &b));

        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            g_0.arc_length().0,
            f64::EPSILON
        ));
        assert_eq!(180.0, Degrees::from(g_0.azi()).0);
    }

    #[test]
    fn test_geodesicarc_90s_0n_50e() {
        let a = LatLong::new(Degrees(-90.0), Degrees(0.0));
        let b = LatLong::new(Degrees(0.0), Degrees(50.0));
        let g_0 = GeodesicSegment::from((&a, &b));

        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            g_0.arc_length().0,
            f64::EPSILON
        ));
        assert_eq!(0.0, Degrees::from(g_0.azi()).0);
    }

    #[test]
    fn test_calculate_atd_and_xtd() {
        // Karney's example
        // Istanbul, Washington and Reyjavik
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let g_0 = GeodesicSegment::from((&istanbul, &washington));

        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));

        // Calculate geodesic along track and across track distances to 1mm precision.
        let precision = Metres(1e-3);

        let (atd, xtd, iterations) = g_0.calculate_atd_and_xtd(&reyjavik, precision);
        assert!(is_within_tolerance(3928788.572, atd.0, precision.0));
        assert!(is_within_tolerance(-1010585.9988368, xtd.0, precision.0));
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);

        // Karney's latitude and longitude from Final result at:
        // https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#8a93
        let position = g_0.lat_long(atd);
        assert!(is_within_tolerance(
            54.92853149711691,
            Degrees::from(position.lat()).0,
            128.0 * f64::EPSILON
        ));
        assert!(is_within_tolerance(
            -21.93729106604878,
            Degrees::from(position.lon()).0,
            2048.0 * f64::EPSILON
        ));

        // Test delta_azimuth at interception, should be PI/2
        let azimuth_1 = g_0.azimuth(atd);
        let g_1 = GeodesicSegment::from((&position, &reyjavik));
        let azimuth_2 = g_1.azi();
        let delta_azimuth = azimuth_2 - azimuth_1;
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            Radians::from(delta_azimuth).0,
            128.0 * f64::EPSILON
        ));

        // opposite geodesic
        let g_0 = GeodesicSegment::from((&washington, &istanbul));
        let (atd, xtd, iterations) = g_0.calculate_atd_and_xtd(&reyjavik, precision);
        assert!(is_within_tolerance(
            g_0.length().0 - 3928788.572,
            atd.0,
            precision.0
        ));
        assert!(is_within_tolerance(1010585.9988368, xtd.0, 1e-3));
        println!("calculate_atd_and_xtd iterations: {:?}", iterations);

        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
        let distance = g_0.shortest_distance(&reyjavik, precision);
        assert!(is_within_tolerance(
            1010585.998836817,
            distance.0,
            precision.0
        ));

        let accra = LatLong::new(Degrees(6.0), Degrees(0.0));
        let distance = g_0.shortest_distance(&accra, precision);
        assert!(is_within_tolerance(
            4891211.398445355,
            distance.0,
            precision.0
        ));

        let chicago = LatLong::new(Degrees(42.0), Degrees(-88.0));
        let distance = g_0.shortest_distance(&chicago, precision);
        assert!(is_within_tolerance(
            989277.1859906457,
            distance.0,
            precision.0
        ));

        let singapore = LatLong::new(Degrees(1.0), Degrees(104.0));
        let distance = g_0.shortest_distance(&singapore, precision);
        assert!(is_within_tolerance(
            8699538.22763653,
            distance.0,
            precision.0
        ));
    }

    #[test]
    fn test_intersection_point_distance() {
        // Karney's example
        // Istanbul, Washington, Reyjavik and Accra
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
        let accra = LatLong::new(Degrees(6.0), Degrees(0.0));

        let g_0 = GeodesicSegment::from((&istanbul, &washington));
        let g_1 = GeodesicSegment::from((&reyjavik, &accra));

        let (d1, d2) = calculate_intersection_distances(&g_0, &g_1, Metres(1e-3));
        assert!(is_within_tolerance(0.5423772210673643, d1.0, 1e-6));
        assert!(is_within_tolerance(0.17504915893226164, d2.0, 1e-6));

        let result = calculate_intersection_point(&g_0, &g_1, Metres(1e-3));
        let lat_lon = result.unwrap();

        assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));
        assert!(is_within_tolerance(
            -14.56385574430775,
            lat_lon.lon().0,
            1e-6
        ));

        // Swap geodesics
        let result = calculate_intersection_point(&g_1, &g_0, Metres(1e-3));
        let lat_lon = result.unwrap();

        assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));
        assert!(is_within_tolerance(
            -14.56385574430775,
            lat_lon.lon().0,
            1e-6
        ));
    }

    #[test]
    fn test_intersection_point_non_intersecting() {
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
        let accra = LatLong::new(Degrees(6.0), Degrees(0.0));

        let g_0 = GeodesicSegment::from((&istanbul, &accra));
        let g_1 = GeodesicSegment::from((&reyjavik, &washington));

        let result = calculate_intersection_point(&g_0, &g_1, Metres(1e-3));
        assert!(result.is_none());
    }

    #[test]
    fn test_intersection_coincident_geodesic_paths() {
        let south_pole_1 = LatLong::new(Degrees(-88.0), Degrees(-180.0));
        let south_pole_2 = LatLong::new(Degrees(-87.0), Degrees(0.0));

        let g_0 = GeodesicSegment::from((&south_pole_1, &south_pole_2));

        let intersection_lengths = calculate_intersection_distances(&g_0, &g_0, Metres(1e-3));
        assert_eq!(g_0.arc_length().half(), intersection_lengths.0);
        assert_eq!(g_0.arc_length().half(), intersection_lengths.1);

        let intersection_point = calculate_intersection_point(&g_0, &g_0, Metres(1e-3));
        let lat_lon = intersection_point.unwrap();
        assert!(is_within_tolerance(-89.5, lat_lon.lat().0, 1e-5));
        assert_eq!(0.0, lat_lon.lon().0);

        let south_pole_3 = LatLong::new(Degrees(-85.0), Degrees(0.0));
        let south_pole_4 = LatLong::new(Degrees(-86.0), Degrees(0.0));
        let g_1 = GeodesicSegment::from((&south_pole_3, &south_pole_4));
        let intersection_point = calculate_intersection_point(&g_0, &g_1, Metres(1e-3));
        assert!(intersection_point.is_none());
    }
}
