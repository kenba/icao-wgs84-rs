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

//! The `intersection` module contains functions for calculating geodesic
//! intersections using vectors.
//!
//! Charles Karney's developed a solution to the
//! [Intersection between two geodesic lines](https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a).
//! using a gnomonic projection.
//!
//! [Baselga and Martinez-Llario(2017)](https://www.researchgate.net/publication/321358300_Intersection_and_point-to-line_solutions_for_geodesics_on_the_ellipsoid)
//! solve geodesic intersection and point-to-line problems by using the
//! correspondence between geodesics on an ellipsoid and great-circles on a
//! unit sphere.
//!
//! Nick Korbey and I [Barker and Korbey(2019)](https://www.researchgate.net/publication/335749834_Geodesic_Geometry)
//! developed Baselga and Martinez-Llario's algorithms by using vectors to
//! solve geodesic intersection and point-to-line problems on the unit
//! sphere.
//!
//! [Karney(2023)](https://arxiv.org/abs/2308.00495) also developed Baselga and
//! Martinez-Llario's algorithms and provides reliable algorithms for solving
//! geodesic intersections on both oblate and prolate ellipsoids.
//!
//! This `intersection` module calculates initial great-circle intersection
//! distances on the unit sphere and then refines the intersection
//! distances using great-circles on the unit sphere at the calculated
//! intersection distances.

use crate::{Angle, GeodesicSegment, Radians, geodesic};
use angle_sc::max;
use unit_sphere::{great_circle, vector};

/// Determine whether two `GeodesicSegment`s are coincident,
/// i.e. they lie on the same geodesic path.
///
/// * `g_0`, `g_1` the `GeodesicSegment`s.
/// * `sin_max_coincident_angle` the sine of the maximum angle between coincident `GeodesicSegment`s.
///
/// returns true if the segments are on the same geodesic path, false otherwise.
#[must_use]
pub fn geodesics_are_coincident(
    g_0: &GeodesicSegment,
    g_1: &GeodesicSegment,
    sin_max_coincident_angle: f64,
) -> bool {
    // calculate the geodesic path between start positions
    let (g_2_azi, _g_2_arc_length, g_2_end_azi, _) = geodesic::aux_sphere_azimuths_length(
        g_0.beta(),
        g_1.beta(),
        g_1.lon() - g_0.lon(),
        Radians(great_circle::MIN_VALUE),
        g_0.ellipsoid(),
    );

    // the segments are coincident if their angle differences to the
    // path between start positions are both within sin_max_coincident_angle
    let delta_azimuth_0_2 = g_2_azi - g_0.azi();
    let delta_azimuth_1_2 = g_2_end_azi - g_1.azi();
    delta_azimuth_0_2.sin().abs().0 < sin_max_coincident_angle
        && delta_azimuth_1_2.sin().abs().0 < sin_max_coincident_angle
}

/// Find the closest intersection distances of two `GeodesicSegment`s.
///
/// * `g_0`, `g_1` the `GeodesicSegment`s.
/// * `use_antipodal_intersection` use the antipodal intersection point
/// * `distance_0`, `distance_1` the initial arc intersection distances from the segment start.
/// * `precision` the precision in `Radians`
/// * `sq_sin_max_coincident_angle` the square of the sine of the maximum angle between
///   coincident `GeodesicSegment`s.
///
/// returns the arc distances from the segment start to the closest intersection point,
/// the relative angle at the intersection point, and the number of iterations required.
#[must_use]
fn find_geodesic_intersection_distances(
    g_0: &GeodesicSegment,
    g_1: &GeodesicSegment,
    use_antipodal_intersection: bool,
    distance_0: Radians,
    distance_1: Radians,
    precision: Radians,
    sq_sin_max_coincident_angle: f64,
) -> (Radians, Radians, u32) {
    const MAX_ITERATIONS: u32 = 10;

    let sq_precision = great_circle::gc2e_distance(precision).powi(2);

    let mut distance_0 = distance_0;
    let mut distance_1 = distance_1;
    let mut iterations = 1;
    while iterations < MAX_ITERATIONS {
        let (point_0, pole_0) = g_0.arc_point_and_pole(distance_0);
        let (point_1, pole_1) = g_1.arc_point_and_pole(distance_1);

        let sq_d = vector::sq_distance(&point_0, &point_1);
        if sq_d < sq_precision {
            break;
        }

        iterations += 1;

        let x = if use_antipodal_intersection {
            vector::intersection::calculate_intersection(
                &pole_1,
                &pole_0,
                sq_sin_max_coincident_angle,
            )
        } else {
            vector::intersection::calculate_intersection(
                &pole_0,
                &pole_1,
                sq_sin_max_coincident_angle,
            )
        };
        match x {
            Some(x) => {
                distance_0 += vector::calculate_great_circle_atd(&point_0, &pole_0, &x);
                distance_1 += vector::calculate_great_circle_atd(&point_1, &pole_1, &x);
            }
            None => break,
        }
    }

    (distance_0, distance_1, iterations)
}

/// Find the closest intersection distances of two `GeodesicSegment`s and the
/// relative angle at the reference point.
///
/// The reference point is the closest intersection point or the centroid
/// if the `GeodesicSegment`s are coincident.
///
/// * `g_0`, `g_1` the `GeodesicSegment`s.
/// * `precision` the precision in `Radians`.
/// * `sin_max_coincident_angle` the sine of the maximum angle between coincident `GeodesicSegment`s.
///
/// returns the arc distances from the segment mid points to the reference point,
/// the relative angle at the reference point, and the number of iterations required.
///
/// # Panics
///
/// The function will panic if the `GeodesicSegment`s are **not** on the same `Ellipsoid`.
#[must_use]
pub fn calculate_arc_reference_distances_and_angle(
    g_0: &GeodesicSegment,
    g_1: &GeodesicSegment,
    precision: Radians,
    sin_max_coincident_angle: f64,
) -> (Radians, Radians, Angle, u32) {
    // The Geodesics MUST be on the same `Ellipsoid`
    assert!(g_0.ellipsoid() == g_1.ellipsoid());

    if g_0 == g_1 {
        return (Radians(0.0), Radians(0.0), Angle::default(), 0);
    }

    let sin_max_angle = max(sin_max_coincident_angle, vector::MIN_SIN_ANGLE);
    let sq_sin_max_coincident_angle = sin_max_angle * sin_max_angle;

    let half_length_0 = g_0.arc_length().half();
    let half_length_1 = g_1.arc_length().half();
    let (mid_point_0, pole_0) = g_0.arc_point_and_pole(half_length_0);
    let (mid_point_1, pole_1) = g_1.arc_point_and_pole(half_length_1);
    let centroid = mid_point_0 + mid_point_1;
    let intersection =
        vector::intersection::calculate_intersection(&pole_0, &pole_1, sq_sin_max_coincident_angle);
    if let Some(c) = intersection
        && !geodesics_are_coincident(g_0, g_1, sin_max_angle)
    {
        // great circles or geodesics interact

        // find the closest intersection
        let use_antipodal_intersection = vector::intersection::use_antipodal_point(&c, &centroid);
        let x = if use_antipodal_intersection { -c } else { c };

        // calculate distances to the closest geodesic intersection from arc start points
        let (distance_0, distance_1, iterations) = find_geodesic_intersection_distances(
            g_0,
            g_1,
            use_antipodal_intersection,
            vector::calculate_great_circle_atd(&g_0.a(), &pole_0, &x),
            vector::calculate_great_circle_atd(&g_1.a(), &pole_1, &x),
            precision,
            sq_sin_max_coincident_angle,
        );

        let angle = g_1.arc_azimuth(Angle::from(distance_1))
            - g_0.arc_azimuth(Angle::from(distance_0)).abs();

        // return distances to the closest geodesic intersection from arc mid points
        (
            distance_0 - half_length_0,
            distance_1 - half_length_1,
            angle.abs(),
            iterations,
        )
    } else {
        // great circles or geodesics are coincident

        // calculate reference distances to the normalised centroid from the mid points
        let c = vector::normalise_centroid(&centroid, &mid_point_0, &pole_0);
        let distance_0 = vector::calculate_great_circle_atd(&mid_point_0, &pole_0, &c);
        let distance_1 = vector::calculate_great_circle_atd(&mid_point_1, &pole_1, &c);

        let angle = g_1.arc_azimuth(Angle::from(distance_1))
            - g_0.arc_azimuth(Angle::from(distance_0)).abs();

        // return distances to the normalized centroid from arc mid points
        (distance_0, distance_1, angle, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeodesicSegment, Metres};
    use angle_sc::{Angle, Degrees, Radians, is_within_tolerance};
    use unit_sphere::LatLong;

    #[test]
    fn test_calculate_arc_reference_distances_and_angle_coincident_great_circles() {
        let latlong_w1 = LatLong::new(Degrees(0.0), Degrees(-1.0));
        let latlong_e1 = LatLong::new(Degrees(0.0), Degrees(1.0));
        let g_0 = GeodesicSegment::from((&latlong_w1, &latlong_e1));

        // 1m precision in Radians on the unit sphere
        let precision = Radians(Metres(1.0).0 / g_0.ellipsoid().a().0);
        let sin_precision = precision.0.sin();

        // same segments
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_0, precision, sin_precision);
        assert_eq!(Radians(0.0), result.0);
        assert_eq!(Radians(0.0), result.1);
        assert_eq!(Degrees(0.0), Degrees::from(result.2));

        // opposite segments and same geodesic paths
        let latlong_w179 = LatLong::new(Degrees(0.0), Degrees(-179.0));
        let latlong_e179 = LatLong::new(Degrees(0.0), Degrees(179.0));
        let g_1 = GeodesicSegment::from((&latlong_e179, &latlong_w179));
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_1, precision, sin_precision);
        assert!(is_within_tolerance(
            -core::f64::consts::FRAC_PI_2,
            result.0.0,
            precision.0
        ));
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            result.1.0,
            precision.0
        ));
        assert_eq!(Degrees(0.0), Degrees::from(result.2));

        // opposite segments and geodesic paths
        let g_2 = GeodesicSegment::from((&latlong_w179, &latlong_e179));
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_2, precision, sin_precision);
        assert!(is_within_tolerance(
            -core::f64::consts::FRAC_PI_2,
            result.0.0,
            precision.0
        ));
        assert!(is_within_tolerance(
            -core::f64::consts::FRAC_PI_2,
            result.1.0,
            precision.0
        ));
        assert_eq!(Degrees(180.0), Degrees::from(result.2));
    }

    #[test]
    fn test_calculate_arc_reference_distances_and_angle_intersecting_great_circles() {
        let latlong_w1 = LatLong::new(Degrees(0.0), Degrees(-1.0));
        let latlong_e1 = LatLong::new(Degrees(0.0), Degrees(1.0));
        let g_0 = GeodesicSegment::from((&latlong_w1, &latlong_e1));

        let latlong_s1 = LatLong::new(Degrees(-1.0), Degrees(0.0));
        let latlong_n1 = LatLong::new(Degrees(1.0), Degrees(0.0));
        let g_1 = GeodesicSegment::from((&latlong_s1, &latlong_n1));

        // 1m precision in Radians on the unit sphere
        let precision = Radians(Metres(1.0).0 / g_0.ellipsoid().a().0);
        let sin_precision = precision.0.sin();

        // intersection, same mid points, acute angle
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_1, precision, sin_precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(0.0, result.1.0, precision.0));
        assert_eq!(Degrees(90.0), Degrees::from(result.2));

        let angle = Degrees(44.81195977064123);

        let latlong_sw1 = LatLong::new(Degrees(-1.0), Degrees(-1.0));
        let latlong_ne1 = LatLong::new(Degrees(1.0), Degrees(1.0));
        let g_2 = GeodesicSegment::from((&latlong_sw1, &latlong_ne1));
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_2, precision, sin_precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(0.0, result.1.0, precision.0));
        assert!(is_within_tolerance(
            angle.0,
            Degrees::from(result.2).0,
            precision.0
        ));

        // intersection, same mid points, obtuse angle
        let g_3 = GeodesicSegment::from((&latlong_ne1, &latlong_sw1));
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_3, precision, sin_precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(0.0, result.1.0, precision.0));
        assert!(is_within_tolerance(
            180.0 - angle.0,
            Degrees::from(result.2).0,
            precision.0
        ));

        // intersection, different mid points, acute angle
        let g_4 = GeodesicSegment::new(
            Angle::default(),
            Angle::default(),
            g_2.arc_azimuth(Angle::from(g_2.arc_length().half())),
            Radians(core::f64::consts::FRAC_PI_2),
            Metres(0.0),
            g_2.ellipsoid(),
        );
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_4, precision, sin_precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(
            -core::f64::consts::FRAC_PI_4,
            result.1.0,
            precision.0
        ));
        assert!(is_within_tolerance(
            angle.0,
            Degrees::from(result.2).0,
            precision.0
        ));

        // intersection, different mid points, obtuse angle
        let g_5 = g_4.reverse();
        let result =
            calculate_arc_reference_distances_and_angle(&g_0, &g_5, precision, sin_precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_4,
            result.1.0,
            precision.0
        ));
        assert!(is_within_tolerance(
            180.0 - angle.0,
            Degrees::from(result.2).0,
            precision.0
        ));
    }
}
