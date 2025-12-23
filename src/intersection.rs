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
use unit_sphere::{great_circle, vector};

/// Calculate the distances along a pair of Geodesics (in Radians) to their
/// closest intersection point given the initial great-circle distances.
/// * `g1`, `g2` the Geodesics.
/// * `sq_precision` the square of the Euclidean precision.
/// * `use_antipodal_intersection` use the antipodal intersection point
/// * `initial_distances` the initial intersection distances in `Radians`.
///
/// returns the distances along a the Geodesics (in Radians), the relative angle
/// at the intersection point and the number of iterations required.
#[must_use]
fn calculate_geodesic_intersection_distances(
    g1: &GeodesicSegment,
    g2: &GeodesicSegment,
    sq_precision: f64,
    use_antipodal_intersection: bool,
    initial_distances: (Radians, Radians),
) -> (Radians, Radians, Angle, u32) {
    const MAX_ITERATIONS: u32 = 10;

    let (mut distance1, mut distance2) = initial_distances;
    let mut iterations = 1;
    while iterations < MAX_ITERATIONS {
        let (pos1, pole1) = g1.arc_point_and_pole(distance1);
        let (pos2, pole2) = g2.arc_point_and_pole(distance2);

        let sq_d = vector::sq_distance(&pos1, &pos2);
        if sq_d < sq_precision {
            break;
        }

        iterations += 1;

        let c = if use_antipodal_intersection {
            vector::intersection::calculate_intersection(&pole2, &pole1, vector::MIN_SQ_NORM)
        } else {
            vector::intersection::calculate_intersection(&pole1, &pole2, vector::MIN_SQ_NORM)
        };
        match c {
            Some(c) => {
                let (delta1, delta2) = vector::intersection::calculate_intersection_distances(
                    &pos1, &pole1, &pos2, &pole2, &c,
                );
                distance1 += delta1;
                distance2 += delta2;
            }
            None => break,
        }
    }

    let angle = g2.arc_azimuth(Angle::from(distance2)) - g1.arc_azimuth(Angle::from(distance1));

    (distance1, distance2, angle, iterations)
}

/// Calculate the distances along a pair of Geodesics (in Radians) to their
/// closest intersection or reference points.
/// * `g1`, `g2` the Geodesics.
/// * `precision` the precision in `Radians`
///
/// returns the distances along the Geodesics to the intersection point or to
/// their closest (reference) points if the Geodesics do not intersect and
/// the number of iterations required.
///
/// # Panics
///
/// The function will panic if the Geodesics are **not** on the same `Ellipsoid`.
#[must_use]
pub fn calculate_sphere_intersection_distances(
    g1: &GeodesicSegment,
    g2: &GeodesicSegment,
    precision: Radians,
) -> (Radians, Radians, Angle, u32) {
    // The Geodesics MUST be on the same `Ellipsoid`
    assert!(g1.ellipsoid() == g2.ellipsoid());

    // Convert precision in `Radians` to the square of Euclidean precision.
    let e_precision = great_circle::gc2e_distance(precision);
    let sq_precision = e_precision * e_precision;

    // Determine whether the geodesics are reciprocal
    let delta_azimuth1_2 = g2.azi() - g1.azi();
    let reciprocal = delta_azimuth1_2.cos().0.is_sign_negative();

    // if the start points are within precision of each other
    let sq_d = vector::sq_distance(&g1.a(), &g2.a());
    if sq_d < sq_precision {
        (Radians(0.0), Radians(0.0), delta_azimuth1_2, 0)
    } else {
        // Calculate geodesic path between start positions
        let (g3_azi, g3_arc_length, g3_end_azi, _) = geodesic::aux_sphere_azimuths_length(
            g1.beta(),
            g2.beta(),
            g2.lon() - g1.lon(),
            Radians(great_circle::MIN_VALUE),
            g1.ellipsoid(),
        );

        // Determine whether the geodesics are coincident
        let delta_azimuth1_3 = g3_azi - g1.azi();
        let delta_azimuth2_3 = g3_end_azi - g2.azi();
        if delta_azimuth1_3.sin().abs().0 < vector::MIN_SIN_ANGLE
            && delta_azimuth2_3.sin().abs().0 < vector::MIN_SIN_ANGLE
        {
            // The geodesics are coincident
            let (a1, pole1) = g1.arc_point_and_pole(Radians(0.0));
            let atd = Radians(
                g3_arc_length
                    .0
                    .copysign(unit_sphere::vector::sin_atd(&a1, &pole1, &g2.a()).0),
            );
            let distances = vector::intersection::calculate_coincident_arc_distances(
                atd,
                reciprocal,
                g1.arc_length(),
                g2.arc_length(),
            );
            let angle = if reciprocal {
                Angle::default().opposite()
            } else {
                Angle::default()
            };
            (distances.0, distances.1, angle, 0)
        } else {
            // Calculate the intersection of the poles at the mid points of the unit
            // sphere great circle arcs
            let half_arc_length1 = g1.arc_length().half();
            let half_arc_length2 = g2.arc_length().half();
            let (a1mid, pole1mid) = g1.arc_point_and_pole(half_arc_length1);
            let (a2mid, pole2mid) = g2.arc_point_and_pole(half_arc_length2);

            // Determine whether the great circles on the unit sphere are coincident
            vector::intersection::calculate_intersection(&pole1mid, &pole2mid, vector::MIN_SQ_NORM)
                .map_or_else(
                    || {
                        // This code should never be executed.
                        let (a1, pole1) = g1.arc_point_and_pole(Radians(0.0));
                        let atd = Radians(
                            g3_arc_length
                                .0
                                .copysign(unit_sphere::vector::sin_atd(&a1, &pole1, &g2.a()).0),
                        );
                        let distances = vector::intersection::calculate_coincident_arc_distances(
                            atd,
                            reciprocal,
                            g1.arc_length(),
                            g2.arc_length(),
                        );
                        let angle = if reciprocal {
                            Angle::default().opposite()
                        } else {
                            Angle::default()
                        };
                        (distances.0, distances.1, angle, 0)
                    },
                    |c| {
                        // find the closest intersection
                        let centroid = 0.5 * (a1mid + a2mid);
                        let use_antipodal_intersection =
                            vector::intersection::use_antipodal_point(&c, &centroid);
                        let c = if use_antipodal_intersection { -c } else { c };
                        // calculate distances along the arcs to the closest intersection
                        let (mut distance1, mut distance2) =
                            vector::intersection::calculate_intersection_distances(
                                &a1mid, &pole1mid, &a2mid, &pole2mid, &c,
                            );
                        // calculate distances from arc start
                        distance1 += half_arc_length1;
                        distance2 += half_arc_length2;
                        calculate_geodesic_intersection_distances(
                            g1,
                            g2,
                            sq_precision,
                            use_antipodal_intersection,
                            (distance1, distance2),
                        )
                    },
                )
        }
    }
}

/// Determine whether two `GeodesicSegment`s are coincident,
/// i.e. they lie on the same geodesic path.
///
/// * `g_0`, `g_1` the `GeodesicSegment`s.
/// * `min_sin_angle` the sine of the minimum angle between coincident `GeodesicSegment`s.
///
/// returns true if the segments are on the same geodesic path, false otherwise.
#[must_use]
pub fn geodesics_are_coincident(
    g_0: &GeodesicSegment,
    g_1: &GeodesicSegment,
    min_sin_angle: f64,
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
    // path between start positions are both within min_sin_angle
    let delta_azimuth_0_2 = g_2_azi - g_0.azi();
    let delta_azimuth_1_2 = g_2_end_azi - g_1.azi();
    delta_azimuth_0_2.sin().abs().0 < min_sin_angle
        && delta_azimuth_1_2.sin().abs().0 < min_sin_angle
}

/// Find the closest intersection distances of two `GeodesicSegment`s.
///
/// * `g_0`, `g_1` the `GeodesicSegment`s.
/// * `use_antipodal_intersection` use the antipodal intersection point
/// * `distance_0`, `distance_1` the initial arc intersection distances from the segment start.
/// * `precision` the precision in `Radians`
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
            vector::intersection::calculate_intersection(&pole_1, &pole_0, vector::MIN_SQ_NORM)
        } else {
            vector::intersection::calculate_intersection(&pole_0, &pole_1, vector::MIN_SQ_NORM)
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
) -> (Radians, Radians, Angle, u32) {
    // The Geodesics MUST be on the same `Ellipsoid`
    assert!(g_0.ellipsoid() == g_1.ellipsoid());

    if g_0 == g_1 {
        return (Radians(0.0), Radians(0.0), Angle::default(), 0);
    }

    let half_length_0 = g_0.arc_length().half();
    let half_length_1 = g_1.arc_length().half();
    let (mid_point_0, pole_0) = g_0.arc_point_and_pole(half_length_0);
    let (mid_point_1, pole_1) = g_1.arc_point_and_pole(half_length_1);
    let centroid = mid_point_0 + mid_point_1;
    let intersection =
        vector::intersection::calculate_intersection(&pole_0, &pole_1, vector::MIN_SQ_NORM);
    if let Some(c) = intersection
        && !geodesics_are_coincident(g_0, g_1, vector::MIN_SIN_ANGLE)
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
        );

        let angle =
            g_1.arc_azimuth(Angle::from(distance_1)) - g_0.arc_azimuth(Angle::from(distance_0));

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

        // determine whether the geodesics oppose each other
        let angle = if pole_0.dot(&pole_1).is_sign_negative() {
            Angle::default().opposite()
        } else {
            Angle::default()
        };

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
    fn test_intersection_same_start_point() {
        // A GeodesicSegment along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(45.0), Degrees(180.0));
        let g1 = GeodesicSegment::from((&a, &b));

        // Convert 1mm precision to Radians
        let precision = Radians(Metres(1e-3).0 / g1.ellipsoid().a().0);

        let c = LatLong::new(Degrees(45.0), Degrees(45.0));
        let g2 = GeodesicSegment::from((&a, &c));

        let (distance1, distance2, angle, iterations) =
            calculate_sphere_intersection_distances(&g1, &g2, precision);
        assert_eq!(0.0, distance1.0);
        assert_eq!(0.0, distance2.0);
        assert!(is_within_tolerance(
            73.67369956843568,
            Degrees::from(angle).0,
            f64::EPSILON
        ));
        assert_eq!(0, iterations);
    }

    #[test]
    fn test_non_intersection_same_geodesic() {
        // A GeodesicSegment along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(50.0), Degrees(0.0));
        let g1 = GeodesicSegment::from((&a, &b));

        // Convert 1mm precision to Radians
        let precision = Radians(Metres(1e-3).0 / g1.ellipsoid().a().0);

        let c = LatLong::new(Degrees(55.0), Degrees(0.0));
        let d = LatLong::new(Degrees(60.0), Degrees(0.0));
        let g2 = GeodesicSegment::from((&c, &d));

        let (distance1, distance2, angle, iterations) =
            calculate_sphere_intersection_distances(&g1, &g2, precision);
        assert_eq!(0.17463328753266863, distance1.0);
        assert_eq!(0.0, distance2.0);
        assert_eq!(0.0, Degrees::from(angle).0);
        assert_eq!(0, iterations);
    }

    #[test]
    fn test_intersection_same_geodesic_split() {
        // A GeodesicSegment along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(1.0), Degrees(0.0));
        let b = LatLong::new(Degrees(-0.998286322222), Degrees(179.296674991667));
        let g = GeodesicSegment::from((&a, &b));

        // Split g into two geodesics
        let half_length = Metres(g.length().0 / 2.0);
        let half_arc_length = g.metres_to_radians(half_length);

        // a geodesic from the start of g to its mid point
        let g1 = GeodesicSegment::new(
            g.beta(),
            g.lon(),
            g.azi(),
            half_arc_length,
            Metres(0.0),
            g.ellipsoid(),
        );
        // a geodesic from the mid point of g to its end
        let half_arc_length_angle = Angle::from(half_arc_length);
        let g2 = GeodesicSegment::new(
            g.arc_beta(half_arc_length_angle),
            g.arc_longitude(half_arc_length),
            g.arc_azimuth(half_arc_length_angle),
            g.arc_length() - half_arc_length,
            Metres(0.0),
            g.ellipsoid(),
        );

        // 1mm precision in Radians on the unit sphere
        let precision = Radians(Metres(1e-3).0 / g.ellipsoid().a().0);

        // geodesics are coincident
        let (distance1, distance2, angle, iterations) =
            calculate_sphere_intersection_distances(&g1, &g2, precision);
        assert!(is_within_tolerance(
            g1.arc_length().0,
            distance1.0,
            f64::EPSILON
        ));
        assert_eq!(0.0, distance2.0);
        assert_eq!(0.0, Degrees::from(angle).0);
        assert_eq!(0, iterations);

        // a geodesic from the mid point of g to another point
        let g3 = GeodesicSegment::new(
            g.arc_beta(half_arc_length_angle),
            g.arc_longitude(half_arc_length),
            g.azi(),
            half_arc_length,
            Metres(0.0),
            g.ellipsoid(),
        );

        // geodesics are NOT coincident
        let (distance1, distance2, angle, iterations) =
            calculate_sphere_intersection_distances(&g1, &g3, precision);
        assert!(is_within_tolerance(
            g1.arc_length().0,
            distance1.0,
            precision.0
        ));
        assert!(is_within_tolerance(0.0, distance2.0, precision.0));
        assert!(is_within_tolerance(
            Degrees::from(g3.azi() - g1.arc_azimuth(half_arc_length_angle)).0,
            Degrees::from(angle).0,
            precision.0
        ));
        assert_eq!(5, iterations);
    }

    #[test]
    fn test_same_geodesic_no_intersection() {
        // A GeodesicSegment along the Equator, heading West
        let a = LatLong::new(Degrees(0.0), Degrees(-4.0));
        let b = LatLong::new(Degrees(0.0), Degrees(0.0));
        let g1 = GeodesicSegment::from((&a, &b));
        let g3 = GeodesicSegment::from((&b, &a));

        // A GeodesicSegment along the Equator, heading East
        let c = LatLong::new(Degrees(0.0), Degrees(0.25));
        let d = LatLong::new(Degrees(0.0), Degrees(04.0));
        let g2 = GeodesicSegment::from((&c, &d));
        let g4 = GeodesicSegment::from((&d, &c));

        // 1m precision in Radians on the unit sphere
        let precision = Radians(Metres(1.0).0 / g1.ellipsoid().a().0);

        let (distance1, distance2, _angle, _iterations) =
            calculate_sphere_intersection_distances(&g1, &g2, precision);
        assert!(is_within_tolerance(
            4.2643_f64.to_radians(),
            distance1.0,
            precision.0
        ));
        assert_eq!(0.0, distance2.0);

        let (distance1, distance2, _angle, _iterations) =
            calculate_sphere_intersection_distances(&g3, &g2, precision);
        assert!(is_within_tolerance(
            -0.25084_f64.to_radians(),
            distance1.0,
            precision.0
        ));
        assert_eq!(-0.00437800174091304, distance1.0);
        assert_eq!(0.0, distance2.0);

        let (distance1, distance2, _angle, _iterations) =
            calculate_sphere_intersection_distances(&g1, &g4, precision);
        assert!(is_within_tolerance(
            4.2643_f64.to_radians(),
            distance1.0,
            precision.0
        ));
        assert_eq!(g4.arc_length().0, distance2.0);

        let (distance1, distance2, _angle, _iterations) =
            calculate_sphere_intersection_distances(&g3, &g4, precision);
        assert_eq!(0.0, distance1.0);
        assert!(is_within_tolerance(
            4.01345_f64.to_radians(),
            distance2.0,
            precision.0
        ));
    }

    #[test]
    fn test_calculate_arc_reference_distances_and_angle_coincident_great_circles() {
        let latlong_w1 = LatLong::new(Degrees(0.0), Degrees(-1.0));
        let latlong_e1 = LatLong::new(Degrees(0.0), Degrees(1.0));
        let g_0 = GeodesicSegment::from((&latlong_w1, &latlong_e1));

        // 1m precision in Radians on the unit sphere
        let precision = Radians(Metres(1.0).0 / g_0.ellipsoid().a().0);

        // same segments
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_0, precision);
        assert_eq!(Radians(0.0), result.0);
        assert_eq!(Radians(0.0), result.1);
        assert_eq!(Degrees(0.0), Degrees::from(result.2));

        // opposite segments and same geodesic paths
        let latlong_w179 = LatLong::new(Degrees(0.0), Degrees(-179.0));
        let latlong_e179 = LatLong::new(Degrees(0.0), Degrees(179.0));
        let g_1 = GeodesicSegment::from((&latlong_e179, &latlong_w179));
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_1, precision);
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
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_2, precision);
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

        // intersection, same mid points, acute angle
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_1, precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(0.0, result.1.0, precision.0));
        assert_eq!(Degrees(90.0), Degrees::from(result.2));

        let angle = Degrees(44.81195977064123);

        let latlong_sw1 = LatLong::new(Degrees(-1.0), Degrees(-1.0));
        let latlong_ne1 = LatLong::new(Degrees(1.0), Degrees(1.0));
        let g_2 = GeodesicSegment::from((&latlong_sw1, &latlong_ne1));
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_2, precision);
        assert!(is_within_tolerance(0.0, result.0.0, precision.0));
        assert!(is_within_tolerance(0.0, result.1.0, precision.0));
        assert!(is_within_tolerance(
            angle.0,
            Degrees::from(result.2).0,
            precision.0
        ));

        // intersection, same mid points, obtuse angle
        let g_3 = GeodesicSegment::from((&latlong_ne1, &latlong_sw1));
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_3, precision);
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
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_4, precision);
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
        let result = calculate_arc_reference_distances_and_angle(&g_0, &g_5, precision);
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
