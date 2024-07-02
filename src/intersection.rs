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

//! The `intersection` module contains functions for calculating geodesic
//! intersections using vectors.
//!
//! Charles Karney's developed a solution to the
//! [Intersection between two geodesic lines](https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a).
//! using a gnomonic projection.
//!
//! [Baselga and Martinez-Llario(2017)](https://www.researchgate.net/publication/321358300_Intersection_and_point-to-line_solutions_for_geodesics_on_the_ellipsoid)
//! solve geodesic intersection and point-to-line problems by using the
//! correspondence between geodesics on an ellipsoid and great-circles on the
//! auxiliary sphere.
//!
//! Nick Korbey and I [Barker and Korbey(2019)](https://www.researchgate.net/publication/335749834_Geodesic_Geometry)
//! developed Baselga and Martinez-Llario's algorithms by using vectors to
//! solve geodesic intersection and point-to-line problems on the auxiliary
//! sphere.
//!
//! [Karney(2023)](https://arxiv.org/abs/2308.00495) also developed Baselga and
//! Martinez-Llario's algorithms and provides reliable algorithms for solving
//! geodesic intersections on both oblate and prolate ellipsoids.
//!
//! This `intersection` module calculates initial great-circle intersection
//! distances on the auxiliary sphere and then refines the intersection
//! distances using great-circles on the auxiliary sphere at the calculated
//! intersection distances.

use crate::{Geodesic, Radians};
use unit_sphere::{great_circle, vector};

/// Calculate the distances along a pair of Geodesics (in Radians) to their
/// closest intersection point given the initial great-circle distances.
/// * `g1`, `g2` the Geodesics.
/// * `sq_precision` the square of the Euclidean precision.
/// * `use_antipodal_intersection` use the antipodal intersection point
/// * `initial_distances` the initial intersection distances in `Radians`.
///
/// returns the distances along a the Geodesics (in Radians) and
/// the number of iterations required.
#[must_use]
pub fn calculate_geodesic_intersection_distances(
    g1: &Geodesic,
    g2: &Geodesic,
    sq_precision: f64,
    use_antipodal_intersection: bool,
    initial_distances: (Radians, Radians),
) -> (Radians, Radians, u32) {
    const MAX_ITERATIONS: u32 = 10;

    let (mut distance1, mut distance2) = initial_distances;
    let mut iterations = 1;
    while iterations < MAX_ITERATIONS {
        let (pos1, pole1) = g1.aux_point_and_pole(distance1);
        let (pos2, pole2) = g2.aux_point_and_pole(distance2);

        let sq_d = vector::sq_distance(&pos1, &pos2);
        if sq_d < sq_precision {
            break;
        }

        iterations += 1;

        let c = if use_antipodal_intersection {
            vector::intersection::calculate_intersection_point(&pole2, &pole1)
        } else {
            vector::intersection::calculate_intersection_point(&pole1, &pole2)
        };
        match c {
            Some(c) => {
                let (delta1, delta2) = vector::intersection::calculate_intersection_distances(
                    &pos1, &pole1, &pos2, &pole2, &c,
                );
                distance1 = distance1 + delta1;
                distance2 = distance2 + delta2;
            }
            None => break,
        }
    }

    (distance1, distance2, iterations)
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
pub fn calculate_aux_intersection_distances(
    g1: &Geodesic,
    g2: &Geodesic,
    precision: Radians,
) -> (Radians, Radians, u32) {
    // The Geodesics MUST be on the same `Ellipsoid`
    assert!(g1.ellipsoid() == g2.ellipsoid());

    // Convert precision in `Radians` to the square of Euclidean precision.
    let e_precision = great_circle::gc2e_distance(precision);
    let sq_precision = e_precision * e_precision;

    // Get the start points and poles
    let (a1, pole1) = g1.aux_point_and_pole(Radians(0.0));
    let (a2, pole2) = g2.aux_point_and_pole(Radians(0.0));

    // if the start points are within precision of each other
    let sq_d = vector::sq_distance(&a1, &a2);
    if sq_d < sq_precision {
        (Radians(0.0), Radians(0.0), 0)
    } else {
        // Determine whether the great circles on the auxiliary sphere are coincident
        vector::intersection::calculate_intersection_point(&pole1, &pole2).map_or_else(
            || {
                let distances = vector::intersection::calculate_coincident_arc_distances(
                    vector::calculate_great_circle_atd(&a1, &pole1, &a2),
                    pole1.dot(&pole2) < 0.0,
                    g1.aux_length(),
                    g2.aux_length(),
                );
                (distances.0, distances.1, 0)
            },
            |c| {
                let centre = vector::normalise(&(g1.mid_point() + g2.mid_point()));
                let use_antipodal_intersection =
                    centre.is_some_and(|x| vector::sq_distance(&c, &x) > 2.0);
                let d = if use_antipodal_intersection { -c } else { c };
                let initial_distances = vector::intersection::calculate_intersection_distances(
                    &a1, &pole1, &a2, &pole2, &d,
                );
                calculate_geodesic_intersection_distances(
                    g1,
                    g2,
                    sq_precision,
                    use_antipodal_intersection,
                    initial_distances,
                )
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Ellipsoid, Geodesic, Metres};
    use unit_sphere::{Degrees, LatLong, Radians};

    #[test]
    fn test_intersection_same_start_point() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // Convert 1mm precision to Radians
        let precision = Radians(Metres(1e-3).0 / wgs84_ellipsoid.a().0);

        // A Geodesic along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(45.0), Degrees(180.0));
        let g1 = Geodesic::from((&a, &b, &wgs84_ellipsoid));

        let c = LatLong::new(Degrees(45.0), Degrees(45.0));
        let g2 = Geodesic::from((&a, &c, &wgs84_ellipsoid));

        let (distance1, distance2, iterations) =
            calculate_aux_intersection_distances(&g1, &g2, precision);
        println!(
            "calculate_aux_intersection_distances iterations: {:?}",
            iterations
        );
        assert_eq!(0.0, distance1.0);
        assert_eq!(0.0, distance2.0);
    }

    #[test]
    fn test_non_intersection_same_geodesic() {
        let wgs84_ellipsoid = Ellipsoid::wgs84();

        // Convert 1mm precision to Radians
        let precision = Radians(Metres(1e-3).0 / wgs84_ellipsoid.a().0);

        // A Geodesic along the Greenwich meridian, over the North pole and down the IDL
        let a = LatLong::new(Degrees(45.0), Degrees(0.0));
        let b = LatLong::new(Degrees(50.0), Degrees(0.0));
        let g1 = Geodesic::from((&a, &b, &wgs84_ellipsoid));

        let c = LatLong::new(Degrees(55.0), Degrees(0.0));
        let d = LatLong::new(Degrees(60.0), Degrees(0.0));
        let g2 = Geodesic::from((&c, &d, &wgs84_ellipsoid));

        let (distance1, distance2, iterations) =
            calculate_aux_intersection_distances(&g1, &g2, precision);
        assert_eq!(0.17463328753266863, distance1.0);
        assert_eq!(0.0, distance2.0);
        assert_eq!(0, iterations);
    }
}
