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

extern crate icao_wgs84;

use angle_sc::{Angle, Degrees, Radians};
use icao_wgs84::{geodesic, Metres, WGS84_ELLIPSOID};
use polars::prelude::*;
use unit_sphere::{great_circle, LatLong};

const FILENAME: &str = "data/rtca_do_283b_geodesic_examples.csv";

/// Calculate the geodesic values for the given start and end point latitudes and longitudes.
///
/// returns start_azimuth, end_azimuth, distance, arc_distance, iterations
fn calculate_geodesic_inverse_values(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
) -> (Angle, Angle, Metres, Radians, u32) {
    let lat1 = Degrees(lat1);
    let a = LatLong::new(lat1, Degrees(lon1));
    let b = LatLong::new(Degrees(lat2), Degrees(lon2));
    let result = geodesic::calculate_azimuths_arc_length(
        &a,
        &b,
        Radians(great_circle::MIN_VALUE),
        &WGS84_ELLIPSOID,
    );

    let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(lat1));
    let result_m = geodesic::convert_radians_to_metres(beta1, result.0, result.1, &WGS84_ELLIPSOID);

    (result.0, result.2, result_m, result.1, result.3)
}

#[test]
#[ignore]
fn test_rtca_do_283b_examples() -> Result<(), Box<dyn std::error::Error>> {
    let lf = LazyCsvReader::new(PlPath::new(FILENAME))
        .finish()?
        .collect()?;
    // println!("{lf}");

    let length = lf.height();
    let mut index_v = Vec::with_capacity(length);
    let mut end_azi_v = Vec::with_capacity(length);
    let mut delta_azimuth_v = Vec::with_capacity(length);
    let mut delta_length_v = Vec::with_capacity(length);
    let mut delta_error_v = Vec::with_capacity(length);

    let mut index = 0;
    lf.column("Departure_Latitude")?
        .f64()?
        .into_no_null_iter()
        .zip(lf.column("Departure_Longitude")?.f64()?.into_no_null_iter())
        .zip(lf.column("Arrival_Latitude")?.f64()?.into_no_null_iter())
        .zip(lf.column("Arrival_Longitude")?.f64()?.into_no_null_iter())
        .zip(lf.column("Departure_Bearing")?.f64()?.into_no_null_iter())
        // .zip(lf.column("Arrival_Bearing")?.f64()?.into_no_null_iter())
        .zip(lf.column("Ranges")?.f64()?.into_no_null_iter())
        .zip(lf.column("Range_Error")?.f64()?.into_no_null_iter())
        .for_each(
            |((((((lat1, lon1), lat2), lon2), azi1), d_metres), range_error)| {
                let (azi, end_azi, distance_m, _arc_length, _iterations) =
                    calculate_geodesic_inverse_values(lat1, lon1, lat2, lon2);

                // Convert azimuths to degrees
                let azi = Degrees::from(azi).0;
                let end_azi = Degrees::from(end_azi).0;

                let delta_azimuth = libm::fabs(azi1 - azi);
                let delta_length_m = libm::fabs(d_metres - distance_m.0);
                let delta_error_m = range_error - delta_length_m;

                // push values into vectors
                index_v.push(index + 1);
                end_azi_v.push(end_azi);
                delta_azimuth_v.push(delta_azimuth);
                delta_length_v.push(delta_length_m);
                delta_error_v.push(delta_error_m);

                index += 1;
            },
        );

    // Create a DataFrame for the output values
    let line_no_column = Column::new("line_no".into(), index_v);
    let end_azi_column = Column::new("end_azi".into(), end_azi_v);
    let delta_azi_column = Column::new("delta_azi".into(), delta_azimuth_v);
    let delta_length_column = Column::new("delta_length".into(), delta_length_v);
    let delta_error_column = Column::new("delta_error".into(), delta_error_v);
    let df = DataFrame::new(vec![
        line_no_column,
        end_azi_column,
        delta_azi_column,
        delta_length_column,
        delta_error_column,
    ])?;

    println!("{df}");

    Ok(())
}
