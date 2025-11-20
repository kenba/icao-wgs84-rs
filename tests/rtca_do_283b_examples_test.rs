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
use csv;
use icao_wgs84::{Metres, WGS84_ELLIPSOID, geodesic};
use unit_sphere::{LatLong, great_circle};

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

const FILENAME: &str = "data/rtca_do_283b_geodesic_examples.csv";

type DataRecord = (f64, f64, f64, f64, f64, f64, f64, f64);

#[test]
#[ignore]
fn test_rtca_do_283b_examples() -> Result<(), Box<dyn std::error::Error>> {
    let mut rdr = csv::Reader::from_path(FILENAME)?;
    println!("index,end_azimuth,delta_azimuth,delta_length_m,delta_error_m");

    let mut index = 0;
    for result in rdr.deserialize::<DataRecord>() {
        let record = result?;
        // println!("{:?}", record);

        let lat1 = record.0;
        let lon1 = record.1;
        let lat2 = record.2;
        let lon2 = record.3;
        let azi1 = record.4;
        let _azi2 = record.5;
        let d_metres = record.6;
        let range_error = record.7;

        let (azi, end_azi, distance_m, _arc_length, _iterations) =
            calculate_geodesic_inverse_values(lat1, lon1, lat2, lon2);

        // Convert azimuths to degrees
        let azi = Degrees::from(azi).0;
        let end_azi = Degrees::from(end_azi).0;

        // Calculate differences
        let delta_azimuth = (azi1 - azi).abs();
        let delta_length_m = (d_metres - distance_m.0).abs();
        let delta_error_m = range_error - delta_length_m;

        // Output values
        println!("{index},{end_azi},{delta_azimuth},{delta_length_m},{delta_error_m}");

        index += 1;
    }

    Ok(())
}
