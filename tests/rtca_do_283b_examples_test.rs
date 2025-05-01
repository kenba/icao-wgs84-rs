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
use itertools::multizip;
use polars::prelude::*;
use unit_sphere::{great_circle, LatLong};

const FILENAME: &str = "data/rtca_do_283b_geodesic_examples.csv";

#[test]
#[ignore]
fn test_rtca_do_283b_examples() -> Result<(), Box<dyn std::error::Error>> {
    // Read file and run tests
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(FILENAME.into()))?
        .finish()?;

    // println!("{:?}", df.head(None));

    // Use multizip to access DataFrame by rows
    let objects = df.take_columns();
    let dep_lats = objects[0].f64()?.iter();
    let dep_lons = objects[1].f64()?.iter();
    let arr_lats = objects[2].f64()?.iter();
    let arr_lons = objects[3].f64()?.iter();
    let dep_azis = objects[4].f64()?.iter();
    let ranges = objects[6].f64()?.iter();
    let range_errors = objects[7].f64()?.iter();
    let combined = multizip((
        dep_lats,
        dep_lons,
        arr_lats,
        arr_lons,
        dep_azis,
        ranges,
        range_errors,
    ));

    println!("line_number,iterations,delta_azimuth,delta_length,delta_error");
    for (index, (lat1, lon1, lat2, lon2, azi1, d_metres, range_error)) in combined.enumerate() {
        let lat1 = Degrees(lat1.unwrap());
        let lon1 = Degrees(lon1.unwrap());
        let lat2 = Degrees(lat2.unwrap());
        let lon2 = Degrees(lon2.unwrap());
        let azi1 = Degrees(azi1.unwrap());
        // let _azi2 = Degrees(azi2.unwrap());
        let d_metres = Metres(d_metres.unwrap());
        let range_error = range_error.unwrap();

        let a = LatLong::new(lat1, lon1);
        let b = LatLong::new(lat2, lon2);
        let result = geodesic::calculate_azimuths_arc_length(
            &a,
            &b,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );

        let delta_azimuth = libm::fabs(azi1.0 - Degrees::from(result.0).0);

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(lat1));
        let result_m =
            geodesic::convert_radians_to_metres(beta1, result.0, result.1, &WGS84_ELLIPSOID);

        let delta_length_m = libm::fabs(d_metres.0 - result_m.0);
        let delta_error_m = range_error - delta_length_m;
        println!(
            "{},{},{},{},{}",
            index + 1,
            result.3,
            delta_azimuth,
            delta_length_m,
            delta_error_m
        );
    }

    Ok(())
}
