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

extern crate icao_wgs84;

use angle_sc::{Angle, Degrees};
use csv::ReaderBuilder;
use icao_wgs84::{geodesic, Metres, WGS84_ELLIPSOID};
use std::path::Path;
use unit_sphere::LatLong;

#[test]
#[ignore]
fn test_rtca_do_283b_examples() {
    let filename = "rtca_do_283b_geodesic_examples.csv";
    let path = Path::new("data");
    let file_path = path.join(filename);
    let mut csv_reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b',')
        .from_path(file_path)
        .expect("Could not read file: rtca_do_283b_geodesic_examples.csv");

    let mut line_number = 1;
    println!("line_number,iterations,delta_azimuth,delta_length,delta_error");
    for result in csv_reader.records() {
        let record = result.unwrap();
        let lat1 = Degrees(record[0].parse::<f64>().unwrap());
        let lon1 = Degrees(record[1].parse::<f64>().unwrap());
        let lat2 = Degrees(record[2].parse::<f64>().unwrap());
        let lon2 = Degrees(record[3].parse::<f64>().unwrap());
        let azi1 = Degrees(record[4].parse::<f64>().unwrap());
        let _azi2 = Degrees(record[5].parse::<f64>().unwrap());
        let d_metres = Metres(record[6].parse::<f64>().unwrap());
        let range_error = record[7].parse::<f64>().unwrap();

        let a = LatLong::new(lat1, lon1);
        let b = LatLong::new(lat2, lon2);
        let result = geodesic::calculate_azimuth_aux_length(&a, &b, &WGS84_ELLIPSOID);

        let delta_azimuth = libm::fabs(azi1.0 - Degrees::from(result.0).0);

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(lat1));
        let result_m =
            geodesic::convert_radians_to_metres(beta1, result.0, result.1, &WGS84_ELLIPSOID);

        let delta_length_m = libm::fabs(d_metres.0 - result_m.0);
        let delta_error_m = range_error - delta_length_m;
        println!(
            "{},{},{},{}",
            line_number, delta_azimuth, delta_length_m, delta_error_m
        );

        line_number += 1;
    }
}
