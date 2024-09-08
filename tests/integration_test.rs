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

// extern crate we're testing, same as any other code would do.
extern crate icao_wgs84;

use angle_sc::{Angle, Degrees};
use csv::ReaderBuilder;
use icao_wgs84::{geodesic, Metres, WGS84_ELLIPSOID};
use std::env;
use std::path::Path;
use unit_sphere::LatLong;

#[test]
#[ignore]
fn test_geodesic_examples() {
    // Read GEODTEST_DIR/GeodTest.dat file and run tests
    let filename = "GeodTest.dat";
    let dir_key = "GEODTEST_DIR";

    let p = env::var(dir_key).expect("Environment variable not found: GEODTEST_DIR");
    let path = Path::new(&p);
    let file_path = path.join(filename);
    let mut csv_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_path(file_path)
        .expect("Could not read file: GeodTest.dat");
    let mut line_number = 1;
    for result in csv_reader.records() {
        let record = result.unwrap();
        // panic!("first record; {:?}", record);

        let lat1 = Degrees(record[0].parse::<f64>().unwrap());
        let lon1 = Degrees(record[1].parse::<f64>().unwrap());
        let azi1 = Degrees(record[2].parse::<f64>().unwrap());
        let lat2 = Degrees(record[3].parse::<f64>().unwrap());
        let lon2 = Degrees(record[4].parse::<f64>().unwrap());
        let _azi2 = Degrees(record[5].parse::<f64>().unwrap());
        let d_metres = Metres(record[6].parse::<f64>().unwrap());
        let d_degrees = Degrees(record[7].parse::<f64>().unwrap());

        // panic!("lon2; {:?}", lon2);
        let a = LatLong::new(lat1, lon1);
        let b = LatLong::new(lat2, lon2);
        let result = geodesic::calculate_azimuth_aux_length(&a, &b, &WGS84_ELLIPSOID);

        let delta_azimuth = libm::fabs(azi1.0 - Degrees::from(result.0).0);
        // reduce tolerance for entries running between or close to vertices
        let azimuth_tolerance = if line_number <= 400000 {
            5.331e-5
        } else {
            0.077
        };
        if azimuth_tolerance < delta_azimuth {
            panic!(
                "azimuth, line: {:?} delta: {:?} azimuth: {:?} delta_long: {:?} ",
                line_number, delta_azimuth, azi1, lon2
            );
        }

        let delta_length = libm::fabs(d_degrees.0.to_radians() - (result.1).0);
        if 3.0e-10 < delta_length {
            panic!(
                "length, line: {:?} delta: {:?} length: {:?} delta_long: {:?} ",
                line_number, delta_length, d_degrees, lon2
            );
        }

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(lat1));
        let result_m =
            geodesic::convert_radians_to_metres(beta1, result.0, result.1, &WGS84_ELLIPSOID);

        let delta_length_m = libm::fabs(d_metres.0 - result_m.0);
        // if a short geodesic, test delta length, not delta length ratio
        if line_number >= 150000 && line_number < 200000 {
            if 9.0e-5 < delta_length_m {
                panic!(
                    "length, line: {:?} delta: {:?} length: {:?} result: {:?} ",
                    line_number, delta_length_m, d_metres, result_m
                );
            }
        } else {
            let delta_length_m_ratio = delta_length_m / d_metres.0;
            if 2.5e-9 < delta_length_m_ratio {
                panic!(
                    "length, line: {:?} delta ratio: {:?} length: {:?} result: {:?} ",
                    line_number, delta_length_m_ratio, d_metres, result_m
                );
            }
        }

        //  random_df = tests_df[:100000]
        //  antipodal_df = tests_df[100000:150000]
        //  short_df = tests_df[150000:200000]
        //  one_pole_df = tests_df[200000:250000]
        //  two_poles_df = tests_df[250000:300000]
        //  near_meridional_df = tests_df[300000:350000]
        //  near_equatorial_df = tests_df[350000:400000]
        //  between_vertices_df = tests_df[400000:450000]
        //  end_by_vertices_df = tests_df[450000:500000]
        line_number += 1;
    }
}
