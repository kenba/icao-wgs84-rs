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

// extern crate we're testing, same as any other code would do.
extern crate icao_wgs84;

use angle_sc::{Angle, Degrees, Radians};
use icao_wgs84::{geodesic, Metres, WGS84_ELLIPSOID};
use itertools::multizip;
use polars::prelude::*;
use std::env;
use std::path::Path;
use unit_sphere::{great_circle, LatLong};

// The location of the file on sourceforge.net
// const FILEPATH: &str = "https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz/download";

#[test]
#[ignore]
fn test_geodesic_examples() -> Result<(), Box<dyn std::error::Error>> {
    // Read GEODTEST_DIR/GeodTest.dat file and run tests
    let filename = "GeodTest.dat";
    let dir_key = "GEODTEST_DIR";

    let p = env::var(dir_key).expect("Environment variable not found: GEODTEST_DIR");
    let path = Path::new(&p);
    let file_path = path.join(filename);

    let df = CsvReadOptions::default()
        .with_has_header(false)
        .with_parse_options(CsvParseOptions::default().with_separator(b' '))
        .try_into_reader_with_file_path(Some(file_path.into()))?
        .finish()?;

    let lines = df.height();

    // println!("{:?}", df.head(None));

    // Use multizip to access DataFrame by rows
    let objects = df.take_columns();
    let dep_lats = objects[0].f64()?.iter();
    // let dep_lons = objects[1].i64()?.iter(); // departure longitudes are all zero
    let dep_azis = objects[2].f64()?.iter();
    let arr_lats = objects[3].f64()?.iter();
    let arr_lons = objects[4].f64()?.iter();
    let arr_azis = objects[5].f64()?.iter();
    let distances = objects[6].f64()?.iter();
    let combined = multizip((
        dep_lats,
        dep_azis,
        arr_lats,
        arr_lons,
        arr_azis,
        distances,
    ));

    let mut invalid_tests = 0;
    let mut iterations = 0;
    for (index, (lat1, azi1, lat2, lon2, azi2, d_metres)) in combined.enumerate() {
        let lat1 = Degrees(lat1.unwrap());
        let lon1 = Degrees(0.0);
        let lat2 = Degrees(lat2.unwrap());
        let lon2 = Degrees(lon2.unwrap());
        let azi1 = Degrees(azi1.unwrap());
        let azi2 = Degrees(azi2.unwrap());
        let d_metres = Metres(d_metres.unwrap());

        // path crosses the equator from North to South, but GeodTest.dat azimuth is northerly
        if (lat1.0 > 0.0) && (lat2.0 < 0.0) && (azi1.0 < 90.0) {
            invalid_tests += 1;
        }

        // panic!("lon2; {:?}", lon2);
        let a = LatLong::new(lat1, lon1);
        let b = LatLong::new(lat2, lon2);
        let result = geodesic::calculate_azimuths_aux_length(
            &a,
            &b,
            Radians(great_circle::MIN_VALUE),
            &WGS84_ELLIPSOID,
        );
        iterations += result.3;

        let azi = Degrees::from(result.0);
        let delta_azimuth = libm::fabs(azi1.0 - azi.0);
        // reduce tolerance for entries running between or close to vertices
        let azimuth_tolerance = if index <= 400000 { 5.331e-5 } else { 2.0e-2 };
        if azimuth_tolerance < delta_azimuth {
            panic!(
                "azimuth, line: {:?} lat1: {:?} delta: {:?} azimuth: {:?} calculated: {:?} delta_long: {:?} ",
                index, lat1, delta_azimuth, azi1, azi, lon2
            );
        }

        // Compare end_azimuths
        let delta_azimuth = libm::fabs(azi2.0 - Degrees::from(result.2).0);
        if azimuth_tolerance < delta_azimuth {
            panic!(
                "end azimuth, line: {:?} delta: {:?} azimuth: {:?} delta_long: {:?} ",
                index, delta_azimuth, azi2, lon2
            );
        }

        let beta1 = WGS84_ELLIPSOID.calculate_parametric_latitude(Angle::from(lat1));
        let result_m =
            geodesic::convert_radians_to_metres(beta1, result.0, result.1, &WGS84_ELLIPSOID);

        let delta_length_m = libm::fabs(d_metres.0 - result_m.0);
        // if a short geodesic, test delta length, not delta length ratio
        if index >= 150000 && index < 200000 {
            if 3.5e-9 < delta_length_m {
                panic!(
                    "length, line: {:?} delta: {:?} length: {:?} result: {:?} ",
                    index, delta_length_m, d_metres, result_m
                );
            }
        } else {
            let delta_length_m_ratio = delta_length_m / d_metres.0;
            if 2.0e-10 < delta_length_m_ratio {
                panic!(
                    "length, line: {:?} delta ratio: {:?} length: {:?} result: {:?} ",
                    index, delta_length_m_ratio, d_metres, result_m
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
    }

    println!("lines,iterations, invalid_tests");
    println!("{},{},{}", lines, iterations, invalid_tests);

    Ok(())
}
