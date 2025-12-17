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
use csv::{ReaderBuilder, WriterBuilder};
use icao_wgs84::{Metres, WGS84_ELLIPSOID, geodesic};
use std::env;
use std::path::Path;
use std::time::Instant;

use unit_sphere::{LatLong, great_circle};

// The location of the file on sourceforge.net
// const FILEPATH: &str = "https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz/download";

//  random_df = tests_df[:100000]
//  antipodal_df = tests_df[100000:150000]
//  short_df = tests_df[150000:200000]
//  one_pole_df = tests_df[200000:250000]
//  two_poles_df = tests_df[250000:300000]
//  near_meridional_df = tests_df[300000:350000]
//  near_equatorial_df = tests_df[350000:400000]
//  between_vertices_df = tests_df[400000:450000]
//  end_by_vertices_df = tests_df[450000:500000]

// The columns of the data file.
const LAT_1: usize = 0;
const LON_1: usize = 1;
const AZI_1: usize = 2;
const LAT_2: usize = 3;
const LON_2: usize = 4;
const AZI_2: usize = 5;
const D_METRES: usize = 6;
const D_DEGREES: usize = 7;
// const M12: usize = 8;
// const AREA: usize = 9;

const OUTPUT_FILE_PATH: &str = "icao_wgs84_data_rust.dat";

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
fn test_geodesic_examples() -> Result<(), Box<dyn std::error::Error>> {
    // Read GEODTEST_DIR/GeodTest.dat file and run tests
    let filename = "GeodTest.dat";
    let dir_key = "GEODTEST_DIR";

    let p = env::var(dir_key).expect("Environment variable not found: GEODTEST_DIR");
    let path = Path::new(&p);
    let file_path = path.join(filename);
    let mut csv_reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_path(file_path)?;

    let mut csv_writer = WriterBuilder::new()
        .has_headers(false)
        .delimiter(b' ')
        .from_path(OUTPUT_FILE_PATH)?;

    let start_time = Instant::now();

    let mut index = 0;
    for result in csv_reader.records() {
        let record = result.unwrap();

        let lat1 = record[LAT_1].parse::<f64>().unwrap();
        let lon1 = record[LON_1].parse::<f64>().unwrap();
        let azi1 = record[AZI_1].parse::<f64>().unwrap();
        let lat2 = record[LAT_2].parse::<f64>().unwrap();
        let lon2 = record[LON_2].parse::<f64>().unwrap();
        let azi2 = record[AZI_2].parse::<f64>().unwrap();
        let d_metres = record[D_METRES].parse::<f64>().unwrap();
        let _d_degrees = record[D_DEGREES].parse::<f64>().unwrap();

        let (azi, end_azi, distance_m, arc_length, iterations) =
            calculate_geodesic_inverse_values(lat1, lon1, lat2, lon2);

        // Convert azimuths to degrees
        let azi = Degrees::from(azi).0;
        let end_azi = Degrees::from(end_azi).0;

        // Compare start azimuths
        let delta_azimuth = (azi1 - azi).abs();
        // reduce tolerance for entries running between or close to vertices
        let azimuth_tolerance = if index <= 400000 { 5.331e-5 } else { 2.0e-2 };
        if azimuth_tolerance < delta_azimuth {
            panic!(
                "azimuth, line: {:?} lat1: {:?} delta: {:?} azimuth: {:?} calculated: {:?} delta_long: {:?} ",
                index, lat1, delta_azimuth, azi1, azi, lon2
            );
        }

        // Compare end_azimuths
        let delta_azimuth = (azi2 - end_azi).abs();
        if azimuth_tolerance < delta_azimuth {
            panic!(
                "end azimuth, line: {:?} delta: {:?} azimuth: {:?} delta_long: {:?} ",
                index, delta_azimuth, azi2, lon2
            );
        }

        // Compare geodesic distances
        let delta_length_m = (d_metres - distance_m.0).abs();
        // if a short geodesic, test delta length, not delta length ratio
        if index >= 150000 && index < 200000 {
            if 3.5e-9 < delta_length_m {
                panic!(
                    "length, line: {:?} delta: {:?} length: {:?} result: {:?} ",
                    index, delta_length_m, d_metres, distance_m
                );
            }
        } else {
            let delta_length_m_ratio = delta_length_m / d_metres;
            if 9.0e-13 < delta_length_m_ratio {
                panic!(
                    "length, line: {:?} delta ratio: {:?} length: {:?} result: {:?} ",
                    index, delta_length_m_ratio, d_metres, distance_m
                );
            }
        }

        // Output the values to a data file in the same format as the input file
        let fields: [String; 9] = [
            record[LAT_1].to_owned(),
            record[LON_1].to_owned(),
            azi.to_string(),
            record[LAT_2].to_owned(),
            record[LON_2].to_owned(),
            end_azi.to_string(),
            distance_m.0.to_string(),
            arc_length.0.to_degrees().to_string(),
            iterations.to_string(),
        ];
        csv_writer.write_record(&fields)?;

        index += 1;
    }

    println!("Elapsed time: {:.3?}", start_time.elapsed());

    csv_writer.flush()?;

    Ok(())
}
