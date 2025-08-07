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
use polars::prelude::*;
use std::env;
use std::path::Path;
use std::time::Instant;

use unit_sphere::{great_circle, LatLong};

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

/// Read the geodesic data file at file_path
fn read_geodesic_data_file(file_path: &str) -> Result<DataFrame, Box<dyn std::error::Error>> {
    // The Schema of the geodesic data file
    let mut schema = Schema::with_capacity(10);
    schema.insert("lat1".into(), DataType::Float64);
    schema.insert("lon1".into(), DataType::Int64);
    schema.insert("azi1".into(), DataType::Float64);
    schema.insert("lat2".into(), DataType::Float64);
    schema.insert("lon2".into(), DataType::Float64);
    schema.insert("azi2".into(), DataType::Float64);
    schema.insert("distance_metres".into(), DataType::Float64);
    schema.insert("arc_length".into(), DataType::Float64);
    schema.insert("m12".into(), DataType::Float64);
    schema.insert("s12".into(), DataType::Float64);

    // Read all of the columns of the geodesic data file
    let lf = LazyCsvReader::new(PlPath::new(file_path))
        .with_has_header(false)
        .with_separator(b' ')
        .with_schema(Some(Arc::new(schema)))
        .finish()?
        .collect()?;

    Ok(lf)
}

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
    let filename = "GeodTest.dat.gz";
    let dir_key = "GEODTEST_DIR";

    let p = env::var(dir_key).expect("Environment variable not found: GEODTEST_DIR");
    let path = Path::new(&p);
    let file_path = path.join(filename);

    let lf = read_geodesic_data_file(file_path.to_str().expect("Bad file_path"))?;
    // println!("{lf}");

    let start_time = Instant::now();

    let length = lf.height();
    let mut azi_v = Vec::with_capacity(length);
    let mut end_azi_v = Vec::with_capacity(length);
    let mut distance_v = Vec::with_capacity(length);
    let mut arc_length_v = Vec::with_capacity(length);
    let mut iterations_v = Vec::with_capacity(length);

    let mut index = 0;
    lf.column("lat1")?
        .f64()?
        .into_no_null_iter()
        .zip(lf.column("lat2")?.f64()?.into_no_null_iter())
        .zip(lf.column("lon2")?.f64()?.into_no_null_iter())
        .zip(lf.column("azi1")?.f64()?.into_no_null_iter())
        .zip(lf.column("azi2")?.f64()?.into_no_null_iter())
        .zip(lf.column("distance_metres")?.f64()?.into_no_null_iter())
        .for_each(|(((((lat1, lat2), lon2), azi1), azi2), d_metres)| {
            let (azi, end_azi, distance_m, arc_length, iterations) =
                calculate_geodesic_inverse_values(lat1, 0.0, lat2, lon2);

            // Convert azimuths to degrees
            let azi = Degrees::from(azi).0;
            let end_azi = Degrees::from(end_azi).0;

            // Compare start azimuths
            let delta_azimuth = libm::fabs(azi1 - azi);
            // reduce tolerance for entries running between or close to vertices
            let azimuth_tolerance = if index <= 400000 { 5.331e-5 } else { 2.0e-2 };
            if azimuth_tolerance < delta_azimuth {
                panic!(
                    "azimuth, line: {:?} lat1: {:?} delta: {:?} azimuth: {:?} calculated: {:?} delta_long: {:?} ",
                    index, lat1, delta_azimuth, azi1, azi, lon2
                );
            }

            // Compare end_azimuths
            let delta_azimuth = libm::fabs(azi2 - end_azi);
            if azimuth_tolerance < delta_azimuth {
                panic!(
                    "end azimuth, line: {:?} delta: {:?} azimuth: {:?} delta_long: {:?} ",
                    index, delta_azimuth, azi2, lon2
                );
            }

            // Compare geodesic distances
            let delta_length_m = libm::fabs(d_metres - distance_m.0);
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

            // push values into vectors
            azi_v.push(azi);
            end_azi_v.push(end_azi);
            distance_v.push(distance_m.0);
            arc_length_v.push(arc_length.0.to_degrees());
            iterations_v.push(iterations);

            index += 1;
        });

    let azi_column = Column::new("azi".into(), azi_v);
    let end_azi_column = Column::new("end_azi`".into(), end_azi_v);
    let distance_column = Column::new("distance_m".into(), distance_v);
    let arc_length_column = Column::new("arc_length".into(), arc_length_v);
    let iterations_column = Column::new("iterations".into(), iterations_v);
    let mut df = DataFrame::new(vec![
        lf.column("lat1")?.clone(),
        lf.column("lon1")?.clone(),
        azi_column,
        lf.column("lat2")?.clone(),
        lf.column("lon2")?.clone(),
        end_azi_column,
        distance_column,
        arc_length_column,
        iterations_column,
    ])?;

    // println!("{df}");
    println!("Elapsed time: {:.3?}", start_time.elapsed());

    // Output the values to a data file in the same format as the input file
    let output_file_name: String = "icao_wgs84_data_rust.dat".to_owned();
    let mut file =
        std::fs::File::create(&output_file_name).expect("Could not write to output file");
    CsvWriter::new(&mut file)
        .include_header(false)
        .with_separator(b' ')
        .finish(&mut df)?;
    println!("Written inverse geodesic values to: {}", &output_file_name);

    Ok(())
}
