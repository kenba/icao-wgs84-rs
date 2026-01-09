# icao-wgs84

[![crates.io](https://img.shields.io/crates/v/icao-wgs84.svg)](https://crates.io/crates/icao-wgs84)
[![docs.io](https://docs.rs/icao-wgs84/badge.svg)](https://docs.rs/icao-wgs84/)
[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
[![Rust](https://github.com/kenba/icao-wgs84-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/kenba/icao-wgs84-rs/actions)
[![codecov](https://codecov.io/gh/kenba/icao-wgs84-rs/graph/badge.svg?token=85TJX5VAHF)](https://codecov.io/gh/kenba/icao-wgs84-rs)

A library for performing geometric calculations on the
[WGS-84](https://via-technology.aero/img/navigation/REF08-Doc9674.pdf)
ellipsoid, see *Figure 1*.

![WGS-84 Ellipsoid](https://upload.wikimedia.org/wikipedia/commons/thumb/3/3e/WGS84_mean_Earth_radius.svg/800px-WGS84_mean_Earth_radius.svg.png)\
*Figure 1 The WGS-84 Ellipsoid (not to scale)  
[Cmglee](https://commons.wikimedia.org/wiki/User:Cmglee), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), via Wikimedia Commons*

[WGS-84](https://via-technology.aero/img/navigation/REF08-Doc9674.pdf)
has become the de facto standard for satellite navigation since its adoption
by the Navstar Global Positioning System ([GPS](https://www.gps.gov/systems/gps/performance/accuracy/)) and the USA making GPS available for civilian use in 1983.

This library uses the WGS-84 primary parameters defined in Table 3-1 of the
[ICAO WGS-84 Implementation Manual](https://via-technology.aero/img/navigation/REF08-Doc9674.pdf).

## Geodesic navigation

The shortest path between two points on the surface of an ellipsoid is a
[geodesic segment](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid).
It is the equivalent of a straight line segment in planar geometry or a
[great circle arc](https://en.wikipedia.org/wiki/Great_circle) on the
surface of a sphere, see *Figure 2*.

<img src="https://via-technology.aero/img/navigation/ellipsoid/sphere_mercator_long_geodesic.png" width="600">

*Figure 2 A geodesic segment (orange) and a great circle  arc (blue)*

This library uses the correspondence between geodesic segments on an ellipsoid
and great-circle arcs on a unit sphere, together with 3D vectors to calculate:

- the length and azimuths of a geodesic segment between two positions;
- the along track and across track distances of a position relative to a geodesic segment;
- and the intersection of two geodesic segments.

See: [geodesic algorithms](https://via-technology.aero/navigation/geodesic-algorithms/).

## Design

The library is based on Charles Karney's [GeographicLib](https://geographiclib.sourceforge.io/) library.

Like `GeographicLib`, it models geodesic segments as great circle arcs on
the surface of a unit sphere. However, it also uses vectors to perform
calculations between geodesic segments.

The `Ellipsoid` class represents an ellipsoid of revolution.  
The static `WGS84_ELLIPSOID` represents the WGS-84 `Ellipsoid` which is used
by the `GeodesicSegment` `From` traits to create `GeodesicSegment`s on the WGS-84 `Ellipsoid`.

The library depends upon the following crates:

- [angle-sc](https://crates.io/crates/angle-sc) - to define `Angle`, `Degrees` and `Radians` and perform trigonometric
calculations;
- [unit-sphere](https://crates.io/crates/unit-sphere) - to define `LatLong` and
perform great-circle and vector calculations.
- [icao_units](https://crates.io/crates/icao-units) - to define `Metres` and
`NauticalMiles` and perform conversions between them.

![Ellipsoid Class Diagram](docs/images/ellipsoid_class_diagram.svg)  
*Figure 3 Class Diagram*

The library is declared [no_std](https://docs.rust-embedded.org/book/intro/no-std.html)
so it can be used in embedded applications.

## Examples

### Calculate geodesic initial azimuths and length

Calculate the distance in Nautical Miles and azimuths (a.k.a bearing) in degrees
between two positions.

```rust
use icao_wgs84::*;

let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
let (azimuth, length, end_azimuth) = calculate_azimuths_and_geodesic_length(&istanbul, &washington, Radians(great_circle::MIN_VALUE), &WGS84_ELLIPSOID);

let distance_nm = NauticalMiles::from(length);
println!("Istanbul-Washington distance: {:?}", distance_nm);

let azimuth_degrees = Degrees::from(azimuth);
println!("Istanbul-Washington initial azimuth: {:?}", azimuth_degrees.0);

let end_azimuth_degrees = Degrees::from(end_azimuth.opposite());
println!("Washington-Istanbul initial azimuth: {:?}", end_azimuth_degrees.0);

```

### Calculate along track and across track distances

Create a `GeodesicSegment` between two positions and then calculate the
along track and across track distances of a third position relative to the `GeodesicSegment`.

The example is based on this reply from C. F. F. Karney :
<https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#8a93>.  
The expected latitude and longitude are from Karney's reply:

> Final result 54.92853149711691 -21.93729106604878

Note: the across track distance (xtd) is negative because Reyjavik is on the
right hand side of the `GeodesicSegment`.  
Across track distances are:

- positive for positions to the left of the `GeodesicSegment`,
- negative for positions to the right of the `GeodesicSegment`
- and zero for positions within the precision of the `GeodesicSegment`.

```rust
use icao_wgs84::*;
use angle_sc::is_within_tolerance;

let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
let g1 = GeodesicSegment::from(&istanbul, &washington);

let distance_nm = NauticalMiles::from(g1.length());
println!("Istanbul-Washington distance: {:?}", distance_nm);

let azimuth_degrees = Degrees::from(g1.azimuth(Metres(0.0)));
println!("Istanbul-Washington initial azimuth: {:?}", azimuth_degrees.0);

let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));

// Calculate geodesic segment along track and across track distances to 1mm precision.
let (atd, xtd, iterations) = g1.calculate_atd_and_xtd(&reyjavik, Metres(1e-3));
assert!(is_within_tolerance(3928788.572, atd.0, 1e-3));
assert!(is_within_tolerance(-1010585.9988368, xtd.0, 1e-3));
println!("calculate_atd_and_xtd iterations: {:?}", iterations);

let position = g1.lat_long(atd);
assert!(is_within_tolerance(
    54.92853149711691,
    Degrees::from(position.lat()).0,
    128.0 * f64::EPSILON
));
assert!(is_within_tolerance(
    -21.93729106604878,
    Degrees::from(position.lon()).0,
    2048.0 * f64::EPSILON
));
```

Also Note: the example uses 1mm precision to match Karney's result.  
In practice, precision should be determined from position accuracy.  
Higher precision requires more iterations and therefore takes longer to
calculate the result.

### Calculate geodesic intersection point

Create two `GeodesicSegment`s, each between two positions and then calculate the
distances from the geodesic segment start points to their intersection point.

The example is based on this reply from C. F. F. Karney :
<https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a>  
The expected latitude and longitude are from Karney's reply:

> Final result 54.7170296089477 -14.56385574430775

Note: Karney's solution requires all 4 positions to be in the same hemisphere
centered at the intersection point.  
This solution does **not** have that requirement.

```rust
use icao_wgs84::*;
use angle_sc::is_within_tolerance;

let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
let accra = LatLong::new(Degrees(6.0), Degrees(0.0));

let g1 = GeodesicSegment::from((&istanbul, &washington));
let g2 = GeodesicSegment::from((&reyjavik, &accra));

// Calculate the intersection point position
let result = calculate_intersection_point(&g1, &g2, Metres(1e-3));

// Get the intersection point position
let lat_lon = result.unwrap();
assert!(is_within_tolerance(54.7170296089477, lat_lon.lat().0, 1e-6));
assert!(is_within_tolerance(-14.56385574430775, lat_lon.lon().0, 1e-6));
```

## Test

The integration test uses Charles Karney's
[Test data for geodesics](https://geographiclib.sourceforge.io/C++/doc/geodesic.html#testgeod) to verify geodesic segment azimuth and distance calculations between
positions on the WGS-84 ellipsoid.
Run the tests using:

```
cargo test -- --ignored
```

## Contribution

If you want to contribute through code or documentation, the [Contributing](CONTRIBUTING.md) guide is the best place to start. If you have any questions, please feel free to ask.
Just please abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

## License

`icao-wgs84` is provided under a MIT license, see [LICENSE](LICENSE).
