# External data sources

Real, published measurements used by this project. Nothing in this directory is
generated or modified; the files are redistributed as obtained.

## Antenna pad coordinates

CASA observatory configuration files, as distributed with the NRAO Common
Astronomy Software Applications package. Obtained 2026-09-18 from the Debian
CASA data mirror <https://github.com/kernsuite-debian/casalite>, path
`data/alma/simmos/`, which redistributes the official CASA `simmos` set.

| File | Array | Stations used | Coordinate system | Dish |
|---|---|---|---|---|
| `alma.all.cfg` | ALMA, complete pad list | 174 twelve-metre pads | local tangent plane (m) | 12 m |
| `aca.all.cfg` | Atacama Compact Array | 22 | local tangent plane (m) | 7 m |
| `vla.a.cfg` | VLA A configuration | 27 | geocentric ITRF (m) | 25 m |

`alma.all.cfg` also lists 18 seven-metre ACA pads and a `MASTER0` survey
reference marker; `thz_opt.arrays.real_arrays.load_cfg` filters both out when
asked for the 12 m array, because mixing dish sizes would apply one shadowing
limit to two different apertures.

### Verification against published values

These are the checks that the files are what they claim to be:

| Quantity | From these files | Published |
|---|---|---|
| ALMA maximum baseline | 16 195 m | ~16 km, ALMA most extended configuration |
| VLA A-configuration maximum baseline | 36 623 m | ~36 km |
| ALMA pad radial extent | 6 to 10 781 m | Chajnantor plateau, ~10 km from array centre |
| ALMA dish diameter | 12 m | 12 m |

The VLA figure is a particularly useful check because those coordinates are
geocentric and had to be rotated into a local frame; recovering the published
36 km baseline confirms the conversion.

## Site parameters

Recorded in `thz_opt.arrays.real_arrays.REAL_SITES`.

| Site | Latitude | Longitude | Elevation |
|---|---|---|---|
| Llano de Chajnantor (ALMA, ACA) | -23.0229 deg | -67.7551 deg | 5058.7 m |
| Plains of San Agustin (VLA) | +34.0784 deg | -107.6184 deg | 2124 m |

## Atmospheric measurements

Tabulated in `src/thz_opt/constraints/atmosphere_data.py` with their sources.

| Quantity | Value | Source |
|---|---|---|
| Path RMS at 1 km, 120 s, no WVR correction | 200 um | ALMA Memo 624 |
| Path RMS at 1 km, below median PWV (1.24 mm) | 115 um | ALMA Memo 624 |
| Path RMS at 1 km, WVR corrected, wind < 10 m/s | 70 um | ALMA Memo 624 |
| Path RMS at 1 km, WVR corrected, wind > 10 m/s | 140 um | ALMA Memo 624 |
| Structure-function exponents, no WVR | 0.65 (b < 1 km), 0.22 (b > 1 km) | Matsushita et al. 2017, via Memo 624 |
| Structure-function exponents, WVR corrected | 0.60, 0.29 | Matsushita et al. 2017, via Memo 624 |
| Monthly PWV percentiles (10/25/50/75) | table | ALMA Memo 624 Table 3 |
| Year-round median PWV | ~1.0 mm | Cortes et al. 2020 |
| January / August median PWV | 2.56 / 0.72 mm | Cortes et al. 2020 |

- L. T. Maud et al., *Updates to ALMA site properties: using the ESO-Allegro
  phase RMS database*, ALMA Memo 624, 2023, <https://arxiv.org/abs/2304.08318>
  (over 17 000 observations, 2015-2019).
- F. Cortes et al., *Twenty years of precipitable water vapor measurements in
  the Chajnantor area*, A&A 640, A126, 2020,
  <https://doi.org/10.1051/0004-6361/202037784>.

### Verification

The implementation is checked against the memo's own published scaling factors
between its summary baselines:

| Scaling | Memo | This implementation |
|---|---|---|
| 500 to 1000 m | 1.51 | 1.516 |
| 1000 to 5000 m | 1.59 | 1.595 |
| 5000 to 10 000 m | 1.22 | 1.223 |

## What is still not measured

The phase values are medians over thousands of observations at a 120 s
timescale, so they characterise the site rather than any particular night, and
the source memo notes its sample omits the very worst conditions, in which no
observation was attempted. Elevation is not normalised: values are as measured
across a 21-88 degree range with a median of 58 degrees.
