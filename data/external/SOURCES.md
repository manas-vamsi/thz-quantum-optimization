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

## What is still not measured

Atmospheric parameters remain assumed rather than site-measured:

- **Precipitable water vapour.** Typical published values for a high, dry site
  are used, not a measured time series for a specific location and season.
- **Phase structure function.** The exponents follow Kolmogorov theory and the
  scale constant is a typical value from the millimetre-interferometry
  literature; the residual factor after phase referencing is a free parameter.

These are the remaining synthetic inputs, and the manuscript says so.
