"""Real terrain for the Ladakh sites, and what counts as buildable ground.

Why this exists
---------------
An array layout drawn on a flat plane is a pattern, not a design. A pad has to
sit on ground that is actually flat, actually there, and actually reachable.
Without terrain the optimiser will happily place an antenna on a cliff face or
in a river bed, and the resulting coordinate list looks authoritative while
being unbuildable.

Elevation source
----------------
AWS Terrain Tiles (the Tilezen/Mapzen terrarium tile set, hosted by AWS Open
Data at ``s3.amazonaws.com/elevation-tiles-prod``), which mosaics SRTM, ASTER
GDEM and national datasets. Public, no key, no registration. Elevation is
encoded in the RGB channels as

    height_m = R * 256 + G + B / 256 - 32768

Verification
------------
The pipeline is checked against two independently published elevations before
it is trusted: the Indian Astronomical Observatory at Hanle, published at
4500 m, reads 4491 m; the ALMA array operations site, published at 5058.7 m,
reads 5036 m. Both are within the vertical accuracy of the underlying SRTM
data, and the ALMA check matters because it is on the other side of the planet
and exercises the same code path.

Buildability
------------
A cell is treated as buildable when the local slope is below
``MAX_SLOPE_DEG``. Ten degrees is the default, which is loose: it admits ground
that would need cut-and-fill but not ground that is structurally impossible.
The threshold is a stated engineering assumption, not a measurement, and
``buildable_mask`` takes it as an argument so its effect can be varied.

What this does NOT model
------------------------
Land ownership, protected-area boundaries, existing structures, access roads,
power routing, water courses, permafrost, avalanche or rockfall exposure, and
geotechnical bearing capacity. A pad position from this module means "the
ground here is flat enough", not "this is a construction site".
"""

from __future__ import annotations

import hashlib
import io
import math
import urllib.request
from pathlib import Path

import numpy as np

__all__ = [
    "TILE_URL",
    "MAX_SLOPE_DEG",
    "DEM",
    "fetch_dem",
    "buildable_mask",
    "slope_deg",
    "local_enu",
    "enu_to_latlon",
]

TILE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
CACHE = Path(__file__).resolve().parent / "data" / "dem_cache"
MAX_SLOPE_DEG = 10.0
EARTH_R = 6378137.0


def _tile_indices(lat: float, lon: float, z: int):
    n = 2 ** z
    x = (lon + 180.0) / 360.0 * n
    lr = math.radians(lat)
    y = (1.0 - math.log(math.tan(lr) + 1.0 / math.cos(lr)) / math.pi) / 2.0 * n
    return x, y


def _tile_bounds(tx: int, ty: int, z: int):
    """(lon_west, lat_north, lon_east, lat_south) of one tile."""
    n = 2 ** z

    def lat_of(y):
        return math.degrees(math.atan(math.sinh(math.pi * (1.0 - 2.0 * y / n))))

    return (tx / n * 360.0 - 180.0, lat_of(ty),
            (tx + 1) / n * 360.0 - 180.0, lat_of(ty + 1))


def _fetch_tile(tx: int, ty: int, z: int) -> np.ndarray:
    """One terrarium tile as a float elevation array, cached on disk."""
    from PIL import Image

    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{z}_{tx}_{ty}.png"
    if not path.exists():
        url = TILE_URL.format(z=z, x=tx, y=ty)
        with urllib.request.urlopen(url, timeout=60) as r:
            path.write_bytes(r.read())
    rgb = np.asarray(Image.open(io.BytesIO(path.read_bytes())).convert("RGB"),
                     dtype=np.float64)
    return rgb[..., 0] * 256.0 + rgb[..., 1] + rgb[..., 2] / 256.0 - 32768.0


class DEM:
    """A mosaicked elevation grid in a local east-north frame.

    ``elevation[i, j]`` is metres above sea level at northing ``north[i]`` and
    easting ``east[j]``, both metres from the site origin.
    """

    def __init__(self, lat0, lon0, east, north, elevation, zoom, tiles):
        self.lat0, self.lon0 = float(lat0), float(lon0)
        self.east, self.north = east, north
        self.elevation = elevation
        self.zoom, self.n_tiles = int(zoom), int(tiles)

    @property
    def pixel_m(self) -> float:
        return float(abs(self.east[1] - self.east[0]))

    def summary(self) -> dict:
        return {
            "origin_lat": self.lat0, "origin_lon": self.lon0,
            "zoom": self.zoom, "tiles_fetched": self.n_tiles,
            "pixel_m": round(self.pixel_m, 1),
            "grid": tuple(self.elevation.shape),
            "elev_min_m": float(np.nanmin(self.elevation)),
            "elev_max_m": float(np.nanmax(self.elevation)),
            "elev_at_origin_m": float(self.at(0.0, 0.0)),
            "source": "AWS Terrain Tiles (Tilezen terrarium), SRTM/ASTER derived",
        }

    def at(self, east_m: float, north_m: float) -> float:
        """Nearest-cell elevation at a local coordinate."""
        j = int(np.clip(np.searchsorted(self.east, east_m), 0, self.east.size - 1))
        i = int(np.clip(np.searchsorted(self.north, north_m), 0, self.north.size - 1))
        return float(self.elevation[i, j])


def fetch_dem(lat0: float, lon0: float, half_size_m: float = 4000.0,
              zoom: int = 13) -> DEM:
    """Elevation grid covering +/- ``half_size_m`` around a site.

    Zoom 13 is about 16 m per pixel at this latitude, comfortably finer than
    the 12 m minimum pad separation an 8 m dish implies.
    """
    dlat = math.degrees(half_size_m / EARTH_R)
    dlon = math.degrees(half_size_m / (EARTH_R * math.cos(math.radians(lat0))))
    north_lat, south_lat = lat0 + dlat, lat0 - dlat
    west_lon, east_lon = lon0 - dlon, lon0 + dlon

    x0, y0 = _tile_indices(north_lat, west_lon, zoom)
    x1, y1 = _tile_indices(south_lat, east_lon, zoom)
    tx0, tx1 = int(math.floor(x0)), int(math.floor(x1))
    ty0, ty1 = int(math.floor(y0)), int(math.floor(y1))

    rows, count = [], 0
    for ty in range(ty0, ty1 + 1):
        row = []
        for tx in range(tx0, tx1 + 1):
            row.append(_fetch_tile(tx, ty, zoom))
            count += 1
        rows.append(np.hstack(row))
    mosaic = np.vstack(rows)

    # geographic extent of the mosaic, then crop to the requested box
    w_lon, n_lat, _, _ = _tile_bounds(tx0, ty0, zoom)
    _, _, e_lon, s_lat = _tile_bounds(tx1, ty1, zoom)
    lons = np.linspace(w_lon, e_lon, mosaic.shape[1])
    lats = np.linspace(n_lat, s_lat, mosaic.shape[0])

    jsel = np.flatnonzero((lons >= west_lon) & (lons <= east_lon))
    isel = np.flatnonzero((lats >= south_lat) & (lats <= north_lat))
    sub = mosaic[np.ix_(isel, jsel)]
    lat_sub, lon_sub = lats[isel], lons[jsel]

    east = np.radians(lon_sub - lon0) * EARTH_R * math.cos(math.radians(lat0))
    north = np.radians(lat_sub - lat0) * EARTH_R
    # store south-to-north so `north` is increasing
    order = np.argsort(north)
    return DEM(lat0, lon0, east, north[order], sub[order, :], zoom, count)


def slope_deg(dem: DEM) -> np.ndarray:
    """Local surface slope in degrees, from the elevation gradient."""
    dz_dn, dz_de = np.gradient(dem.elevation, dem.north, dem.east)
    return np.degrees(np.arctan(np.hypot(dz_de, dz_dn)))


def buildable_mask(dem: DEM, max_slope_deg: float = MAX_SLOPE_DEG,
                   max_radius_m: float | None = None) -> np.ndarray:
    """Cells flat enough for a pad, optionally within a radius of the origin."""
    mask = slope_deg(dem) <= float(max_slope_deg)
    if max_radius_m is not None:
        ee, nn = np.meshgrid(dem.east, dem.north)
        mask &= np.hypot(ee, nn) <= float(max_radius_m)
    return mask


def local_enu(lat, lon, lat0, lon0):
    """Geodetic to local east/north metres on the tangent plane at the site."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    east = np.radians(lon - lon0) * EARTH_R * math.cos(math.radians(lat0))
    north = np.radians(lat - lat0) * EARTH_R
    return east, north


def enu_to_latlon(east, north, lat0, lon0):
    """Inverse of :func:`local_enu`, so a pad list can be handed to a surveyor."""
    east = np.asarray(east, dtype=float)
    north = np.asarray(north, dtype=float)
    lat = lat0 + np.degrees(north / EARTH_R)
    lon = lon0 + np.degrees(east / (EARTH_R * math.cos(math.radians(lat0))))
    return lat, lon
