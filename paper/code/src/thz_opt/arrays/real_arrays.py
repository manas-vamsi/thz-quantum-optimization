"""Real, published antenna pad positions.

Every layout elsewhere in this package is generated from a formula. This module
reads *measured* pad coordinates for existing arrays, so the selection problem
can be posed on real geometry instead of synthetic points.

The primary dataset is the complete ALMA pad list -- 174 twelve-metre pads on the
Chajnantor plateau at 5000 m, from which each ALMA configuration cycle selects a
subset. That is precisely the problem this project formulates: choose ``N``
antennas from ``M`` prepared pads. Using it removes the largest synthetic
assumption in the study.

Format
------
CASA observatory configuration files (``.cfg``): comment lines beginning ``#``,
then whitespace-separated ``x y z diameter padname``. The header declares the
coordinate system:

``coordsys=LOC``
    local tangent plane in metres, east/north/up. ALMA and ACA use this, so the
    ``(x, y)`` columns are directly the coordinates this package expects.

``coordsys=XYZ``
    geocentric ITRF metres. Converted to a local east/north frame about the
    array centroid by :func:`itrf_to_local`.

Provenance is recorded in ``data/external/SOURCES.md``; the files are
redistributed with CASA under its licence and are not modified here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["ArrayConfig", "load_cfg", "itrf_to_local", "available_arrays", "REAL_SITES"]

# Published site parameters for the arrays shipped with this project.
# Latitudes are geodetic; they set the Earth-rotation geometry.
REAL_SITES = {
    "alma": {
        "name": "ALMA / Atacama Large Millimeter-submillimeter Array",
        "site": "Llano de Chajnantor, Chile",
        "latitude_deg": -23.0229,
        "longitude_deg": -67.7551,
        "elevation_m": 5058.7,
        "dish_diameter_m": 12.0,
        "n_pads": 174,
        "reference": "CASA observatory configuration alma.all.cfg",
    },
    "aca": {
        "name": "ACA / Atacama Compact Array (Morita Array)",
        "site": "Llano de Chajnantor, Chile",
        "latitude_deg": -23.0229,
        "longitude_deg": -67.7551,
        "elevation_m": 5058.7,
        "dish_diameter_m": 7.0,
        "n_pads": 22,
        "reference": "CASA observatory configuration aca.all.cfg",
    },
    "vla": {
        "name": "Karl G. Jansky Very Large Array",
        "site": "Plains of San Agustin, New Mexico, USA",
        "latitude_deg": 34.0784,
        "longitude_deg": -107.6184,
        "elevation_m": 2124.0,
        "dish_diameter_m": 25.0,
        "n_pads": 27,
        "reference": "CASA observatory configuration vla.a.cfg",
    },
}


@dataclass(frozen=True)
class ArrayConfig:
    """A real array configuration: pad coordinates, dish sizes and names."""

    xy: np.ndarray            # (n, 2) local east/north metres
    z: np.ndarray             # (n,) local up metres, kept for reference
    diameters: np.ndarray     # (n,) metres
    names: list               # pad identifiers as published
    coordsys: str
    source: str

    def __len__(self) -> int:
        return int(self.xy.shape[0])

    def summary(self) -> dict:
        r = np.hypot(self.xy[:, 0], self.xy[:, 1])
        d = np.hypot(self.xy[:, None, 0] - self.xy[None, :, 0],
                     self.xy[:, None, 1] - self.xy[None, :, 1])
        off = d[~np.eye(len(self), dtype=bool)]
        return {
            "source": self.source,
            "n_pads": len(self),
            "coordsys": self.coordsys,
            "dish_diameter_m": float(np.median(self.diameters)),
            "radius_min_m": float(r.min()),
            "radius_max_m": float(r.max()),
            "separation_min_m": float(off.min()),
            "baseline_max_m": float(off.max()),
        }


def itrf_to_local(xyz: np.ndarray) -> np.ndarray:
    """Geocentric ITRF metres to a local east/north/up frame at the centroid.

    The reference point is the centroid of the supplied stations, so the result
    is a local tangent plane for that array -- adequate for array design, where
    only relative positions matter.
    """
    p = np.asarray(xyz, dtype=float).reshape(-1, 3)
    c = p.mean(axis=0)
    lon = np.arctan2(c[1], c[0])
    lat = np.arctan2(c[2], np.hypot(c[0], c[1]))

    sl, cl = np.sin(lon), np.cos(lon)
    sb, cb = np.sin(lat), np.cos(lat)
    rot = np.array([
        [-sl, cl, 0.0],                 # east
        [-sb * cl, -sb * sl, cb],       # north
        [cb * cl, cb * sl, sb],         # up
    ])
    return (p - c) @ rot.T


def load_cfg(path: str | Path, dish_diameter_m: float | None = None,
             drop_markers: bool = True) -> ArrayConfig:
    """Read a CASA ``.cfg`` observatory configuration file.

    ``dish_diameter_m`` keeps only pads for that dish size. ``alma.all.cfg``
    lists both the 175 twelve-metre pads and the 18 seven-metre ACA pads, and
    mixing them would silently model an array of two different dish sizes with
    one shadowing limit; pass 12.0 for the main array.

    ``drop_markers`` removes non-pad rows -- ``alma.all.cfg`` ends with a
    ``MASTER0`` survey reference point, which is not a station.
    """
    p = Path(path)
    text = p.read_text(encoding="utf-8", errors="replace")

    coordsys = "LOC"
    rows, names = [], []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if "coordsys" in stripped.lower():
                coordsys = stripped.split("=", 1)[1].strip().split()[0].upper()
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            vals = [float(v) for v in parts[:4]]
        except ValueError:
            continue
        rows.append(vals)
        names.append(parts[4] if len(parts) > 4 else f"P{len(names):03d}")

    if not rows:
        raise ValueError(f"no station rows parsed from {p}")

    arr = np.asarray(rows, dtype=float)
    keep = np.ones(arr.shape[0], dtype=bool)
    if drop_markers:
        keep &= np.array([not n.upper().startswith("MASTER") for n in names])
    if dish_diameter_m is not None:
        keep &= np.isclose(arr[:, 3], float(dish_diameter_m))
    arr = arr[keep]
    names = [n for n, k in zip(names, keep) if k]
    if arr.shape[0] == 0:
        raise ValueError(f"no stations left in {p} after filtering")

    xyz, diam = arr[:, :3], arr[:, 3]

    if coordsys.startswith("XYZ") or coordsys.startswith("GEO"):
        local = itrf_to_local(xyz)
    else:
        local = xyz

    return ArrayConfig(xy=local[:, :2].copy(), z=local[:, 2].copy(),
                       diameters=diam, names=names, coordsys=coordsys,
                       source=p.name)


def available_arrays(directory: str | Path) -> dict:
    """Map array key to its configuration file in ``directory``."""
    d = Path(directory)
    found = {}
    for key, pattern in (("alma", "alma.all.cfg"), ("aca", "aca.all.cfg"),
                         ("vla", "vla.a.cfg")):
        f = d / pattern
        if f.exists():
            found[key] = f
    return found
