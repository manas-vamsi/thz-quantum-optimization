"""Earth-rotation aperture synthesis.

Adopted convention
------------------
This module uses the standard interferometry convention of Thompson, Moran &
Swenson, *Interferometry and Synthesis in Radio Astronomy*.  It is stated here
explicitly because the sign conventions differ between texts; the unit tests in
``tests/test_uv.py`` pin the behaviour so that a future change of convention
cannot pass silently.

Step 1 -- local ENU baseline ``(e, n, h)`` (East, North, Up, metres) at site
latitude ``phi`` to equatorial components::

    X = -n * sin(phi) + h * cos(phi)     (towards H = 0 on the equator)
    Y =  e                               (towards H = -6h on the equator)
    Z =  n * cos(phi) + h * sin(phi)     (towards the celestial pole)

Step 2 -- projection towards a source at hour angle ``H`` and declination
``dec``::

    u =           sin(H) * X +           cos(H) * Y
    v = -sin(dec)*cos(H) * X + sin(dec)*sin(H) * Y + cos(dec) * Z
    w =  cos(dec)*cos(H) * X - cos(dec)*sin(H) * Y + sin(dec) * Z

and ``(u, v, w)`` are divided by the wavelength to obtain wavelengths.

Consistency with the snapshot model
-----------------------------------
For a coplanar array (``h = 0``) observed at the zenith (``dec = phi``,
``H = 0``) this reduces exactly to ``u = e / lambda``, ``v = n / lambda``, which
is the snapshot model of :mod:`thz_opt.interferometry.uv`.  That identity is
asserted in the tests.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .baselines import compute_baselines


@dataclass(frozen=True)
class ObservationConfig:
    """Parameters of a tracked observation.

    Angles are in degrees on input (readable in YAML/JSON) and converted
    internally.  ``hour_angle_start/end`` are in hours; one hour of hour angle
    is 15 degrees of Earth rotation.
    """

    frequency_hz: float = 3.0e11          # 300 GHz
    declination_deg: float = -30.0
    latitude_deg: float = -23.0           # illustrative high-site latitude
    hour_angle_start_h: float = -2.0
    hour_angle_end_h: float = 2.0
    n_times: int = 41

    @property
    def wavelength(self) -> float:
        from .uv import wavelength_from_frequency

        return wavelength_from_frequency(self.frequency_hz)

    @property
    def hour_angles_rad(self) -> np.ndarray:
        return np.deg2rad(
            15.0 * np.linspace(self.hour_angle_start_h, self.hour_angle_end_h, self.n_times)
        )

    def as_dict(self) -> dict:
        d = asdict(self)
        d["wavelength_m"] = self.wavelength
        return d


def enu_to_xyz(b_enu: np.ndarray, latitude_rad: float) -> np.ndarray:
    """Local ENU baselines ``(n, 3)`` -> equatorial ``(X, Y, Z)`` in metres."""
    b = np.asarray(b_enu, dtype=float).reshape(-1, 3)
    e, n, h = b[:, 0], b[:, 1], b[:, 2]
    sp, cp = np.sin(latitude_rad), np.cos(latitude_rad)
    return np.column_stack((-n * sp + h * cp, e, n * cp + h * sp))


def project_uvw(
    b_xyz: np.ndarray, hour_angles_rad: np.ndarray, declination_rad: float, wavelength: float
) -> np.ndarray:
    """Project equatorial baselines onto ``(u, v, w)`` in wavelengths.

    Returns an array of shape ``(n_baselines * n_times, 3)``, ordered baseline-
    major (all times of baseline 0, then baseline 1, ...).
    """
    b = np.asarray(b_xyz, dtype=float).reshape(-1, 3)
    H = np.atleast_1d(np.asarray(hour_angles_rad, dtype=float))
    sd, cd = np.sin(declination_rad), np.cos(declination_rad)
    sH, cH = np.sin(H)[None, :], np.cos(H)[None, :]
    X, Y, Z = b[:, 0:1], b[:, 1:2], b[:, 2:3]

    u = sH * X + cH * Y
    v = -sd * cH * X + sd * sH * Y + cd * Z
    w = cd * cH * X - cd * sH * Y + sd * Z
    return np.stack((u, v, w), axis=-1).reshape(-1, 3) / wavelength


def layout_to_uv_tracks(
    xy: np.ndarray, config: ObservationConfig, include_conjugate: bool = True,
    heights: np.ndarray | None = None
) -> np.ndarray:
    """Earth-rotation ``(u, v)`` tracks of a layout, in wavelengths.

    ``heights`` gives each station's elevation in metres. Passing it makes the
    projection fully three-dimensional: height differences enter the equatorial
    baseline through :func:`enu_to_xyz`, which already carries the ``Up`` term,
    and so change ``(u, v)`` as well as ``w``.

    Omitting it keeps the coplanar approximation ``Up = 0``, which is what the
    controlled studies in this repository use so that layouts differ only in
    their horizontal geometry. The approximation is safe while height
    differences are small against baseline length, and stops being safe when
    they are not -- a real site can put 200 m of relief across a 3 km array.
    """
    b2d, i_idx, j_idx = compute_baselines(xy)
    if heights is None:
        up = np.zeros(len(b2d))
    else:
        h = np.asarray(heights, dtype=float).ravel()
        if h.shape[0] != np.asarray(xy).shape[0]:
            raise ValueError("heights must have one entry per station")
        up = h[j_idx] - h[i_idx]
    b_enu = np.column_stack((b2d, up))
    b_xyz = enu_to_xyz(b_enu, np.deg2rad(config.latitude_deg))
    uvw = project_uvw(
        b_xyz, config.hour_angles_rad, np.deg2rad(config.declination_deg), config.wavelength
    )
    uv = uvw[:, :2]
    if include_conjugate:
        uv = np.vstack((uv, -uv))
    return uv


def elevation_deg(config: ObservationConfig) -> np.ndarray:
    """Source elevation at each sampled hour angle, in degrees.

    Provided so an experiment can check that the requested hour-angle range
    keeps the source above the horizon; no elevation cut is applied silently.
    """
    phi = np.deg2rad(config.latitude_deg)
    dec = np.deg2rad(config.declination_deg)
    H = config.hour_angles_rad
    sin_el = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(H)
    return np.rad2deg(np.arcsin(np.clip(sin_el, -1.0, 1.0)))
