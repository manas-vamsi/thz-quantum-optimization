"""Sky brightness models on the image grid conjugate to a UV grid.

The image grid is fixed by the UV grid, not chosen: a ``UVGrid`` with
``n_cells`` cells of size ``cell_size`` wavelengths is the discrete Fourier
partner of an ``n_cells`` square image with pixels of
``1 / (n_cells * cell_size)`` radians. Building the sky on any other grid would
require interpolation in the transform and would blur the comparison the whole
package exists to make.

Three models are provided, in increasing order of how much they can embarrass
an array:

* **point sources** -- the classic test. Every array can find a single point;
  differences show up in the sidelobes around it and in the ability to separate
  a close pair.
* **a Gaussian** -- smooth and extended, so it punishes missing short spacings.
  An array with no short baselines resolves it out and reports too little flux.
* **a disc with an annular gap** -- the HL Tau-like model of
  :mod:`thz_opt.science.disk_gap`, which needs both short and long spacings at
  once and is the case the science-derived objective was built for.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "SkyModel",
    "pixel_scale_arcsec",
    "image_coordinates",
    "point_sources",
    "gaussian_source",
    "disc_with_gap_image",
]

ARCSEC = np.pi / (180.0 * 3600.0)


@dataclass(frozen=True)
class SkyModel:
    """A brightness image and the grid it lives on."""

    image: np.ndarray
    pixel_rad: float
    name: str = "sky"

    @property
    def n(self) -> int:
        return int(self.image.shape[0])

    @property
    def pixel_arcsec(self) -> float:
        return float(self.pixel_rad / ARCSEC)

    @property
    def total_flux(self) -> float:
        return float(self.image.sum())

    def summary(self) -> dict:
        return {
            "name": self.name,
            "grid": self.n,
            "pixel_arcsec": round(self.pixel_arcsec, 5),
            "field_of_view_arcsec": round(self.n * self.pixel_arcsec, 3),
            "total_flux": round(self.total_flux, 6),
            "peak": float(self.image.max()),
        }


def pixel_scale_arcsec(grid) -> float:
    """Image pixel size implied by a UV grid, in arcseconds."""
    return float(np.rad2deg(1.0 / (grid.n_cells * grid.cell_size)) * 3600.0)


def image_coordinates(grid):
    """``(l, m)`` in radians for every pixel, matching the FFT centring."""
    n = grid.n_cells
    pix = 1.0 / (n * grid.cell_size)
    axis = (np.arange(n) - n // 2) * pix
    return np.meshgrid(axis, axis, indexing="ij"), pix


def point_sources(grid, offsets_arcsec=((0.0, 0.0),), fluxes=(1.0,),
                  name="point sources") -> SkyModel:
    """Delta functions at given offsets from the field centre.

    Offsets are snapped to the nearest pixel, which is exact for the transform
    and is what a gridded comparison can represent.
    """
    (ll, mm), pix = image_coordinates(grid)
    n = grid.n_cells
    img = np.zeros((n, n), dtype=float)
    fluxes = np.atleast_1d(np.asarray(fluxes, dtype=float))
    for (dl, dm), f in zip(offsets_arcsec, fluxes):
        i = int(round(dl * ARCSEC / pix)) + n // 2
        j = int(round(dm * ARCSEC / pix)) + n // 2
        if not (0 <= i < n and 0 <= j < n):
            raise ValueError(f"source at ({dl}, {dm}) arcsec falls outside the field")
        img[i, j] += float(f)
    return SkyModel(img, pix, name)


def gaussian_source(grid, fwhm_arcsec: float, flux: float = 1.0,
                    name=None) -> SkyModel:
    """A circular Gaussian, the simplest source that needs short spacings."""
    (ll, mm), pix = image_coordinates(grid)
    sigma = fwhm_arcsec * ARCSEC / 2.3548
    img = np.exp(-0.5 * (ll ** 2 + mm ** 2) / sigma ** 2)
    img *= flux / img.sum()
    return SkyModel(img, pix, name or f"Gaussian, {fwhm_arcsec} arcsec FWHM")


def disc_with_gap_image(grid, model=None, flux: float = 1.0) -> SkyModel:
    """The protoplanetary disc of :mod:`thz_opt.science.disk_gap`, as an image.

    Uses the same radial profile the Fisher weighting was derived from, so the
    science case and the imaging test describe one source rather than two.
    """
    from ..science.disk_gap import HL_TAU, brightness_profile

    model = model or HL_TAU
    (ll, mm), pix = image_coordinates(grid)
    r_rad = np.hypot(ll, mm)
    r_au = r_rad / ARCSEC * model.distance_pc

    img = np.zeros_like(r_au)
    inside = (r_au >= model.r_in_au) & (r_au <= model.r_out_au)
    if inside.any():
        img[inside] = brightness_profile(r_au[inside], model)
    total = img.sum()
    if total > 0:
        img *= flux / total
    return SkyModel(img, pix, f"disc with gap ({model.name})")
