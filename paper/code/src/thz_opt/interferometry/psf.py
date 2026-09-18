"""Dirty beam (point-spread function) from a sampled UV grid.

Definition used here
--------------------
Let ``S(u, v)`` be the sampling function on the square UV grid of
:class:`thz_opt.interferometry.uv.UVGrid`.  Two weightings are supported:

``"natural"``   S = number of samples in the cell (dense regions dominate)
``"uniform"``   S = 1 where the cell is occupied, 0 elsewhere

The dirty beam is the inverse discrete Fourier transform

    PSF(l, m) = F^-1[ S(u, v) ]

evaluated with :func:`numpy.fft.ifft2` and centred with ``ifftshift`` /
``fftshift`` so that zero spatial frequency sits at the array centre and the
beam peak lands at the centre of the image.  The result is real up to numerical
noise when the sampling is Hermitian (which it is whenever conjugate points are
included); the imaginary part is discarded and its magnitude is reported so the
assumption can be checked rather than assumed.

The PSF is normalised to unit peak.  Pixel scale: with ``n_cells`` cells of
size ``du`` wavelengths, the image spans ``1 / du`` radians with a pixel of
``1 / (n_cells * du)`` radians.

Assumptions and limitations
---------------------------
* Flat-sky, ``w = 0`` (no w-projection); valid for small fields and short
  baselines relative to the array size.
* The grid is a hard binning of the samples -- no convolutional gridding
  kernel, no tapering.  Absolute sidelobe numbers therefore depend on the cell
  size and must only be compared between configurations sharing one grid.
"""

from __future__ import annotations

import numpy as np

from .uv import UVGrid, uv_occupancy


def sampling_function(uv: np.ndarray, grid: UVGrid, weighting: str = "uniform") -> np.ndarray:
    """Sampling function ``S(u, v)`` on the grid, indexed ``[iu, iv]``."""
    occ, _ = uv_occupancy(uv, grid)
    if weighting == "natural":
        return occ.astype(float)
    if weighting == "uniform":
        return (occ > 0).astype(float)
    raise ValueError("weighting must be 'natural' or 'uniform'")


def dirty_beam(uv: np.ndarray, grid: UVGrid, weighting: str = "uniform"):
    """Return ``(psf, info)``.

    ``psf`` is the real part of the centred inverse FFT of the sampling
    function, normalised to a peak of 1 and oriented so that axis 0 is ``l``
    (conjugate to ``u``) and axis 1 is ``m`` (conjugate to ``v``).  ``info``
    records the imaginary residual and the pixel scale.
    """
    s = sampling_function(uv, grid, weighting)
    if s.sum() == 0:
        raise ValueError("sampling function is empty -- no samples fell on the grid")

    img = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(s)))
    peak = np.abs(img).max()
    imag_frac = float(np.abs(img.imag).max() / peak)
    psf = img.real / img.real.max()

    pixel_rad = 1.0 / (grid.n_cells * grid.cell_size)
    info = {
        "weighting": weighting,
        "imag_residual_fraction": imag_frac,
        "pixel_scale_rad": pixel_rad,
        "pixel_scale_arcsec": float(np.rad2deg(pixel_rad) * 3600.0),
        "field_of_view_rad": pixel_rad * grid.n_cells,
        "n_occupied_cells": int((s > 0).sum()),
    }
    return psf, info
