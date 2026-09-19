"""Predict visibilities, sample them on an array's coverage, make a dirty image.

Conventions
-----------
The forward transform is the exact inverse of the one
:func:`thz_opt.interferometry.psf.dirty_beam` applies, so that sampling every
cell returns the input sky to machine precision. That round trip is the first
test in the suite, because every fidelity number downstream is meaningless if
the transform pair is inconsistent.

    V(u,v) = fftshift( fft2( ifftshift( I(l,m) ) ) )
    I(l,m) = fftshift( ifft2( ifftshift( V(u,v) ) ) )

What this models and what it does not
-------------------------------------
Modelled: the array's true UV sampling (including Earth rotation and station
heights), gridded onto cells; per-visibility thermal noise; the resulting dirty
image and dirty beam.

Not modelled: the primary beam, bandwidth and time smearing, non-coplanar *w*
projection during imaging, calibration and pointing errors, and atmospheric
phase noise as a time series. Thermal noise here is white and Gaussian on the
gridded visibilities, which is the standard idealisation and is optimistic.
"""

from __future__ import annotations

import numpy as np

from ..interferometry.uv import uv_occupancy

__all__ = ["predict_visibilities", "thermal_noise_sigma", "observe"]


def predict_visibilities(sky) -> np.ndarray:
    """``V(u, v)`` of a sky image, on the grid conjugate to it."""
    img = np.asarray(sky.image if hasattr(sky, "image") else sky, dtype=float)
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))


def thermal_noise_sigma(sky, snr: float, occupancy: np.ndarray) -> float:
    """Per-cell noise giving a requested peak signal-to-noise in the dirty map.

    Defining noise by the achieved image SNR rather than by a system
    temperature keeps the comparison between arrays fair: every array is given
    the same sensitivity per visibility and the same total observing effort, so
    differences in the reconstruction come from geometry alone.
    """
    n_used = int(np.count_nonzero(occupancy))
    if n_used == 0 or snr <= 0:
        return 0.0
    peak = float(np.abs(np.asarray(sky.image)).max())
    # the dirty-image noise averages over the occupied cells
    return float(peak * np.sqrt(n_used) / (snr * occupancy.size))


def observe(sky, uv, grid, snr: float | None = None, seed: int = 0,
            weighting: str = "uniform"):
    """Observe ``sky`` with the UV coverage ``uv`` and return the dirty image.

    Returns ``(dirty, psf, info)``. ``dirty`` is in the same brightness units
    as the input sky, scaled so that a fully sampled grid reproduces it.
    """
    occ, n_outside = uv_occupancy(np.asarray(uv).reshape(-1, 2), grid)
    if not occ.any():
        raise ValueError("this array samples no cell of the grid")

    if weighting == "uniform":
        s = (occ > 0).astype(float)
    elif weighting == "natural":
        s = occ.astype(float)
    else:
        raise ValueError("weighting must be 'uniform' or 'natural'")

    v_true = predict_visibilities(sky)
    v_obs = v_true * s

    sigma = 0.0
    if snr:
        sigma = thermal_noise_sigma(sky, snr, s)
        rng = np.random.default_rng(seed)
        mask = s > 0
        noise = (rng.normal(0.0, sigma, size=s.shape)
                 + 1j * rng.normal(0.0, sigma, size=s.shape))
        # a real sky has Hermitian visibilities; symmetrise so the dirty image
        # stays real rather than acquiring a spurious imaginary part
        noise = 0.5 * (noise + np.conj(noise[::-1, ::-1]))
        v_obs = v_obs + noise * mask

    dirty = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(v_obs))).real
    beam = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(s))).real
    peak = beam.max()
    if peak <= 0:
        raise ValueError("degenerate sampling function")

    info = {
        "occupied_cells": int(np.count_nonzero(s)),
        "fill_fraction": float(np.count_nonzero(s) / s.size),
        "samples_outside_grid": int(n_outside),
        "noise_sigma": sigma,
        "weighting": weighting,
        "beam_peak_before_normalisation": float(peak),
    }
    # normalise so the beam peaks at one and the map is in sky units
    return dirty / peak, beam / peak, info
