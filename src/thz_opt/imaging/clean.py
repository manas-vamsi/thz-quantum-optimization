"""Hogbom CLEAN, and restoration with a fitted Gaussian beam.

CLEAN is what makes the comparison meaningful. A dirty image is the sky
convolved with the array's point-spread function, so ranking arrays by dirty
images mostly ranks their sidelobes. Deconvolution is the step where an array
with poor UV coverage actually fails: the algorithm cannot recover structure
the array never sampled, and the error shows up in the residual and in the
restored image rather than as a cosmetic artefact.

Hogbom (1974) is the original and the simplest: repeatedly find the brightest
pixel in the residual, record a fraction of it as a component, and subtract a
scaled point-spread function centred there. It is chosen here over Clark or
Cotton-Schwab because it is transparent, has one tuning parameter that matters,
and is sufficient for comparing arrays observed identically -- no claim is made
that it is the best deconvolver available.

The restoring beam is a Gaussian fitted to the main lobe of the dirty beam,
which is standard practice: it discards the sidelobe structure that CLEAN has
already accounted for, leaving a resolution element the components can be
convolved with.
"""

from __future__ import annotations

import numpy as np

__all__ = ["hogbom_clean", "fit_restoring_beam", "restore"]


def hogbom_clean(dirty: np.ndarray, psf: np.ndarray, gain: float = 0.1,
                 n_iter: int = 2000, threshold: float | None = None):
    """Deconvolve ``dirty`` by ``psf``.

    Returns ``(components, residual, info)``. ``threshold`` defaults to three
    times the robust noise estimate of the residual, so the loop stops when it
    is subtracting noise rather than signal.
    """
    res = np.array(dirty, dtype=float, copy=True)
    comps = np.zeros_like(res)
    n = res.shape[0]
    c = n // 2
    if psf.shape != dirty.shape:
        raise ValueError("psf and dirty image must share a grid")

    if threshold is None:
        mad = np.median(np.abs(res - np.median(res)))
        threshold = 3.0 * 1.4826 * mad if mad > 0 else 1e-12

    used = 0
    for used in range(1, n_iter + 1):
        idx = int(np.argmax(np.abs(res)))
        i, j = divmod(idx, n)
        peak = res[i, j]
        if abs(peak) < threshold:
            used -= 1
            break
        amp = gain * peak
        comps[i, j] += amp
        # subtract the psf centred on (i, j), with edge clipping
        di, dj = i - c, j - c
        si = slice(max(0, di), min(n, n + di))
        sj = slice(max(0, dj), min(n, n + dj))
        pi = slice(max(0, -di), min(n, n - di))
        pj = slice(max(0, -dj), min(n, n - dj))
        res[si, sj] -= amp * psf[pi, pj]

    info = {
        "iterations": int(used),
        "gain": gain,
        "threshold": float(threshold),
        "converged": used < n_iter,
        "residual_rms": float(np.sqrt(np.mean(res ** 2))),
        "component_flux": float(comps.sum()),
    }
    return comps, res, info


def fit_restoring_beam(psf: np.ndarray) -> tuple[np.ndarray, dict]:
    """A circular Gaussian matched to the main lobe of the dirty beam.

    The width is taken from the half-power radius of the main lobe, measured
    outward from the peak until the beam first falls below one half. That is
    robust to the sidelobe structure a least-squares fit would be pulled by.
    """
    n = psf.shape[0]
    c = n // 2
    profile = psf[c, c:]
    below = np.flatnonzero(profile < 0.5)
    half_width_px = float(below[0]) if below.size else 1.0
    sigma = max(half_width_px / 1.1774, 0.5)      # half-width at half max

    axis = np.arange(n) - c
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    beam = np.exp(-0.5 * (xx ** 2 + yy ** 2) / sigma ** 2)
    return beam, {"sigma_px": sigma, "fwhm_px": 2.3548 * sigma}


def restore(components: np.ndarray, residual: np.ndarray, psf: np.ndarray):
    """Convolve the components with a clean beam and add the residual back.

    Adding the residual is not cosmetic: it is what keeps flux that CLEAN never
    converged on in the map, so the restored image can still be compared to the
    truth on total flux.
    """
    beam, binfo = fit_restoring_beam(psf)
    ft = np.fft.fft2(np.fft.ifftshift(beam))
    conv = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(
        np.fft.ifftshift(components)) * ft)).real
    # unit-peak restoring beam conserves the component flux
    conv /= beam.sum()
    return conv + residual, binfo
