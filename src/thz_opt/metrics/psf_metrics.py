"""Scalar diagnostics of a dirty beam.

None of these is a complete description of imaging quality; they are the
quantities most often quoted, computed with explicitly stated definitions.

* **Peak sidelobe level** -- the largest absolute PSF value outside the main
  lobe, where the main lobe is the region around the peak bounded by the first
  minimum of the azimuthally averaged profile.
* **Sidelobe RMS** -- root-mean-square of the PSF outside the main lobe.
* **FWHM** -- full width at half maximum of the azimuthally averaged profile,
  in pixels and (if a pixel scale is supplied) in arcseconds.  For a strongly
  elliptical beam an azimuthal average is a poor summary, so the ellipticity is
  reported alongside it.
"""

from __future__ import annotations

import numpy as np


def _radial_average(psf: np.ndarray):
    n0, n1 = psf.shape
    y = np.arange(n0) - n0 // 2
    x = np.arange(n1) - n1 // 2
    ir = np.hypot(y[:, None], x[None, :]).astype(int)
    nbin = ir.max() + 1
    total = np.bincount(ir.ravel(), weights=psf.ravel(), minlength=nbin)
    count = np.bincount(ir.ravel(), minlength=nbin)
    return np.arange(nbin), total / np.maximum(count, 1)


def main_lobe_radius(psf: np.ndarray) -> int:
    """Radius in pixels of the first minimum of the azimuthally averaged PSF."""
    r, prof = _radial_average(psf)
    rising = np.flatnonzero(np.diff(prof) > 0)
    return int(r[rising[0]]) if rising.size else int(r[-1])


def fwhm_pixels(psf: np.ndarray) -> float:
    """FWHM of the azimuthally averaged profile, in pixels (linear interp)."""
    r, prof = _radial_average(psf)
    half = 0.5 * prof[0]
    below = np.flatnonzero(prof < half)
    if below.size == 0:
        return float("nan")
    k = int(below[0])
    if k == 0:
        return 0.0
    r0, r1 = float(r[k - 1]), float(r[k])
    p0, p1 = float(prof[k - 1]), float(prof[k])
    return 2.0 * (r0 + (p0 - half) * (r1 - r0) / (p0 - p1))


def beam_ellipticity(psf: np.ndarray, level: float = 0.5) -> float:
    """``1 - b/a`` of the pixel set above ``level * peak``, from second moments.

    0 is circular; values approaching 1 mean a strongly elongated beam.
    """
    m = psf >= level * psf.max()
    if m.sum() < 3:
        return 0.0
    ys, xs = np.nonzero(m)
    w = psf[m]
    ys = ys - np.average(ys, weights=w)
    xs = xs - np.average(xs, weights=w)
    cov = np.cov(np.vstack((xs, ys)), aweights=w)
    ev = np.sqrt(np.maximum(np.linalg.eigvalsh(cov), 0.0))
    a, b = float(ev[1]), float(ev[0])
    return 1.0 - b / a if a > 0 else 0.0


def psf_metrics(psf: np.ndarray, pixel_scale_arcsec: float | None = None) -> dict:
    """Scalar PSF diagnostics; see the module docstring for definitions."""
    p = np.asarray(psf, dtype=float)
    n0, n1 = p.shape
    y = np.arange(n0) - n0 // 2
    x = np.arange(n1) - n1 // 2
    r = np.hypot(y[:, None], x[None, :])

    r_main = main_lobe_radius(p)
    side = p[r > r_main]
    fwhm_px = fwhm_pixels(p)

    out = {
        "peak": float(p.max()),
        "main_lobe_radius_px": int(r_main),
        "peak_sidelobe_level": float(np.abs(side).max()) if side.size else float("nan"),
        "sidelobe_rms": float(np.sqrt((side ** 2).mean())) if side.size else float("nan"),
        "min_value": float(p.min()),
        "fwhm_px": float(fwhm_px),
        "ellipticity": beam_ellipticity(p),
    }
    if pixel_scale_arcsec is not None:
        out["fwhm_arcsec"] = float(fwhm_px * pixel_scale_arcsec)
        out["pixel_scale_arcsec"] = float(pixel_scale_arcsec)
    return out
