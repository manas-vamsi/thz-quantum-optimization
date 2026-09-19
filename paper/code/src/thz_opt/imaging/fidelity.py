"""How wrong is the reconstruction? Several answers, because one is not enough.

No single scalar captures image quality, and picking one would repeat the
mistake this package exists to correct. Four are reported, each sensitive to a
different failure:

``rms_error``
    Root-mean-square difference from the truth, after both are convolved to the
    same resolution. The general-purpose number.

``dynamic_range``
    Peak of the restored image over the RMS of the residual. This is the figure
    radio astronomers quote, and it is dominated by how well sidelobes were
    deconvolved.

``flux_recovery``
    Recovered flux over true flux. Catches the specific failure of an array
    with no short spacings: it resolves extended emission out and quietly
    reports a fraction of the source. A cell-count metric cannot see this at
    all.

``fidelity``
    Median of ``|truth| / |truth - image|`` over pixels above a cutoff, the
    conventional image-fidelity measure. Robust to a handful of bad pixels in a
    way the RMS is not.

The truth is always smoothed by the same restoring beam before comparison.
Comparing a restored image against an unsmoothed model would charge every array
for resolution it never claimed to have.
"""

from __future__ import annotations

import numpy as np

__all__ = ["fidelity_metrics", "smooth_to_beam"]


def smooth_to_beam(image: np.ndarray, psf: np.ndarray) -> np.ndarray:
    """Convolve a model with the restoring beam fitted to ``psf``."""
    from .clean import fit_restoring_beam

    beam, _ = fit_restoring_beam(psf)
    ft = np.fft.fft2(np.fft.ifftshift(beam))
    out = np.fft.fftshift(np.fft.ifft2(
        np.fft.fft2(np.fft.ifftshift(image)) * ft)).real
    return out / beam.sum()


def fidelity_metrics(restored: np.ndarray, truth, psf: np.ndarray,
                     residual: np.ndarray | None = None,
                     cutoff_fraction: float = 0.05) -> dict:
    """Compare a restored image with the model it came from."""
    model = np.asarray(truth.image if hasattr(truth, "image") else truth,
                       dtype=float)
    smooth = smooth_to_beam(model, psf)
    diff = restored - smooth

    peak = float(np.abs(smooth).max())
    mask = np.abs(smooth) > cutoff_fraction * peak
    denom = np.abs(diff[mask])
    denom = np.where(denom > 0, denom, np.finfo(float).tiny)

    res_rms = (float(np.sqrt(np.mean(residual ** 2)))
               if residual is not None else float(np.sqrt(np.mean(diff ** 2))))

    return {
        "rms_error": float(np.sqrt(np.mean(diff ** 2))),
        "rms_error_relative": float(np.sqrt(np.mean(diff ** 2)) / peak),
        "max_error_relative": float(np.abs(diff).max() / peak),
        "dynamic_range": float(np.abs(restored).max() / res_rms)
        if res_rms > 0 else float("inf"),
        "flux_recovery": float(restored.sum() / model.sum())
        if model.sum() != 0 else float("nan"),
        "fidelity": float(np.median(np.abs(smooth[mask]) / denom))
        if mask.any() else float("nan"),
        "pixels_compared": int(mask.sum()),
    }
