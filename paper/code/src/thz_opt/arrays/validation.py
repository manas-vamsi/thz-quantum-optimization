"""Sanity checks shared by all layout generators."""

from __future__ import annotations

import numpy as np


def validate_layout(xy: np.ndarray) -> np.ndarray:
    """Return ``xy`` as a validated float ``(n, 2)`` array.

    Raises on wrong shape, NaN or infinite values.  Coincident antennas are not
    an error here (a QUBO candidate set may legitimately contain pads that are
    never co-selected) -- use
    :func:`thz_opt.constraints.separation.check_minimum_separation` for that.
    """
    arr = np.asarray(xy, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"layout must have shape (n, 2), got {arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError("layout contains NaN or infinite coordinates")
    return arr


def layout_extent(xy: np.ndarray) -> dict:
    """Radial summary of a layout, in the same units as ``xy`` (metres)."""
    arr = validate_layout(xy)
    r = np.hypot(arr[:, 0], arr[:, 1])
    return {
        "n": int(arr.shape[0]),
        "r_min": float(r.min()),
        "r_max": float(r.max()),
        "r_mean": float(r.mean()),
        "centroid_x": float(arr[:, 0].mean()),
        "centroid_y": float(arr[:, 1].mean()),
    }


def recentre(xy: np.ndarray) -> np.ndarray:
    """Translate a layout so its centroid is at the origin.

    Baselines are translation invariant, so this changes no UV result; it only
    affects plots and radial summaries.
    """
    arr = validate_layout(xy)
    return arr - arr.mean(axis=0, keepdims=True)
