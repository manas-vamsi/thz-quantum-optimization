"""Precipitable-water-vapour (PWV) cost -- EXPLICIT PLACEHOLDER.

Nothing in this module is measured atmospheric data.  No real site statistics
are shipped with this repository, and none are invented.  What is provided is
the *interface*: a per-pad scalar cost that enters the QUBO as a linear term,

    H_pwv = lambda_pwv * sum_i c_i x_i

together with two synthetic cost models whose only purpose is to exercise that
interface in tests and toy examples:

``constant``   every pad costs the same (the term becomes a constant offset
               under a fixed-N selection constraint and changes nothing)
``altitude``   cost decreasing with a supplied pad altitude, following the
               qualitative expectation that higher, drier pads have lower PWV

Both are labelled ``synthetic=True`` in their metadata.  Replacing them with
real site data is a prerequisite before any PWV-related claim can be made, and
this is listed as an open limitation in the README.
"""

from __future__ import annotations

import numpy as np


def constant_pwv_cost(n_pads: int, value: float = 1.0) -> np.ndarray:
    """Uniform synthetic cost vector ``c_i = value``."""
    return np.full(int(n_pads), float(value))


def altitude_pwv_cost(
    altitude_m: np.ndarray, scale_height_m: float = 2000.0, c0: float = 1.0
) -> np.ndarray:
    """Synthetic cost ``c_i = c0 * exp(-h_i / H)``.

    An exponential fall-off with altitude is the textbook first-order shape of
    a water-vapour column, but the constants here are illustrative only; this
    function must not be used to claim a PWV value for a real site.
    """
    h = np.asarray(altitude_m, dtype=float)
    if scale_height_m <= 0:
        raise ValueError("scale_height_m must be positive")
    return float(c0) * np.exp(-h / scale_height_m)


def normalise_cost(c: np.ndarray) -> np.ndarray:
    """Rescale a cost vector onto ``[0, 1]`` so its weight is comparable with
    the other QUBO terms (see :mod:`thz_opt.qubo.coefficients`)."""
    arr = np.asarray(c, dtype=float)
    lo, hi = float(arr.min()), float(arr.max())
    if hi - lo <= 0:
        return np.zeros_like(arr)
    return (arr - lo) / (hi - lo)


def pwv_metadata(model: str) -> dict:
    """Provenance record attached to any result that used a PWV cost."""
    return {"pwv_model": model, "synthetic": True,
            "note": "illustrative cost model, not measured atmospheric data"}
