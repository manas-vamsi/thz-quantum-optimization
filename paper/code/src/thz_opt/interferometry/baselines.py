"""Antenna-pair baselines.

Convention (fixed for the whole repository)
-------------------------------------------
For antenna positions ``r_i = (x_i, y_i)`` in local ENU metres (x East,
y North), the baseline of the ordered pair ``(i, j)`` with ``i < j`` is

    b_ij = r_j - r_i = (x_j - x_i,  y_j - y_i)

Only ``i < j`` pairs are enumerated, giving exactly ``N(N-1)/2`` baselines.
The reversed pair ``(j, i)`` carries ``-b_ij``; it is physically sampled and is
added explicitly by :func:`thz_opt.interferometry.uv.baselines_to_uv` when
``include_conjugate=True``.  A negative ``u`` or ``v`` is a perfectly ordinary
sample, never a sign of poor coverage.
"""

from __future__ import annotations

import numpy as np

from ..arrays.validation import validate_layout


def pair_indices(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Row/column index arrays of all ``i < j`` pairs, in row-major order."""
    i, j = np.triu_indices(n, k=1)
    return i, j


def n_baselines(n: int) -> int:
    """``N(N-1)/2``."""
    return n * (n - 1) // 2


def compute_baselines(xy: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(b, i, j)`` for every unique pair.

    ``b`` has shape ``(N(N-1)/2, 2)`` in metres, ``b[k] = xy[j[k]] - xy[i[k]]``.
    """
    arr = validate_layout(xy)
    i, j = pair_indices(arr.shape[0])
    return arr[j] - arr[i], i, j


def baseline_lengths(xy: np.ndarray) -> np.ndarray:
    """Euclidean baseline lengths ``|b_ij|`` in metres, ``i < j``."""
    b, _, _ = compute_baselines(xy)
    return np.hypot(b[:, 0], b[:, 1])


def baseline_matrix(xy: np.ndarray) -> np.ndarray:
    """Full symmetric ``(n, n)`` distance matrix in metres (zero diagonal)."""
    arr = validate_layout(xy)
    d = arr[:, None, :] - arr[None, :, :]
    return np.hypot(d[..., 0], d[..., 1])
