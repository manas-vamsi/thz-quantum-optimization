"""UV-coverage metrics.

Every metric here is defined mathematically and none of them is "the" UV
coverage score.  They measure different things and can disagree; an experiment
must report several.  In particular the widely quoted "percentage UV coverage"
is *filled-cell fraction on a particular grid over a particular extent* and is
meaningless without both.

Notation
--------
``G`` is the set of grid cells considered, ``|G|`` its size.  ``n_c`` is the
number of samples in cell ``c``, ``O = {c : n_c > 0}`` the occupied cells.
"""

from __future__ import annotations

import numpy as np

from ..interferometry.uv import UVGrid, grid_indices


def occupied_fraction(occ: np.ndarray, mask: np.ndarray | None = None) -> float:
    """``|O| / |G|``.

    ``mask`` restricts ``G`` to a subset of cells (e.g. an annulus between the
    shortest and longest baseline); without it, ``G`` is the full square grid
    including the corners, which no circularly bounded array can ever fill.
    Always report which ``G`` was used.
    """
    o = np.asarray(occ)
    if mask is None:
        return float((o > 0).sum() / o.size)
    m = np.asarray(mask, dtype=bool)
    denom = int(m.sum())
    if denom == 0:
        return 0.0
    return float(((o > 0) & m).sum() / denom)


def unique_cells(occ: np.ndarray) -> int:
    """``|O|`` -- the number of distinct cells carrying at least one sample."""
    return int((np.asarray(occ) > 0).sum())


def annulus_mask(grid: UVGrid, r_inner: float, r_outer: float) -> np.ndarray:
    """Boolean ``(n_cells, n_cells)`` mask of cells whose centre radius lies in
    ``[r_inner, r_outer]`` wavelengths.  Indexed ``[iu, iv]`` like the grid."""
    k = np.arange(grid.n_cells) - grid.n_cells / 2.0
    cu = k * grid.cell_size
    r = np.hypot(cu[:, None], cu[None, :])
    return (r >= r_inner) & (r <= r_outer)


def radial_profile(uv: np.ndarray, n_bins: int = 30, r_max: float | None = None):
    """Histogram of ``|(u, v)|`` -- the radial UV sample distribution.

    Returns ``(bin_centres, counts, density_per_area)`` where the third array
    divides the count by the annulus area ``pi (r_out^2 - r_in^2)`` so that a
    sample distribution uniform per unit UV area gives a flat curve.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    r = np.hypot(arr[:, 0], arr[:, 1])
    hi = float(r.max()) if r_max is None else float(r_max)
    edges = np.linspace(0.0, max(hi, 1e-12), n_bins + 1)
    counts, _ = np.histogram(r, bins=edges)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    return 0.5 * (edges[1:] + edges[:-1]), counts, counts / area


def angular_profile(uv: np.ndarray, n_bins: int = 36):
    """Histogram of baseline position angle folded onto ``[0, pi)``.

    Folding is correct because ``(u, v)`` and ``(-u, -v)`` are the same physical
    baseline; without folding every array looks trivially symmetric.
    Returns ``(bin_centres_rad, counts)``.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    ang = np.mod(np.arctan2(arr[:, 1], arr[:, 0]), np.pi)
    edges = np.linspace(0.0, np.pi, n_bins + 1)
    counts, _ = np.histogram(ang, bins=edges)
    return 0.5 * (edges[1:] + edges[:-1]), counts


def angular_uniformity(uv: np.ndarray, n_bins: int = 36) -> float:
    """Normalised entropy of the folded position-angle histogram, in ``[0, 1]``.

    ``H = -sum p log p / log(n_bins)``; 1 means every angular bin holds the same
    number of samples, 0 means all samples share one bin.  Empty bins contribute
    zero.  This is a *uniformity* measure, not an imaging-quality measure.
    """
    _, counts = angular_profile(uv, n_bins)
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum() / np.log(n_bins))


def density_uniformity(occ: np.ndarray) -> float:
    """Normalised entropy of the occupancy distribution over *occupied* cells.

    1 means every occupied cell holds the same number of samples (no
    redundancy imbalance), lower values mean the samples pile up in few cells.
    Cells outside ``O`` are excluded, so this is independent of the grid extent
    and complements :func:`occupied_fraction`, which depends on it strongly.
    """
    o = np.asarray(occ, dtype=float)
    vals = o[o > 0]
    if vals.size <= 1:
        return 0.0 if vals.size == 0 else 1.0
    p = vals / vals.sum()
    return float(-(p * np.log(p)).sum() / np.log(vals.size))


def coverage_summary(uv: np.ndarray, occ: np.ndarray, grid: UVGrid) -> dict:
    """Bundle of the metrics above, with the grid definition attached.

    ``occupied_fraction_full`` uses the whole square grid; ``occupied_fraction_
    annulus`` restricts ``G`` to the annulus actually reachable by the array,
    between the shortest and longest sampled UV radius.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    r = np.hypot(arr[:, 0], arr[:, 1])
    mask = annulus_mask(grid, float(r.min()), float(r.max()))
    _, _, inside = grid_indices(arr, grid)
    return {
        "n_samples": int(arr.shape[0]),
        "n_samples_outside_grid": int((~inside).sum()),
        "unique_cells": unique_cells(occ),
        "occupied_fraction_full": occupied_fraction(occ),
        "occupied_fraction_annulus": occupied_fraction(occ, mask),
        "annulus_cells": int(mask.sum()),
        "density_uniformity": density_uniformity(occ),
        "angular_uniformity": angular_uniformity(arr),
        "uv_radius_min_lambda": float(r.min()),
        "uv_radius_max_lambda": float(r.max()),
        **grid.as_dict(),
    }
