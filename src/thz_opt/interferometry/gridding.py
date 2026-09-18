"""UV gridding kernels, and why the hard-edged version is not enough.

Everything measured so far bins each UV sample into exactly one cell. That is
simple and exactly reproducible, but it has a real artifact: a sample sitting a
hair either side of a cell boundary lands in a different cell, so
coverage-style metrics jump discontinuously as an antenna moves by centimetres.
Two configurations that are physically indistinguishable can score differently,
and an optimiser will happily exploit that.

Two standard fixes, both implemented here:

``gaussian_gridded_occupancy``
    Spread each sample over neighbouring cells with a Gaussian kernel, the way
    real imaging pipelines grid visibilities (they use prolate spheroidal
    functions; a Gaussian is the readable approximation). Cell occupancies
    become non-integer.

``shifted_grid_average``
    Evaluate a metric on several sub-cell-shifted grids and average. This costs
    nothing conceptually and directly measures how much of a result is grid
    artifact.

**The formulation is unaffected.** With a kernel, ``a_ck`` is a float instead
of an integer, but ``n_c = sum_k a_ck y_k`` is still *linear* in the baseline
variables, so ``sum_c n_c^2`` is still exactly quadratic. Softening the grid
buys robustness and costs nothing in the QUBO -- worth stating explicitly,
because the opposite would have been a good reason not to do it.
"""

from __future__ import annotations

import numpy as np

from .uv import UVGrid

__all__ = [
    "gaussian_gridded_occupancy",
    "gaussian_cell_weights",
    "shifted_grid_average",
    "grid_sensitivity",
]


def gaussian_cell_weights(
    uv: np.ndarray,
    grid: UVGrid,
    kernel_cells: float = 0.7,
    truncate: float = 3.0,
    offset: tuple = (0.0, 0.0),
) -> dict:
    """Per-cell weights for each sample under a Gaussian gridding kernel.

    Returns ``{flat_cell_id: total_weight}``. Each sample contributes total
    weight 1 (up to truncation), spread over the cells within ``truncate``
    kernel widths. ``kernel_cells`` is the Gaussian sigma in units of the cell
    size; 0 reproduces nearest-cell binning. ``offset`` shifts the grid origin
    by that fraction of a cell, for :func:`shifted_grid_average`.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    n = grid.n_cells
    du = grid.cell_size

    # continuous cell coordinates, cell centres at integer positions
    fu = (arr[:, 0] + grid.uv_max) / du - offset[0]
    fv = (arr[:, 1] + grid.uv_max) / du - offset[1]

    if kernel_cells <= 0:
        iu = np.floor(fu + 0.5).astype(int)
        iv = np.floor(fv + 0.5).astype(int)
        keep = (iu >= 0) & (iu < n) & (iv >= 0) & (iv < n)
        out: dict = {}
        for a, b in zip(iu[keep], iv[keep]):
            key = int(a) * n + int(b)
            out[key] = out.get(key, 0.0) + 1.0
        return out

    radius = int(np.ceil(truncate * kernel_cells))
    steps = np.arange(-radius, radius + 1)
    base_u = np.floor(fu + 0.5).astype(int)
    base_v = np.floor(fv + 0.5).astype(int)

    out = {}
    two_s2 = 2.0 * kernel_cells * kernel_cells
    for su in steps:
        cu = base_u + su
        wu = np.exp(-((cu - fu) ** 2) / two_s2)
        for sv in steps:
            cv = base_v + sv
            w = wu * np.exp(-((cv - fv) ** 2) / two_s2)
            keep = (cu >= 0) & (cu < n) & (cv >= 0) & (cv < n) & (w > 1e-6)
            if not keep.any():
                continue
            ids = cu[keep] * n + cv[keep]
            ws = w[keep]
            for key, val in zip(ids.tolist(), ws.tolist()):
                out[key] = out.get(key, 0.0) + val

    # normalise so every sample carries unit total weight
    total_per_sample = 0.0
    if out:
        total_per_sample = sum(out.values()) / max(arr.shape[0], 1)
    if total_per_sample > 0:
        scale = 1.0 / total_per_sample
        out = {k: v * scale for k, v in out.items()}
    return out


def gaussian_gridded_occupancy(
    uv: np.ndarray, grid: UVGrid, kernel_cells: float = 0.7, truncate: float = 3.0
) -> np.ndarray:
    """Occupancy grid with a Gaussian gridding kernel, indexed ``[iu, iv]``.

    Values are floats. With ``kernel_cells = 0`` this reduces exactly to
    :func:`thz_opt.interferometry.uv.uv_occupancy`.
    """
    w = gaussian_cell_weights(uv, grid, kernel_cells, truncate)
    occ = np.zeros((grid.n_cells, grid.n_cells))
    if not w:
        return occ
    ids = np.fromiter(w.keys(), dtype=np.int64, count=len(w))
    vals = np.fromiter(w.values(), dtype=float, count=len(w))
    occ[ids // grid.n_cells, ids % grid.n_cells] = vals
    return occ


def shifted_grid_average(metric, uv: np.ndarray, grid: UVGrid,
                         n_shifts: int = 4, kernel_cells: float = 0.0) -> dict:
    """Average a metric over sub-cell grid shifts, and report the spread.

    ``metric(occupancy) -> float`` is evaluated on ``n_shifts**2`` grids offset
    by fractions of a cell. The returned ``relative_spread`` is the honest
    measure of how much of a single-grid number is artifact: a metric that
    changes by 10 % under a half-cell shift should not be quoted to three
    significant figures.
    """
    vals = []
    for a in range(n_shifts):
        for b in range(n_shifts):
            off = (a / n_shifts, b / n_shifts)
            w = gaussian_cell_weights(uv, grid, kernel_cells, offset=off)
            occ = np.zeros((grid.n_cells, grid.n_cells))
            if w:
                ids = np.fromiter(w.keys(), dtype=np.int64, count=len(w))
                vv = np.fromiter(w.values(), dtype=float, count=len(w))
                occ[ids // grid.n_cells, ids % grid.n_cells] = vv
            vals.append(float(metric(occ)))
    v = np.asarray(vals)
    mean = float(v.mean())
    return {
        "mean": mean,
        "std": float(v.std()),
        "min": float(v.min()),
        "max": float(v.max()),
        "relative_spread": float(v.std() / abs(mean)) if mean else float("nan"),
        "n_grids": int(v.size),
    }


def grid_sensitivity(uv: np.ndarray, grid: UVGrid, n_shifts: int = 4) -> dict:
    """How grid-dependent are the two headline metrics for this sampling?"""
    cells = shifted_grid_average(lambda o: float((o > 0).sum()), uv, grid, n_shifts)
    energy = shifted_grid_average(lambda o: float((o ** 2).sum()), uv, grid, n_shifts)
    return {
        "unique_cells": cells,
        "sidelobe_energy": energy,
        "n_cells_per_axis": grid.n_cells,
    }
