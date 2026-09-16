"""Baselines -> (u, v) coordinates and UV-plane discretisation.

Snapshot model
--------------
For a coplanar array observing at the zenith (source at the local zenith, i.e.
declination equal to the site latitude at hour angle zero), the projected
baseline is the physical baseline and

    u = dx / lambda,    v = dy / lambda          [wavelengths]

This is the *simplified snapshot* used for the first experiments.  A general
pointing requires the projection in
:mod:`thz_opt.interferometry.earth_rotation`; the two are deliberately kept in
separate modules so that snapshot results are never silently contaminated by a
rotation model.

UV grid
-------
The grid is square, centred on the origin, spanning ``[-uv_max, +uv_max]`` in
both axes with ``n_cells`` bins per axis, so the cell size is
``2 * uv_max / n_cells`` wavelengths.  Samples outside the extent are dropped
and counted separately -- they are never wrapped or clipped into edge cells.

Cells are *centred* on ``u = (k - n_cells/2) * cell_size`` (``k = 0 ...
n_cells-1``), i.e. one cell is centred exactly on the origin.  This matters:
with an even ``n_cells`` it makes the binning exactly symmetric under
``(u, v) -> (-u, -v)``, so a Hermitian sample set produces a Hermitian
occupancy grid and hence a real dirty beam.  Binning on cell *edges* instead
introduces a half-cell offset and a spurious imaginary part in the PSF.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .baselines import compute_baselines

C_LIGHT = 299_792_458.0  # m / s


def wavelength_from_frequency(freq_hz: float) -> float:
    """lambda = c / nu, in metres."""
    if freq_hz <= 0:
        raise ValueError("frequency must be positive")
    return C_LIGHT / freq_hz


def baselines_to_uv(
    b: np.ndarray, wavelength: float, include_conjugate: bool = True
) -> np.ndarray:
    """Convert baseline vectors in metres to ``(u, v)`` in wavelengths.

    With ``include_conjugate`` the returned array holds both ``(u, v)`` and
    ``(-u, -v)`` for every baseline, i.e. twice as many rows.
    """
    if wavelength <= 0:
        raise ValueError("wavelength must be positive")
    uv = np.asarray(b, dtype=float).reshape(-1, 2) / wavelength
    if include_conjugate:
        uv = np.vstack((uv, -uv))
    return uv


def layout_to_uv(
    xy: np.ndarray, wavelength: float, include_conjugate: bool = True
) -> np.ndarray:
    """Snapshot ``(u, v)`` samples of a layout, in wavelengths."""
    b, _, _ = compute_baselines(xy)
    return baselines_to_uv(b, wavelength, include_conjugate)


@dataclass(frozen=True)
class UVGrid:
    """Square UV grid definition.

    ``n_cells`` bins per axis over ``[-uv_max, uv_max]``; ``cell_size`` is
    derived.  Frozen so a grid can be recorded verbatim in experiment metadata.
    """

    uv_max: float
    n_cells: int

    def __post_init__(self) -> None:
        if self.uv_max <= 0:
            raise ValueError("uv_max must be positive")
        if self.n_cells < 1:
            raise ValueError("n_cells must be >= 1")

    @property
    def cell_size(self) -> float:
        return 2.0 * self.uv_max / self.n_cells

    @property
    def edges(self) -> np.ndarray:
        return np.linspace(-self.uv_max, self.uv_max, self.n_cells + 1)

    @property
    def total_cells(self) -> int:
        return self.n_cells * self.n_cells

    def as_dict(self) -> dict:
        return {
            "uv_max_lambda": self.uv_max,
            "n_cells": self.n_cells,
            "cell_size_lambda": self.cell_size,
            "total_cells": self.total_cells,
        }

    @classmethod
    def from_uv(cls, uv: np.ndarray, n_cells: int, pad: float = 1.02) -> "UVGrid":
        """Grid that just contains ``uv`` (``pad`` > 1 leaves a small margin)."""
        m = float(np.max(np.abs(np.asarray(uv, dtype=float)))) * pad
        return cls(uv_max=max(m, 1e-12), n_cells=n_cells)


def grid_indices(uv: np.ndarray, grid: UVGrid):
    """Map samples to cell indices.

    Returns ``(iu, iv, inside)``: integer indices along u and v for the samples
    with ``inside == True``, plus the boolean mask itself.  Index ``k``
    corresponds to the cell centred on ``(k - n_cells/2) * cell_size``.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    inside = (np.abs(arr[:, 0]) <= grid.uv_max) & (np.abs(arr[:, 1]) <= grid.uv_max)
    sel = arr[inside]
    iu = np.clip(
        np.floor((sel[:, 0] + grid.uv_max) / grid.cell_size + 0.5).astype(int),
        0,
        grid.n_cells - 1,
    )
    iv = np.clip(
        np.floor((sel[:, 1] + grid.uv_max) / grid.cell_size + 0.5).astype(int),
        0,
        grid.n_cells - 1,
    )
    return iu, iv, inside


def uv_occupancy(uv: np.ndarray, grid: UVGrid):
    """Integer occupancy matrix and the number of samples outside the extent.

    ``occ[a, b]`` counts samples whose u-index is ``a`` and v-index is ``b``,
    i.e. the array is indexed ``[iu, iv]``.  Plot ``occ.T`` with
    ``origin="lower"`` so that u runs along the horizontal axis.
    """
    iu, iv, inside = grid_indices(uv, grid)
    occ = np.zeros((grid.n_cells, grid.n_cells), dtype=int)
    np.add.at(occ, (iu, iv), 1)
    return occ, int((~inside).sum())


def cell_ids(uv: np.ndarray, grid: UVGrid) -> np.ndarray:
    """Flat cell id ``iu * n_cells + iv`` for each in-extent sample."""
    iu, iv, _ = grid_indices(uv, grid)
    return iu * grid.n_cells + iv
