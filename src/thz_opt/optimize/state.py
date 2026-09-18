"""Incremental UV occupancy for a pad selection.

Every optimiser in this package needs to ask "what happens to the UV coverage
if I move this one antenna?" thousands of times.  Recomputing the whole
occupancy grid each time is what makes naive search unusable, so this class
maintains the occupancy incrementally.

Adding pad ``p`` to a selection ``S`` activates exactly ``|S|`` new baselines,
those between ``p`` and each already-selected pad.  Updating the cell counts
therefore costs ``O(|S| x cells-per-baseline)``, not ``O(|S|^2 x ...)``.

Two quantities are tracked, because they are the two competing objectives
identified in ``docs/objective_function.md``:

``unique_cells``
    ``|{c : n_c > 0}|`` -- how much of the UV plane is touched at all.

``sum_sq``
    ``sum_c n_c^2`` -- by Parseval this is the dirty beam's total energy, so
    minimising it minimises PSF sidelobes.

Both are exact at every step: :meth:`UVState.recompute` re-derives them from
scratch and the tests assert the incremental and from-scratch values agree
after long random walks.
"""

from __future__ import annotations

import numpy as np

from ..interferometry.baselines import pair_indices

__all__ = ["UVState", "build_pair_tables"]


def build_pair_tables(multiplicities: list):
    """Convert ``{cell: count}`` dicts into flat arrays for fast updates.

    Returns ``(cells, counts)``: lists of ``int`` / ``int`` arrays, one entry
    per candidate pair, in ``i < j`` order.
    """
    cells, counts = [], []
    for d in multiplicities:
        if d:
            ks = np.fromiter(d.keys(), dtype=np.int64, count=len(d))
            vs = np.fromiter(d.values(), dtype=np.int64, count=len(d))
            order = np.argsort(ks)
            cells.append(ks[order])
            counts.append(vs[order])
        else:
            cells.append(np.empty(0, dtype=np.int64))
            counts.append(np.empty(0, dtype=np.int64))
    return cells, counts


class UVState:
    """Mutable UV occupancy of a pad selection, with O(|S|) updates.

    Parameters
    ----------
    cells, counts : per-pair cell ids and sample multiplicities, from
        :func:`build_pair_tables`.
    n_pads : number of candidate pads.
    n_cells_total : size of the flat cell-id space (``grid.n_cells ** 2``).
    """

    def __init__(self, cells: list, counts: list, n_pads: int, n_cells_total: int,
                 cell_weights: np.ndarray | None = None):
        self.n_pads = int(n_pads)
        self._cells = cells
        self._counts = counts
        i_idx, j_idx = pair_indices(self.n_pads)
        self._pair_of = np.full((self.n_pads, self.n_pads), -1, dtype=np.int64)
        self._pair_of[i_idx, j_idx] = np.arange(i_idx.size)
        self._pair_of[j_idx, i_idx] = np.arange(i_idx.size)

        dtype = float if any(c.dtype.kind == "f" for c in counts if c.size) else np.int64
        self._n = np.zeros(int(n_cells_total), dtype=dtype)
        self.selected = np.zeros(self.n_pads, dtype=bool)
        self.unique_cells = 0
        self.sum_sq = 0.0

        # Optional per-cell weighting.  A science case that cares about one part
        # of the UV plane (outer for compact sources, inner for extended
        # emission) expresses that as a weight per cell; the incremental update
        # already touches only the affected cells, so tracking the weighted
        # quantities costs nothing extra.
        self._w = None if cell_weights is None else np.asarray(cell_weights, dtype=float)
        if self._w is not None and self._w.shape[0] != self._n.shape[0]:
            raise ValueError("cell_weights must have one entry per grid cell")
        self.weighted_cells = 0.0
        self.weighted_sum_sq = 0.0

    # -- bookkeeping ------------------------------------------------------

    def _apply(self, k: int, sign: int) -> None:
        c, a = self._cells[k], self._counts[k]
        if c.size == 0:
            return
        old = self._n[c]
        new = old + sign * a
        self.sum_sq += float(np.sum(new * new - old * old))
        self.unique_cells += int(np.count_nonzero(new) - np.count_nonzero(old))
        if self._w is not None:
            w = self._w[c]
            self.weighted_sum_sq += float(np.sum(w * (new * new - old * old)))
            self.weighted_cells += float(np.sum(w * ((new != 0).astype(float)
                                                     - (old != 0).astype(float))))
        self._n[c] = new

    def add_pad(self, p: int) -> None:
        """Select pad ``p`` and activate its baselines to the current set."""
        if self.selected[p]:
            return
        for q in np.flatnonzero(self.selected):
            self._apply(int(self._pair_of[p, q]), +1)
        self.selected[p] = True

    def remove_pad(self, p: int) -> None:
        """Deselect pad ``p`` and deactivate its baselines."""
        if not self.selected[p]:
            return
        self.selected[p] = False
        for q in np.flatnonzero(self.selected):
            self._apply(int(self._pair_of[p, q]), -1)

    def set_selection(self, selection) -> None:
        """Reset to an arbitrary selection (boolean mask or index list)."""
        sel = np.zeros(self.n_pads, dtype=bool)
        idx = np.asarray(selection)
        if idx.size:
            sel[idx.astype(bool) if idx.dtype == bool else idx.astype(int)] = True
        for p in np.flatnonzero(self.selected & ~sel):
            self.remove_pad(int(p))
        for p in np.flatnonzero(sel & ~self.selected):
            self.add_pad(int(p))

    # -- queries ----------------------------------------------------------

    @property
    def n_selected(self) -> int:
        return int(self.selected.sum())

    def selection_indices(self) -> np.ndarray:
        return np.flatnonzero(self.selected)

    def copy_selection(self) -> np.ndarray:
        return self.selected.copy()

    def recompute(self):
        """Recompute both quantities from scratch (validation reference)."""
        n = np.zeros_like(self._n)
        sel = np.flatnonzero(self.selected)
        for a in range(sel.size):
            for b in range(a + 1, sel.size):
                k = int(self._pair_of[sel[a], sel[b]])
                c, m = self._cells[k], self._counts[k]
                if c.size:
                    n[c] += m
        return int(np.count_nonzero(n)), float(np.sum(n * n))

    def recompute_weighted(self):
        """From-scratch weighted quantities (validation reference)."""
        if self._w is None:
            return 0.0, 0.0
        n = np.zeros_like(self._n)
        sel = np.flatnonzero(self.selected)
        for a in range(sel.size):
            for b in range(a + 1, sel.size):
                k = int(self._pair_of[sel[a], sel[b]])
                c, m = self._cells[k], self._counts[k]
                if c.size:
                    n[c] += m
        occupied = (n != 0).astype(float)
        return float(np.sum(self._w * occupied)), float(np.sum(self._w * n * n))

    def occupancy(self) -> np.ndarray:
        """Flat occupancy vector (copy)."""
        return self._n.copy()
