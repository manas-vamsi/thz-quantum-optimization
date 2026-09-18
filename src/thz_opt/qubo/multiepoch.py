"""Multi-epoch (reconfigurable) array selection as a QUBO.

This is the formulation the golden-spiral paper was reaching for. That paper
proposed *moving* antennas during an observation to fill the UV plane, but the
mechanism it described -- rigidly rotating the whole array -- cannot help,
because a rotation is an isometry and leaves every baseline length unchanged.
What *does* help is occupying a **different set of pads** in each epoch, so the
union of the epochs samples more of the UV plane than any single configuration
could.

That is a selection problem again, one per epoch, coupled by two things:

1. **Shared coverage.** The objective sees the union over all epochs, so a
   baseline used in epoch 1 changes what is worth using in epoch 2.
2. **Movement cost.** Physically relocating an antenna costs time, money and
   recalibration, so a configuration that changes little between epochs is
   preferable at equal image quality.

Variables
---------
``x[i,t]``  pad ``i`` occupied in epoch ``t``     -- ``M*T`` variables
``y[k,t]``  baseline ``k`` active in epoch ``t``  -- ``P*T`` variables

Terms
-----
* per-epoch cardinality: ``lambda_sel * (sum_i x[i,t] - N)^2`` for each ``t``
* Rosenberg: ``y[k,t] = x[i,t] x[j,t]`` for each ``t``
* objective over the pooled sampling function
  ``n_c = sum_{k,t} a_ck y[k,t]`` -- still exactly quadratic, just over a
  larger index set
* movement cost, using the exact binary identity

      |x - x'| = x + x' - 2 x x'

  so ``lambda_move * sum_{i,t} |x[i,t+1] - x[i,t]|`` is quadratic with no
  auxiliary variables. Note it counts each *change* (a move out plus a move in
  counts as two), which is the natural measure when antennas are transported
  individually.

Cost, stated honestly
---------------------
Variables grow linearly in ``T``, which is fine. The coupling count is the risk:
the pooled objective couples ``(k,t)`` to ``(l,s)`` across epochs, so the dense
worst case is ``(P*T)^2``. In practice the overlap graph is sparse -- two
baselines interact only if their UV tracks share a cell -- and
:func:`multiepoch_size_report` measures the real number before anything is
built, so the model is never promised at a size it cannot reach.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..interferometry.baselines import pair_indices

__all__ = [
    "MultiEpochSpec",
    "multiepoch_size_report",
    "build_multiepoch_qubo",
    "complete_multiepoch_auxiliary",
    "movement_count",
]


@dataclass(frozen=True)
class MultiEpochSpec:
    """Problem size and penalty weights for the multi-epoch model."""

    n_pads: int
    n_epochs: int
    n_select: int
    lambda_select: float
    lambda_rosenberg: float
    lambda_move: float = 0.0

    @property
    def n_pairs(self) -> int:
        return self.n_pads * (self.n_pads - 1) // 2

    @property
    def n_variables(self) -> int:
        return self.n_epochs * (self.n_pads + self.n_pairs)

    def x(self, i: int, t: int) -> int:
        return t * (self.n_pads + self.n_pairs) + i

    def y(self, k: int, t: int) -> int:
        return t * (self.n_pads + self.n_pairs) + self.n_pads + k

    def as_dict(self) -> dict:
        return {
            "n_pads": self.n_pads, "n_epochs": self.n_epochs,
            "n_select": self.n_select, "n_pairs": self.n_pairs,
            "n_variables": self.n_variables,
            "lambda_select": self.lambda_select,
            "lambda_rosenberg": self.lambda_rosenberg,
            "lambda_move": self.lambda_move,
        }


def multiepoch_size_report(multiplicities_per_epoch: list, n_pads: int,
                           n_epochs: int) -> dict:
    """Measure the real coupling count before committing to building the model.

    ``multiplicities_per_epoch[t][k]`` is ``{cell: count}`` for pair ``k`` in
    epoch ``t`` (epochs may share tracks, or differ if their hour angles do).
    """
    n_pairs = n_pads * (n_pads - 1) // 2
    n_units = n_pairs * n_epochs

    by_cell: dict = {}
    for t in range(n_epochs):
        for k, d in enumerate(multiplicities_per_epoch[t]):
            for c in d:
                by_cell.setdefault(c, set()).add((k, t))

    couplings = set()
    for owners in by_cell.values():
        if len(owners) < 2:
            continue
        ow = sorted(owners)
        for a in range(len(ow)):
            for b in range(a + 1, len(ow)):
                couplings.add((ow[a], ow[b]))

    dense = n_units * (n_units - 1) // 2
    return {
        "n_epochs": n_epochs,
        "n_pads": n_pads,
        "n_pairs_per_epoch": n_pairs,
        "n_baseline_variables": n_units,
        "n_variables_total": n_epochs * (n_pads + n_pairs),
        "n_objective_couplings": len(couplings),
        "n_objective_couplings_dense": dense,
        "fill_fraction": len(couplings) / dense if dense else 0.0,
        "n_cells_touched": len(by_cell),
    }


def build_multiepoch_qubo(multiplicities_per_epoch: list, spec: MultiEpochSpec,
                          max_variables: int = 200000):
    """Assemble the multi-epoch QUBO minimising pooled PSF sidelobe energy.

    The objective is ``sum_c n_c^2`` with ``n_c`` pooled over epochs -- the same
    Parseval-backed quantity as the single-epoch model, so the two are directly
    comparable. Returns ``(Q, offset, meta)`` upper triangular, in the
    convention of :mod:`thz_opt.qubo.objective`.
    """
    if spec.n_variables > max_variables:
        raise ValueError(f"model needs {spec.n_variables} variables "
                         f"(> max_variables={max_variables})")

    i_idx, j_idx = pair_indices(spec.n_pads)
    V = spec.n_variables
    Q = np.zeros((V, V))
    offset = 0.0

    def add(u: int, v: int, val: float) -> None:
        a, b = (u, v) if u <= v else (v, u)
        Q[a, b] += val

    # per-epoch cardinality and Rosenberg
    for t in range(spec.n_epochs):
        for i in range(spec.n_pads):
            add(spec.x(i, t), spec.x(i, t), spec.lambda_select * (1.0 - 2.0 * spec.n_select))
            for j in range(i + 1, spec.n_pads):
                add(spec.x(i, t), spec.x(j, t), 2.0 * spec.lambda_select)
        offset += spec.lambda_select * spec.n_select ** 2

        for k in range(spec.n_pairs):
            xi, xj, yk = spec.x(int(i_idx[k]), t), spec.x(int(j_idx[k]), t), spec.y(k, t)
            add(xi, xj, spec.lambda_rosenberg)
            add(xi, yk, -2.0 * spec.lambda_rosenberg)
            add(xj, yk, -2.0 * spec.lambda_rosenberg)
            add(yk, yk, 3.0 * spec.lambda_rosenberg)

    # pooled sidelobe energy over all (pair, epoch) units
    by_cell: dict = {}
    for t in range(spec.n_epochs):
        for k, d in enumerate(multiplicities_per_epoch[t]):
            for c, v in d.items():
                by_cell.setdefault(c, []).append((spec.y(k, t), float(v)))

    for owners in by_cell.values():
        for a in range(len(owners)):
            va, wa = owners[a]
            add(va, va, wa * wa)
            for b in range(a + 1, len(owners)):
                vb, wb = owners[b]
                add(va, vb, 2.0 * wa * wb)

    # movement cost: |x[i,t+1] - x[i,t]| = x[i,t+1] + x[i,t] - 2 x[i,t+1] x[i,t]
    n_move_terms = 0
    if spec.lambda_move:
        for t in range(spec.n_epochs - 1):
            for i in range(spec.n_pads):
                a, b = spec.x(i, t), spec.x(i, t + 1)
                add(a, a, spec.lambda_move)
                add(b, b, spec.lambda_move)
                add(a, b, -2.0 * spec.lambda_move)
                n_move_terms += 1

    meta = {
        **spec.as_dict(),
        "n_quadratic_terms": int(np.count_nonzero(np.triu(Q, 1))),
        "n_movement_terms": n_move_terms,
        "n_cells_touched": len(by_cell),
        "objective": "pooled_psf_sidelobe_energy",
    }
    return Q, float(offset), meta


def complete_multiepoch_auxiliary(x_matrix: np.ndarray, spec: MultiEpochSpec) -> np.ndarray:
    """Full variable vector from an ``(n_epochs, n_pads)`` selection matrix."""
    xm = np.asarray(x_matrix, dtype=int)
    if xm.shape != (spec.n_epochs, spec.n_pads):
        raise ValueError(f"expected shape {(spec.n_epochs, spec.n_pads)}, got {xm.shape}")
    i_idx, j_idx = pair_indices(spec.n_pads)
    v = np.zeros(spec.n_variables, dtype=int)
    for t in range(spec.n_epochs):
        for i in range(spec.n_pads):
            v[spec.x(i, t)] = xm[t, i]
        prod = xm[t][i_idx] * xm[t][j_idx]
        for k in range(spec.n_pairs):
            v[spec.y(k, t)] = int(prod[k])
    return v


def movement_count(x_matrix: np.ndarray) -> int:
    """Total number of pad-occupancy changes between consecutive epochs."""
    xm = np.asarray(x_matrix, dtype=int)
    return int(np.abs(np.diff(xm, axis=0)).sum())
