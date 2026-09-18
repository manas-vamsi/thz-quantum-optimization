"""Pareto fronts -- reporting the trade-off instead of guessing a weight.

The objectives in this project genuinely conflict: coverage wants many distinct
UV cells, sidelobe energy wants the samples spread evenly, coherence wants short
baselines, resolution wants long ones. Combining them as
``alpha*A + beta*B + eta*C`` forces a choice of weights that nobody can justify
from first principles, and the choice silently determines the answer.

A Pareto front avoids the question. A layout is **dominated** if some other
layout is at least as good on every objective and strictly better on one. The
non-dominated set is the complete menu of sensible compromises; a scientist
picks a point on it using mission priorities, which is a decision they are
qualified to make and an optimiser is not.

What this module does not do: claim the front is complete. Sweeping scalar
weights finds only the *convex* part of the front -- points in a non-convex
dent are unreachable by any weighting. That limitation is stated in
:func:`sweep_weights` and is why the front is reported as "found", not "the".
"""

from __future__ import annotations

import numpy as np

__all__ = ["dominates", "pareto_indices", "sweep_weights"]


def dominates(a, b, tol: float = 0.0) -> bool:
    """Does objective vector ``a`` dominate ``b``? (all objectives minimised)"""
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return bool(np.all(x <= y + tol) and np.any(x < y - tol))


def pareto_indices(points, tol: float = 0.0) -> np.ndarray:
    """Indices of the non-dominated rows of ``points`` (all minimised).

    ``O(n^2)``, which is ample for the tens-to-hundreds of points a weight
    sweep produces.
    """
    p = np.asarray(points, dtype=float)
    if p.ndim != 2:
        raise ValueError("points must be 2-D: one row per candidate, one column per objective")
    keep = []
    for i in range(p.shape[0]):
        if not any(dominates(p[j], p[i], tol) for j in range(p.shape[0]) if j != i):
            keep.append(i)
    return np.asarray(keep, dtype=int)


def sweep_weights(run_one, weight_grid, objective_names) -> dict:
    """Run an optimiser across a grid of scalarisation weights.

    ``run_one(weights) -> (objective_vector, payload)``; ``weight_grid`` is an
    iterable of weight tuples. Returns the full set of runs, the indices of the
    non-dominated ones, and the front itself.

    **Caveat, stated because it is easy to forget:** scalarised sweeps recover
    only the convex hull of the Pareto front. If the true front has a concave
    region, no choice of positive weights will find points inside it. For a
    complete front one needs epsilon-constraint scanning (optimise one
    objective subject to bounds on the others), which is a straightforward but
    separate piece of work.
    """
    objectives, payloads, weights = [], [], []
    for w in weight_grid:
        vec, payload = run_one(w)
        objectives.append(np.asarray(vec, dtype=float))
        payloads.append(payload)
        weights.append(tuple(float(x) for x in np.atleast_1d(w)))

    obj = np.vstack(objectives) if objectives else np.empty((0, len(objective_names)))
    front = pareto_indices(obj) if obj.size else np.empty(0, dtype=int)
    return {
        "objective_names": list(objective_names),
        "weights": weights,
        "objectives": obj,
        "payloads": payloads,
        "front_indices": front,
        "n_runs": int(obj.shape[0]),
        "n_nondominated": int(front.size),
        "note": "scalarised sweep recovers only the convex part of the front",
    }
