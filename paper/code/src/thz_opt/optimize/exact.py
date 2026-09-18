"""Provably optimal pad selection by mixed-integer linear programming.

Why this matters more than it looks
-----------------------------------
Every result in the array-configuration literature is the output of a
heuristic -- simulated annealing, a genetic algorithm, a greedy pass -- and is
reported without any statement of how far from optimal it is. That is not a
criticism of those papers; the problem is NP-hard (it is densest-k-subgraph,
see ``docs/PROJECT_BRIEF.md`` §5.5) and exact methods do not scale.

But at moderate ``M`` they *do* run, and then the answer is not "good" but
**optimal, with proof**. That converts every other number in the project from
"our method scored X" into "our method scored X against a known ceiling of Y" --
which is the only honest way to evaluate a quantum solver later.

The linearisation
-----------------
Maximising unique UV-cell coverage is naturally a *linear* integer program,
even though it is quartic as a QUBO:

    maximise    sum_c z_c

    subject to  y_k <= x_i,  y_k <= x_j          for pair k = (i,j)
                z_c <= sum_{k : c in C_k} y_k    for every cell c
                sum_i x_i = N
                x_i + x_j <= 1                   for shadowed pairs
                x in {0,1}^M,  y in {0,1}^P,  0 <= z_c <= 1

The first pair of constraints makes ``y_k <= x_i x_j``; because coverage is
being maximised there is no incentive to set ``y_k`` below that bound, so
equality holds at the optimum without needing ``y_k >= x_i + x_j - 1``. The
same argument lets ``z_c`` stay continuous: it is pushed to
``min(1, sum y_k)`` automatically, which is integral whenever ``y`` is. Both
simplifications shrink the model substantially, and the tests verify the
result against brute-force enumeration.

Solved with HiGHS via :func:`scipy.optimize.milp`, so no extra dependency.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from ..interferometry.baselines import pair_indices

__all__ = ["max_coverage_milp"]


def max_coverage_milp(
    cell_sets: list,
    n_pads: int,
    n_select: int,
    forbidden_pairs=None,
    time_limit: float | None = 60.0,
    mip_rel_gap: float = 0.0,
):
    """Maximise the number of distinct UV cells covered, exactly.

    Parameters
    ----------
    cell_sets : per-pair arrays of UV cell ids (unique), in ``i < j`` order.
    n_pads, n_select : ``M`` and ``N``.
    forbidden_pairs : ``(M, M)`` boolean matrix or list of ``(i, j)`` pairs that
        may not both be selected (shadowing).
    time_limit : seconds; ``None`` for no limit. On timeout the best incumbent
        is returned with ``proven_optimal=False``.
    mip_rel_gap : stop when the relative gap is below this. 0 means prove
        optimality.

    Returns
    -------
    dict with the selection, the optimal coverage, the solver status, and
    ``proven_optimal`` -- which is the whole point of using this.
    """
    i_idx, j_idx = pair_indices(n_pads)
    n_pairs = len(cell_sets)
    if i_idx.size != n_pairs:
        raise ValueError(f"cell_sets has {n_pairs} entries, expected {i_idx.size}")

    cells = sorted({int(c) for cs in cell_sets for c in np.asarray(cs).tolist()})
    cell_index = {c: t for t, c in enumerate(cells)}
    n_cells = len(cells)
    if n_cells == 0:
        raise ValueError("no UV cells are covered by any candidate pair")

    # variable layout: [ x (M) | y (P) | z (C) ]
    off_y, off_z = n_pads, n_pads + n_pairs
    n_var = n_pads + n_pairs + n_cells

    rows, cols, data, lb, ub = [], [], [], [], []
    r = 0

    def add_row(entries, lo, hi):
        nonlocal r
        for col, val in entries:
            rows.append(r)
            cols.append(col)
            data.append(val)
        lb.append(lo)
        ub.append(hi)
        r += 1

    # y_k - x_i <= 0  and  y_k - x_j <= 0
    for k in range(n_pairs):
        add_row([(off_y + k, 1.0), (int(i_idx[k]), -1.0)], -np.inf, 0.0)
        add_row([(off_y + k, 1.0), (int(j_idx[k]), -1.0)], -np.inf, 0.0)

    # z_c - sum_k y_k <= 0
    cell_to_pairs: dict = {}
    for k, cs in enumerate(cell_sets):
        for c in np.asarray(cs).tolist():
            cell_to_pairs.setdefault(int(c), []).append(k)
    for c, ks in cell_to_pairs.items():
        entries = [(off_z + cell_index[c], 1.0)] + [(off_y + k, -1.0) for k in ks]
        add_row(entries, -np.inf, 0.0)

    # sum_i x_i == N
    add_row([(i, 1.0) for i in range(n_pads)], float(n_select), float(n_select))

    # shadowing: x_i + x_j <= 1
    n_forbidden = 0
    if forbidden_pairs is not None:
        arr = np.asarray(forbidden_pairs)
        pairs = (np.argwhere(np.triu(arr.astype(bool), 1))
                 if arr.shape == (n_pads, n_pads) else np.asarray(forbidden_pairs))
        for i, j in pairs:
            add_row([(int(i), 1.0), (int(j), 1.0)], -np.inf, 1.0)
            n_forbidden += 1

    A = coo_matrix((data, (rows, cols)), shape=(r, n_var))
    constraints = LinearConstraint(A.tocsr(), np.array(lb), np.array(ub))

    c_obj = np.zeros(n_var)
    c_obj[off_z:] = -1.0                      # maximise sum z  ->  minimise -sum z

    integrality = np.zeros(n_var)
    integrality[: off_z] = 1                  # x and y binary; z can stay continuous
    bounds = Bounds(np.zeros(n_var), np.ones(n_var))

    options = {"mip_rel_gap": mip_rel_gap}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)

    res = milp(c=c_obj, constraints=constraints, integrality=integrality,
               bounds=bounds, options=options)

    out = {
        "status": int(res.status),
        "message": str(res.message),
        "proven_optimal": bool(res.status == 0),
        "n_variables": n_var,
        "n_constraints": r,
        "n_x": n_pads,
        "n_y": n_pairs,
        "n_z": n_cells,
        "n_forbidden_pairs": n_forbidden,
    }
    if res.x is None:
        out.update({"selection": None, "coverage": None})
        return out

    x = np.asarray(res.x[:n_pads])
    selection = x > 0.5
    out["selection"] = np.flatnonzero(selection).tolist()
    out["coverage"] = int(round(-res.fun))
    out["mip_gap"] = float(getattr(res, "mip_gap", np.nan))

    # Even when optimality is not proven, the dual bound is a *certificate of a
    # ceiling*: no selection whatsoever can cover more than this. That is worth
    # reporting on its own -- a heuristic within 1 % of the dual bound is within
    # 1 % of optimal, proven, without the solver ever closing the gap.
    dual = getattr(res, "mip_dual_bound", None)
    out["upper_bound"] = float(-dual) if dual is not None else None
    if out["upper_bound"] is not None and out["coverage"] is not None:
        ub = out["upper_bound"]
        out["incumbent_gap_percent"] = (100.0 * (ub - out["coverage"]) / ub) if ub > 0 else 0.0
    return out
