"""From UV physics to pairwise QUBO coefficients.

This module is where the central research question of the repository lives:
how much of a *global* UV-coverage objective can honestly be pushed into
*pairwise* coefficients ``w_ij`` that a QUBO can hold.

The offline/online split
------------------------
Everything here is computed ONCE, offline, from the candidate pad geometry:
baselines, their UV tracks, and which grid cells they touch.  The optimiser
never runs an FFT or a forward simulation; it only sees the resulting numbers.

The exact objective and why it is not quadratic
-----------------------------------------------
Let ``C(i, j)`` be the set of UV cells touched by the baseline between pads
``i`` and ``j``.  For a selection ``x`` the set of covered cells is

    Cov(x) = union over {i<j : x_i = x_j = 1} of C(i, j)

and the natural coverage objective is ``|Cov(x)|``.  Writing ``y_ij = x_i x_j``
and ``z_c = OR over pairs covering c of y_ij``, the OR is what breaks
quadraticity: ``z_c = 1 - prod (1 - y_ij)`` expands into products of arbitrarily
many ``y`` terms, i.e. into monomials of degree up to ``2 * (number of pairs
covering c)`` in ``x``.  So ``|Cov(x)|`` is in general a higher-order
pseudo-Boolean function, NOT a quadratic form -- unless every cell is touched by
at most one candidate pair, in which case ``|Cov(x)|`` equals
``sum_{i<j} |C(i,j)| x_i x_j`` exactly.  That special case is checked by
:func:`coverage_is_exactly_quadratic`.

Three pairwise weightings are therefore offered, each a different compromise:

``"count"``
    ``w_ij = |C(i, j)|``.  Equals the exact coverage only in the no-overlap
    case; otherwise it systematically over-counts, because a cell touched by
    two selected baselines is paid for twice.  It is an upper bound:
    ``sum w_ij y_ij >= |Cov(x)|``.

``"shared"``
    ``w_ij = sum over c in C(i,j) of 1 / m_c`` where ``m_c`` is the number of
    *candidate* pairs that touch cell ``c``.  Credit for a contested cell is
    split among all pairs that could supply it.  This is exact when no cell is
    contested; in general it is neither an upper nor a lower bound, but it
    removes the systematic double counting of ``"count"`` and tracks
    ``|Cov(x)|`` far more closely (measured in experiment 05).

``"radial"``
    ``"shared"`` with each cell additionally weighted by
    ``1 + r_c / r_ref``, so that filling the sparsely sampled outer UV plane is
    worth more than adding another short baseline.  This is an explicit *design
    preference*, not a physical result.

A fourth route -- auxiliary binary variables ``z_c`` with the standard
OR-linearisation, giving an exact but much larger QUBO -- is implemented in
:mod:`thz_opt.qubo.objective` as :func:`~thz_opt.qubo.objective.build_qubo_with_cell_variables`
and its cost is discussed in ``docs/qubo_mapping.md``.
"""

from __future__ import annotations

import numpy as np

from ..interferometry.baselines import compute_baselines, pair_indices
from ..interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from ..interferometry.uv import UVGrid, cell_ids

WEIGHTINGS = ("count", "shared", "radial")


def baseline_cell_sets(
    pads: np.ndarray,
    wavelength: float,
    grid: UVGrid,
    observation: ObservationConfig | None = None,
    include_conjugate: bool = True,
) -> list:
    """Cells touched by each candidate pair, in ``i < j`` order.

    Returns a list of ``N(N-1)/2`` arrays of unique flat cell ids.  With
    ``observation=None`` the snapshot model is used (each pair touches one cell,
    or two with its conjugate); with an :class:`ObservationConfig` the pair
    touches every cell along its Earth-rotation track.
    """
    pads = np.asarray(pads, dtype=float)
    n = pads.shape[0]
    i_idx, j_idx = pair_indices(n)
    out = []

    if observation is None:
        b, _, _ = compute_baselines(pads)
        for k in range(b.shape[0]):
            uv = b[k : k + 1] / wavelength
            if include_conjugate:
                uv = np.vstack((uv, -uv))
            out.append(np.unique(cell_ids(uv, grid)))
    else:
        for a, bidx in zip(i_idx, j_idx):
            uv = layout_to_uv_tracks(pads[[a, bidx]], observation, include_conjugate)
            out.append(np.unique(cell_ids(uv, grid)))
    return out


def cell_claim_counts(cell_sets: list) -> dict:
    """``m_c`` -- how many candidate pairs touch each cell."""
    counts: dict = {}
    for cells in cell_sets:
        for c in cells.tolist():
            counts[c] = counts.get(c, 0) + 1
    return counts


def coverage_is_exactly_quadratic(cell_sets: list) -> bool:
    """``True`` iff no UV cell is touched by more than one candidate pair.

    In that -- and only that -- case ``|Cov(x)|`` is exactly
    ``sum_{i<j} |C(i,j)| x_i x_j``, so the ``"count"`` weighting is not an
    approximation.
    """
    return all(v <= 1 for v in cell_claim_counts(cell_sets).values())


def cell_radius(cell: int, grid: UVGrid) -> float:
    """UV radius in wavelengths of the centre of a flat cell id."""
    iu, iv = divmod(int(cell), grid.n_cells)
    cu = (iu - grid.n_cells / 2.0) * grid.cell_size
    cv = (iv - grid.n_cells / 2.0) * grid.cell_size
    return float(np.hypot(cu, cv))


def pairwise_uv_weights(
    cell_sets: list,
    grid: UVGrid,
    weighting: str = "shared",
    r_ref: float | None = None,
) -> np.ndarray:
    """Pairwise UV utilities ``w_ij >= 0``, one per ``i < j`` pair.

    See the module docstring for the definition and the bias of each
    ``weighting``.  ``r_ref`` (wavelengths) sets the radial weighting scale for
    ``"radial"``; ``None`` uses ``grid.uv_max / 4``.
    """
    if weighting not in WEIGHTINGS:
        raise ValueError("weighting must be one of " + str(WEIGHTINGS))

    if weighting == "count":
        return np.array([float(len(c)) for c in cell_sets])

    claims = cell_claim_counts(cell_sets)
    ref = grid.uv_max / 4.0 if r_ref is None else float(r_ref)

    w = np.empty(len(cell_sets))
    for k, cells in enumerate(cell_sets):
        total = 0.0
        for c in cells.tolist():
            share = 1.0 / claims[c]
            if weighting == "radial":
                share *= 1.0 + cell_radius(c, grid) / ref
            total += share
        w[k] = total
    return w


def exact_cells_covered(selection: np.ndarray, cell_sets: list, n_pads: int) -> int:
    """``|Cov(x)|`` -- the exact number of distinct UV cells covered.

    This is the non-quadratic ground truth every surrogate is measured against.
    """
    sel = np.asarray(selection, dtype=bool)
    i_idx, j_idx = pair_indices(n_pads)
    covered = set()
    for k, (a, b) in enumerate(zip(i_idx, j_idx)):
        if sel[a] and sel[b]:
            covered.update(cell_sets[k].tolist())
    return len(covered)


def surrogate_value(selection: np.ndarray, w: np.ndarray, n_pads: int) -> float:
    """``sum_{i<j} w_ij x_i x_j`` -- the quadratic surrogate of coverage."""
    sel = np.asarray(selection, dtype=float)
    i_idx, j_idx = pair_indices(n_pads)
    return float(np.sum(np.asarray(w, dtype=float) * sel[i_idx] * sel[j_idx]))


def suggest_penalty_weights(
    w_uv: np.ndarray,
    n_select: int,
    shadow_pairs_count: int = 0,
    safety: float = 2.0,
) -> dict:
    """A documented, non-arbitrary rule for the penalty weights.

    The objective term ``H_uv = -sum w_ij x_i x_j`` can gain at most

        U_max = sum of the ``n_select*(n_select-1)/2`` largest ``w_ij``

    from any selection of ``n_select`` pads.  A constraint penalty is *safe*
    (cannot be bought off by UV gain) when violating it by one unit costs more
    than ``U_max``.  Adding one antenna beyond ``N`` changes ``H_select`` by at
    least ``lambda_select`` and adding one shadowed pair changes ``H_shadow`` by
    ``lambda_shadow``, so

        lambda_select = lambda_shadow = safety * U_max,   safety > 1

    These are *sufficient* bounds and deliberately conservative: smaller values
    often work and give a better-conditioned energy landscape, which is why the
    rule and its inputs are returned alongside the numbers rather than hidden.
    """
    w = np.sort(np.asarray(w_uv, dtype=float))[::-1]
    k = max(n_select * (n_select - 1) // 2, 1)
    u_max = float(w[:k].sum())
    lam = safety * max(u_max, 1e-12)
    return {
        "lambda_select": lam,
        "lambda_shadow": lam,
        "lambda_phase": lam,
        "lambda_pwv": 0.0,
        "u_max": u_max,
        "safety_factor": safety,
        "rule": "lambda = safety * (sum of the N(N-1)/2 largest w_ij)",
        "n_shadow_pairs": int(shadow_pairs_count),
    }
