"""Ready-made objectives for the optimisers.

Each factory returns a callable ``score(state) -> float`` that is **minimised**,
so rewards are negated. They are deliberately thin wrappers around the
quantities :class:`thz_opt.optimize.state.UVState` already tracks exactly, so
that a search evaluates the *true* objective, not a surrogate -- the surrogate
question of ``docs/qubo_mapping.md`` arises only when the problem is handed to
a QUBO solver, and keeping the two separate is what makes the comparison fair.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "max_unique_cells",
    "min_sidelobe_energy",
    "weighted_objective",
    "coherence_weighted_cells",
]


def max_unique_cells():
    """Maximise distinct UV cells touched (returns ``-unique_cells``)."""
    def score(state) -> float:
        return -float(state.unique_cells)
    score.name = "max_unique_cells"
    return score


def min_sidelobe_energy(normalise: bool = True):
    """Minimise ``sum_c n_c^2``, i.e. dirty-beam energy at fixed peak.

    With ``normalise`` the value is divided by ``(sum_c n_c)^2``, which is the
    quantity that is actually proportional to the unit-peak beam energy; without
    it, selections with more samples are unfairly penalised. Since every
    selection here has the same number of baselines, the two agree up to a
    constant -- the flag matters only if the number of active baselines varies.
    """
    def score(state) -> float:
        if not normalise:
            return float(state.sum_sq)
        total = float(state.occupancy().sum())
        return float(state.sum_sq) / (total * total) if total > 0 else np.inf
    score.name = "min_sidelobe_energy"
    return score


def weighted_objective(w_cells: float = 1.0, w_sidelobe: float = 0.0,
                       cell_scale: float = 1.0, sidelobe_scale: float = 1.0):
    """Linear combination of the two competing objectives.

    ``score = -w_cells * cells/cell_scale + w_sidelobe * sum_sq/sidelobe_scale``

    The scales exist because the two quantities differ by orders of magnitude;
    pass the values measured for a reference configuration so the weights mean
    what they look like. This is the explicit version of the trade-off in
    ``docs/objective_function.md`` eq. (5.4) -- there is no principled default,
    which is exactly why it is a parameter.
    """
    def score(state) -> float:
        return (-w_cells * state.unique_cells / cell_scale
                + w_sidelobe * state.sum_sq / sidelobe_scale)
    score.name = f"weighted(cells={w_cells},sidelobe={w_sidelobe})"
    return score


def coherence_weighted_cells(pair_gamma: np.ndarray, pair_cells: list,
                             n_pads: int, threshold: float = 0.0):
    """Maximise coverage counting only baselines that survive decorrelation.

    ``pair_gamma[k]`` is the coherence factor of pair ``k`` from
    :func:`thz_opt.constraints.coherence.decorrelation_factor`. Cells are
    credited in proportion to the coherence of the baseline that samples them,
    so a long baseline that the atmosphere destroys contributes little even
    though it technically touches new UV cells.

    This is the THz-specific objective: the same UV coverage is worth
    different amounts depending on whether the phase survives. It cannot use
    the incremental state (the weighting is per-baseline, not per-cell), so it
    is slower -- use it for final scoring or small problems.
    """
    from ..interferometry.baselines import pair_indices

    gamma = np.asarray(pair_gamma, dtype=float)
    i_idx, j_idx = pair_indices(n_pads)

    def score(state) -> float:
        sel = state.selected
        best: dict = {}
        for k in range(len(pair_cells)):
            if not (sel[i_idx[k]] and sel[j_idx[k]]):
                continue
            g = gamma[k]
            if g <= threshold:
                continue
            for c in pair_cells[k].tolist():
                if g > best.get(c, 0.0):
                    best[c] = g
        return -float(sum(best.values()))

    score.name = "coherence_weighted_cells"
    return score
