"""Mandatory checks before any solver is allowed near the problem.

Three independent things are verified:

1. :func:`validate_qubo` -- the matrix form of the objective reproduces the
   term-by-term definition for *every* configuration, not just the optimum.
2. :func:`check_convention_roundtrip` -- the upper-triangular convention
   ``sum_i Q_ii x_i + sum_{i<j} Q_ij x_i x_j`` and the symmetric convention
   ``x^T A x`` give identical energies (the factor-of-two trap).
3. :func:`compare_surrogate_to_exact` -- how well a pairwise UV surrogate
   ranks selections against the exact, non-quadratic ``|Cov(x)|``.  This one is
   a *measurement*, not a pass/fail test: the surrogate is known to be an
   approximation and the point is to quantify by how much.
"""

from __future__ import annotations

import numpy as np

from .coefficients import exact_cells_covered, surrogate_value
from .exhaustive import enumerate_energies, feasible_configurations
from .objective import QUBOTerms, build_qubo, qubo_energy, symmetric_energy, to_symmetric


def validate_qubo(terms: QUBOTerms, tol: float = 1e-9, feasible_only: bool = False) -> dict:
    """Assert ``direct_objective == qubo_energy`` over all configurations."""
    Q, offset = build_qubo(terms)
    configs, direct, quad = enumerate_energies(terms, Q, offset, feasible_only=feasible_only)
    diff = np.abs(direct - quad)
    k = int(np.argmax(diff))
    return {
        "n_configurations": int(configs.shape[0]),
        "max_abs_difference": float(diff.max()),
        "within_tolerance": bool(diff.max() <= tol),
        "tolerance": tol,
        "worst_config": "".join(str(int(b)) for b in configs[k]),
        "worst_direct": float(direct[k]),
        "worst_qubo": float(quad[k]),
    }


def check_convention_roundtrip(Q: np.ndarray, offset: float = 0.0, n_trials: int = 64,
                               seed: int = 0, tol: float = 1e-9) -> dict:
    """Upper-triangular and symmetric conventions must agree on random inputs."""
    rng = np.random.default_rng(seed)
    n = Q.shape[0]
    A = to_symmetric(Q)
    worst = 0.0
    for _ in range(n_trials):
        x = rng.integers(0, 2, n)
        worst = max(worst, abs(qubo_energy(x, Q, offset) - symmetric_energy(x, A, offset)))
    return {"max_abs_difference": float(worst), "within_tolerance": bool(worst <= tol),
            "n_trials": n_trials}


def compare_surrogate_to_exact(cell_sets: list, w: np.ndarray, n_pads: int,
                               n_select: int) -> dict:
    """Quantify the pairwise surrogate against exact unique-cell coverage.

    Over all ``C(M, N)`` feasible selections, computes

    * the Pearson and Spearman correlation between ``sum w_ij x_i x_j`` and
      ``|Cov(x)|``,
    * whether the surrogate's argmax is also an exact argmax,
    * the coverage shortfall of the surrogate's best selection relative to the
      true best, in cells and in per cent,
    * the mean relative over-count ``(surrogate - exact) / exact``.

    A high correlation is evidence that the surrogate is usable as an optimiser
    objective; it is NOT evidence that the surrogate equals the coverage.
    """
    configs = feasible_configurations(n_pads, n_select)
    exact = np.array([exact_cells_covered(c, cell_sets, n_pads) for c in configs], dtype=float)
    surr = np.array([surrogate_value(c, w, n_pads) for c in configs], dtype=float)

    degenerate_exact = exact.std() == 0.0
    degenerate_surrogate = surr.std() == 0.0
    degenerate = degenerate_exact or degenerate_surrogate

    def _corr(a, b):
        """Correlation, or NaN when either side is constant.

        A constant side means every feasible selection scores the same, so the
        surrogate is neither right nor wrong -- reporting 1.0 there (which the
        rank transform would do) would be an artefact.
        """
        if a.std() == 0.0 or b.std() == 0.0:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    def _spearman(a, b):
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        return _corr(ra, rb) if not degenerate else float("nan")

    best_exact = float(exact.max())
    k_surr = int(np.argmax(surr))
    cov_at_surr_best = float(exact[k_surr])
    rel = (surr - exact) / np.maximum(exact, 1.0)

    return {
        "n_configurations": int(configs.shape[0]),
        "degenerate": bool(degenerate),
        "degenerate_exact": bool(degenerate_exact),
        "degenerate_surrogate": bool(degenerate_surrogate),
        "pearson_r": _corr(surr, exact),
        "spearman_rho": _spearman(surr, exact),
        "surrogate_argmax_is_exact_argmax": bool(cov_at_surr_best >= best_exact - 1e-12),
        "best_exact_cells": best_exact,
        "cells_at_surrogate_optimum": cov_at_surr_best,
        "coverage_shortfall_cells": best_exact - cov_at_surr_best,
        "coverage_shortfall_percent": 100.0 * (best_exact - cov_at_surr_best) / best_exact
        if best_exact > 0 else 0.0,
        "mean_relative_overcount": float(rel.mean()),
        "max_relative_overcount": float(rel.max()),
    }
