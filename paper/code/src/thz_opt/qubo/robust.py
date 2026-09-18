"""Robustness over scenarios -- the part that is quadratic.

Weather changes. So does the residual after phase referencing, and the
declination you actually get scheduled. A layout that is optimal for one
atmosphere and poor for the rest is not the useful answer.

The **mean** over a finite scenario set is a weighted sum of QUBOs, so it stays
exactly quadratic and costs nothing:

    E_theta[H] = sum_s p_s H(x; theta_s)

The **variance** is not. ``Var[H] = E[H^2] - E[H]^2`` and ``H`` is already
quadratic in the baseline variables, so ``H^2`` is quartic in ``y`` (degree 8
in ``x``). Mean-variance robustness is therefore not a QUBO, and no amount of
auxiliary variables makes it cheap. What is cheap and honest is to optimise the
mean and *report* the spread afterwards, which :func:`scenario_spread` does.

``ponytail: mean only. CVaR needs an epigraph variable plus binary-encoded
slacks per scenario -- add it if a solver run shows the mean-optimal layout has
an unacceptable tail.``
"""

from __future__ import annotations

import numpy as np

from .baseline_qubo import BaselineTerms, objective_value

__all__ = ["scenario_mean", "scenario_spread"]


def _normalised(weights, n: int) -> np.ndarray:
    if weights is None:
        return np.full(n, 1.0 / n)
    w = np.asarray(weights, dtype=float)
    if w.shape[0] != n or np.any(w < 0) or w.sum() <= 0:
        raise ValueError("weights must be one non-negative entry per scenario, not all zero")
    return w / w.sum()


def scenario_mean(terms_list: list, weights=None) -> BaselineTerms:
    """Probability-weighted mean of several scenarios' objectives.

    Every scenario must cover the same pads and pairs; the scenarios differ in
    their coefficients (different coherence, declination, weather), not in the
    decision variables.
    """
    if not terms_list:
        raise ValueError("need at least one scenario")
    p = _normalised(weights, len(terms_list))
    n_pads, n_pairs = terms_list[0].n_pads, terms_list[0].n_pairs
    if any(t.n_pads != n_pads or t.n_pairs != n_pairs for t in terms_list):
        raise ValueError("all scenarios must share the same pad and pair count")

    lin = np.zeros(n_pairs)
    quad: dict = {}
    offset = 0.0
    for ps, t in zip(p, terms_list):
        lin += ps * np.asarray(t.linear, dtype=float)
        offset += ps * t.offset
        for kl, c in t.quadratic.items():
            quad[kl] = quad.get(kl, 0.0) + ps * c

    return BaselineTerms(n_pads=n_pads, linear=lin, quadratic=quad, offset=offset,
                         name=f"scenario_mean({len(terms_list)})")


def scenario_spread(y: np.ndarray, terms_list: list, weights=None) -> dict:
    """How a given selection performs across every scenario.

    Reported after optimisation, because the spread is what the mean hides: two
    layouts with the same expected objective can differ completely in their
    worst case.
    """
    p = _normalised(weights, len(terms_list))
    vals = np.array([objective_value(y, t) for t in terms_list], dtype=float)
    mean = float(vals @ p)
    var = float(((vals - mean) ** 2) @ p)
    return {
        "n_scenarios": int(vals.size),
        "mean": mean,
        "std": float(np.sqrt(var)),
        "min": float(vals.min()),
        "max": float(vals.max()),
        "worst_case": float(vals.max()),          # objectives here are minimised
        "values": vals.tolist(),
    }
