"""Brute-force enumeration for small candidate sets.

For ``M`` pads there are ``2**M`` binary configurations; this module enumerates
them so that, for toy sizes, the exact optimum is known and no solver -- quantum
or classical -- has to be trusted.  ``M = 20`` is already a million
configurations, so a guard is applied.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

from .objective import QUBOTerms, direct_objective, qubo_energy


def all_configurations(n_pads: int, max_pads: int = 22) -> np.ndarray:
    """Every binary vector of length ``n_pads``, as a ``(2**n_pads, n_pads)`` array."""
    if n_pads > max_pads:
        raise ValueError(f"refusing to enumerate 2**{n_pads} configurations")
    idx = np.arange(2 ** n_pads, dtype=np.uint64)
    bits = ((idx[:, None] >> np.arange(n_pads, dtype=np.uint64)[None, :]) & np.uint64(1))
    return bits.astype(np.int8)


def feasible_configurations(n_pads: int, n_select: int) -> np.ndarray:
    """Only the ``C(M, N)`` vectors with exactly ``n_select`` ones."""
    rows = []
    for combo in combinations(range(n_pads), n_select):
        v = np.zeros(n_pads, dtype=np.int8)
        v[list(combo)] = 1
        rows.append(v)
    return np.array(rows, dtype=np.int8)


def enumerate_energies(terms: QUBOTerms, Q: np.ndarray, offset: float, feasible_only: bool = False):
    """Energies of every configuration under both implementations.

    Returns ``(configs, direct, qubo)`` where ``direct`` follows the term
    definitions and ``qubo`` uses the matrix form.  They must agree; that is
    what :func:`thz_opt.qubo.validation.validate_qubo` checks.
    """
    configs = (
        feasible_configurations(terms.n_pads, terms.n_select)
        if feasible_only
        else all_configurations(terms.n_pads)
    )
    direct = np.array([direct_objective(c, terms) for c in configs])
    quad = np.array([qubo_energy(c, Q, offset) for c in configs])
    return configs, direct, quad


def best_configuration(configs: np.ndarray, energies: np.ndarray):
    """``(config, energy, index)`` of the lowest-energy configuration."""
    k = int(np.argmin(energies))
    return configs[k], float(energies[k]), k


def rank_table(configs: np.ndarray, energies: np.ndarray, top: int = 10):
    """The ``top`` lowest-energy configurations as ``(bitstring, energy)`` rows."""
    order = np.argsort(energies)[:top]
    return [("".join(str(int(b)) for b in configs[k]), float(energies[k])) for k in order]
