"""The QUBO objective: direct form, matrix form, and an exact (large) variant.

Convention (used by every function in this package)
---------------------------------------------------
The QUBO is stored as an **upper-triangular** matrix ``Q`` of shape ``(M, M)``
plus a scalar ``offset``, and its energy is

    E(x) = sum_i Q_ii x_i  +  sum_{i<j} Q_ij x_i x_j  +  offset          (1)

i.e. each unordered pair appears exactly ONCE, with the whole coefficient on
``Q_ij``, ``i < j``.  The alternative convention

    E(x) = x^T A x + offset                                              (2)

puts half of each pair coefficient on ``A_ij`` and half on ``A_ji``.  The two
are related by ``A = (Q + Q^T)/2`` off the diagonal and ``A_ii = Q_ii`` on it,
which :func:`to_symmetric` performs; :func:`to_upper_triangular` inverts it.
Getting this factor of two wrong is the classic QUBO bug, so
:func:`thz_opt.qubo.validation.check_convention_roundtrip` asserts that (1) and
(2) agree for random ``x``.

Terms
-----
::

    H_total = H_select + H_shadow + H_pwv + H_uv + H_phase

    H_select = lambda_select * (sum_i x_i - N)^2
             = lambda_select * [ (1 - 2N) sum_i x_i + 2 sum_{i<j} x_i x_j + N^2 ]

    H_shadow = lambda_shadow * sum_{i<j, d_ij < d_min} x_i x_j

    H_pwv    = lambda_pwv * sum_i c_i x_i          (synthetic cost, see constraints.pwv)

    H_uv     = - sum_{i<j} w_ij x_i x_j            (w_ij >= 0, from qubo.coefficients)

    H_phase  = lambda_phase * sum_{i<j} P_ij x_i x_j   (connectivity surrogate)

The expansion of ``H_select`` uses ``x_i^2 = x_i`` for binary variables, which
is what turns the square into the linear coefficient ``(1 - 2N)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..interferometry.baselines import pair_indices


@dataclass
class QUBOTerms:
    """Everything needed to build the objective for one candidate pad set.

    ``w_uv`` is indexed by ``i < j`` pair order (see
    :func:`thz_opt.interferometry.baselines.pair_indices`); ``shadow`` and
    ``phase`` are full ``(M, M)`` matrices of which only the upper triangle is
    read.
    """

    n_pads: int
    n_select: int
    w_uv: np.ndarray | None = None
    shadow: np.ndarray | None = None
    phase: np.ndarray | None = None
    pwv_cost: np.ndarray | None = None
    lambda_select: float = 1.0
    lambda_shadow: float = 1.0
    lambda_phase: float = 0.0
    lambda_pwv: float = 0.0
    metadata: dict = field(default_factory=dict)


def _pair_matrix(values: np.ndarray, n: int) -> np.ndarray:
    """Scatter a pair-ordered vector into an upper-triangular ``(n, n)`` matrix."""
    i, j = pair_indices(n)
    m = np.zeros((n, n))
    m[i, j] = np.asarray(values, dtype=float)
    return m


def build_qubo(terms: QUBOTerms):
    """Return ``(Q, offset)`` in the upper-triangular convention of eq. (1)."""
    n, k = terms.n_pads, terms.n_select
    Q = np.zeros((n, n))

    # H_select
    lam = terms.lambda_select
    np.fill_diagonal(Q, Q.diagonal() + lam * (1.0 - 2.0 * k))
    iu = np.triu_indices(n, k=1)
    Q[iu] += 2.0 * lam
    offset = lam * k * k

    # H_shadow
    if terms.shadow is not None:
        s = np.asarray(terms.shadow, dtype=float)
        Q[iu] += terms.lambda_shadow * s[iu]

    # H_phase
    if terms.phase is not None:
        p = np.asarray(terms.phase, dtype=float)
        Q[iu] += terms.lambda_phase * p[iu]

    # H_uv (reward, hence negative)
    if terms.w_uv is not None:
        Q[iu] -= _pair_matrix(terms.w_uv, n)[iu]

    # H_pwv (linear)
    if terms.pwv_cost is not None:
        Q[np.diag_indices(n)] += terms.lambda_pwv * np.asarray(terms.pwv_cost, dtype=float)

    return Q, float(offset)


def qubo_energy(x: np.ndarray, Q: np.ndarray, offset: float = 0.0) -> float:
    """Energy of eq. (1) for an upper-triangular ``Q`` and a binary ``x``.

    For binary variables ``x_i^2 = x_i``, so ``x^T triu(Q) x`` already gives
    ``sum_i Q_ii x_i + sum_{i<j} Q_ij x_i x_j``.
    """
    v = np.asarray(x, dtype=float)
    if not np.all((v == 0) | (v == 1)):
        raise ValueError("qubo_energy expects a binary vector")
    return float(v @ np.triu(np.asarray(Q, dtype=float)) @ v + offset)


def to_symmetric(Q: np.ndarray) -> np.ndarray:
    """Upper-triangular ``Q`` -> symmetric ``A`` with ``E = x^T A x`` (eq. 2)."""
    q = np.asarray(Q, dtype=float)
    d = np.diag(np.diag(q))
    off = q - d
    return d + 0.5 * (off + off.T)


def to_upper_triangular(A: np.ndarray) -> np.ndarray:
    """Symmetric ``A`` -> the upper-triangular ``Q`` of eq. (1)."""
    a = np.asarray(A, dtype=float)
    d = np.diag(np.diag(a))
    off = a - d
    return d + np.triu(off + off.T, k=1)


def symmetric_energy(x: np.ndarray, A: np.ndarray, offset: float = 0.0) -> float:
    """Energy of eq. (2)."""
    v = np.asarray(x, dtype=float)
    return float(v @ np.asarray(A, dtype=float) @ v + offset)


def to_qubo_dict(Q: np.ndarray, tol: float = 0.0) -> dict:
    """``{(i, j): coefficient}`` in the convention of eq. (1).

    This is the form accepted by ``dimod.BinaryQuadraticModel.from_qubo`` and by
    D-Wave samplers, which use exactly the same one-entry-per-pair convention.
    """
    q = np.asarray(Q, dtype=float)
    out = {}
    for i in range(q.shape[0]):
        for j in range(i, q.shape[1]):
            if abs(q[i, j]) > tol:
                out[(i, j)] = float(q[i, j])
    return out


# --------------------------------------------------------------------------
# Direct (definition-following) objective -- the reference implementation the
# matrix form is validated against.
# --------------------------------------------------------------------------

def direct_terms(x: np.ndarray, terms: QUBOTerms) -> dict:
    """Each ``H`` term computed straight from its mathematical definition.

    Deliberately written with loops and explicit sums rather than matrix
    algebra, so that it cannot share a bug with :func:`build_qubo`.
    """
    v = np.asarray(x, dtype=int)
    n = terms.n_pads
    sel = np.flatnonzero(v)

    h_select = terms.lambda_select * (int(v.sum()) - terms.n_select) ** 2

    h_shadow = 0.0
    if terms.shadow is not None:
        for a in range(len(sel)):
            for b in range(a + 1, len(sel)):
                h_shadow += terms.lambda_shadow * float(terms.shadow[sel[a], sel[b]])

    h_phase = 0.0
    if terms.phase is not None:
        for a in range(len(sel)):
            for b in range(a + 1, len(sel)):
                h_phase += terms.lambda_phase * float(terms.phase[sel[a], sel[b]])

    h_uv = 0.0
    if terms.w_uv is not None:
        i_idx, j_idx = pair_indices(n)
        for k in range(len(i_idx)):
            if v[i_idx[k]] and v[j_idx[k]]:
                h_uv -= float(terms.w_uv[k])

    h_pwv = 0.0
    if terms.pwv_cost is not None:
        h_pwv = terms.lambda_pwv * float(np.sum(np.asarray(terms.pwv_cost)[sel]))

    total = h_select + h_shadow + h_phase + h_uv + h_pwv
    return {
        "H_select": float(h_select),
        "H_shadow": float(h_shadow),
        "H_phase": float(h_phase),
        "H_uv": float(h_uv),
        "H_pwv": float(h_pwv),
        "H_total": float(total),
    }


def direct_objective(x: np.ndarray, terms: QUBOTerms) -> float:
    """``H_total`` from the definitions (see :func:`direct_terms`)."""
    return direct_terms(x, terms)["H_total"]


# --------------------------------------------------------------------------
# Exact coverage QUBO with auxiliary variables (approach B of docs/qubo_mapping.md)
# --------------------------------------------------------------------------

def build_qubo_with_cell_variables(
    cell_sets: list,
    n_pads: int,
    n_select: int,
    lambda_select: float,
    lambda_aux: float,
    max_variables: int = 4096,
):
    """Exact ``-|Cov(x)|`` QUBO using auxiliary pair and cell variables.

    Variable blocks, in index order:

    1. ``x_i``      pad selection, ``M`` variables
    2. ``y_k``      pair activation, one per candidate pair with a non-empty
       cell set, enforced to equal ``x_i x_j`` by the Rosenberg penalty
       ``lambda_aux * (x_i x_j - 2 x_i y - 2 x_j y + 3 y)``, which is zero iff
       ``y = x_i x_j`` and positive otherwise
    3. ``z_c``      cell coverage, one per cell touched by at least one pair
    4. ``s_c,b``    binary-encoded slack enforcing ``z_c <= sum_k y_k`` through
       ``lambda_aux * (sum_k y_k - z_c - s_c)^2``

    The reward is ``-sum_c z_c``, so the ground state covers as many distinct
    cells as the selection constraint allows -- exactly, with no pairwise
    approximation.  The price is the variable count, which is reported in the
    returned metadata and is why this construction is only practical for toy
    problems; ``max_variables`` guards against building something enormous by
    accident.

    Returns ``(Q, offset, meta)`` in the convention of eq. (1).
    """
    i_idx, j_idx = pair_indices(n_pads)
    pairs = [(int(i_idx[k]), int(j_idx[k]), cell_sets[k]) for k in range(len(cell_sets))
             if len(cell_sets[k]) > 0]

    cells = sorted({int(c) for _, _, cs in pairs for c in cs.tolist()})
    cell_to_pairs = {c: [] for c in cells}
    for p, (_, _, cs) in enumerate(pairs):
        for c in cs.tolist():
            cell_to_pairs[int(c)].append(p)

    n_x = n_pads
    n_y = len(pairs)
    n_z = len(cells)
    slack_bits = [max(int(np.ceil(np.log2(len(cell_to_pairs[c]) + 1))), 1) for c in cells]
    n_s = int(sum(slack_bits))
    total_vars = n_x + n_y + n_z + n_s
    if total_vars > max_variables:
        raise ValueError(
            f"exact construction needs {total_vars} variables (> max_variables={max_variables})"
        )

    off_y = n_x
    off_z = off_y + n_y
    off_s = off_z + n_z

    Q = np.zeros((total_vars, total_vars))
    offset = 0.0

    def add(i: int, j: int, val: float) -> None:
        a, b = (i, j) if i <= j else (j, i)
        Q[a, b] += val

    # 1. selection constraint on x
    for i in range(n_x):
        add(i, i, lambda_select * (1.0 - 2.0 * n_select))
        for j in range(i + 1, n_x):
            add(i, j, 2.0 * lambda_select)
    offset += lambda_select * n_select * n_select

    # 2. y_k = x_i x_j  (Rosenberg)
    for p, (i, j, _) in enumerate(pairs):
        y = off_y + p
        add(i, j, lambda_aux)
        add(i, y, -2.0 * lambda_aux)
        add(j, y, -2.0 * lambda_aux)
        add(y, y, 3.0 * lambda_aux)

    # 3. reward + 4. slack equality  (sum_k y_k - z_c - s_c)^2
    s_cursor = off_s
    for ci, c in enumerate(cells):
        z = off_z + ci
        add(z, z, -1.0)  # reward one unit per covered cell

        bits = slack_bits[ci]
        coeffs = [(off_y + p, 1.0) for p in cell_to_pairs[c]]
        coeffs.append((z, -1.0))
        coeffs.extend((s_cursor + b, -(2.0 ** b)) for b in range(bits))
        s_cursor += bits

        for a, (va, ca) in enumerate(coeffs):
            add(va, va, lambda_aux * ca * ca)
            for vb, cb in coeffs[a + 1 :]:
                add(va, vb, 2.0 * lambda_aux * ca * cb)

    meta = {
        "n_variables": total_vars,
        "n_x": n_x,
        "n_y": n_y,
        "n_z": n_z,
        "n_slack": n_s,
        "n_cells": n_z,
        "lambda_select": lambda_select,
        "lambda_aux": lambda_aux,
        "index_blocks": {"x": (0, n_x), "y": (off_y, off_y + n_y),
                         "z": (off_z, off_z + n_z), "s": (off_s, total_vars)},
        "pairs": [(i, j) for i, j, _ in pairs],
        "cells": cells,
        "slack_bits": slack_bits,
        "cell_to_pairs": cell_to_pairs,
    }
    return Q, float(offset), meta


def complete_auxiliary(x: np.ndarray, meta: dict) -> np.ndarray:
    """Optimal auxiliary assignment for a given ``x`` in the exact construction.

    Sets ``y_k = x_i x_j``, ``z_c = OR_k y_k`` and the slack to
    ``sum_k y_k - z_c``, which is the unique assignment satisfying every
    penalty; the resulting full vector therefore has energy
    ``H_select(x) - |Cov(x)|``.
    """
    v = np.asarray(x, dtype=int)
    y = np.array([v[i] * v[j] for i, j in meta["pairs"]], dtype=int)
    z, s_bits = [], []
    for c, bits in zip(meta["cells"], meta["slack_bits"]):
        total = int(sum(y[p] for p in meta["cell_to_pairs"][c]))
        zc = 1 if total > 0 else 0
        z.append(zc)
        slack = total - zc
        s_bits.extend((slack >> b) & 1 for b in range(bits))
    return np.concatenate([v, y, np.array(z, dtype=int), np.array(s_bits, dtype=int)])
