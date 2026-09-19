"""Tests for actually solving the QUBO with a sampler.

These matter more than they look. Every earlier QUBO test checks that the
*formulation* is right -- that the energy of a known assignment equals the
objective. None of them checks that a solver handed the QUBO comes back with a
feasible answer, which is a different claim: constraints are enforced by
penalty, so an insufficient penalty weight produces a QUBO whose ground state
is infeasible, and every energy reported from it is meaningless.
"""

import itertools

import numpy as np
import pytest

dimod = pytest.importorskip("dimod")
samplers = pytest.importorskip("dwave.samplers")

from thz_opt.interferometry.baselines import pair_indices
from thz_opt.qubo.baseline_qubo import (
    BaselineTerms,
    bonferroni_coverage_terms,
    build_baseline_qubo,
    complete_baseline_auxiliary,
    objective_value,
    selection_to_y,
)
from thz_opt.qubo.objective import qubo_energy


def _toy_cell_sets(m_pads, rng, n_cells=40):
    """One small random set of UV cells per candidate pair."""
    p = m_pads * (m_pads - 1) // 2
    return [np.unique(rng.integers(0, n_cells, size=rng.integers(3, 9)))
            for _ in range(p)]


def _bqm(Q, offset):
    b = dimod.BinaryQuadraticModel(Q, "BINARY")
    b.offset = offset
    return b


def _decode(sample, m_pads):
    v = np.array([sample[i] for i in range(len(sample))], dtype=int)
    return v[:m_pads], v[m_pads:]


def test_dimod_energy_agrees_with_our_convention():
    """Our upper-triangular Q and dimod's BQM must mean the same thing.

    A silent factor of two here would make every solver result wrong while
    every formulation test still passed.
    """
    rng = np.random.default_rng(0)
    m = 6
    terms = bonferroni_coverage_terms(_toy_cell_sets(m, rng), m)
    Q, off, _ = build_baseline_qubo(terms, n_select=3)
    bqm = _bqm(Q, off)
    for _ in range(20):
        x = np.zeros(m, dtype=int)
        x[rng.choice(m, size=3, replace=False)] = 1
        full = complete_baseline_auxiliary(x, m)
        assert bqm.energy({i: int(b) for i, b in enumerate(full)}) == pytest.approx(
            qubo_energy(full, Q, off))


def test_sampler_returns_feasible_assignments():
    """The penalty floors must hold against an actual solver, not in theory.

    ``selection_penalty_floor`` and ``rosenberg_penalty_floor`` are derived
    bounds; this is the empirical check that a sampler cannot find a cheaper
    infeasible state.
    """
    rng = np.random.default_rng(1)
    m, n = 8, 4
    terms = bonferroni_coverage_terms(_toy_cell_sets(m, rng), m)
    Q, off, _ = build_baseline_qubo(terms, n_select=n)
    ss = samplers.SimulatedAnnealingSampler().sample(
        _bqm(Q, off), num_reads=200, seed=7)
    i_idx, j_idx = pair_indices(m)
    for sample in ss.samples():
        x, y = _decode(sample, m)
        assert x.sum() == n
        assert np.array_equal(y, x[i_idx] * x[j_idx])


def test_ground_state_is_the_true_qubo_optimum():
    """Exhaustive check that the QUBO's minimum sits on the right selection."""
    rng = np.random.default_rng(2)
    m, n = 8, 4
    cell_sets = _toy_cell_sets(m, rng)
    terms = bonferroni_coverage_terms(cell_sets, m)
    Q, off, _ = build_baseline_qubo(terms, n_select=n)

    best_x, best_obj = None, np.inf
    for combo in itertools.combinations(range(m), n):
        x = np.zeros(m, dtype=int)
        x[list(combo)] = 1
        obj = objective_value(selection_to_y(x, m), terms)
        if obj < best_obj:
            best_x, best_obj = x, obj

    ss = samplers.SimulatedAnnealingSampler().sample(
        _bqm(Q, off), num_reads=500, seed=11)
    x_found, _ = _decode(ss.first.sample, m)
    assert objective_value(selection_to_y(x_found, m), terms) == pytest.approx(best_obj)
    # the energy of the ground state equals the objective, penalties all zero
    assert ss.first.energy == pytest.approx(
        qubo_energy(complete_baseline_auxiliary(best_x, m), Q, off))


def test_an_insufficient_penalty_actually_breaks_feasibility():
    """The floors are not decorative: below them the ground state is infeasible.

    Without this, a passing feasibility test proves nothing -- it could pass
    for any weight at all.
    """
    rng = np.random.default_rng(3)
    m, n = 8, 4
    terms = bonferroni_coverage_terms(_toy_cell_sets(m, rng), m)
    Q, off, meta = build_baseline_qubo(terms, n_select=n)
    weak = build_baseline_qubo(terms, n_select=n,
                               lambda_select=meta["lambda_select"],
                               lambda_rosenberg=1e-3)[0]
    ss = samplers.SimulatedAnnealingSampler().sample(
        _bqm(weak, 0.0), num_reads=100, seed=5)
    i_idx, j_idx = pair_indices(m)
    x, y = _decode(ss.first.sample, m)
    assert not np.array_equal(y, x[i_idx] * x[j_idx])


def test_variable_count_is_quadratic_in_pads():
    """The scaling cost of the baseline-variable trick, recorded as a fact.

    This is what decides whether a given instance fits on hardware, so it is
    asserted rather than left implicit.
    """
    rng = np.random.default_rng(4)
    for m in (8, 12, 20):
        terms = bonferroni_coverage_terms(_toy_cell_sets(m, rng), m)
        _, _, meta = build_baseline_qubo(terms, n_select=4)
        assert meta["n_variables"] == m + m * (m - 1) // 2
