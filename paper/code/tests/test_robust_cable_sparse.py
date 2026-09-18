"""Tests for scenario-mean robustness, cable cost, and sparse QUBO output."""

import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.arrays.random_array import random_layout
from thz_opt.constraints.cable import cable_summary, mst_length, star_cable_cost
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid
from thz_opt.qubo.baseline_qubo import (
    BaselineTerms,
    objective_value,
    sidelobe_energy_terms,
    track_cell_multiplicities,
)
from thz_opt.qubo.multiepoch import (
    MultiEpochSpec,
    build_multiepoch_qubo,
    complete_multiepoch_auxiliary,
)
from thz_opt.qubo.objective import qubo_energy
from thz_opt.qubo.robust import scenario_mean, scenario_spread

LAM = 1e-3
OBS = [ObservationConfig(frequency_hz=3e11, declination_deg=d,
                         latitude_deg=-23.0, n_times=5)
       for d in (-50.0, -30.0, -10.0)]


# ------------------------------------------------------------------ robust

def _terms(n_pairs, lin, quad, n_pads=6):
    return BaselineTerms(n_pads=n_pads, linear=np.asarray(lin, dtype=float),
                         quadratic=dict(quad))


def test_scenario_mean_averages_coefficients():
    a = _terms(3, [1.0, 2.0, 3.0], {(0, 1): 4.0})
    b = _terms(3, [3.0, 2.0, 1.0], {(0, 1): 0.0, (1, 2): 8.0})
    m = scenario_mean([a, b])
    assert np.allclose(m.linear, [2.0, 2.0, 2.0])
    assert m.quadratic[(0, 1)] == pytest.approx(2.0)
    assert m.quadratic[(1, 2)] == pytest.approx(4.0)


def test_scenario_mean_respects_weights():
    a = _terms(2, [0.0, 0.0], {})
    b = _terms(2, [10.0, 10.0], {})
    m = scenario_mean([a, b], weights=[3.0, 1.0])
    assert np.allclose(m.linear, [2.5, 2.5])


def test_scenario_mean_equals_the_mean_objective():
    """The whole reason the mean is usable: it is itself a quadratic form."""
    pads = golden_spiral_layout(7, r_min=20.0, r_max=500.0)
    grid = UVGrid.from_uv(layout_to_uv_tracks(pads, OBS[1]), 16)
    terms = [sidelobe_energy_terms(track_cell_multiplicities(pads, LAM, grid, o), 7)
             for o in OBS]
    mean_terms = scenario_mean(terms)

    rng = np.random.default_rng(0)
    for _ in range(10):
        y = rng.integers(0, 2, terms[0].n_pairs)
        direct = float(np.mean([objective_value(y, t) for t in terms]))
        assert objective_value(y, mean_terms) == pytest.approx(direct)


def test_scenario_spread_reports_the_tail_the_mean_hides():
    a = _terms(2, [0.0, 0.0], {})
    b = _terms(2, [10.0, 0.0], {})
    y = np.array([1, 0])
    s = scenario_spread(y, [a, b])
    assert s["mean"] == pytest.approx(5.0)
    assert s["min"] == pytest.approx(0.0) and s["max"] == pytest.approx(10.0)
    assert s["worst_case"] == pytest.approx(10.0)
    assert s["std"] == pytest.approx(5.0)


def test_scenario_validation():
    a = _terms(2, [0.0, 0.0], {})
    with pytest.raises(ValueError):
        scenario_mean([])
    with pytest.raises(ValueError):
        scenario_mean([a, _terms(3, [0.0, 0.0, 0.0], {})])
    with pytest.raises(ValueError):
        scenario_mean([a, a], weights=[0.0, 0.0])


# ------------------------------------------------------------------- cable

def test_star_cost_is_distance_to_the_hub():
    pads = np.array([[0.0, 0.0], [3.0, 4.0], [-6.0, 8.0]])
    c = star_cable_cost(pads, hub=np.array([0.0, 0.0]))
    assert np.allclose(c, [0.0, 5.0, 10.0])


def test_star_hub_defaults_to_the_centroid():
    pads = random_layout(12, 30.0, 1000.0, seed=5)
    c = star_cable_cost(pads)
    assert np.allclose(c, np.linalg.norm(pads - pads.mean(axis=0), axis=1))


def test_star_never_undercharges_relative_to_the_mst():
    """A star is one spanning tree, so it is an upper bound on the minimum."""
    pads = random_layout(20, 30.0, 1000.0, seed=6)
    for sel in (np.arange(20), np.arange(0, 20, 2), np.array([1, 4, 9])):
        s = cable_summary(pads, sel)
        assert s["star_length_m"] >= s["mst_length_m"] - 1e-9
        assert s["star_over_mst"] >= 1.0 - 1e-9


def test_mst_of_two_collinear_pads_is_their_span():
    pads = np.array([[0.0, 0.0], [10.0, 0.0]])
    # hub at the centroid (5, 0): the tree is hub-left plus hub-right
    assert mst_length(pads) == pytest.approx(10.0)


# ------------------------------------------------------------------ sparse

def test_sparse_and_dense_multiepoch_agree_exactly():
    M, T, N = 9, 2, 4
    pads = random_layout(M, 30.0, 800.0, seed=7)
    grid = UVGrid.from_uv(layout_to_uv_tracks(pads, OBS[1]), 16)
    mult = [track_cell_multiplicities(pads, LAM, grid, o) for o in OBS[:T]]
    spec = MultiEpochSpec(M, T, N, 1e5, 1e4, lambda_move=50.0)

    Qd, offd, md = build_multiepoch_qubo(mult, spec)
    Qs, offs, ms = build_multiepoch_qubo(mult, spec, sparse=True)
    assert offd == pytest.approx(offs)
    assert ms["sparse"] is True and md["sparse"] is False

    rng = np.random.default_rng(1)
    for _ in range(10):
        xm = np.zeros((T, M), dtype=int)
        for t in range(T):
            xm[t, rng.choice(M, N, replace=False)] = 1
        v = complete_multiepoch_auxiliary(xm, spec)
        e_sparse = sum(c * v[a] * v[b] for (a, b), c in Qs.items()) + offs
        assert qubo_energy(v, Qd, offd) == pytest.approx(e_sparse)


def test_sparse_bypasses_the_dense_variable_ceiling():
    """The dense guard exists to stop an accidental multi-GB allocation; the
    sparse path is how a realistic size is actually reached."""
    M, T = 20, 2
    pads = random_layout(M, 30.0, 800.0, seed=8)
    grid = UVGrid.from_uv(layout_to_uv_tracks(pads, OBS[1]), 12)
    mult = [track_cell_multiplicities(pads, LAM, grid, o) for o in OBS[:T]]
    spec = MultiEpochSpec(M, T, 6, 1e5, 1e4)

    with pytest.raises(ValueError, match="sparse=True"):
        build_multiepoch_qubo(mult, spec, max_variables=10)

    Q, _, meta = build_multiepoch_qubo(mult, spec, max_variables=10, sparse=True)
    assert isinstance(Q, dict)
    assert meta["n_variables"] == T * (M + M * (M - 1) // 2)
    assert len(Q) < meta["n_variables"] ** 2  # genuinely sparse
