import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.interferometry.baselines import compute_baselines
from thz_opt.interferometry.earth_rotation import ObservationConfig
from thz_opt.interferometry.uv import UVGrid, layout_to_uv
from thz_opt.qubo.baseline_qubo import (
    baseline_cell_overlaps,
    bonferroni_bound,
    bonferroni_coverage_terms,
    build_baseline_qubo,
    complete_baseline_auxiliary,
    cornwell_terms,
    density_match_terms,
    objective_value,
    rosenberg_penalty_floor,
    selection_to_y,
)
from thz_opt.qubo.coefficients import baseline_cell_sets, exact_cells_covered
from thz_opt.qubo.exhaustive import feasible_configurations
from thz_opt.qubo.objective import qubo_energy

LAM = 1e-3
M, N = 8, 4


@pytest.fixture(scope="module")
def toy():
    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    uv = layout_to_uv(pads, LAM)
    grid = UVGrid.from_uv(uv, 16)
    cell_sets = baseline_cell_sets(pads, LAM, grid)
    b, _, _ = compute_baselines(pads)
    tracks = [np.vstack((b[k : k + 1] / LAM, -b[k : k + 1] / LAM)) for k in range(b.shape[0])]
    return pads, grid, cell_sets, tracks


def test_selection_to_y_matches_products(toy):
    x = np.zeros(M, dtype=int)
    x[[0, 2, 3, 6]] = 1
    y = selection_to_y(x, M)
    i, j = np.triu_indices(M, 1)
    assert np.array_equal(y, x[i] * x[j])
    assert y.sum() == N * (N - 1) // 2


def test_overlaps_match_brute_force_intersections(toy):
    _, _, cs, _ = toy
    ov = baseline_cell_overlaps(cs)
    for (k, l), v in list(ov.items())[:40]:
        assert v == len(set(cs[k].tolist()) & set(cs[l].tolist()))
    # pairs absent from the dict must genuinely share nothing
    assert len(set(cs[0].tolist()) & set(cs[1].tolist())) == ov.get((0, 1), 0)


def test_bonferroni_is_a_lower_bound_on_exact_coverage(toy):
    _, _, cs, _ = toy
    ov = baseline_cell_overlaps(cs)
    for x in feasible_configurations(M, N):
        y = selection_to_y(x, M)
        assert bonferroni_bound(y, cs, ov) <= exact_cells_covered(x, cs, M) + 1e-9


def test_bonferroni_is_exact_without_triple_coverage(toy):
    """With snapshot sampling on this grid no cell is covered three times, so
    the second-order bound is not a bound at all -- it is the answer."""
    _, _, cs, _ = toy
    ov = baseline_cell_overlaps(cs)
    for x in feasible_configurations(M, N):
        y = selection_to_y(x, M)
        assert bonferroni_bound(y, cs, ov) == pytest.approx(exact_cells_covered(x, cs, M))


def test_qubo_energy_equals_minus_the_bound(toy):
    _, _, cs, _ = toy
    terms = bonferroni_coverage_terms(cs, M)
    Q, off, meta = build_baseline_qubo(terms, N)
    assert meta["n_variables"] == M + M * (M - 1) // 2
    for x in feasible_configurations(M, N):
        full = complete_baseline_auxiliary(x, M)
        assert qubo_energy(full, Q, off) == pytest.approx(-bonferroni_bound(selection_to_y(x, M), cs))


def test_rosenberg_penalty_makes_cheating_unprofitable(toy):
    _, _, cs, _ = toy
    terms = bonferroni_coverage_terms(cs, M)
    Q, off, meta = build_baseline_qubo(terms, N)
    rng = np.random.default_rng(0)
    for x in feasible_configurations(M, N)[:15]:
        full = complete_baseline_auxiliary(x, M)
        e = qubo_energy(full, Q, off)
        for _ in range(100):
            other = full.copy()
            other[M:] = rng.integers(0, 2, meta["n_y"])
            assert qubo_energy(other, Q, off) >= e - 1e-9


def test_infeasible_selection_cannot_win(toy):
    _, _, cs, _ = toy
    terms = bonferroni_coverage_terms(cs, M)
    Q, off, meta = build_baseline_qubo(terms, N)
    best_feasible = min(
        qubo_energy(complete_baseline_auxiliary(x, M), Q, off)
        for x in feasible_configurations(M, N)
    )
    for k in (N - 1, N + 1, M):
        for x in feasible_configurations(M, k)[:20]:
            assert qubo_energy(complete_baseline_auxiliary(x, M), Q, off) > best_feasible


def test_cornwell_terms_reproduce_the_direct_energy(toy):
    _, _, _, tracks = toy
    terms = cornwell_terms(tracks, M)
    allv = np.vstack(tracks)
    R = float(np.hypot(allv[:, 0], allv[:, 1]).max())
    floor = 1e-3 * R
    for x in feasible_configurations(M, N)[:10]:
        y = selection_to_y(x, M)
        active = np.vstack([tracks[k] for k in range(len(tracks)) if y[k]])
        d = np.hypot(active[:, None, 0] - active[None, :, 0],
                     active[:, None, 1] - active[None, :, 1])
        iu = np.triu_indices(active.shape[0], 1)
        direct = R ** 2 * float(np.sum(1.0 / np.maximum(d[iu], floor) ** 2))
        assert objective_value(y, terms) == pytest.approx(direct, rel=1e-9)


def test_cornwell_cutoff_sparsifies(toy):
    _, _, _, tracks = toy
    dense = cornwell_terms(tracks, M)
    sparse = cornwell_terms(tracks, M, coupling_cutoff=1e5)
    assert sparse.sparsity()["n_couplings"] < dense.sparsity()["n_couplings"]


def test_density_terms_reproduce_the_direct_chi2(toy):
    _, _, _, tracks = toy
    allv = np.vstack(tracks)
    edges = np.linspace(0.0, float(np.hypot(allv[:, 0], allv[:, 1]).max()), 12)
    target = np.full(11, 4.0)
    terms = density_match_terms(tracks, M, target, edges)
    for x in feasible_configurations(M, N)[:10]:
        y = selection_to_y(x, M)
        active = np.vstack([tracks[k] for k in range(len(tracks)) if y[k]])
        h, _ = np.histogram(np.hypot(active[:, 0], active[:, 1]), bins=edges)
        assert objective_value(y, terms) == pytest.approx(float(np.sum((h - target) ** 2)))


def test_penalty_floor_scales_with_the_objective(toy):
    _, _, cs, _ = toy
    terms = bonferroni_coverage_terms(cs, M)
    lam = rosenberg_penalty_floor(terms, safety=2.0)
    assert lam >= 2.0 * max(len(c) for c in cs)


def test_bonferroni_beats_the_pairwise_surrogate_under_earth_rotation():
    """The measurement that justifies the whole baseline-variable formulation."""
    from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
    from thz_opt.qubo.coefficients import pairwise_uv_weights, surrogate_value

    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    obs = ObservationConfig(frequency_hz=3e11, declination_deg=-30.0,
                            latitude_deg=-23.0, n_times=11)
    grid = UVGrid.from_uv(layout_to_uv_tracks(pads, obs), 16)
    cs = baseline_cell_sets(pads, LAM, grid, observation=obs)
    ov = baseline_cell_overlaps(cs)
    w = pairwise_uv_weights(cs, grid, "shared")

    configs = feasible_configurations(M, N)
    exact = np.array([exact_cells_covered(x, cs, M) for x in configs], dtype=float)
    bonf = np.array([bonferroni_bound(selection_to_y(x, M), cs, ov) for x in configs])
    surr = np.array([surrogate_value(x, w, M) for x in configs])

    def rho(a, b):
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        return float(np.corrcoef(ra, rb)[0, 1])

    assert np.all(bonf <= exact + 1e-9)
    assert rho(bonf, exact) >= rho(surr, exact)
    assert exact[int(np.argmax(bonf))] >= exact.max() - 1e-9
