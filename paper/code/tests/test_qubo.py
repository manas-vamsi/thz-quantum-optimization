import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.interferometry.uv import UVGrid, layout_to_uv
from thz_opt.qubo.coefficients import (
    baseline_cell_sets,
    cell_claim_counts,
    coverage_is_exactly_quadratic,
    exact_cells_covered,
    pairwise_uv_weights,
    suggest_penalty_weights,
    surrogate_value,
)
from thz_opt.qubo.exhaustive import (
    all_configurations,
    best_configuration,
    enumerate_energies,
    feasible_configurations,
)
from thz_opt.qubo.objective import (
    QUBOTerms,
    build_qubo,
    build_qubo_with_cell_variables,
    direct_terms,
    qubo_energy,
    symmetric_energy,
    to_qubo_dict,
    to_symmetric,
    complete_auxiliary,
    to_upper_triangular,
)
from thz_opt.qubo.validation import (
    check_convention_roundtrip,
    compare_surrogate_to_exact,
    validate_qubo,
)

LAM = 1e-3
M, N_SEL = 8, 4


@pytest.fixture(scope="module")
def toy():
    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    uv = layout_to_uv(pads, LAM)
    grid = UVGrid.from_uv(uv, 16)  # coarse on purpose: cells are contested
    cell_sets = baseline_cell_sets(pads, LAM, grid)
    w = pairwise_uv_weights(cell_sets, grid, "shared")
    return pads, grid, cell_sets, w


def _terms(pads, w, lam=10.0):
    return QUBOTerms(
        n_pads=M,
        n_select=N_SEL,
        w_uv=w,
        shadow=shadow_pairs(pads, SeparationConfig(6.0, 1.5)).astype(float),
        lambda_select=lam,
        lambda_shadow=lam,
    )


def test_select_penalty_is_zero_exactly_on_the_constraint():
    t = QUBOTerms(n_pads=M, n_select=N_SEL, lambda_select=3.0)
    for x in feasible_configurations(M, N_SEL)[:5]:
        assert direct_terms(x, t)["H_select"] == 0.0
    off = np.zeros(M, dtype=int)
    off[:5] = 1
    assert direct_terms(off, t)["H_select"] == pytest.approx(3.0)


def test_qubo_matches_direct_objective_everywhere(toy):
    pads, _, _, w = toy
    r = validate_qubo(_terms(pads, w))
    assert r["n_configurations"] == 2 ** M
    assert r["within_tolerance"], r


def test_qubo_matches_direct_with_all_terms_active(toy):
    pads, _, _, w = toy
    t = _terms(pads, w)
    t.phase = np.triu(np.ones((M, M)), 1)
    t.lambda_phase = 0.7
    t.pwv_cost = np.linspace(0.0, 1.0, M)
    t.lambda_pwv = 1.3
    assert validate_qubo(t)["within_tolerance"]


def test_conventions_agree(toy):
    pads, _, _, w = toy
    Q, off = build_qubo(_terms(pads, w))
    assert check_convention_roundtrip(Q, off)["within_tolerance"]
    A = to_symmetric(Q)
    assert np.allclose(A, A.T)
    assert np.allclose(to_upper_triangular(A), np.triu(Q))
    x = np.zeros(M, dtype=int)
    x[[0, 2, 5, 7]] = 1
    assert qubo_energy(x, Q, off) == pytest.approx(symmetric_energy(x, A, off))


def test_qubo_dict_roundtrip(toy):
    pads, _, _, w = toy
    Q, off = build_qubo(_terms(pads, w))
    d = to_qubo_dict(Q)
    x = np.zeros(M, dtype=int)
    x[[1, 3, 4, 6]] = 1
    e = sum(v * x[i] * x[j] for (i, j), v in d.items()) + off
    assert e == pytest.approx(qubo_energy(x, Q, off))


def test_optimum_respects_the_selection_constraint(toy):
    pads, _, _, w = toy
    t = _terms(pads, w, lam=suggest_penalty_weights(w, N_SEL)["lambda_select"])
    Q, off = build_qubo(t)
    configs, direct, quad = enumerate_energies(t, Q, off)
    best, _, _ = best_configuration(configs, quad)
    assert best.sum() == N_SEL
    assert np.allclose(direct, quad)


def test_penalty_rule_is_sufficient(toy):
    """With the suggested weights no infeasible configuration can win."""
    pads, _, _, w = toy
    p = suggest_penalty_weights(w, N_SEL)
    t = _terms(pads, w, lam=p["lambda_select"])
    Q, off = build_qubo(t)
    configs, _, quad = enumerate_energies(t, Q, off)
    feas = configs.sum(axis=1) == N_SEL
    assert quad[feas].min() < quad[~feas].min()


def test_exact_coverage_and_surrogate_relationship(toy):
    _, grid, cell_sets, w = toy
    claims = cell_claim_counts(cell_sets)
    assert max(claims.values()) > 1, "toy grid should have contested cells"
    assert not coverage_is_exactly_quadratic(cell_sets)

    w_count = pairwise_uv_weights(cell_sets, grid, "count")
    for x in feasible_configurations(M, N_SEL):
        exact = exact_cells_covered(x, cell_sets, M)
        assert surrogate_value(x, w_count, M) >= exact - 1e-9  # count is an upper bound
        assert exact > 0


def test_count_weighting_is_exact_when_no_cell_is_contested():
    """Fine grid -> every baseline owns its cells -> the QUBO form is exact."""
    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    uv = layout_to_uv(pads, LAM)
    grid = UVGrid.from_uv(uv, 256)
    cs = baseline_cell_sets(pads, LAM, grid)
    assert coverage_is_exactly_quadratic(cs)
    w = pairwise_uv_weights(cs, grid, "count")
    for x in feasible_configurations(M, N_SEL)[:20]:
        assert surrogate_value(x, w, M) == pytest.approx(exact_cells_covered(x, cs, M))


def test_surrogate_comparison_report(toy):
    _, _, cell_sets, w = toy
    r = compare_surrogate_to_exact(cell_sets, w, M, N_SEL)
    assert r["n_configurations"] == 70
    assert -1.0 <= r["spearman_rho"] <= 1.0
    assert r["coverage_shortfall_cells"] >= 0.0


def test_radial_weighting_favours_long_baselines(toy):
    _, grid, cell_sets, _ = toy
    w_shared = pairwise_uv_weights(cell_sets, grid, "shared")
    w_radial = pairwise_uv_weights(cell_sets, grid, "radial")
    assert np.all(w_radial >= w_shared - 1e-12)


def test_exact_auxiliary_qubo_reproduces_coverage():
    """The auxiliary-variable construction is exact where the pairwise one is not."""
    pads = golden_spiral_layout(5, r_min=20.0, r_max=300.0)
    uv = layout_to_uv(pads, LAM)
    grid = UVGrid.from_uv(uv, 8)
    cs = baseline_cell_sets(pads, LAM, grid)
    Q, off, meta = build_qubo_with_cell_variables(cs, 5, 3, lambda_select=50.0, lambda_aux=50.0)
    assert meta["n_variables"] > 5 and meta["n_z"] == len(meta["cells"])

    rng = np.random.default_rng(0)
    for x in feasible_configurations(5, 3):
        full = complete_auxiliary(x, meta)
        e = qubo_energy(full, Q, off)
        # the intended completion reproduces the exact coverage reward
        assert e == pytest.approx(-exact_cells_covered(x, cs, 5))
        # and no other auxiliary assignment can do better
        for _ in range(50):
            other = full.copy()
            other[5:] = rng.integers(0, 2, meta["n_variables"] - 5)
            assert qubo_energy(other, Q, off) >= e - 1e-9


def test_enumeration_helpers():
    assert all_configurations(4).shape == (16, 4)
    assert feasible_configurations(6, 2).shape == (15, 6)
    assert np.all(feasible_configurations(6, 2).sum(axis=1) == 2)
    with pytest.raises(ValueError):
        all_configurations(30)


def test_qubo_energy_rejects_non_binary():
    with pytest.raises(ValueError):
        qubo_energy(np.array([0.5, 1.0]), np.zeros((2, 2)))
