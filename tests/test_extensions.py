"""Tests for multi-frequency, gridding, science cases, Pareto and multi-epoch."""

import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.gridding import (
    gaussian_gridded_occupancy,
    grid_sensitivity,
    shifted_grid_average,
)
from thz_opt.interferometry.multifrequency import (
    FrequencyBand,
    channel_coherence_weights,
    multi_frequency_uv,
)
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.optimize.pareto import dominates, pareto_indices, sweep_weights
from thz_opt.optimize.science_cases import (
    SCIENCE_CASES,
    cell_radii,
    compact_source_weights,
    extended_emission_weights,
    radial_weights,
    science_case_objective,
)
from thz_opt.optimize.state import UVState, build_pair_tables
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities
from thz_opt.qubo.multiepoch import (
    MultiEpochSpec,
    build_multiepoch_qubo,
    complete_multiepoch_auxiliary,
    movement_count,
    multiepoch_size_report,
)
from thz_opt.qubo.objective import qubo_energy

LAM = 1e-3
OBS = ObservationConfig(frequency_hz=3e11, declination_deg=-30.0,
                        latitude_deg=-23.0, n_times=9)


# ---------------------------------------------------------------- frequency

def test_zero_bandwidth_reduces_to_single_frequency():
    xy = golden_spiral_layout(8, r_min=20.0, r_max=500.0)
    band = FrequencyBand(3e11, 0.0, 1)
    assert np.allclose(multi_frequency_uv(xy, OBS, band), layout_to_uv_tracks(xy, OBS))
    assert band.uv_radius_ratio == pytest.approx(1.0)


def test_band_scales_uv_radius_as_frequency():
    band = FrequencyBand(3e11, 0.5, 3)
    f = band.frequencies
    assert f[0] == pytest.approx(2.25e11) and f[-1] == pytest.approx(3.75e11)
    # a fixed baseline sits at radius proportional to frequency
    assert band.uv_radius_ratio == pytest.approx(f[-1] / f[0])


def test_multi_frequency_adds_coverage():
    xy = golden_spiral_layout(10, r_min=20.0, r_max=800.0)
    grid = UVGrid(uv_max=2.5e6, n_cells=96)
    mono = uv_occupancy(multi_frequency_uv(xy, OBS, FrequencyBand(3e11, 0.0, 1)), grid)[0]
    wide = uv_occupancy(multi_frequency_uv(xy, OBS, FrequencyBand(3e11, 0.3, 5)), grid)[0]
    assert int((wide > 0).sum()) > int((mono > 0).sum())


def test_high_frequency_end_decorrelates_first():
    band = FrequencyBand(3e11, 0.4, 5)
    g = channel_coherence_weights(600.0, band)
    assert g.shape == (5,)
    assert np.all(np.diff(g) <= 0)  # coherence falls with increasing frequency


def test_band_validation():
    with pytest.raises(ValueError):
        FrequencyBand(-1.0, 0.1, 3)
    with pytest.raises(ValueError):
        FrequencyBand(3e11, 0.1, 0)


# ----------------------------------------------------------------- gridding

def test_zero_kernel_reproduces_hard_binning():
    xy = golden_spiral_layout(10, r_min=20.0, r_max=800.0)
    uv = layout_to_uv_tracks(xy, OBS)
    grid = UVGrid(uv_max=2.5e6, n_cells=64)
    hard, _ = uv_occupancy(uv, grid)
    assert np.allclose(hard, gaussian_gridded_occupancy(uv, grid, kernel_cells=0.0))


def test_gaussian_kernel_conserves_total_weight():
    xy = golden_spiral_layout(10, r_min=20.0, r_max=800.0)
    uv = layout_to_uv_tracks(xy, OBS)
    grid = UVGrid(uv_max=2.5e6, n_cells=64)
    occ = gaussian_gridded_occupancy(uv, grid, kernel_cells=0.7)
    assert occ.sum() == pytest.approx(uv.shape[0], rel=1e-6)
    assert np.any((occ > 0) & (occ < 1))  # genuinely spread, not just binned


def test_shifted_grid_average_reports_spread():
    xy = golden_spiral_layout(10, r_min=20.0, r_max=800.0)
    uv = layout_to_uv_tracks(xy, OBS)
    grid = UVGrid(uv_max=2.5e6, n_cells=64)
    out = shifted_grid_average(lambda o: float((o > 0).sum()), uv, grid, n_shifts=3)
    assert out["n_grids"] == 9
    assert out["min"] <= out["mean"] <= out["max"]
    assert out["relative_spread"] >= 0.0
    s = grid_sensitivity(uv, grid, n_shifts=2)
    assert "unique_cells" in s and "sidelobe_energy" in s


# ------------------------------------------------------------ science cases

def test_radial_weight_profiles_point_the_right_way():
    grid = UVGrid(uv_max=1.0e6, n_cells=32)
    r = cell_radii(grid)
    outer = compact_source_weights(grid)
    inner = extended_emission_weights(grid)
    far, near = r > np.percentile(r, 80), r < np.percentile(r, 20)
    assert outer[far].mean() > outer[near].mean()
    assert inner[near].mean() > inner[far].mean()
    for profile in ("flat", "outer", "inner", "annulus"):
        assert radial_weights(grid, profile).mean() == pytest.approx(1.0)
    with pytest.raises(ValueError):
        radial_weights(grid, "nonsense")


def test_weighted_state_matches_recomputation():
    xy = golden_spiral_layout(9, r_min=20.0, r_max=600.0)
    grid = UVGrid.from_uv(layout_to_uv_tracks(xy, OBS), 16)
    cells, counts = build_pair_tables(track_cell_multiplicities(xy, LAM, grid, OBS))
    w = compact_source_weights(grid)
    st = UVState(cells, counts, 9, grid.n_cells ** 2, cell_weights=w)

    rng = np.random.default_rng(0)
    for _ in range(20):
        st.set_selection(rng.choice(9, 4, replace=False))
        wc, wss = st.recompute_weighted()
        assert st.weighted_cells == pytest.approx(wc)
        assert st.weighted_sum_sq == pytest.approx(wss)


def test_science_case_objective_is_a_minimisation():
    xy = golden_spiral_layout(9, r_min=20.0, r_max=600.0)
    grid = UVGrid.from_uv(layout_to_uv_tracks(xy, OBS), 16)
    cells, counts = build_pair_tables(track_cell_multiplicities(xy, LAM, grid, OBS))
    st = UVState(cells, counts, 9, grid.n_cells ** 2,
                 cell_weights=extended_emission_weights(grid))
    score = science_case_objective("extended_emission")
    st.set_selection([0, 1, 2])
    few = score(st)
    st.set_selection([0, 1, 2, 3, 4, 5])
    many = score(st)
    assert many < few  # more antennas cover more weighted cells, so lower score
    with pytest.raises(ValueError):
        science_case_objective("not_a_case")
    assert set(SCIENCE_CASES) == {"compact_source", "extended_emission"}


# ------------------------------------------------------------------- pareto

def test_domination_rules():
    assert dominates([1.0, 1.0], [2.0, 2.0])
    assert dominates([1.0, 2.0], [1.0, 3.0])
    assert not dominates([1.0, 3.0], [2.0, 2.0])
    assert not dominates([1.0, 1.0], [1.0, 1.0])  # equal does not dominate


def test_pareto_front_selection():
    pts = np.array([[1.0, 5.0], [2.0, 3.0], [5.0, 1.0], [3.0, 4.0], [6.0, 6.0]])
    front = pareto_indices(pts)
    # [2,3] beats [3,4] on both objectives, and [6,6] is beaten by everything
    assert set(front.tolist()) == {0, 1, 2}


def test_sweep_weights_collects_a_front():
    def run_one(w):
        a, b = float(w[0]), float(w[1])
        return (a, b), {"weights": (a, b)}

    out = sweep_weights(run_one, [(1.0, 5.0), (2.0, 3.0), (6.0, 6.0)], ["A", "B"])
    assert out["n_runs"] == 3
    assert out["n_nondominated"] == 2
    assert "convex" in out["note"]


# -------------------------------------------------------------- multi-epoch

def test_multiepoch_qubo_reproduces_the_direct_objective():
    M, N, T = 6, 3, 2
    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    grid = UVGrid.from_uv(layout_to_uv_tracks(pads, OBS), 16)
    epochs = [ObservationConfig(frequency_hz=3e11, declination_deg=-30.0,
                                latitude_deg=-23.0, hour_angle_start_h=a,
                                hour_angle_end_h=b, n_times=5)
              for a, b in ((-2.0, 0.0), (0.0, 2.0))]
    mult = [track_cell_multiplicities(pads, LAM, grid, o) for o in epochs]
    spec = MultiEpochSpec(M, T, N, lambda_select=1e5, lambda_rosenberg=1e4,
                          lambda_move=50.0)
    Q, off, meta = build_multiepoch_qubo(mult, spec)
    assert meta["n_variables"] == T * (M + M * (M - 1) // 2)

    i_idx, j_idx = np.triu_indices(M, 1)
    rng = np.random.default_rng(1)
    for _ in range(15):
        xm = np.zeros((T, M), dtype=int)
        for t in range(T):
            xm[t, rng.choice(M, N, replace=False)] = 1
        n: dict = {}
        for t in range(T):
            for k, d in enumerate(mult[t]):
                if xm[t, i_idx[k]] and xm[t, j_idx[k]]:
                    for c, v in d.items():
                        n[c] = n.get(c, 0) + v
        direct = sum(v * v for v in n.values()) + spec.lambda_move * movement_count(xm)
        assert qubo_energy(complete_multiepoch_auxiliary(xm, spec), Q, off) == pytest.approx(direct)


def test_movement_identity_and_count():
    """|x - y| = x + y - 2xy for binaries, which is what makes the term quadratic."""
    for x in (0, 1):
        for y in (0, 1):
            assert abs(x - y) == x + y - 2 * x * y
    xm = np.array([[1, 0, 1], [1, 1, 0]])
    assert movement_count(xm) == 2


def test_size_report_matches_the_built_model():
    M, T = 5, 2
    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    grid = UVGrid.from_uv(layout_to_uv_tracks(pads, OBS), 12)
    mult = [track_cell_multiplicities(pads, LAM, grid, OBS) for _ in range(T)]
    rep = multiepoch_size_report(mult, M, T)
    assert rep["n_variables_total"] == T * (M + M * (M - 1) // 2)
    assert 0.0 <= rep["fill_fraction"] <= 1.0
    spec = MultiEpochSpec(M, T, 3, 1e4, 1e3)
    _, _, meta = build_multiepoch_qubo(mult, spec)
    assert meta["n_variables"] == rep["n_variables_total"]


def test_multiepoch_refuses_to_build_something_enormous():
    mult = [[{}] * (50 * 49 // 2)]
    spec = MultiEpochSpec(50, 1, 10, 1.0, 1.0)
    with pytest.raises(ValueError):
        build_multiepoch_qubo(mult, spec, max_variables=100)
