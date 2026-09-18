import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.constraints.phase import PhaseConfig, pairwise_phase_penalty, unsupported_stations
from thz_opt.constraints.pwv import altitude_pwv_cost, constant_pwv_cost, normalise_cost
from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation, shadow_pairs
from thz_opt.interferometry.baselines import baseline_lengths
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, layout_to_uv, uv_occupancy
from thz_opt.metrics.baseline_distribution import (
    baseline_summary,
    empirical_cdf,
    length_histogram,
    log_spacing_uniformity,
    powerlaw_fit,
)
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.metrics.redundancy import multiplicity_histogram, redundancy_metrics
from thz_opt.metrics.uv_coverage import (
    angular_uniformity,
    annulus_mask,
    coverage_summary,
    density_uniformity,
    occupied_fraction,
    radial_profile,
    unique_cells,
)

LAM = 1e-3


@pytest.fixture(scope="module")
def sampled():
    xy = golden_spiral_layout(20, r_min=20.0, r_max=1000.0)
    uv = layout_to_uv(xy, LAM)
    grid = UVGrid.from_uv(uv, 64)
    occ, _ = uv_occupancy(uv, grid)
    return xy, uv, grid, occ


def test_occupied_fraction_bounds(sampled):
    _, _, _, occ = sampled
    f = occupied_fraction(occ)
    assert 0.0 <= f <= 1.0
    assert unique_cells(occ) == int((occ > 0).sum())


def test_occupied_fraction_with_mask(sampled):
    _, _, grid, occ = sampled
    mask = annulus_mask(grid, 0.0, grid.uv_max)
    f_full, f_mask = occupied_fraction(occ), occupied_fraction(occ, mask)
    assert f_mask >= f_full  # the mask removes unreachable corner cells


def test_redundancy_accounting(sampled):
    _, uv, _, occ = sampled
    r = redundancy_metrics(occ)
    assert r["n_samples"] == occ.sum()
    assert r["unique_cells"] + r["redundant_samples"] == r["n_samples"]
    m, counts = multiplicity_histogram(occ)
    assert (counts * m).sum() == r["n_samples"]


def test_fully_redundant_grid():
    occ = np.zeros((4, 4), dtype=int)
    occ[1, 1] = 10
    r = redundancy_metrics(occ)
    assert r["unique_cells"] == 1
    assert r["redundancy_fraction"] == pytest.approx(0.9)
    assert density_uniformity(occ) == 1.0  # single occupied cell, by convention


def test_density_uniformity_is_one_for_flat_occupancy():
    assert density_uniformity(np.full((3, 3), 5)) == pytest.approx(1.0)


def test_angular_uniformity_extremes():
    ring = np.stack(
        [np.cos(np.linspace(0, np.pi, 36, endpoint=False)),
         np.sin(np.linspace(0, np.pi, 36, endpoint=False))], axis=1
    )
    assert angular_uniformity(ring, n_bins=36) == pytest.approx(1.0)
    spoke = np.repeat(np.array([[1.0, 0.0]]), 36, axis=0)
    assert angular_uniformity(spoke, n_bins=36) == pytest.approx(0.0)


def test_radial_profile_conserves_samples(sampled):
    _, uv, _, _ = sampled
    _, counts, density = radial_profile(uv, n_bins=20)
    assert counts.sum() == uv.shape[0]
    assert np.all(np.isfinite(density))


def test_coverage_summary_has_no_nans(sampled):
    _, uv, grid, occ = sampled
    s = coverage_summary(uv, occ, grid)
    assert all(np.isfinite(v) for v in s.values() if isinstance(v, float))
    assert s["n_samples"] == uv.shape[0]
    assert 0.0 <= s["occupied_fraction_full"] <= 1.0


def test_baseline_histogram_and_cdf(sampled):
    xy, _, _, _ = sampled
    d = baseline_lengths(xy)
    _, counts, _ = length_histogram(d, n_bins=15)
    assert counts.sum() == d.size
    x, cdf = empirical_cdf(d)
    assert np.all(np.diff(x) >= 0) and cdf[-1] == pytest.approx(1.0)


def test_powerlaw_fit_recovers_a_known_exponent():
    """Counts per logarithmic bin must scale as d**alpha.

    A log bin at d has width proportional to d, so a sample uniform in d has
    N(d) proportional to d, i.e. alpha = 1.  A sample uniform in log10(d) is
    flat per log bin, i.e. alpha = 0.
    """
    rng = np.random.default_rng(0)
    uniform_in_d = rng.uniform(1.0, 1000.0, 200000)
    fit = powerlaw_fit(uniform_in_d, n_bins=20)
    assert fit["alpha"] == pytest.approx(1.0, abs=0.2)
    assert fit["r_squared"] > 0.9

    uniform_in_log = 10.0 ** rng.uniform(0.0, 3.0, 200000)
    flat = powerlaw_fit(uniform_in_log, n_bins=20)
    assert flat["alpha"] == pytest.approx(0.0, abs=0.2)
    assert log_spacing_uniformity(uniform_in_log, 20) > 0.99


def test_baseline_summary_keys(sampled):
    xy, _, _, _ = sampled
    s = baseline_summary(baseline_lengths(xy))
    assert s["n_baselines"] == 190
    assert s["d_min_m"] <= s["d_median_m"] <= s["d_max_m"]
    assert s["dynamic_range"] == pytest.approx(s["d_max_m"] / s["d_min_m"])


def test_psf_metrics_on_a_known_beam(sampled):
    _, uv, grid, _ = sampled
    psf, info = dirty_beam(uv, grid)
    m = psf_metrics(psf, info["pixel_scale_arcsec"])
    assert m["peak"] == pytest.approx(1.0)
    assert 0.0 <= m["peak_sidelobe_level"] <= 1.0
    assert m["fwhm_px"] > 0
    assert 0.0 <= m["ellipticity"] < 1.0


def test_separation_check_flags_violations():
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [100.0, 0.0]])
    cfg = SeparationConfig(dish_diameter_m=6.0, minimum_separation_factor=1.5)
    r = check_minimum_separation(xy, cfg)
    assert r["d_min_m"] == pytest.approx(9.0)
    assert r["n_violations"] == 1 and r["violations"][0][:2] == (0, 1)
    assert not r["valid"]
    m = shadow_pairs(xy, cfg)
    assert m[0, 1] and not m[0, 2] and not m.diagonal().any()


def test_separation_passes_for_a_spread_layout():
    xy = golden_spiral_layout(15, r_min=50.0, r_max=1000.0)
    assert check_minimum_separation(xy)["valid"]


def test_pwv_costs_are_synthetic_but_well_formed():
    assert np.allclose(constant_pwv_cost(5, 2.0), 2.0)
    c = altitude_pwv_cost(np.array([0.0, 2000.0]))
    assert c[1] < c[0]
    n = normalise_cost(c)
    assert n.min() == 0.0 and n.max() == 1.0
    assert np.allclose(normalise_cost(np.ones(4)), 0.0)


def test_phase_connectivity_surrogate():
    """A chain of pads is fully supported; an isolated distant pad is not."""
    cfg = PhaseConfig(coherence_length_m=100.0, core_radius_m=80.0)
    chain = np.array([[0.0, 0.0], [90.0, 0.0], [180.0, 0.0], [5000.0, 0.0]])
    r = unsupported_stations(chain, config=cfg, centre=np.array([0.0, 0.0]))
    assert r["n_unsupported"] == 1 and r["unsupported_indices"] == [3]
    p = pairwise_phase_penalty(chain, cfg)
    assert p[0, 3] == 1.0 and p[0, 1] == 0.0
