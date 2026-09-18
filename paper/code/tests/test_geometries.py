import numpy as np
import pytest

from thz_opt.arrays.geometries import (
    gaussian_layout,
    hierarchical_layout,
    multi_arm_spiral,
    phyllotaxis_layout,
    reuleaux_base,
    reuleaux_layout,
)
from thz_opt.arrays.golden_spiral import PHI
from thz_opt.interferometry.baselines import baseline_lengths, baseline_matrix
from thz_opt.interferometry.uv import layout_to_uv
from thz_opt.metrics.density_matching import (
    cornwell_energy,
    density_match_chi2,
    fit_gaussian_sigma,
    radial_density_profile,
)

EXTENT = dict(r_min=20.0, r_max=1000.0)


def test_multi_arm_spiral_shape_and_extent():
    xy = multi_arm_spiral(30, n_arms=3, **EXTENT)
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert xy.shape == (30, 2)
    assert r.min() == pytest.approx(EXTENT["r_min"])
    assert r.max() == pytest.approx(EXTENT["r_max"])


def test_multi_arm_spiral_arms_are_rotated_copies():
    a = multi_arm_spiral(30, n_arms=3, **EXTENT)
    r = np.hypot(a[:, 0], a[:, 1])
    # arms are dealt in equal blocks, so the radial sequences must coincide
    assert np.allclose(r[:10], r[10:20])
    assert np.allclose(r[:10], r[20:])


def test_multi_arm_spiral_growth_is_configurable():
    """growth is the radius multiplier per turn of the unnormalised spiral."""
    xy = multi_arm_spiral(5, n_arms=1, growth=3.0, n_turns=2.0, normalize=False)
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert r[-1] / r[0] == pytest.approx(9.0)  # two turns at x3 per turn
    with pytest.raises(ValueError):
        multi_arm_spiral(5, growth=0.9)


def test_reuleaux_has_constant_width():
    """Defining property: every pair of boundary points is at most `width`
    apart, and the diameter is attained in every direction."""
    p = reuleaux_base(width=100.0, n_points=360)
    d = baseline_matrix(p)
    assert d.max() == pytest.approx(100.0, rel=1e-3)
    # the support function is constant: width measured along many directions
    for ang in np.linspace(0, np.pi, 17):
        u = np.array([np.cos(ang), np.sin(ang)])
        proj = p @ u
        # 360 sampled boundary points, so the extremum in a given direction can
        # fall between samples: a few parts in 1e3 is discretisation, not shape
        assert proj.max() - proj.min() == pytest.approx(100.0, rel=5e-3)


def test_reuleaux_layout_extent_and_rings():
    xy = reuleaux_layout(21, n_rings=3, **EXTENT)
    assert xy.shape == (21, 2)
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert r.max() <= EXTENT["r_max"] * 1.001


def test_hierarchical_layout_size_and_selfsimilarity():
    base = reuleaux_base(1.0, 3)
    xy = hierarchical_layout(base, levels=3, scale=0.35, r_max=1000.0)
    assert xy.shape == (27, 2)
    assert np.hypot(xy[:, 0], xy[:, 1]).max() == pytest.approx(1000.0)
    # a hierarchy of repeated patterns is deliberately redundant: the baseline
    # set contains repeated lengths at every scale
    d = np.round(baseline_lengths(xy), 6)
    assert len(np.unique(d)) < d.size


def test_gaussian_layout_gives_gaussian_uv():
    """Difference of two Gaussians is Gaussian with sqrt(2) times the width."""
    xy = gaussian_layout(600, r_min=0.0, r_max=1e9, seed=4, sigma_fraction=1e-7)
    sigma_a = fit_gaussian_sigma(xy)
    uv = layout_to_uv(xy, 1.0)
    sigma_uv = fit_gaussian_sigma(uv)
    assert sigma_uv / sigma_a == pytest.approx(np.sqrt(2.0), rel=0.05)


def test_gaussian_layout_respects_extent_when_truncated():
    xy = gaussian_layout(200, seed=2, **EXTENT)
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert r.min() >= EXTENT["r_min"] and r.max() <= EXTENT["r_max"]


def test_phyllotaxis_reduces_to_vogel_by_default():
    from thz_opt.arrays.fibonacci import fibonacci_layout

    a = phyllotaxis_layout(30, **EXTENT)
    b = fibonacci_layout(30, variant="golden_angle", **EXTENT)
    # same angular rule (1/PHI^2 turns == the golden angle) and radial exponent
    assert np.allclose(np.hypot(*a.T), np.hypot(*b.T))
    assert 1.0 / PHI ** 2 == pytest.approx(1.0 - 1.0 / PHI)


def test_phyllotaxis_alpha_changes_the_layout():
    a = phyllotaxis_layout(30, alpha=1.0 / PHI ** 2, **EXTENT)
    b = phyllotaxis_layout(30, alpha=np.sqrt(2) - 1.0, **EXTENT)
    assert not np.allclose(a, b)


def test_rational_alpha_collapses_position_angles():
    """A rational alpha = p/q repeats angles every q antennas; an irrational
    one never does.  This is the whole Fibonacci-vs-golden distinction."""
    def n_directions(xy):
        r = np.hypot(xy[:, 0], xy[:, 1])
        return len(np.unique(np.round(xy / r[:, None], 6), axis=0))

    assert n_directions(phyllotaxis_layout(24, alpha=1.0 / 3.0, **EXTENT)) == 3
    assert n_directions(phyllotaxis_layout(24, alpha=5.0 / 8.0, **EXTENT)) == 8
    assert n_directions(phyllotaxis_layout(24, alpha=1.0 / PHI ** 2, **EXTENT)) == 24


def test_radial_density_profile_normalisation():
    uv = layout_to_uv(gaussian_layout(100, seed=5, **EXTENT), 1e-3)
    centres, density, edges = radial_density_profile(uv, 25)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    assert float((density * area).sum()) == pytest.approx(1.0)
    assert centres.shape == density.shape


def test_density_match_chi2_is_small_for_a_gaussian_array():
    """A Gaussian antenna distribution must match a Gaussian UV target closely,
    and a ring array must not."""
    g = density_match_chi2(layout_to_uv(gaussian_layout(200, seed=6, **EXTENT), 1e-3))
    from thz_opt.arrays.random_array import uniform_ring_layout

    ring = density_match_chi2(layout_to_uv(uniform_ring_layout(200, n_rings=1, **EXTENT), 1e-3))
    assert g["normalised_chi2"] < ring["normalised_chi2"]


def test_cornwell_energy_prefers_the_spread_configuration():
    """Clustered UV points must score worse (higher energy) than spread ones."""
    spread = np.column_stack((np.linspace(-1, 1, 50), np.zeros(50)))
    clustered = np.column_stack((np.concatenate([np.linspace(-1, -0.9, 25),
                                                 np.linspace(0.9, 1, 25)]), np.zeros(50)))
    assert cornwell_energy(clustered)["cornwell_F"] > cornwell_energy(spread)["cornwell_F"]


def test_cornwell_energy_floor_reported():
    dup = np.zeros((10, 2))
    dup[:, 0] = np.repeat([1.0, 2.0], 5)
    out = cornwell_energy(dup)
    assert out["pairs_at_floor"] > 0
    assert np.isfinite(out["cornwell_F"])
