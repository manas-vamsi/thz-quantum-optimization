"""Tests for the factorial layout construction and the variance decomposition.

The ANOVA is checked against cases whose answer is known by construction --
data with only a row effect, only a column effect, or only interaction --
because a decomposition that silently attributes variance to the wrong factor
would produce a confident and completely wrong conclusion about which design
choice matters.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))

from thz_opt.arrays.factorial import (  # noqa: E402
    ANGULAR_LAWS,
    GOLDEN_ANGLE,
    RADIAL_LAWS,
    angular_law,
    layout_from_laws,
    radial_law,
)


# --------------------------------------------------------------------------
# radial laws
# --------------------------------------------------------------------------

def test_every_radial_law_spans_the_unit_interval():
    for name in RADIAL_LAWS:
        r = radial_law(name, 20)
        assert r.shape == (20,)
        assert r.min() == pytest.approx(0.0)
        assert r.max() == pytest.approx(1.0)
        assert np.all(np.diff(r) >= -1e-12), f"{name} is not monotonic"


def test_sqrt_law_is_uniform_in_area():
    """The defining property: equal area between successive radii.

    This is the null hypothesis the other laws are measured against, so it has
    to actually be what it claims.
    """
    r = radial_law("sqrt", 200)
    area = r ** 2
    steps = np.diff(area)
    assert steps.std() / steps.mean() < 1e-9


def test_linear_law_is_uniform_in_radius():
    r = radial_law("linear", 100)
    steps = np.diff(r)
    assert steps.std() / steps.mean() < 1e-9


def test_power_law_exponent_is_respected():
    t = np.linspace(0, 1, 50)
    for gamma in (0.5, 2.0, 3.0):
        r = radial_law("power", 50, gamma=gamma)
        assert r == pytest.approx(t ** gamma)


def test_exponential_law_concentrates_outward():
    """A log-spiral radial law puts most antennas at large radius."""
    r = radial_law("exponential", 101, b=3.0)
    assert np.median(r) < 0.5


def test_radial_law_rejects_unknown_names_and_tiny_arrays():
    with pytest.raises(ValueError):
        radial_law("quadratic-ish", 10)
    with pytest.raises(ValueError):
        radial_law("linear", 1)


# --------------------------------------------------------------------------
# angular laws
# --------------------------------------------------------------------------

def test_golden_angle_has_the_right_value():
    """2 pi / phi^2, about 137.507 degrees."""
    assert np.rad2deg(GOLDEN_ANGLE) == pytest.approx(137.50776, abs=1e-4)


def test_golden_angle_never_repeats_a_bearing():
    th = angular_law("golden", 200)
    d = np.abs(th[:, None] - th[None, :])
    off = d[~np.eye(len(th), dtype=bool)]
    assert off.min() > 1e-3


def test_rational_law_is_deliberately_degenerate():
    """It revisits only q bearings, which is the point of including it."""
    th = angular_law("rational", 60, p=2, q=5)
    unique = np.unique(np.round(th, 9))
    assert len(unique) == 5


def test_multi_arm_law_uses_exactly_its_arms():
    th = angular_law("multi_arm", 60, arms=4)
    assert len(np.unique(np.round(th, 9))) == 4


def test_uniform_law_spreads_evenly_over_the_circle():
    n = 36
    th = np.sort(angular_law("uniform", n))
    gaps = np.diff(th)
    assert gaps.std() < 1e-9
    assert gaps.mean() == pytest.approx(2 * np.pi / n, rel=1e-6)


def test_random_law_is_reproducible_from_a_seed():
    a = angular_law("random", 30, rng=np.random.default_rng(3))
    b = angular_law("random", 30, rng=np.random.default_rng(3))
    assert np.array_equal(a, b)


def test_angular_law_rejects_unknown_names():
    with pytest.raises(ValueError):
        angular_law("fibonacci-ish", 10)


# --------------------------------------------------------------------------
# assembled layouts
# --------------------------------------------------------------------------

def test_layout_respects_the_requested_radial_extent():
    for radial in RADIAL_LAWS:
        for angular in ANGULAR_LAWS:
            xy = layout_from_laws(radial, angular, 24, 30.0, 2000.0, seed=1)
            r = np.hypot(xy[:, 0], xy[:, 1])
            assert r.min() == pytest.approx(30.0, rel=1e-6)
            assert r.max() == pytest.approx(2000.0, rel=1e-6)


def test_every_combination_is_constructible():
    for radial in RADIAL_LAWS:
        for angular in ANGULAR_LAWS:
            xy = layout_from_laws(radial, angular, 16, 20.0, 1000.0, seed=0)
            assert xy.shape == (16, 2)
            assert np.all(np.isfinite(xy))


def test_jitter_perturbs_but_does_not_relocate():
    base = layout_from_laws("sqrt", "golden", 24, 30.0, 2000.0, seed=5)
    jit = layout_from_laws("sqrt", "golden", 24, 30.0, 2000.0,
                           jitter_m=10.0, seed=5)
    shift = np.hypot(*(jit - base).T)
    assert shift.max() > 0.0
    assert shift.mean() < 60.0      # a few sigma of a 10 m jitter


def test_zero_jitter_is_deterministic():
    a = layout_from_laws("linear", "uniform", 20, 10.0, 500.0, seed=2)
    b = layout_from_laws("linear", "uniform", 20, 10.0, 500.0, seed=2)
    assert np.array_equal(a, b)


# --------------------------------------------------------------------------
# the variance decomposition
# --------------------------------------------------------------------------

def _anova():
    mod = __import__("16_factorial_layout")
    return mod.two_way_anova


def _build(cell_values, n_rep=4, noise=0.0, seed=0):
    """Replicated observations from a table of true cell means."""
    rng = np.random.default_rng(seed)
    vals, ri, ci = [], [], []
    cell_values = np.asarray(cell_values, dtype=float)
    for i in range(cell_values.shape[0]):
        for j in range(cell_values.shape[1]):
            for _ in range(n_rep):
                vals.append(cell_values[i, j] + rng.normal(0.0, noise))
                ri.append(i)
                ci.append(j)
    return vals, ri, ci, cell_values.shape


def test_a_pure_row_effect_is_attributed_to_rows():
    anova = _anova()
    table = np.tile(np.array([[1.0], [5.0], [9.0]]), (1, 4))
    vals, ri, ci, shape = _build(table, noise=1e-9)
    a = anova(vals, ri, ci, *shape)
    assert a["eta2_radial"] == pytest.approx(1.0, abs=1e-6)
    assert a["eta2_angular"] == pytest.approx(0.0, abs=1e-6)
    assert a["eta2_interaction"] == pytest.approx(0.0, abs=1e-6)


def test_a_pure_column_effect_is_attributed_to_columns():
    anova = _anova()
    table = np.tile(np.array([[2.0, 4.0, 8.0, 16.0]]), (3, 1))
    vals, ri, ci, shape = _build(table, noise=1e-9)
    a = anova(vals, ri, ci, *shape)
    assert a["eta2_angular"] == pytest.approx(1.0, abs=1e-6)
    assert a["eta2_radial"] == pytest.approx(0.0, abs=1e-6)


def test_a_pure_interaction_is_not_mistaken_for_a_main_effect():
    """The case that would most damage the conclusion if it were wrong."""
    anova = _anova()
    table = np.array([[1.0, -1.0], [-1.0, 1.0]])     # zero row and col means
    vals, ri, ci, shape = _build(table, noise=1e-9)
    a = anova(vals, ri, ci, *shape)
    assert a["eta2_interaction"] == pytest.approx(1.0, abs=1e-6)
    assert a["eta2_radial"] == pytest.approx(0.0, abs=1e-6)
    assert a["eta2_angular"] == pytest.approx(0.0, abs=1e-6)


def test_the_components_sum_to_the_whole():
    anova = _anova()
    rng = np.random.default_rng(11)
    table = rng.normal(0.0, 3.0, size=(5, 5))
    vals, ri, ci, shape = _build(table, n_rep=5, noise=0.7, seed=4)
    a = anova(vals, ri, ci, *shape)
    total = (a["eta2_radial"] + a["eta2_angular"]
             + a["eta2_interaction"] + a["eta2_residual"])
    assert total == pytest.approx(1.0, abs=1e-9)


def test_pure_noise_lands_in_the_residual():
    anova = _anova()
    table = np.zeros((4, 4))
    vals, ri, ci, shape = _build(table, n_rep=8, noise=1.0, seed=9)
    a = anova(vals, ri, ci, *shape)
    assert a["eta2_residual"] > 0.5
    assert a["f_radial"] < 10.0        # no real effect to find
