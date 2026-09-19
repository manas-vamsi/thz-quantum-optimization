import numpy as np
import pytest
from scipy.special import j1

from thz_opt.science.disk_gap import (
    ARCSEC,
    HL_TAU,
    PARAMETERS,
    DiskModel,
    brightness_profile,
    depth_error_bar,
    fisher_information,
    fisher_weight_profile,
    gap_derivative_profile,
    hankel_transform,
    resolution_requirement,
    uv_cell_weights,
    visibility_profile,
)
from thz_opt.interferometry.uv import UVGrid

LAM = 1e-3


def test_hankel_transform_matches_the_analytic_uniform_disc():
    """A uniform disc of angular radius R has V(q) = pi R^2 * 2 J1(x)/x.

    This is the only closed form available here, so it is the one check that
    the numerical transform -- on which every weight in this module depends --
    is correct rather than merely plausible.
    """
    m = DiskModel(r_in_au=1e-6, r_out_au=100.0, gap_depth=0.0, gamma=0.0)
    r = np.linspace(m.r_in_au, m.r_out_au, 20000)
    q = np.array([1e4, 5e4, 1e5, 2e5])
    numeric = hankel_transform(np.ones_like(r), r, q, m)
    x = 2.0 * np.pi * q * m.theta_out
    analytic = np.pi * m.theta_out ** 2 * 2.0 * j1(x) / x
    assert numeric == pytest.approx(analytic, rel=1e-6)


def test_flux_conserving_derivative_integrates_to_zero():
    """The defining property: deepening the gap redistributes, not removes.

    If this fails the weight collapses onto the shortest baseline, because
    total flux dominates the information.
    """
    r = np.linspace(HL_TAU.r_in_au, HL_TAU.r_out_au, 4096)
    for par in PARAMETERS:
        d = gap_derivative_profile(r, HL_TAU, flux_conserving=True, parameter=par)
        scale = np.trapezoid(np.abs(d) * r, r)
        assert abs(np.trapezoid(d * r, r)) < 1e-10 * scale


def test_unconstrained_derivative_is_the_degenerate_case():
    """Documents why flux conservation is needed, rather than asserting it."""
    r = np.linspace(HL_TAU.r_in_au, HL_TAU.r_out_au, 4096)
    d = gap_derivative_profile(r, HL_TAU, flux_conserving=False)
    assert np.trapezoid(d * r, r) < 0            # flux is removed

    q = np.logspace(3.5, 7.0, 400)
    loose = fisher_weight_profile(q, flux_conserving=False)
    tight = fisher_weight_profile(q, flux_conserving=True)
    # the degenerate weight peaks at the shortest baseline; the correct one
    # peaks at a resolved scale
    assert q[int(np.argmax(loose))] < q[int(np.argmax(tight))] / 10


def test_gap_radius_demands_longer_baselines_than_gap_depth():
    """The result that makes 'optimise for this source' under-specified."""
    q = np.logspace(3.5, 7.0, 900)
    peaks = {p: q[int(np.argmax(fisher_weight_profile(q, parameter=p)))]
             for p in PARAMETERS}
    assert peaks["radius"] > 3.0 * peaks["depth"]


def test_visibility_is_real_and_peaks_at_zero_spacing():
    q = np.array([0.0, 1e4, 1e5, 1e6])
    v = visibility_profile(q)
    assert np.all(np.isfinite(v))
    assert v[0] == pytest.approx(np.abs(v).max())


def test_fisher_information_is_additive_and_error_bar_shrinks():
    """More samples cannot hurt: F is a sum of non-negative terms."""
    q1 = np.logspace(4.5, 5.5, 50)
    q2 = np.logspace(4.5, 5.5, 200)
    assert fisher_information(q2) > fisher_information(q1)
    assert depth_error_bar(q2) < depth_error_bar(q1)
    assert depth_error_bar(np.array([])) == float("inf")


def test_error_bar_scales_with_noise():
    """sigma(p) is proportional to the per-visibility noise."""
    q = np.logspace(4.5, 5.5, 100)
    assert depth_error_bar(q, HL_TAU, 2.0) == pytest.approx(
        2.0 * depth_error_bar(q, HL_TAU, 1.0))


def test_resolution_requirement_is_self_consistent():
    req = resolution_requirement(HL_TAU, LAM)
    assert req["b_max_required_m"] > req["b_min_allowed_m"]
    # a gap 22x smaller than the disc needs a much longer baseline
    assert req["b_max_required_m"] / req["b_min_allowed_m"] > 10.0
    # both limits scale linearly with wavelength
    req2 = resolution_requirement(HL_TAU, 2 * LAM)
    assert req2["b_max_required_m"] == pytest.approx(2 * req["b_max_required_m"])


def test_angular_conversion_uses_the_au_per_pc_identity():
    """1 au at 1 pc is 1 arcsec, by definition of the parsec."""
    m = DiskModel(distance_pc=1.0)
    assert m.angular(1.0) == pytest.approx(ARCSEC)


def test_uv_cell_weights_match_the_grid():
    grid = UVGrid(uv_max=5e5, n_cells=16)
    w = uv_cell_weights(grid, HL_TAU)
    assert w.shape == (16 * 16,)
    assert np.all(w >= 0.0) and np.all(np.isfinite(w))
    assert w.max() > 0.0


def test_unknown_parameter_is_rejected():
    r = np.linspace(1.0, 120.0, 100)
    with pytest.raises(ValueError):
        gap_derivative_profile(r, HL_TAU, parameter="eccentricity")


def test_brightness_profile_has_a_gap_where_the_model_says():
    r = np.linspace(HL_TAU.r_in_au, HL_TAU.r_out_au, 4000)
    prof = brightness_profile(r)
    near = np.abs(r - HL_TAU.gap_radius_au) < 1.0
    far = np.abs(r - HL_TAU.gap_radius_au - 25.0) < 1.0
    assert prof[near].mean() < prof[far].mean()
    assert np.all(prof >= 0.0)
