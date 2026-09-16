import numpy as np
import pytest

from thz_opt.arrays.fibonacci import GOLDEN_ANGLE, fibonacci_layout, fibonacci_sequence
from thz_opt.arrays.golden_spiral import PHI, golden_spiral_layout
from thz_opt.arrays.random_array import random_layout, uniform_ring_layout
from thz_opt.arrays.validation import layout_extent, recentre, validate_layout

EXTENT = dict(r_min=10.0, r_max=1000.0)


def test_phi_value():
    assert PHI == pytest.approx(1.6180339887, abs=1e-9)
    assert PHI**2 == pytest.approx(PHI + 1.0)


def test_golden_angle():
    assert np.rad2deg(GOLDEN_ANGLE) == pytest.approx(137.507764, abs=1e-5)


def test_fibonacci_sequence_recurrence():
    f = fibonacci_sequence(12)
    assert list(f[:8]) == [1, 1, 2, 3, 5, 8, 13, 21]
    assert np.allclose(f[2:], f[1:-1] + f[:-2])


def test_fibonacci_ratio_converges_to_phi():
    f = fibonacci_sequence(40)
    assert f[-1] / f[-2] == pytest.approx(PHI, abs=1e-8)


def test_golden_spiral_growth_per_turn():
    """Before normalisation the radius must grow by exactly PHI per 2*pi."""
    xy = golden_spiral_layout(9, n_turns=4.0, normalize=False, a=1.0)
    r = np.hypot(xy[:, 0], xy[:, 1])
    # 9 points over 4 turns -> every second point is one full turn apart
    assert np.allclose(r[2:] / r[:-2], PHI, rtol=1e-10)


@pytest.mark.parametrize("variant", ["radial", "golden_angle", "rational_angle"])
def test_layouts_respect_radial_extent(variant):
    xy = fibonacci_layout(25, variant=variant, **EXTENT)
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert r.min() == pytest.approx(EXTENT["r_min"])
    assert r.max() == pytest.approx(EXTENT["r_max"])


def test_golden_spiral_respects_radial_extent():
    xy = golden_spiral_layout(25, **EXTENT)
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert r.min() == pytest.approx(EXTENT["r_min"])
    assert r.max() == pytest.approx(EXTENT["r_max"])


def test_shapes_and_finiteness():
    for xy in (
        golden_spiral_layout(17, **EXTENT),
        fibonacci_layout(17, **EXTENT),
        random_layout(17, seed=3, **EXTENT),
        uniform_ring_layout(17, **EXTENT),
    ):
        assert xy.shape == (17, 2)
        assert np.all(np.isfinite(xy))
        validate_layout(xy)


def test_random_layout_is_seeded():
    assert np.array_equal(random_layout(12, seed=7), random_layout(12, seed=7))
    assert not np.array_equal(random_layout(12, seed=7), random_layout(12, seed=8))


def test_rational_angle_is_periodic_but_golden_angle_is_not():
    """The qualitative Fibonacci/golden distinction: rational steps close."""
    k = 8
    step_rational = 2 * np.pi * fibonacci_sequence(k + 1)[k - 1] / fibonacci_sequence(k + 1)[k]
    f_k = int(fibonacci_sequence(k + 1)[k])
    assert np.mod(f_k * step_rational, 2 * np.pi) == pytest.approx(0.0, abs=1e-9)
    assert np.mod(f_k * GOLDEN_ANGLE, 2 * np.pi) > 1e-3


def test_orientation_is_a_rigid_rotation():
    a = golden_spiral_layout(20, **EXTENT)
    b = golden_spiral_layout(20, orientation=0.7, **EXTENT)
    assert np.allclose(np.hypot(*a.T), np.hypot(*b.T))


def test_recentre_and_extent():
    xy = random_layout(30, seed=1, **EXTENT)
    c = recentre(xy)
    assert np.allclose(c.mean(axis=0), 0.0, atol=1e-9)
    e = layout_extent(xy)
    assert e["n"] == 30 and e["r_max"] <= EXTENT["r_max"] + 1e-9


def test_validation_rejects_bad_input():
    with pytest.raises(ValueError):
        validate_layout(np.zeros((4, 3)))
    with pytest.raises(ValueError):
        validate_layout(np.array([[0.0, np.nan]]))
    with pytest.raises(ValueError):
        golden_spiral_layout(10, r_min=100.0, r_max=50.0)
    with pytest.raises(ValueError):
        fibonacci_layout(10, variant="not_a_variant")
