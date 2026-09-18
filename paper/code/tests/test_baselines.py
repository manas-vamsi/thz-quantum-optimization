import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.interferometry.baselines import (
    baseline_lengths,
    baseline_matrix,
    compute_baselines,
    n_baselines,
    pair_indices,
)

SQUARE = np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])


def test_baseline_count():
    for n in (2, 5, 10, 50):
        xy = golden_spiral_layout(n)
        assert compute_baselines(xy)[0].shape == (n_baselines(n), 2)
        assert n_baselines(n) == n * (n - 1) // 2


def test_pairs_are_unique_and_ordered():
    i, j = pair_indices(6)
    assert np.all(i < j)
    assert len({(a, b) for a, b in zip(i.tolist(), j.tolist())}) == n_baselines(6)


def test_baseline_convention():
    b, i, j = compute_baselines(SQUARE)
    assert np.allclose(b[0], SQUARE[1] - SQUARE[0])
    assert (i[0], j[0]) == (0, 1)


def test_known_square_lengths():
    d = np.sort(baseline_lengths(SQUARE))
    assert np.allclose(d, [10, 10, 10, 10, np.sqrt(200), np.sqrt(200)])


def test_translation_invariance():
    xy = golden_spiral_layout(12)
    a = compute_baselines(xy)[0]
    b = compute_baselines(xy + np.array([123.0, -456.0]))[0]
    assert np.allclose(a, b)


def test_rotation_preserves_lengths():
    xy = golden_spiral_layout(12)
    t = 0.9
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    assert np.allclose(baseline_lengths(xy), baseline_lengths(xy @ rot.T))


def test_baseline_matrix_is_symmetric_with_zero_diagonal():
    m = baseline_matrix(golden_spiral_layout(8))
    assert np.allclose(m, m.T)
    assert np.allclose(np.diag(m), 0.0)


def test_matrix_matches_pair_list():
    xy = golden_spiral_layout(9)
    i, j = pair_indices(9)
    assert np.allclose(baseline_matrix(xy)[i, j], baseline_lengths(xy))


def test_no_nans():
    assert np.all(np.isfinite(baseline_lengths(golden_spiral_layout(40))))


def test_single_antenna_has_no_baselines():
    assert compute_baselines(np.zeros((1, 2)))[0].shape == (0, 2)


def test_rejects_bad_layout():
    with pytest.raises(ValueError):
        compute_baselines(np.zeros(5))
