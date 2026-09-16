import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.interferometry.baselines import compute_baselines, n_baselines
from thz_opt.interferometry.earth_rotation import (
    ObservationConfig,
    elevation_deg,
    enu_to_xyz,
    layout_to_uv_tracks,
)
from thz_opt.interferometry.psf import dirty_beam, sampling_function
from thz_opt.interferometry.uv import (
    C_LIGHT,
    UVGrid,
    baselines_to_uv,
    cell_ids,
    grid_indices,
    layout_to_uv,
    uv_occupancy,
    wavelength_from_frequency,
)

LAM = 1e-3  # 1 mm, illustrative only


def test_wavelength_conversion():
    assert wavelength_from_frequency(3e11) == pytest.approx(C_LIGHT / 3e11)
    assert wavelength_from_frequency(C_LIGHT) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        wavelength_from_frequency(0.0)


def test_uv_scaling_is_baseline_over_wavelength():
    b = np.array([[100.0, 0.0]])
    uv = baselines_to_uv(b, LAM, include_conjugate=False)
    assert uv[0, 0] == pytest.approx(100.0 / LAM)


def test_conjugate_points_are_added_in_pairs():
    xy = golden_spiral_layout(7)
    uv = layout_to_uv(xy, LAM)
    assert uv.shape[0] == 2 * n_baselines(7)
    half = uv.shape[0] // 2
    assert np.allclose(uv[:half], -uv[half:])
    assert np.allclose(uv.sum(axis=0), 0.0, atol=1e-6)


def test_grid_geometry():
    g = UVGrid(uv_max=100.0, n_cells=10)
    assert g.cell_size == pytest.approx(20.0)
    assert g.total_cells == 100
    assert g.edges[0] == -100.0 and g.edges[-1] == 100.0
    with pytest.raises(ValueError):
        UVGrid(uv_max=-1.0, n_cells=4)


def test_origin_sits_in_a_single_cell_and_binning_is_symmetric():
    g = UVGrid(uv_max=100.0, n_cells=10)
    iu, iv, _ = grid_indices(np.array([[0.0, 0.0]]), g)
    assert (iu[0], iv[0]) == (5, 5)
    pts = np.array([[37.0, -12.0]])
    a = grid_indices(pts, g)[0][0], grid_indices(pts, g)[1][0]
    b = grid_indices(-pts, g)[0][0], grid_indices(-pts, g)[1][0]
    assert (a[0] + b[0], a[1] + b[1]) == (10, 10)  # mirror about the centre cell


def test_occupancy_conserves_samples():
    xy = golden_spiral_layout(15)
    uv = layout_to_uv(xy, LAM)
    g = UVGrid.from_uv(uv, 64)
    occ, outside = uv_occupancy(uv, g)
    assert occ.sum() + outside == uv.shape[0]
    assert outside == 0


def test_samples_outside_extent_are_dropped_not_clipped():
    g = UVGrid(uv_max=1.0, n_cells=4)
    uv = np.array([[0.5, 0.5], [10.0, 0.0]])
    occ, outside = uv_occupancy(uv, g)
    assert outside == 1 and occ.sum() == 1


def test_cell_ids_match_grid_indices():
    xy = golden_spiral_layout(10)
    uv = layout_to_uv(xy, LAM)
    g = UVGrid.from_uv(uv, 32)
    iu, iv, _ = grid_indices(uv, g)
    assert np.array_equal(cell_ids(uv, g), iu * g.n_cells + iv)


def test_earth_rotation_reduces_to_snapshot_at_zenith():
    """dec == latitude and H == 0 must reproduce u = dx/lambda, v = dy/lambda."""
    xy = golden_spiral_layout(8)
    cfg = ObservationConfig(
        frequency_hz=C_LIGHT / LAM,
        declination_deg=-23.0,
        latitude_deg=-23.0,
        hour_angle_start_h=0.0,
        hour_angle_end_h=0.0,
        n_times=1,
    )
    a = layout_to_uv_tracks(xy, cfg)
    b = layout_to_uv(xy, LAM)
    assert np.allclose(a, b, rtol=1e-9, atol=1e-3)
    assert elevation_deg(cfg)[0] == pytest.approx(90.0)


def test_enu_to_xyz_pole_component():
    """A due-North baseline at latitude phi has Z = n cos(phi)."""
    phi = np.deg2rad(30.0)
    xyz = enu_to_xyz(np.array([[0.0, 100.0, 0.0]]), phi)
    assert xyz[0, 2] == pytest.approx(100.0 * np.cos(phi))
    assert xyz[0, 0] == pytest.approx(-100.0 * np.sin(phi))


def test_earth_rotation_track_count():
    xy = golden_spiral_layout(6)
    cfg = ObservationConfig(n_times=9)
    uv = layout_to_uv_tracks(xy, cfg)
    assert uv.shape == (2 * n_baselines(6) * 9, 2)


def test_earth_rotation_at_the_pole_traces_circles():
    """For dec = 90 deg the UV tracks are exact circles of constant radius."""
    xy = golden_spiral_layout(5)
    cfg = ObservationConfig(declination_deg=90.0, latitude_deg=30.0, n_times=32)
    uv = layout_to_uv_tracks(xy, cfg, include_conjugate=False)
    r = np.hypot(uv[:, 0], uv[:, 1]).reshape(n_baselines(5), 32)
    assert np.allclose(r.std(axis=1) / r.mean(axis=1), 0.0, atol=1e-9)


def test_psf_is_real_and_peaks_at_the_centre():
    xy = golden_spiral_layout(12)
    uv = layout_to_uv(xy, LAM)
    g = UVGrid.from_uv(uv, 64)
    psf, info = dirty_beam(uv, g)
    assert info["imag_residual_fraction"] < 1e-10
    assert psf.max() == pytest.approx(1.0)
    assert np.unravel_index(np.argmax(psf), psf.shape) == (32, 32)


def test_sampling_function_weightings():
    xy = golden_spiral_layout(10)
    uv = layout_to_uv(xy, LAM)
    g = UVGrid.from_uv(uv, 32)
    s_u = sampling_function(uv, g, "uniform")
    s_n = sampling_function(uv, g, "natural")
    assert set(np.unique(s_u)).issubset({0.0, 1.0})
    assert s_n.sum() == uv.shape[0]
    with pytest.raises(ValueError):
        sampling_function(uv, g, "robust")
