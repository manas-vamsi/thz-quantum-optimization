"""Tests for end-to-end imaging.

The transform pair is checked first and hardest. Every fidelity number in this
package is a difference between two images, so an inconsistency between the
forward and inverse transforms would not announce itself -- it would quietly
bias every array comparison by the same unknown factor.
"""

import numpy as np
import pytest

from thz_opt.imaging import (
    disc_with_gap_image,
    fidelity_metrics,
    gaussian_source,
    hogbom_clean,
    observe,
    point_sources,
    predict_visibilities,
    restore,
)
from thz_opt.imaging.clean import fit_restoring_beam
from thz_opt.imaging.sky import image_coordinates
from thz_opt.interferometry.uv import UVGrid

GRID = UVGrid(uv_max=5e5, n_cells=64)


def _all_cells(grid):
    """Every cell centre: cell k is centred on ``(k - n/2) * cell_size``."""
    c = (np.arange(grid.n_cells) - grid.n_cells // 2) * grid.cell_size
    uu, vv = np.meshgrid(c, c, indexing="ij")
    return np.column_stack([uu.ravel(), vv.ravel()])


# --------------------------------------------------------------------------
# the transform pair
# --------------------------------------------------------------------------

def test_forward_and_inverse_transforms_are_exact_inverses():
    sky = point_sources(GRID, [(0.0, 0.0), (2.0, -1.5)], [1.0, 0.4])
    v = predict_visibilities(sky)
    back = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(v))).real
    assert np.abs(back - sky.image).max() < 1e-12


def test_complete_sampling_returns_the_sky_and_a_delta_beam():
    """If every cell is observed the array is perfect and must vanish.

    This is the end-to-end version of the previous test: it exercises
    ``observe`` itself, including the gridding and the normalisation, not just
    the raw FFT calls.
    """
    sky = point_sources(GRID, [(0.0, 0.0), (1.0, 2.0)], [1.0, 0.7])
    dirty, beam, info = observe(sky, _all_cells(GRID), GRID)
    assert info["occupied_cells"] == GRID.n_cells ** 2
    assert np.abs(dirty - sky.image).max() < 1e-12
    assert beam.max() == pytest.approx(1.0)
    assert np.count_nonzero(beam > 1e-9) == 1        # a true delta


def test_a_point_source_at_the_centre_gives_a_flat_visibility_amplitude():
    """The Fourier pair every interferometry text starts from."""
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    v = predict_visibilities(sky)
    assert np.allclose(np.abs(v), 1.0)


def test_image_grid_is_the_fourier_partner_of_the_uv_grid():
    (_, _), pix = image_coordinates(GRID)
    assert pix == pytest.approx(1.0 / (GRID.n_cells * GRID.cell_size))


# --------------------------------------------------------------------------
# sky models
# --------------------------------------------------------------------------

def test_sky_models_carry_the_flux_they_are_given():
    g = gaussian_source(GRID, fwhm_arcsec=2.0, flux=3.0)
    assert g.total_flux == pytest.approx(3.0, rel=1e-6)
    d = disc_with_gap_image(GRID, flux=2.0)
    assert d.total_flux == pytest.approx(2.0, rel=1e-6)
    assert np.all(d.image >= 0.0)


def test_a_source_outside_the_field_is_rejected_not_wrapped():
    """Silently aliasing a source to the other edge would be worse than failing."""
    with pytest.raises(ValueError):
        point_sources(GRID, [(1e6, 0.0)], [1.0])


def test_the_disc_image_shows_the_gap():
    from thz_opt.science.disk_gap import HL_TAU

    sky = disc_with_gap_image(GRID)
    (ll, mm), pix = image_coordinates(GRID)
    r_au = np.hypot(ll, mm) / (np.pi / (180 * 3600)) * HL_TAU.distance_pc
    at_gap = np.abs(r_au - HL_TAU.gap_radius_au) < 3.0
    outside = np.abs(r_au - HL_TAU.gap_radius_au - 20.0) < 3.0
    if at_gap.any() and outside.any():
        assert sky.image[at_gap].mean() < sky.image[outside].mean()


# --------------------------------------------------------------------------
# CLEAN
# --------------------------------------------------------------------------

def test_clean_recovers_a_point_source_from_a_real_array():
    """With adequate coverage, a point source should come back almost exactly."""
    rng = np.random.default_rng(0)
    n = GRID.n_cells
    cells = _all_cells(GRID)
    keep = rng.random(len(cells)) < 0.3          # 30 % sparse sampling
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    dirty, beam, _ = observe(sky, cells[keep], GRID)
    comps, res, info = hogbom_clean(dirty, beam, gain=0.1, n_iter=2000)
    restored, _ = restore(comps, res, beam)
    m = fidelity_metrics(restored, sky, beam, res)
    assert m["flux_recovery"] == pytest.approx(1.0, abs=0.15)
    assert m["rms_error_relative"] < 0.05
    assert info["converged"]


def test_clean_on_a_delta_beam_is_the_identity():
    """Nothing to deconvolve means nothing should change."""
    sky = point_sources(GRID, [(0.0, 0.0), (3.0, 1.0)], [1.0, 0.5])
    dirty, beam, _ = observe(sky, _all_cells(GRID), GRID)
    comps, res, _ = hogbom_clean(dirty, beam, gain=1.0, n_iter=50)
    assert np.abs(comps - sky.image).max() < 1e-6
    assert np.abs(res).max() < 1e-6


def test_clean_conserves_flux_between_components_and_residual():
    rng = np.random.default_rng(3)
    cells = _all_cells(GRID)
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    dirty, beam, _ = observe(sky, cells[rng.random(len(cells)) < 0.4], GRID)
    comps, res, _ = hogbom_clean(dirty, beam, gain=0.1, n_iter=500)
    # every subtraction removes amp * psf, and sum(psf) is preserved, so the
    # residual plus what the components account for must equal the dirty map
    assert res.sum() + comps.sum() * beam.sum() == pytest.approx(dirty.sum(), abs=1e-6)


def test_restoring_beam_is_narrower_than_the_field_and_positive():
    rng = np.random.default_rng(5)
    cells = _all_cells(GRID)
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    _, beam, _ = observe(sky, cells[rng.random(len(cells)) < 0.25], GRID)
    clean_beam, info = fit_restoring_beam(beam)
    assert 0.5 <= info["sigma_px"] < GRID.n_cells / 4
    assert clean_beam.max() == pytest.approx(1.0)
    # a Gaussian underflows to exactly zero in the far corners; what matters is
    # that it is non-negative everywhere and peaked at the centre
    assert np.all(clean_beam >= 0)
    c = GRID.n_cells // 2
    assert clean_beam[c, c] == clean_beam.max()


def test_clean_stops_on_the_noise_threshold():
    """The loop must not keep going once it is subtracting noise."""
    rng = np.random.default_rng(7)
    cells = _all_cells(GRID)
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    dirty, beam, _ = observe(sky, cells[rng.random(len(cells)) < 0.3], GRID,
                             snr=20, seed=1)
    _, _, info = hogbom_clean(dirty, beam, gain=0.1, n_iter=100000)
    assert info["converged"]
    assert info["iterations"] < 100000


# --------------------------------------------------------------------------
# what the fidelity metrics are for
# --------------------------------------------------------------------------

def test_missing_short_spacings_lose_extended_flux():
    """The failure a UV-cell count cannot see.

    An array with a hole at the centre of the UV plane resolves out smooth
    emission and reports a fraction of the source, while still occupying plenty
    of cells. This is the reason image fidelity is measured at all.
    """
    cells = _all_cells(GRID)
    q = np.hypot(cells[:, 0], cells[:, 1])
    inner_cut = q > 0.25 * GRID.uv_max          # no short baselines
    full = q <= 1e30

    sky = gaussian_source(GRID, fwhm_arcsec=3.0)
    recovered = {}
    for tag, mask in (("with short spacings", full), ("without", inner_cut)):
        dirty, beam, _ = observe(sky, cells[mask], GRID)
        comps, res, _ = hogbom_clean(dirty, beam, gain=0.1, n_iter=3000)
        restored, _ = restore(comps, res, beam)
        recovered[tag] = fidelity_metrics(restored, sky, beam, res)["flux_recovery"]

    assert recovered["with short spacings"] > 0.9
    assert recovered["without"] < 0.6
    assert recovered["without"] < recovered["with short spacings"]


def test_a_perfect_reconstruction_scores_perfectly():
    """Complete sampling, cleaned and restored, must reproduce the model.

    The comparison has to be made against the *restored* image, not the dirty
    one. ``fidelity_metrics`` smooths the model by the restoring beam, because
    charging an array for resolution it never claimed would be unfair; handing
    it an unsmoothed dirty map therefore compares two different resolutions and
    reports a spurious 6 % error even for a perfect array.
    """
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    dirty, beam, _ = observe(sky, _all_cells(GRID), GRID)
    comps, res, _ = hogbom_clean(dirty, beam, gain=1.0, n_iter=20)
    restored, _ = restore(comps, res, beam)
    m = fidelity_metrics(restored, sky, beam, res)
    assert m["rms_error_relative"] < 1e-6
    assert m["flux_recovery"] == pytest.approx(1.0, abs=1e-6)


def test_noise_degrades_dynamic_range_monotonically():
    rng = np.random.default_rng(11)
    cells = _all_cells(GRID)
    mask = rng.random(len(cells)) < 0.3
    sky = point_sources(GRID, [(0.0, 0.0)], [1.0])
    drs = []
    for snr in (10, 100, 1000):
        dirty, beam, _ = observe(sky, cells[mask], GRID, snr=snr, seed=2)
        comps, res, _ = hogbom_clean(dirty, beam, gain=0.1, n_iter=2000)
        restored, _ = restore(comps, res, beam)
        drs.append(fidelity_metrics(restored, sky, beam, res)["dynamic_range"])
    assert drs[0] < drs[1] < drs[2]


# --------------------------------------------------------------------------
# instrumental effects
# --------------------------------------------------------------------------

from thz_opt.imaging import (  # noqa: E402
    apply_primary_beam,
    primary_beam,
    primary_beam_fwhm_arcsec,
    radiometer_noise,
    sefd,
    smear_image,
    system_temperature,
)

ARCSEC = np.pi / (180.0 * 3600.0)


def test_primary_beam_width_follows_lambda_over_d():
    lam = 299_792_458.0 / 230e9
    w8 = primary_beam_fwhm_arcsec(8.0, lam)
    assert w8 == pytest.approx(1.13 * lam / 8.0 / ARCSEC)
    # twice the dish, half the beam
    assert primary_beam_fwhm_arcsec(16.0, lam) == pytest.approx(w8 / 2)
    # twice the wavelength, twice the beam
    assert primary_beam_fwhm_arcsec(8.0, 2 * lam) == pytest.approx(2 * w8)


def test_primary_beam_is_unity_at_centre_and_half_at_the_fwhm():
    lam = 299_792_458.0 / 230e9
    a = primary_beam(GRID, dish_m=8.0, wavelength_m=lam)
    c = GRID.n_cells // 2
    assert a[c, c] == pytest.approx(1.0)

    fwhm = primary_beam_fwhm_arcsec(8.0, lam) * ARCSEC
    pix = 1.0 / (GRID.n_cells * GRID.cell_size)
    k = int(round(0.5 * fwhm / pix))
    if 0 < k < c:
        assert a[c + k, c] == pytest.approx(0.5, abs=0.02)


def test_primary_beam_attenuates_an_off_axis_source():
    lam = 299_792_458.0 / 230e9
    pix_arcsec = 1.0 / (GRID.n_cells * GRID.cell_size) / ARCSEC
    sky = point_sources(GRID, [(0.0, 0.0), (20 * pix_arcsec, 0.0)], [1.0, 1.0])
    attenuated = apply_primary_beam(sky, 8.0, lam, GRID)
    c = GRID.n_cells // 2
    assert attenuated.image[c, c] == pytest.approx(1.0)
    assert attenuated.image[c + 20, c] < 1.0
    assert attenuated.total_flux < sky.total_flux


def test_zero_smearing_is_exactly_the_identity():
    """The no-smear limit must not leak flux to interpolation.

    An earlier polar-resampling implementation lost up to 24 % of a point
    source here while applying no actual smearing, so this is asserted exactly
    rather than approximately.
    """
    sky = point_sources(GRID, [(0.0, 0.0), (3.0, -2.0)], [1.0, 0.5])
    out = smear_image(sky.image, fractional_bandwidth=0.0, integration_s=0.0)
    assert np.array_equal(out, sky.image)


def test_smearing_grows_with_distance_from_the_phase_centre():
    """The defining property: smearing defines a field of view."""
    pix_arcsec = 1.0 / (GRID.n_cells * GRID.cell_size) / ARCSEC
    radii = (0, 8, 16, 24)
    sky = point_sources(GRID, [(r * pix_arcsec, 0.0) for r in radii],
                        [1.0] * len(radii))
    out = smear_image(sky.image, fractional_bandwidth=0.08)
    c = GRID.n_cells // 2
    peaks = [out[c + r, c] for r in radii]
    assert peaks[0] == pytest.approx(1.0, abs=0.05)   # centre barely touched
    assert peaks == sorted(peaks, reverse=True)       # monotonic falloff
    assert peaks[-1] < 0.6 * peaks[0]


def test_smearing_conserves_total_flux():
    """Smearing redistributes brightness; it does not destroy it."""
    sky = gaussian_source(GRID, fwhm_arcsec=2.0, flux=1.0)
    for bw, t in ((0.05, 0.0), (0.0, 600.0), (0.05, 600.0)):
        out = smear_image(sky.image, fractional_bandwidth=bw, integration_s=t,
                          declination_deg=30.0)
        assert out.sum() == pytest.approx(1.0, rel=0.03)


def test_time_smearing_vanishes_at_the_pole():
    """A source at the celestial pole has no tangential smear: cos(dec) = 0."""
    sky = point_sources(GRID, [(5.0, 0.0)], [1.0])
    at_pole = smear_image(sky.image, integration_s=600.0, declination_deg=90.0)
    assert np.array_equal(at_pole, sky.image)


# --------------------------------------------------------------------------
# noise on a physical scale
# --------------------------------------------------------------------------

def test_system_temperature_rises_with_opacity_and_airmass():
    clear = system_temperature(0.05, elevation_deg=90.0)
    thick = system_temperature(0.40, elevation_deg=90.0)
    low = system_temperature(0.05, elevation_deg=30.0)
    assert thick["t_sys_k"] > clear["t_sys_k"]
    assert low["t_sys_k"] > clear["t_sys_k"]
    assert clear["transmission"] > thick["transmission"]
    assert low["airmass"] == pytest.approx(2.0)


def test_zero_opacity_leaves_only_receiver_and_cmb():
    t = system_temperature(0.0, elevation_deg=90.0, t_receiver_k=60.0)
    assert t["transmission"] == pytest.approx(1.0)
    assert t["t_sys_k"] == pytest.approx(60.0 + 2.73, abs=0.01)


def test_invalid_elevation_is_rejected():
    with pytest.raises(ValueError):
        system_temperature(0.1, elevation_deg=0.0)


def test_sefd_scales_as_tsys_over_area():
    a = sefd(100.0, 8.0)
    assert sefd(200.0, 8.0) == pytest.approx(2 * a)       # twice Tsys
    assert sefd(100.0, 16.0) == pytest.approx(a / 4)      # twice the diameter


def test_radiometer_noise_follows_the_root_law():
    s = sefd(120.0, 8.0)
    base = radiometer_noise(s, 16, 8e9, 3600.0)["sigma_jy"]
    longer = radiometer_noise(s, 16, 8e9, 14400.0)["sigma_jy"]
    assert longer == pytest.approx(base / 2)              # 4x time, half noise
    wider = radiometer_noise(s, 16, 32e9, 3600.0)["sigma_jy"]
    assert wider == pytest.approx(base / 2)               # 4x bandwidth
    assert radiometer_noise(s, 16, 8e9, 3600.0)["n_baselines"] == 120


def test_an_interferometer_needs_two_antennas():
    with pytest.raises(ValueError):
        radiometer_noise(1000.0, 1, 8e9, 60.0)
