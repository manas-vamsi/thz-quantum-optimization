import numpy as np
import pytest

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.constraints.coherence import (
    CoherenceConfig,
    coherence_penalty_matrix,
    coherence_summary,
    decorrelation_factor,
    path_rms,
    phase_rms,
)
from thz_opt.interferometry.baselines import pair_indices
from thz_opt.qubo.objective import QUBOTerms, build_qubo
from thz_opt.qubo.validation import validate_qubo

LAM = 1e-3


def test_path_rms_matches_the_reference_point():
    """By definition sigma_path = sigma_1km at a 1 km baseline."""
    cfg = CoherenceConfig(sigma_1km_m=1.2e-3, l_3d_m=1000.0)
    assert path_rms(1000.0, cfg) == pytest.approx(1.2e-3)


def test_path_rms_exponents_in_each_regime():
    """The measured ALMA exponents, not the idealised Kolmogorov ones."""
    cfg = CoherenceConfig.from_measured("wvr_corrected")
    assert path_rms(200.0, cfg) / path_rms(100.0, cfg) == pytest.approx(2 ** 0.60)
    assert path_rms(4000.0, cfg) / path_rms(2000.0, cfg) == pytest.approx(2 ** 0.29)

    uncorrected = CoherenceConfig.from_measured("uncorrected")
    assert path_rms(200.0, uncorrected) / path_rms(100.0, uncorrected) == pytest.approx(2 ** 0.65)
    assert path_rms(4000.0, uncorrected) / path_rms(2000.0, uncorrected) == pytest.approx(2 ** 0.22)

    # saturation is available but switched off by default, since none is
    # measured within ALMA's 16 km extent
    sat = CoherenceConfig(l_3d_m=1000.0, l_out_m=6000.0)
    assert path_rms(20000.0, sat) == pytest.approx(path_rms(6000.0, sat))


def test_model_reproduces_the_published_scaling_factors():
    """ALMA Memo 624 states the factors used to scale between its summary
    baselines. Recovering them is the check that this implementation matches
    the measured structure function rather than merely resembling it."""
    cfg = CoherenceConfig.from_measured("wvr_corrected")
    for short, long_, published in ((500.0, 1000.0, 1.51),
                                    (1000.0, 5000.0, 1.59),
                                    (5000.0, 10000.0, 1.22)):
        assert path_rms(long_, cfg) / path_rms(short, cfg) == pytest.approx(
            published, abs=0.02)


def test_measured_anchor_values():
    """Each condition must reproduce its published 1 km path RMS."""
    for condition, micron in (("wvr_corrected", 70.0), ("wvr_corrected_windy", 140.0),
                              ("uncorrected", 200.0), ("uncorrected_low_pwv", 115.0)):
        cfg = CoherenceConfig.from_measured(condition)
        assert path_rms(1000.0, cfg) * 1e6 == pytest.approx(micron)


def test_measured_exponents_are_shallower_than_kolmogorov():
    """A real consequence: textbook constants overstate coherence loss."""
    measured = CoherenceConfig.from_measured("wvr_corrected")
    textbook = CoherenceConfig(sigma_1km_m=measured.sigma_1km_m,
                               alpha_3d=5 / 6, alpha_2d=1 / 3)
    assert path_rms(10000.0, textbook) > path_rms(10000.0, measured)


def test_path_rms_is_continuous_at_the_breakpoints():
    cfg = CoherenceConfig(l_3d_m=1000.0, l_out_m=6000.0)
    for b in (1000.0, 6000.0):
        lo = path_rms(b * (1 - 1e-6), cfg)
        hi = path_rms(b * (1 + 1e-6), cfg)
        assert lo == pytest.approx(hi, rel=1e-4)


def test_path_rms_is_monotonic():
    b = np.logspace(0, 5, 200)
    p = path_rms(b)
    assert np.all(np.diff(p) >= -1e-15)


def test_phase_rms_scales_inversely_with_wavelength():
    cfg = CoherenceConfig()
    assert phase_rms(500.0, 1e-3, cfg) / phase_rms(500.0, 3e-3, cfg) == pytest.approx(3.0)
    with pytest.raises(ValueError):
        phase_rms(500.0, 0.0, cfg)


def test_decorrelation_bounds_and_monotonicity():
    g = decorrelation_factor(np.logspace(0, 4, 100), LAM)
    assert np.all((g > 0) & (g <= 1.0))
    assert np.all(np.diff(g) <= 1e-15)
    # a zero-length baseline is perfectly coherent
    assert decorrelation_factor(0.0, LAM) == pytest.approx(1.0)


def test_wvr_correction_improves_measured_coherence():
    """The measured effect of water-vapour-radiometer correction, at 300 GHz."""
    b = 5000.0
    raw = decorrelation_factor(b, LAM, CoherenceConfig.from_measured("uncorrected"))
    wvr = decorrelation_factor(b, LAM, CoherenceConfig.from_measured("wvr_corrected"))
    windy = decorrelation_factor(b, LAM, CoherenceConfig.from_measured("wvr_corrected_windy"))
    assert wvr > windy > raw
    assert wvr > 0.5          # a 5 km baseline stays usable after correction


def test_longer_wavelength_is_more_forgiving():
    cfg = CoherenceConfig()
    assert decorrelation_factor(1000.0, 3e-3, cfg) > decorrelation_factor(1000.0, 1e-3, cfg)


def test_penalty_matrix_shape_and_range():
    xy = golden_spiral_layout(12, r_min=20.0, r_max=1000.0)
    p = coherence_penalty_matrix(xy, LAM)
    assert p.shape == (12, 12)
    assert np.allclose(p, p.T)
    assert np.allclose(np.diag(p), 0.0)
    i, j = pair_indices(12)
    assert np.all((p[i, j] >= 0.0) & (p[i, j] <= 1.0))


def test_penalty_grows_with_baseline_length():
    xy = np.array([[0.0, 0.0], [50.0, 0.0], [2000.0, 0.0]])
    p = coherence_penalty_matrix(xy, LAM)
    assert p[0, 2] > p[0, 1]


def test_summary_reports_the_constants_as_measured():
    s = coherence_summary(golden_spiral_layout(10, r_min=20.0, r_max=1000.0), LAM)
    assert s["constants_measured"] is True
    assert "ALMA Memo 624" in s["source"]
    assert 0.0 <= s["coherence_mean"] <= 1.0
    assert 0.0 <= s["fraction_below_0p5"] <= 1.0


def test_unknown_condition_is_rejected():
    from thz_opt.constraints.atmosphere_data import phase_conditions

    with pytest.raises(ValueError):
        phase_conditions("perfect_weather")


def test_coherence_term_is_exactly_quadratic_in_the_qubo():
    """The whole point of this model: it drops straight into Q_ij."""
    M, N = 8, 4
    xy = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    terms = QUBOTerms(
        n_pads=M,
        n_select=N,
        phase=coherence_penalty_matrix(xy, LAM),
        lambda_select=10.0,
        lambda_phase=5.0,
    )
    result = validate_qubo(terms)
    assert result["within_tolerance"], result
    Q, _ = build_qubo(terms)
    assert Q.shape == (M, M)  # no auxiliary variables needed
