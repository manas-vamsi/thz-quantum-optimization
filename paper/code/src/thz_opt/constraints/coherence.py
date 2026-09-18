"""Atmospheric phase coherence -- the physical replacement for the placeholder.

:mod:`thz_opt.constraints.phase` implements a graph "stepping-stone" surrogate
and is explicitly labelled a placeholder.  This module implements the standard
millimetre/submillimetre treatment instead, which has two advantages: it is
derived from measured atmospheric behaviour, and it is **exactly quadratic in
the pad-selection variables** (the graph condition was not).

The model
---------
Tropospheric water vapour produces an excess path length whose RMS follows a
broken power law in baseline length (Carilli & Holdaway 1999, Radio Science 34,
817):

    sigma_path(b) = kappa * sigma_1km * (b / 1 km) ** alpha

The constants used here are **measured at the ALMA site**, from over 17 000
observations analysed in ALMA Memo 624 (Maud et al. 2023):

    sigma_1km = 70 um   (WVR corrected, wind < 10 m/s, 120 s timescale)
              = 140 um  (WVR corrected, wind > 10 m/s)
              = 200 um  (no WVR correction, median of all observations)
              = 115 um  (no WVR, below the 1.24 mm median PWV)

    alpha = 0.60 below a 1 km baseline, 0.29 above   (WVR corrected)
          = 0.65 below,                  0.22 above  (uncorrected)

The exponents are from the spatial structure function of Matsushita et al.
(2017) as adopted in that memo. They are notably **shallower than the idealised
Kolmogorov values** of 5/6 and 1/3: using the textbook numbers overstates how
fast coherence degrades with baseline length. No saturation is seen within
ALMA's 16 km extent, so the outer-scale branch is disabled by default.

``kappa`` remains available as a multiplier for exploring hypothetical extra
correction, but it is 1.0 by default because the correction is already in the
measured numbers.

The RMS phase is ``sigma_phi = 2 pi sigma_path / lambda``, and for Gaussian
phase errors the measured visibility amplitude is reduced by the standard
decorrelation factor

    gamma = exp(-sigma_phi**2 / 2)

``gamma`` in (0, 1] is the fraction of a baseline's coherence that survives.
Because it depends only on the pair, the penalty

    H_phase = lambda_ph * sum_{i<j} (1 - gamma_ij) x_i x_j

is exactly quadratic -- one entry per ``Q_ij``, no auxiliary variables.

The model form is standard; the constants are measured at one site (ALMA) and
should be replaced for any other. What remains unmeasured is temporal: these are
medians over thousands of observations at a 120 s timescale, so they describe
the site rather than a particular night, and the source memo notes its sample
omits the very worst conditions, in which no observation was attempted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..interferometry.baselines import baseline_matrix

__all__ = [
    "CoherenceConfig",
    "path_rms",
    "phase_rms",
    "decorrelation_factor",
    "coherence_penalty_matrix",
    "coherence_summary",
]


@dataclass(frozen=True)
class CoherenceConfig:
    """Parameters of the phase structure function.

    The defaults are the **measured** ALMA values, not textbook constants:
    a median path RMS of 70 um on a 1 km baseline over 120 s with
    water-vapour-radiometer correction applied, and structure-function
    exponents of 0.60 below the 1 km break and 0.29 above it. Source and the
    other measured observing conditions are in
    :mod:`thz_opt.constraints.atmosphere_data`; build a config for any of them
    with :meth:`from_measured`.

    Two things are worth noting about the measured exponents. They are
    shallower than the idealised Kolmogorov values of 5/6 and 1/3, so an
    analysis using the textbook numbers overestimates how quickly coherence is
    lost with baseline length. And the measured data show no saturation within
    ALMA's 16 km extent, so ``l_out_m`` defaults beyond it and the outer-scale
    branch is inactive unless deliberately enabled.
    """

    sigma_1km_m: float = 70.0e-6     # RMS path on a 1 km baseline, metres
    alpha_3d: float = 0.60           # measured exponent below the break
    alpha_2d: float = 0.29           # measured exponent above the break
    l_3d_m: float = 1000.0           # measured break scale
    l_out_m: float = 20000.0         # no saturation measured within ALMA's extent
    kappa: float = 1.0               # extra residual factor, 1.0 = as measured

    #: set when the constants came from :meth:`from_measured`
    condition: str = "wvr_corrected"
    measured: bool = True

    @classmethod
    def from_measured(cls, condition: str = "wvr_corrected", kappa: float = 1.0):
        """Config from the measured ALMA phase conditions.

        ``condition`` selects one of the observing regimes reported in ALMA
        Memo 624 -- ``uncorrected``, ``uncorrected_low_pwv``, ``wvr_corrected``
        or ``wvr_corrected_windy``. ``kappa`` scales the measured path RMS and
        exists only to explore hypothetical further correction; leave it at 1.0
        to use the site as measured.
        """
        from .atmosphere_data import BREAK_SCALE_M, phase_conditions

        c = phase_conditions(condition)
        return cls(sigma_1km_m=c["sigma_1km_m"], alpha_3d=c["alpha_short"],
                   alpha_2d=c["alpha_long"], l_3d_m=BREAK_SCALE_M,
                   l_out_m=20000.0, kappa=kappa, condition=condition,
                   measured=True)

    def as_dict(self) -> dict:
        return {
            "sigma_1km_m": self.sigma_1km_m,
            "sigma_1km_um": self.sigma_1km_m * 1e6,
            "alpha_3d": self.alpha_3d,
            "alpha_2d": self.alpha_2d,
            "l_3d_m": self.l_3d_m,
            "l_out_m": self.l_out_m,
            "kappa": self.kappa,
            "condition": self.condition,
            "model": "broken power law with measured exponents "
                     "+ Gaussian decorrelation",
            "constants_measured": self.measured,
            "source": "ALMA Memo 624 (Maud et al. 2023), arXiv:2304.08318",
        }


def path_rms(b: np.ndarray, config: CoherenceConfig | None = None) -> np.ndarray:
    """RMS excess path length in metres for baseline lengths ``b`` (metres).

    The three regimes are joined continuously: the 2D branch is scaled so it
    meets the 3D branch at ``l_3d_m``, and the saturated branch holds the value
    reached at ``l_out_m``.
    """
    cfg = config or CoherenceConfig()
    d = np.asarray(b, dtype=float)
    km = np.maximum(d, 0.0) / 1000.0
    l3 = cfg.l_3d_m / 1000.0
    lo = cfg.l_out_m / 1000.0

    inner = km ** cfg.alpha_3d
    # continuity at l_3d: outer branch matched to the inner value there
    scale_2d = (l3 ** cfg.alpha_3d) / (l3 ** cfg.alpha_2d) if l3 > 0 else 1.0
    outer = scale_2d * km ** cfg.alpha_2d
    saturated = scale_2d * lo ** cfg.alpha_2d

    shape = np.where(km <= l3, inner, np.where(km <= lo, outer, saturated))
    return cfg.kappa * cfg.sigma_1km_m * shape


def phase_rms(b: np.ndarray, wavelength: float,
              config: CoherenceConfig | None = None) -> np.ndarray:
    """RMS phase in radians, ``2 pi sigma_path / lambda``."""
    if wavelength <= 0:
        raise ValueError("wavelength must be positive")
    return 2.0 * np.pi * path_rms(b, config) / wavelength


def decorrelation_factor(b: np.ndarray, wavelength: float,
                         config: CoherenceConfig | None = None) -> np.ndarray:
    """``gamma = exp(-sigma_phi^2 / 2)``, the surviving coherence fraction.

    1 means fully coherent, 0 means the baseline measures nothing.  Monotonically
    decreasing in baseline length, which is the physical content: long baselines
    are harder at high frequency.
    """
    s = phase_rms(b, wavelength, config)
    return np.exp(-0.5 * s ** 2)


def coherence_penalty_matrix(xy: np.ndarray, wavelength: float,
                             config: CoherenceConfig | None = None) -> np.ndarray:
    """``(n, n)`` matrix of ``1 - gamma_ij`` -- the ``H_phase`` coefficients.

    Feed straight into the ``phase`` slot of
    :class:`thz_opt.qubo.objective.QUBOTerms`; the diagonal is zeroed since a
    pad has no baseline with itself.
    """
    d = baseline_matrix(xy)
    p = 1.0 - decorrelation_factor(d, wavelength, config)
    np.fill_diagonal(p, 0.0)
    return p


def coherence_summary(xy: np.ndarray, wavelength: float,
                      config: CoherenceConfig | None = None) -> dict:
    """Diagnostics of how hard the atmosphere makes this array at this lambda."""
    cfg = config or CoherenceConfig()
    from ..interferometry.baselines import baseline_lengths

    d = baseline_lengths(xy)
    g = decorrelation_factor(d, wavelength, cfg)
    return {
        "wavelength_m": wavelength,
        "n_baselines": int(d.size),
        "coherence_mean": float(g.mean()),
        "coherence_min": float(g.min()),
        "coherence_max": float(g.max()),
        "fraction_below_0p5": float((g < 0.5).mean()),
        "path_rms_at_dmax_m": float(path_rms(d.max(), cfg)),
        "phase_rms_at_dmax_rad": float(phase_rms(d.max(), wavelength, cfg)),
        **cfg.as_dict(),
    }
