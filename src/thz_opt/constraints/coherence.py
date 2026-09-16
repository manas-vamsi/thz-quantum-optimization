"""Atmospheric phase coherence -- the physical replacement for the placeholder.

:mod:`thz_opt.constraints.phase` implements a graph "stepping-stone" surrogate
and is explicitly labelled a placeholder.  This module implements the standard
millimetre/submillimetre treatment instead, which has two advantages: it is
derived from measured atmospheric behaviour, and it is **exactly quadratic in
the pad-selection variables** (the graph condition was not).

The model
---------
Tropospheric water vapour follows Kolmogorov turbulence, so the RMS excess path
length on a baseline of length ``b`` follows a broken power law (Carilli &
Holdaway 1999, Radio Science 34, 817):

    sigma_path(b) = kappa * sigma_1km * (b / 1 km) ** alpha

    alpha = 5/6   for b <~ L_3D      (3D turbulence, L_3D ~ 0.5-2 km)
    alpha = 1/3   for L_3D <~ b <~ L_out   (2D turbulence)
    saturated     for b >~ L_out     (outer scale, ~5-10 km)

with the power law made continuous at each breakpoint.  ``kappa <= 1`` is the
residual factor left after phase referencing or water-vapour radiometry; at
ALMA this residual is large at high frequency (>200 um path on 10 km baselines),
so ``kappa`` is not a small number and should be measured, not assumed.

The RMS phase is ``sigma_phi = 2 pi sigma_path / lambda``, and for Gaussian
phase errors the measured visibility amplitude is reduced by the standard
decorrelation factor

    gamma = exp(-sigma_phi**2 / 2)

``gamma`` in (0, 1] is the fraction of a baseline's coherence that survives.
Because it depends only on the pair, the penalty

    H_phase = lambda_ph * sum_{i<j} (1 - gamma_ij) x_i x_j

is exactly quadratic -- one entry per ``Q_ij``, no auxiliary variables.

Everything here is a *model*.  The constants are typical values for a good
high site, not measurements of any particular site, and they are exposed as
parameters for that reason.
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

    Defaults are order-of-magnitude values for a good high-altitude site; they
    are not a measurement of any specific site.
    """

    sigma_1km_m: float = 1.0e-3      # RMS path on a 1 km baseline, metres
    alpha_3d: float = 5.0 / 6.0      # 3D Kolmogorov exponent
    alpha_2d: float = 1.0 / 3.0      # 2D exponent above the layer thickness
    l_3d_m: float = 1000.0           # 3D -> 2D breakpoint
    l_out_m: float = 6000.0          # outer scale; saturates beyond this
    kappa: float = 1.0               # residual after phase referencing / WVR

    def as_dict(self) -> dict:
        return {
            "sigma_1km_m": self.sigma_1km_m,
            "alpha_3d": self.alpha_3d,
            "alpha_2d": self.alpha_2d,
            "l_3d_m": self.l_3d_m,
            "l_out_m": self.l_out_m,
            "kappa": self.kappa,
            "model": "Kolmogorov broken power law + Gaussian decorrelation",
            "constants_measured": False,
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
