"""Multi-frequency synthesis -- extra UV coverage for free.

A physical baseline ``b`` maps to ``(u, v) = b / lambda``, so observing the same
baseline at several frequencies places it at several *radii* in the UV plane,
along the same position angle.  A real receiver observes a band, not a single
frequency, so this coverage exists whether or not the design model accounts for
it.  For terahertz receivers the fractional bandwidth is large enough that
ignoring it materially understates what an array actually samples.

The important structural point for this project: multi-frequency synthesis adds
**no decision variables at all**.  Each candidate pair simply contributes more
UV samples, so its per-cell occupancy becomes

    a_ck = sum_f a_ckf

and every objective built on ``a_ck`` -- sidelobe energy, coverage, density
matching -- stays exactly as quadratic as it was.  The QUBO grows in the values
of its coefficients, not in its size.

Caveats stated plainly
----------------------
* Combining channels this way assumes the sky brightness does not change across
  the band.  That is the standard multi-frequency-synthesis assumption and it is
  wrong for sources with strong spectral structure; real pipelines fit a
  spectral index per pixel (MFS with Taylor terms).  For array *design* the
  approximation is the usual one, but it is an approximation.
* Decorrelation is frequency dependent: ``sigma_phi = 2 pi sigma_path / lambda``,
  so the top of the band is always less coherent than the bottom.  When channels
  are weighted by coherence, use
  :func:`channel_coherence_weights` rather than treating the band as uniform.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .earth_rotation import ObservationConfig, layout_to_uv_tracks
from .uv import C_LIGHT

__all__ = [
    "FrequencyBand",
    "multi_frequency_uv",
    "channel_coherence_weights",
]


@dataclass(frozen=True)
class FrequencyBand:
    """A receiver band sampled by ``n_channels`` discrete frequencies.

    ``fractional_bandwidth`` is ``delta_nu / nu_centre``; 0 collapses to a
    single frequency and reproduces the monochromatic model exactly.
    """

    centre_hz: float = 3.0e11
    fractional_bandwidth: float = 0.25
    n_channels: int = 5

    def __post_init__(self) -> None:
        if self.centre_hz <= 0:
            raise ValueError("centre_hz must be positive")
        if not 0.0 <= self.fractional_bandwidth < 2.0:
            raise ValueError("fractional_bandwidth must be in [0, 2)")
        if self.n_channels < 1:
            raise ValueError("n_channels must be >= 1")

    @property
    def frequencies(self) -> np.ndarray:
        """Channel centre frequencies, linearly spaced across the band."""
        if self.n_channels == 1 or self.fractional_bandwidth == 0.0:
            return np.array([self.centre_hz])
        half = 0.5 * self.fractional_bandwidth * self.centre_hz
        return np.linspace(self.centre_hz - half, self.centre_hz + half, self.n_channels)

    @property
    def wavelengths(self) -> np.ndarray:
        return C_LIGHT / self.frequencies

    @property
    def uv_radius_ratio(self) -> float:
        """Ratio of the largest to smallest UV radius a fixed baseline spans.

        This is the radial smearing multi-frequency synthesis buys: a baseline
        that would be one point becomes a radial streak of this length ratio.
        """
        f = self.frequencies
        return float(f.max() / f.min())

    def as_dict(self) -> dict:
        return {
            "centre_hz": self.centre_hz,
            "fractional_bandwidth": self.fractional_bandwidth,
            "n_channels": self.n_channels,
            "frequencies_hz": self.frequencies.tolist(),
            "uv_radius_ratio": self.uv_radius_ratio,
        }


def multi_frequency_uv(
    xy: np.ndarray,
    config: ObservationConfig,
    band: FrequencyBand | None = None,
    include_conjugate: bool = True,
) -> np.ndarray:
    """Earth-rotation UV samples stacked over every channel of the band.

    ``band=None`` falls back to the single frequency in ``config``, so callers
    can switch multi-frequency on and off with one argument.
    """
    if band is None:
        return layout_to_uv_tracks(xy, config, include_conjugate)

    chunks = []
    for nu in band.frequencies:
        cfg = ObservationConfig(
            frequency_hz=float(nu),
            declination_deg=config.declination_deg,
            latitude_deg=config.latitude_deg,
            hour_angle_start_h=config.hour_angle_start_h,
            hour_angle_end_h=config.hour_angle_end_h,
            n_times=config.n_times,
        )
        chunks.append(layout_to_uv_tracks(xy, cfg, include_conjugate))
    return np.vstack(chunks)


def channel_coherence_weights(
    baseline_length_m: float | np.ndarray,
    band: FrequencyBand,
    coherence_config=None,
) -> np.ndarray:
    """Decorrelation factor per channel for a given baseline length.

    Shape ``(n_channels,)`` for a scalar length, or ``(n_baselines,
    n_channels)`` for an array.  Because ``sigma_phi`` scales as ``1/lambda``,
    the high-frequency end of the band decorrelates first -- so a wide band does
    not add coverage uniformly, it adds well-measured short-radius samples and
    poorly-measured long-radius ones.
    """
    from ..constraints.coherence import decorrelation_factor

    b = np.atleast_1d(np.asarray(baseline_length_m, dtype=float))
    lam = band.wavelengths
    out = np.empty((b.size, lam.size))
    for c, wl in enumerate(lam):
        out[:, c] = decorrelation_factor(b, float(wl), coherence_config)
    return out[0] if np.isscalar(baseline_length_m) or out.shape[0] == 1 else out
