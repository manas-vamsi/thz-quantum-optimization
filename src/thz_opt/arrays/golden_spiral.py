"""Model 1 -- golden-ratio logarithmic spiral antenna layout.

Mathematical definition (see docs/mathematical_formulation.md Sec. 3)::

    r(theta) = a * PHI ** (theta / (2 * pi)),   PHI = (1 + sqrt(5)) / 2

so the radius grows by exactly one factor of PHI per full turn.  ``N`` antennas
are placed at angles sampled uniformly in ``theta`` over ``n_turns`` turns,
starting at ``theta_0``.

The raw spiral radius depends on ``a``; because every comparison in this
repository must use the same physical radial extent, the generated radii are
affinely rescaled onto ``[r_min, r_max]`` (``normalize=True``, the default).
This is an IMPLEMENTATION CHOICE, not a result from the source paper: it keeps
the *shape* of the logarithmic spiral (the relative radial spacing) while
pinning the two endpoints.
"""

from __future__ import annotations

import numpy as np

PHI = (1.0 + np.sqrt(5.0)) / 2.0


def golden_spiral_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    n_turns: float = 3.0,
    theta_0: float = 0.0,
    orientation: float = 0.0,
    normalize: bool = True,
    a: float = 1.0,
) -> np.ndarray:
    """Return an ``(n, 2)`` array of antenna positions in metres.

    Parameters
    ----------
    n : number of antennas.
    r_min, r_max : target radial extent in metres (used when ``normalize``).
    n_turns : number of full 2*pi turns spanned by the sampled angles.
    theta_0 : initial spiral angle in radians.
    orientation : rigid rotation of the whole layout in radians.
    normalize : affinely rescale radii onto ``[r_min, r_max]``.
    a : scale constant of the unnormalised spiral (ignored when ``normalize``).
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if r_max <= r_min:
        raise ValueError("r_max must exceed r_min")

    theta = theta_0 + np.linspace(0.0, 2.0 * np.pi * n_turns, n)
    r = a * PHI ** ((theta - theta_0) / (2.0 * np.pi))

    if normalize:
        r = _rescale(r, r_min, r_max)

    ang = theta + orientation
    return np.column_stack((r * np.cos(ang), r * np.sin(ang)))


def _rescale(r: np.ndarray, r_min: float, r_max: float) -> np.ndarray:
    """Affine map of ``r`` onto ``[r_min, r_max]`` (constant input -> r_min)."""
    lo, hi = float(np.min(r)), float(np.max(r))
    if hi - lo <= 0.0:
        return np.full_like(r, r_min, dtype=float)
    return r_min + (r_max - r_min) * (r - lo) / (hi - lo)
