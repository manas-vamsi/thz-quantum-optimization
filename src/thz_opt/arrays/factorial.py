"""Layouts built from an explicit radial law crossed with an angular law.

Why the existing comparison cannot answer its own question
----------------------------------------------------------
The layout shootout scores twelve named families -- golden spiral, Vogel
phyllotaxis, Reuleaux, hierarchical, Gaussian, random and so on -- against each
other at matched N and extent. It is a fair comparison of those twelve objects,
and it is silent about *why* any of them wins, because they differ in two ways
at once. The golden spiral changes both how radius grows with index and how
angle advances; so does phyllotaxis; so does a hierarchical array. When one of
them scores higher, nothing in the experiment says whether the radial law, the
angular law, or their combination is responsible.

That matters practically. "Use the golden angle" and "spread antennas uniformly
in area" are different pieces of advice, and a reader cannot extract either one
from a ranking of composite objects.

This module separates the two. A layout is

    r_k = R_max * f_radial(k / (N - 1))
    theta_k = f_angular(k, N)

so any radial law can be crossed with any angular law. Scoring the full grid
turns "which family is best" into two answerable questions: how much of the
variation is explained by the radial law, how much by the angular law, and how
much only by their interaction.

The laws
--------
Radial, all mapping ``t`` in [0, 1] onto a radius fraction:

``linear``
    ``r = t``. Uniform in radius, so the antenna surface density falls as 1/r
    and the array is centrally concentrated.
``sqrt``
    ``r = sqrt(t)``. Uniform in *area*, the natural null hypothesis.
``power``
    ``r = t**gamma``, the general case the first two are members of.
``exponential``
    ``r ~ exp(b t)``, rescaled to [0, 1]. This is the radial behaviour of a
    logarithmic spiral, isolated from its angular behaviour.
``gaussian``
    The inverse Gaussian CDF, which realises the Gaussian UV density Boone
    (2001) recommends.

Angular, mapping index ``k`` onto an angle:

``golden``
    ``k`` times the golden angle, 2 pi / phi^2. Never repeats.
``uniform``
    ``2 pi k / N``. Every antenna at a different bearing, evenly spaced.
``rational``
    ``2 pi k p / q``. Deliberately degenerate: it revisits only ``q`` distinct
    bearings, so it isolates what non-repetition is worth.
``random``
    Uniform random bearings, the null hypothesis.
``multi_arm``
    ``m`` straight arms, the configuration most existing arrays actually use.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "RADIAL_LAWS",
    "ANGULAR_LAWS",
    "radial_law",
    "angular_law",
    "layout_from_laws",
]

PHI = (1.0 + np.sqrt(5.0)) / 2.0
GOLDEN_ANGLE = 2.0 * np.pi / PHI ** 2

RADIAL_LAWS = ("linear", "sqrt", "power", "exponential", "gaussian")
ANGULAR_LAWS = ("golden", "uniform", "rational", "random", "multi_arm")


def radial_law(name: str, n: int, gamma: float = 2.0, b: float = 2.5,
               sigma: float = 0.45) -> np.ndarray:
    """Radius fractions in [0, 1] for ``n`` antennas under one radial law."""
    if n < 2:
        raise ValueError("need at least two antennas")
    t = np.linspace(0.0, 1.0, n)

    if name == "linear":
        r = t
    elif name == "sqrt":
        r = np.sqrt(t)
    elif name == "power":
        r = t ** gamma
    elif name == "exponential":
        r = (np.exp(b * t) - 1.0) / (np.exp(b) - 1.0)
    elif name == "gaussian":
        # inverse CDF of a half-Gaussian, giving a Gaussian radial density
        from scipy.special import erfinv
        q = np.linspace(1e-3, 1.0 - 1e-3, n)
        r = sigma * np.sqrt(2.0) * erfinv(q)
        r = r / r.max()
    else:
        raise ValueError(f"unknown radial law {name!r}; choose from {RADIAL_LAWS}")

    lo, hi = float(r.min()), float(r.max())
    return (r - lo) / (hi - lo) if hi > lo else r


def angular_law(name: str, n: int, p: int = 2, q: int = 5, arms: int = 3,
                rng=None) -> np.ndarray:
    """Bearings in radians for ``n`` antennas under one angular law."""
    k = np.arange(n, dtype=float)

    if name == "golden":
        return (k * GOLDEN_ANGLE) % (2.0 * np.pi)
    if name == "uniform":
        return (2.0 * np.pi * k / n) % (2.0 * np.pi)
    if name == "rational":
        return (2.0 * np.pi * k * p / q) % (2.0 * np.pi)
    if name == "random":
        rng = np.random.default_rng() if rng is None else rng
        return rng.uniform(0.0, 2.0 * np.pi, size=n)
    if name == "multi_arm":
        return (2.0 * np.pi * (k % arms) / arms) % (2.0 * np.pi)
    raise ValueError(f"unknown angular law {name!r}; choose from {ANGULAR_LAWS}")


def layout_from_laws(radial: str, angular: str, n: int, r_min: float,
                     r_max: float, jitter_m: float = 0.0, seed: int | None = None,
                     **kwargs) -> np.ndarray:
    """An ``(n, 2)`` layout from one radial law crossed with one angular law.

    ``jitter_m`` adds isotropic Gaussian noise to each position. It exists so a
    factorial experiment has replication: without within-cell variation there
    is no error term, and no way to say whether a difference between laws is
    larger than the noise.
    """
    rng = np.random.default_rng(seed)
    frac = radial_law(radial, n, **{k: v for k, v in kwargs.items()
                                    if k in ("gamma", "b", "sigma")})
    theta = angular_law(angular, n, rng=rng,
                        **{k: v for k, v in kwargs.items()
                           if k in ("p", "q", "arms")})
    r = r_min + (r_max - r_min) * frac
    xy = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    if jitter_m > 0:
        xy = xy + rng.normal(0.0, jitter_m, size=xy.shape)
    return xy
