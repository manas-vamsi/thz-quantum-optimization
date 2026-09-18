"""Models 3 and 4 -- reference layouts used as controls.

These exist so that any golden-vs-Fibonacci difference can be read against a
neutral baseline.  Both are normalised to the same ``[r_min, r_max]`` extent as
the spiral models.
"""

from __future__ import annotations

import numpy as np


def random_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    seed: int = 0,
    area_uniform: bool = True,
) -> np.ndarray:
    """Pseudo-random layout inside the annulus ``r_min <= r <= r_max``.

    ``area_uniform`` draws radii with density proportional to ``r`` so points
    are uniform per unit *area*; otherwise radii are uniform in ``r`` (which
    concentrates points towards the centre in area terms).
    """
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0.0, 2.0 * np.pi, n)
    if area_uniform:
        r = np.sqrt(rng.uniform(r_min**2, r_max**2, n))
    else:
        r = rng.uniform(r_min, r_max, n)
    return np.column_stack((r * np.cos(theta), r * np.sin(theta)))


def uniform_ring_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    n_rings: int = 3,
    theta_0: float = 0.0,
) -> np.ndarray:
    """Quasi-uniform concentric-ring layout (Model 4).

    Antennas are split as evenly as possible over ``n_rings`` rings whose radii
    are linearly spaced on ``[r_min, r_max]``.  Consecutive rings are offset in
    angle by half a ring spacing to avoid radial spokes.
    """
    if n_rings < 1:
        raise ValueError("n_rings must be >= 1")
    n_rings = min(n_rings, n)
    radii = np.linspace(r_min, r_max, n_rings) if n_rings > 1 else np.array([r_max])
    counts = np.full(n_rings, n // n_rings, dtype=int)
    counts[: n % n_rings] += 1

    pts = []
    for k, (rad, cnt) in enumerate(zip(radii, counts)):
        if cnt == 0:
            continue
        ang = theta_0 + np.arange(cnt) * (2.0 * np.pi / cnt) + k * np.pi / max(cnt, 1)
        pts.append(np.column_stack((rad * np.cos(ang), rad * np.sin(ang))))
    return np.vstack(pts)
