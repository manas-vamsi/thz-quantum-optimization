"""Array geometries from the interferometry literature.

The golden spiral and the Fibonacci constructions in the sibling modules are
*proposals*.  This module implements the geometries that the radio-array design
literature actually converged on, so that a proposal can be measured against
them rather than against nothing.

Implemented families
--------------------
``multi_arm_spiral``
    Logarithmic spiral with ``n_arms`` arms and a **free growth rate**.  Two
    points matter here.  First, real arrays (ALMA's "zoom spiral" of Conway
    1998/2000, the VLA wye, SKA-Mid) use *several* arms; a single arm gives
    poor angular coverage because all baselines cluster around a few position
    angles.  Second, the growth-per-turn is exposed as a parameter, so the
    claim that the golden ratio specifically is the right growth rate becomes
    testable instead of assumed.

``reuleaux_layout``
    Antennas on a Reuleaux triangle -- the curve of constant width with the
    fewest sides.  Keto (1997, ApJ 475, 843) showed that curves of constant
    width give the most uniform snapshot UV sampling and the lowest sidelobes
    of any closed-curve array; the Submillimeter Array's pads are laid out as
    four nested Reuleaux triangles because of that result.

``hierarchical_layout``
    H-array of Keto (2012, JAI): repeat one small well-behaved sub-array at
    several scales.  The antenna positions are all sums
    ``s^0 p_{i0} + s^1 p_{i1} + ...`` over a base pattern ``p``.  The resulting
    baseline set is self-similar, which gives a smooth multi-scale UV
    distribution without running a numerical optimiser.

``gaussian_layout``
    Antenna positions drawn from an isotropic 2D Gaussian.  This is the
    cheapest way to hit the target that Boone (2002, A&A 386, 1160) argues
    for: a Gaussian radial UV density with uniform azimuth, which Fourier
    transforms to a Gaussian beam, i.e. a beam with no sidelobe rings at all.
    The construction is exact, not approximate -- the difference of two
    independent ``N(0, sigma^2)`` vectors is ``N(0, 2 sigma^2)``, so Gaussian
    antennas give an exactly Gaussian baseline (UV) distribution with
    ``sigma_uv = sqrt(2) * sigma_antenna``.

``phyllotaxis_layout``
    Two-parameter generalisation of the sunflower construction:
    ``theta_n = 2 pi alpha n`` and ``r_n proportional to n**p``.  Vogel's
    "Fibonacci" sunflower is the single point ``alpha = 1/PHI^2``, ``p = 1/2``.
    Exposing ``alpha`` is what makes "is the golden ratio special?" an
    experiment: any badly-approximable irrational should behave almost
    identically, and rationals should be visibly worse.

All generators return ``(n, 2)`` metre coordinates and, where meaningful,
rescale onto ``[r_min, r_max]`` so comparisons hold the physical extent fixed.
"""

from __future__ import annotations

import numpy as np

from .golden_spiral import PHI, _rescale

__all__ = [
    "multi_arm_spiral",
    "reuleaux_layout",
    "hierarchical_layout",
    "gaussian_layout",
    "phyllotaxis_layout",
    "reuleaux_base",
]


def multi_arm_spiral(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    n_arms: int = 3,
    growth: float | None = None,
    n_turns: float = 1.0,
    theta_0: float = 0.0,
    orientation: float = 0.0,
    normalize: bool = True,
) -> np.ndarray:
    """Logarithmic spiral with ``n_arms`` arms.

    ``growth`` is the factor by which the radius multiplies per full turn;
    ``None`` uses the golden ratio, which makes this a multi-arm version of
    Model 1.  Antennas are dealt round-robin to the arms, each arm being the
    same spiral rotated by ``2 pi k / n_arms``.

    With ``normalize`` the radii are rescaled onto ``[r_min, r_max]`` exactly
    as in :func:`~thz_opt.arrays.golden_spiral.golden_spiral_layout`.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n_arms < 1:
        raise ValueError("n_arms must be >= 1")
    g = PHI if growth is None else float(growth)
    if g <= 1.0:
        raise ValueError("growth must exceed 1 (radius must increase per turn)")

    per_arm = np.full(n_arms, n // n_arms, dtype=int)
    per_arm[: n % n_arms] += 1

    pts = []
    for k, cnt in enumerate(per_arm):
        if cnt == 0:
            continue
        t = np.linspace(0.0, 2.0 * np.pi * n_turns, cnt) if cnt > 1 else np.zeros(1)
        r = g ** (t / (2.0 * np.pi))
        if normalize:
            r = _rescale(r, r_min, r_max) if cnt > 1 else np.full(1, r_max)
        ang = theta_0 + t + 2.0 * np.pi * k / n_arms + orientation
        pts.append(np.column_stack((r * np.cos(ang), r * np.sin(ang))))
    return np.vstack(pts)


def reuleaux_base(width: float = 1.0, n_points: int = 3, theta_0: float = 0.0) -> np.ndarray:
    """``n_points`` equally spaced (in arc length) on a Reuleaux triangle.

    The triangle has vertices of an equilateral triangle of side ``width``; the
    boundary is the three 60-degree arcs of radius ``width`` centred on the
    opposite vertices.  The result is centred on the centroid, so the maximum
    radius is ``width / sqrt(3)`` and -- the defining property -- the maximum
    separation between any two boundary points is exactly ``width``.
    """
    if n_points < 1:
        raise ValueError("n_points must be >= 1")
    w = float(width)
    ang = theta_0 + np.array([np.pi / 2.0, np.pi / 2.0 + 2.0 * np.pi / 3.0,
                              np.pi / 2.0 + 4.0 * np.pi / 3.0])
    verts = (w / np.sqrt(3.0)) * np.column_stack((np.cos(ang), np.sin(ang)))

    # arc k is centred on vertex k and runs between the other two vertices
    s = (np.arange(n_points) + 0.5) / n_points  # fractional position around the perimeter
    arc = np.minimum((s * 3.0).astype(int), 2)
    frac = s * 3.0 - arc                        # position within the arc, in [0, 1)

    pts = np.empty((n_points, 2))
    for k in range(3):
        sel = arc == k
        if not sel.any():
            continue
        c = verts[k]
        a = verts[(k + 1) % 3] - c
        start = np.arctan2(a[1], a[0])
        phi = start + frac[sel] * (np.pi / 3.0)
        pts[sel] = c + w * np.column_stack((np.cos(phi), np.sin(phi)))
    return pts - pts.mean(axis=0, keepdims=True)


def reuleaux_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    n_rings: int = 1,
    ring_ratio: float = 0.5,
    theta_0: float = 0.0,
    jitter: float = 0.0,
    seed: int = 0,
) -> np.ndarray:
    """Antennas on one or several nested Reuleaux triangles.

    ``n_rings > 1`` reproduces the SMA-style nested arrangement: successive
    triangles are scaled by ``ring_ratio`` and rotated by 60 degrees so their
    baselines interleave rather than repeat.  ``r_max`` sets the outer
    triangle's circumradius, so the longest baseline is
    ``sqrt(3) * r_max`` (the constant width).

    ``jitter`` perturbs each antenna along the curve by up to that fraction of
    the local spacing; Keto (1997) notes that a small perturbation off the
    equally spaced positions improves UV uniformity.  ``r_min`` is unused for a
    single ring and is accepted only so that every generator shares a
    signature.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    n_rings = max(1, min(int(n_rings), n))
    rng = np.random.default_rng(seed)

    counts = np.full(n_rings, n // n_rings, dtype=int)
    counts[: n % n_rings] += 1

    pts = []
    for k, cnt in enumerate(counts):
        if cnt == 0:
            continue
        width = np.sqrt(3.0) * r_max * (ring_ratio ** k)
        offset = theta_0 + k * np.pi / 3.0
        p = reuleaux_base(width=width, n_points=cnt, theta_0=offset)
        if jitter > 0.0:
            spacing = 2.0 * np.pi * width / max(cnt, 1) / 3.0
            p = p + rng.uniform(-jitter, jitter, p.shape) * spacing
        pts.append(p)
    return np.vstack(pts)


def hierarchical_layout(
    base: np.ndarray,
    levels: int = 2,
    scale: float = 0.35,
    r_max: float = 1000.0,
    normalize: bool = True,
) -> np.ndarray:
    """H-array (Keto 2012): every sum of scaled copies of ``base``.

    With a base pattern of ``b`` points and ``levels`` levels the array has
    exactly ``b ** levels`` antennas, at positions

        p_{i0} + s * p_{i1} + s^2 * p_{i2} + ...

    ``scale`` (``s``) controls the trade-off the paper is about: small ``s``
    separates the scales, giving a spiky multi-scale UV distribution with high
    resolution; ``s`` near ``1/b`` blends them into a smoother, lower-sidelobe
    distribution.  ``normalize`` rescales the whole array so its maximum radius
    is ``r_max``.
    """
    b = np.asarray(base, dtype=float).reshape(-1, 2)
    if levels < 1:
        raise ValueError("levels must be >= 1")

    pts = b.copy()
    for k in range(1, levels):
        pts = (pts[:, None, :] + (scale ** k) * b[None, :, :]).reshape(-1, 2)
    pts = pts - pts.mean(axis=0, keepdims=True)

    if normalize:
        rr = np.hypot(pts[:, 0], pts[:, 1]).max()
        if rr > 0:
            pts = pts * (r_max / rr)
    return pts


def gaussian_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    seed: int = 0,
    sigma_fraction: float = 0.4,
    truncate: bool = True,
) -> np.ndarray:
    """Antennas drawn from an isotropic 2D Gaussian of width
    ``sigma = sigma_fraction * r_max``.

    Rationale (Boone 2002): a Gaussian radial UV density with uniform azimuth
    transforms to a Gaussian synthesised beam, which has no sidelobe rings.
    Because the difference of two independent Gaussians is Gaussian, drawing
    the *antennas* from a Gaussian gives that UV density exactly, with
    ``sigma_uv = sqrt(2) sigma_antenna`` -- no optimisation required.

    ``truncate`` resamples any point outside ``[r_min, r_max]`` so the array
    respects the same extent as the other generators.  Truncation distorts the
    Gaussian tail, which is the price of a bounded site; with the default
    ``sigma_fraction`` about 4 % of draws are affected.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    rng = np.random.default_rng(seed)
    sigma = sigma_fraction * r_max

    pts = np.empty((n, 2))
    filled = 0
    while filled < n:
        cand = rng.normal(0.0, sigma, size=(max(n - filled, 1) * 2, 2))
        if truncate:
            rad = np.hypot(cand[:, 0], cand[:, 1])
            cand = cand[(rad >= r_min) & (rad <= r_max)]
        take = min(len(cand), n - filled)
        pts[filled : filled + take] = cand[:take]
        filled += take
    return pts


def phyllotaxis_layout(
    n: int,
    r_min: float = 10.0,
    r_max: float = 1000.0,
    alpha: float | None = None,
    radial_exponent: float = 0.5,
    theta_0: float = 0.0,
) -> np.ndarray:
    """Generalised sunflower: ``theta_n = 2 pi alpha n``, ``r_n ~ n**p``.

    ``alpha = None`` uses ``1/PHI**2 = 0.381966...``, i.e. the golden angle,
    which makes this identical to the Vogel construction.  Any other value is
    allowed, and that is the point: the golden angle's virtue is that ``PHI``
    is the *most badly approximable* irrational, so the angular sequence avoids
    near-repeats for as long as possible.  If that is what matters for array
    design, then other badly approximable irrationals (``sqrt(2) - 1``,
    ``sqrt(3) - 1``) should perform almost identically and rationals should be
    clearly worse -- which is a testable prediction, not a belief.

    ``radial_exponent`` 0.5 spreads antennas uniformly per unit area; larger
    values push them outwards, smaller values concentrate them centrally.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    a = (1.0 / PHI ** 2) if alpha is None else float(alpha)
    idx = np.arange(n, dtype=float)
    r = _rescale((idx + 0.5) ** float(radial_exponent), r_min, r_max) if n > 1 else np.array([r_max])
    theta = theta_0 + 2.0 * np.pi * a * idx
    return np.column_stack((r * np.cos(theta), r * np.sin(theta)))
