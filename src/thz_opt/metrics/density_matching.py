"""Objectives the array-design literature actually optimises.

The metrics in :mod:`thz_opt.metrics.uv_coverage` describe a UV distribution.
The two here *score* it against a design goal, which is what an optimiser
needs.  Both are standard; neither is invented for this repository.

1. **Target radial UV density** (Boone 2002, A&A 386, 1160).  The synthesised
   beam is the Fourier transform of the UV sampling density, so the density's
   *shape* sets the beam's shape.  A Gaussian radial density with uniform
   azimuth transforms to a Gaussian beam -- no sidelobe rings at all -- while a
   uniformly filled disc transforms to an Airy-like ``J1(r)/r`` beam whose
   first sidelobe is ~13 dB down and whose rings decay slowly.  "Fill as many
   cells as possible" is therefore not the same goal as "image cleanly", and
   the two can point in opposite directions.

2. **UV-point repulsion energy** (Cornwell 1988, IEEE AP-36, 1165; used for
   pad selection by Karastergiou, Neri & Gurwell 2006, ApJS 164, 552):

       F = R^2 * sum_{i<j} 1 / |V_i - V_j|^2

   over UV samples ``V``, with ``R`` the maximum baseline.  Treating the UV
   points as mutually repelling charges and minimising ``F`` spreads them as
   evenly as the geometry allows.  This is the closest thing the field has to a
   standard scalar objective for configuration design, and -- important for
   this project -- it is a sum over *pairs of baselines*, which is exactly the
   structure a QUBO over baseline-activation variables can hold (see
   :mod:`thz_opt.qubo.baseline_qubo`).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "radial_density_profile",
    "gaussian_target_profile",
    "density_match_chi2",
    "cornwell_energy",
    "fit_gaussian_sigma",
]


def radial_density_profile(uv: np.ndarray, n_bins: int = 30, r_max: float | None = None):
    """Normalised radial UV *density* (samples per unit UV area).

    Returns ``(bin_centres, density, edges)`` where ``density`` integrates to 1
    over the plane, so arrays with different sample counts are comparable.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    r = np.hypot(arr[:, 0], arr[:, 1])
    hi = float(r.max()) if r_max is None else float(r_max)
    edges = np.linspace(0.0, max(hi, 1e-12), n_bins + 1)
    counts, _ = np.histogram(r, bins=edges)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    density = counts / area
    total = float((density * area).sum())
    if total > 0:
        density = density / total
    return 0.5 * (edges[1:] + edges[:-1]), density, edges


def gaussian_target_profile(centres: np.ndarray, sigma: float) -> np.ndarray:
    """Isotropic Gaussian UV density ``exp(-r^2 / 2 sigma^2) / (2 pi sigma^2)``.

    Normalised the same way as :func:`radial_density_profile`, so the two can
    be differenced directly.
    """
    c = np.asarray(centres, dtype=float)
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    return np.exp(-(c ** 2) / (2.0 * sigma ** 2)) / (2.0 * np.pi * sigma ** 2)


def fit_gaussian_sigma(uv: np.ndarray) -> float:
    """Maximum-likelihood ``sigma`` of an isotropic 2D Gaussian fitted to ``uv``.

    For an isotropic 2D Gaussian, ``E[r^2] = 2 sigma^2``.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    return float(np.sqrt(np.mean(arr[:, 0] ** 2 + arr[:, 1] ** 2) / 2.0))


def density_match_chi2(
    uv: np.ndarray,
    sigma: float | None = None,
    n_bins: int = 30,
    r_max: float | None = None,
    weight_by_area: bool = True,
) -> dict:
    """Mismatch between the measured radial UV density and a Gaussian target.

    ``sigma`` defaults to the maximum-likelihood fit to the samples, which
    measures *shape* mismatch only; pass an explicit ``sigma`` to score against
    a design target.  ``weight_by_area`` weights each annulus by its area so
    that the statistic is an approximation to
    ``integral (rho - rho_target)^2 dA`` rather than a bin-count sum.

    Lower is better.  The value is dimensionless and normalised by the target's
    own squared norm, so it is comparable between arrays of different size.
    """
    centres, density, edges = radial_density_profile(uv, n_bins, r_max)
    s = fit_gaussian_sigma(uv) if sigma is None else float(sigma)
    target = gaussian_target_profile(centres, s)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2) if weight_by_area else np.ones_like(centres)

    num = float(np.sum(area * (density - target) ** 2))
    den = float(np.sum(area * target ** 2))
    return {
        "sigma_lambda": s,
        "chi2": num,
        "normalised_chi2": num / den if den > 0 else float("nan"),
        "n_bins": int(n_bins),
    }


def cornwell_energy(
    uv: np.ndarray,
    max_samples: int = 4000,
    seed: int = 0,
    power: float = 2.0,
    floor_fraction: float = 1e-3,
) -> dict:
    """Cornwell's UV-point repulsion energy, normalised for comparability.

    ``F = R^power * mean_{i<j} |V_i - V_j|^-power``, where ``R`` is the maximum
    UV radius.  The mean (rather than the sum) and the ``R^power`` factor make
    the value dimensionless and independent of the number of samples, so arrays
    with different sample counts can be compared -- Cornwell's original sum
    cannot be.

    Coincident UV samples (exactly redundant baselines) would make the energy
    infinite, so separations are floored at ``floor_fraction * R``.  That floor
    is a stated choice, not physics: it caps the penalty a perfectly redundant
    pair can contribute.  The number of pairs hitting the floor is reported so
    its influence can be judged.

    Cost is ``O(n^2)`` in the number of UV samples.  Above ``max_samples`` the
    samples are subsampled with a fixed seed and ``subsampled`` is set in the
    result; the value is then an unbiased estimate of the mean but not exact.
    """
    arr = np.asarray(uv, dtype=float).reshape(-1, 2)
    n_total = arr.shape[0]
    if n_total < 2:
        return {"cornwell_F": float("nan"), "n_samples_used": n_total, "subsampled": False}

    subsampled = n_total > max_samples
    if subsampled:
        rng = np.random.default_rng(seed)
        arr = arr[rng.choice(n_total, max_samples, replace=False)]

    n = arr.shape[0]
    R = float(np.hypot(arr[:, 0], arr[:, 1]).max())
    floor = max(floor_fraction * R, 1e-12)
    total = 0.0
    n_pairs = 0
    n_floored = 0
    chunk = max(1, int(2e7 // max(n, 1)))
    for start in range(0, n, chunk):
        block = arr[start : start + chunk]
        d = np.hypot(
            block[:, None, 0] - arr[None, :, 0],
            block[:, None, 1] - arr[None, :, 1],
        )
        rows = np.arange(start, start + block.shape[0])
        mask = rows[:, None] < np.arange(n)[None, :]
        dd = d[mask]
        n_floored += int((dd < floor).sum())
        total += float(np.sum(1.0 / np.maximum(dd, floor) ** power))
        n_pairs += int(mask.sum())

    return {
        "cornwell_F": (R ** power) * total / max(n_pairs, 1),
        "n_samples_used": n,
        "subsampled": bool(subsampled),
        "uv_radius_max": R,
        "pairs_at_floor": n_floored,
        "floor_fraction": floor_fraction,
    }
