"""Baseline-length distribution diagnostics.

The project onboarding note asserts that useful UV sampling requires a
logarithmic / power-law spread of baseline lengths bridging the compact core
and the extended stations.  This module provides the tools to *test* that claim
on a given layout rather than assume it: an empirical CDF, log-spaced
histograms, and a least-squares fit of a power law to the baseline-length
distribution together with its goodness of fit.
"""

from __future__ import annotations

import numpy as np


def length_histogram(lengths: np.ndarray, n_bins: int = 30, log_bins: bool = False):
    """``(bin_centres, counts, edges)`` of the baseline-length distribution."""
    d = np.asarray(lengths, dtype=float)
    d = d[d > 0]
    if d.size == 0:
        raise ValueError("no positive baseline lengths")
    if log_bins:
        edges = np.logspace(np.log10(d.min()), np.log10(d.max()), n_bins + 1)
        centres = np.sqrt(edges[1:] * edges[:-1])
    else:
        edges = np.linspace(d.min(), d.max(), n_bins + 1)
        centres = 0.5 * (edges[1:] + edges[:-1])
    counts, _ = np.histogram(d, bins=edges)
    return centres, counts, edges


def empirical_cdf(lengths: np.ndarray):
    """``(sorted_lengths, cdf)`` with ``cdf = (k+1)/n``."""
    d = np.sort(np.asarray(lengths, dtype=float))
    return d, np.arange(1, d.size + 1) / d.size


def powerlaw_fit(lengths: np.ndarray, n_bins: int = 20) -> dict:
    """Least-squares fit of ``log10 N(d) = alpha * log10 d + c`` on log bins.

    Returns the exponent ``alpha``, the intercept, the coefficient of
    determination ``r_squared`` of the fit in log-log space, and the number of
    non-empty bins used.  A near-flat count per *logarithmic* bin
    (``alpha ~ 0``) is what "logarithmic distribution of baseline lengths"
    means operationally; ``r_squared`` says whether a power law describes the
    distribution at all.  Empty bins are dropped, which biases the fit for
    sparse arrays -- ``n_bins_used`` is reported so that can be judged.
    """
    centres, counts, _ = length_histogram(lengths, n_bins=n_bins, log_bins=True)
    keep = counts > 0
    if keep.sum() < 3:
        return {"alpha": float("nan"), "intercept": float("nan"),
                "r_squared": float("nan"), "n_bins_used": int(keep.sum())}
    x = np.log10(centres[keep])
    y = np.log10(counts[keep])
    alpha, c = np.polyfit(x, y, 1)
    resid = y - (alpha * x + c)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid**2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    return {"alpha": float(alpha), "intercept": float(c),
            "r_squared": float(r2), "n_bins_used": int(keep.sum())}


def log_spacing_uniformity(lengths: np.ndarray, n_bins: int = 20) -> float:
    """Normalised entropy of the baseline lengths over log-spaced bins.

    1.0 means the baselines are spread evenly across decades of length -- the
    operational reading of "logarithmic coverage".  0 means they all sit in one
    log bin.  Depends on ``n_bins`` and on the length range, so compare only
    across configurations sharing both.
    """
    _, counts, _ = length_histogram(lengths, n_bins=n_bins, log_bins=True)
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum() / np.log(n_bins))


def baseline_summary(lengths: np.ndarray, n_bins: int = 20) -> dict:
    """Bundle of scalar baseline-length statistics, all in metres."""
    d = np.asarray(lengths, dtype=float)
    return {
        "n_baselines": int(d.size),
        "d_min_m": float(d.min()),
        "d_max_m": float(d.max()),
        "d_mean_m": float(d.mean()),
        "d_median_m": float(np.median(d)),
        "dynamic_range": float(d.max() / d.min()) if d.min() > 0 else float("inf"),
        "log_spacing_uniformity": log_spacing_uniformity(d, n_bins),
        **{f"powerlaw_{k}": v for k, v in powerlaw_fit(d, n_bins).items()},
    }
