"""Redundancy of UV sampling.

A *redundant* sample is one that lands in a cell already occupied by another
sample.  Redundancy is not intrinsically bad -- it buys signal-to-noise and
calibration closure -- so these functions only quantify it; deciding whether it
should be rewarded or penalised is a modelling choice made in the QUBO layer.
"""

from __future__ import annotations

import numpy as np


def redundancy_metrics(occ: np.ndarray) -> dict:
    """Redundancy statistics of an occupancy matrix.

    * ``n_samples``            sum of ``n_c``
    * ``unique_cells``         ``|O|``
    * ``redundant_samples``    ``sum_c (n_c - 1)`` over occupied cells
    * ``redundancy_fraction``  ``redundant_samples / n_samples``
    * ``mean_multiplicity``    ``n_samples / |O|``
    * ``max_multiplicity``     ``max_c n_c``
    """
    o = np.asarray(occ)
    total = int(o.sum())
    uniq = int((o > 0).sum())
    redundant = total - uniq
    return {
        "n_samples": total,
        "unique_cells": uniq,
        "redundant_samples": redundant,
        "redundancy_fraction": float(redundant / total) if total else 0.0,
        "mean_multiplicity": float(total / uniq) if uniq else 0.0,
        "max_multiplicity": int(o.max()) if o.size else 0,
    }


def multiplicity_histogram(occ: np.ndarray):
    """``(multiplicity, n_cells_with_that_multiplicity)`` for occupied cells."""
    o = np.asarray(occ)
    vals = o[o > 0]
    if vals.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    m = np.arange(1, int(vals.max()) + 1)
    return m, np.array([(vals == k).sum() for k in m], dtype=int)
