"""Minimum-separation (shadowing) constraint.

Source-derived statement
------------------------
The project onboarding note states that antenna pads must satisfy a minimum
separation of roughly ``d_min >= 1.5 D``, where ``D`` is the dish diameter, to
avoid mutual shadowing.  The factor is exposed as a parameter with that default
rather than hard-coded, and no tighter claim is made here: 1.5 is what the note
says, and the true value depends on dish design and the elevation range the
array is meant to work at.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..interferometry.baselines import baseline_matrix, pair_indices


@dataclass(frozen=True)
class SeparationConfig:
    """Dish geometry that sets the shadowing limit."""

    dish_diameter_m: float = 6.0
    minimum_separation_factor: float = 1.5

    @property
    def d_min(self) -> float:
        return self.dish_diameter_m * self.minimum_separation_factor

    def as_dict(self) -> dict:
        return {
            "dish_diameter_m": self.dish_diameter_m,
            "minimum_separation_factor": self.minimum_separation_factor,
            "d_min_m": self.d_min,
        }


def check_minimum_separation(xy: np.ndarray, config: SeparationConfig | None = None) -> dict:
    """Report every pair closer than ``d_min``.

    Returns the number of violations, the list of ``(i, j, distance)`` triples,
    the smallest separation present, and a boolean ``valid`` flag.  Nothing is
    modified or removed -- repair is the caller's decision.
    """
    cfg = config or SeparationConfig()
    d = baseline_matrix(xy)
    i, j = pair_indices(d.shape[0])
    dij = d[i, j]
    bad = dij < cfg.d_min
    return {
        "valid": bool(not bad.any()),
        "n_violations": int(bad.sum()),
        "violations": [(int(a), int(b), float(c)) for a, b, c in zip(i[bad], j[bad], dij[bad])],
        "min_separation_m": float(dij.min()) if dij.size else float("inf"),
        **cfg.as_dict(),
    }


def shadow_pairs(xy: np.ndarray, config: SeparationConfig | None = None) -> np.ndarray:
    """Boolean ``(n, n)`` matrix, ``True`` where ``i != j`` and ``d_ij < d_min``.

    This is the matrix consumed by the shadow term of the QUBO.
    """
    cfg = config or SeparationConfig()
    d = baseline_matrix(xy)
    mask = d < cfg.d_min
    np.fill_diagonal(mask, False)
    return mask
