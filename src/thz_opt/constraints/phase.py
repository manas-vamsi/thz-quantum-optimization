"""Phase / "stepping-stone" connectivity -- RESEARCH PLACEHOLDER.

This module deliberately does NOT model atmospheric phase physics.  There is no
justified coherence model available to the project yet, so implementing one
would amount to inventing numbers.  What is implemented instead is a purely
graph-theoretic abstraction of the idea in the onboarding note -- that a long
baseline is only calibratable if a chain of intermediate stations connects it
back towards the compact core:

* Nodes are the selected pads.
* An edge joins two selected pads whose separation is at most
  ``coherence_length_m`` (a free parameter, not a measured coherence length).
* The "core" is the set of selected pads within ``core_radius_m`` of a
  reference point, which defaults to the centroid of all candidate pads but can
  be supplied explicitly (the centroid is a poor reference for a layout with a
  few very remote outriggers).
* A selected pad is *supported* if the graph contains a path from it to any
  core pad.

The penalty is the number of unsupported pads.  This is a connectivity
surrogate and is labelled as such everywhere it is used; replacing it with a
physically derived coherence criterion is an explicit open task.

Note for the QUBO layer: path existence is a global property of the selected
set and is NOT a quadratic function of the selection variables.  The pairwise
surrogate offered here (:func:`pairwise_phase_penalty`) penalises *individual
long-baseline pairs that have no common intermediate pad*, which is quadratic
but strictly weaker than the path condition.  The difference is discussed in
``docs/qubo_mapping.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..interferometry.baselines import baseline_matrix


@dataclass(frozen=True)
class PhaseConfig:
    """Free parameters of the connectivity surrogate (not physical constants)."""

    coherence_length_m: float = 300.0
    core_radius_m: float = 150.0

    def as_dict(self) -> dict:
        return {
            "coherence_length_m": self.coherence_length_m,
            "core_radius_m": self.core_radius_m,
            "model": "graph connectivity surrogate",
            "physically_justified": False,
        }


def connectivity_graph(xy: np.ndarray, config: PhaseConfig | None = None) -> np.ndarray:
    """Boolean adjacency matrix: ``True`` where ``d_ij <= coherence_length_m``."""
    cfg = config or PhaseConfig()
    adj = baseline_matrix(xy) <= cfg.coherence_length_m
    np.fill_diagonal(adj, False)
    return adj


def core_nodes(
    xy: np.ndarray, config: PhaseConfig | None = None, centre: np.ndarray | None = None
) -> np.ndarray:
    """Boolean mask of pads within ``core_radius_m`` of the reference point.

    ``centre`` defaults to the centroid of ``xy``.
    """
    cfg = config or PhaseConfig()
    arr = np.asarray(xy, dtype=float)
    ref = arr.mean(axis=0) if centre is None else np.asarray(centre, dtype=float)
    return np.linalg.norm(arr - ref, axis=1) <= cfg.core_radius_m


def unsupported_stations(
    xy: np.ndarray,
    selection: np.ndarray | None = None,
    config: PhaseConfig | None = None,
    centre: np.ndarray | None = None,
) -> dict:
    """Count selected pads with no connectivity path back to the core.

    ``selection`` is a boolean/0-1 vector over the pads; ``None`` selects all.
    The traversal is a plain breadth-first search restricted to selected pads.
    """
    cfg = config or PhaseConfig()
    arr = np.asarray(xy, dtype=float)
    sel = np.ones(arr.shape[0], dtype=bool) if selection is None else np.asarray(selection, dtype=bool)
    adj = connectivity_graph(arr, cfg) & sel[:, None] & sel[None, :]
    core = core_nodes(arr, cfg, centre) & sel

    reached = core.copy()
    if core.any():
        frontier = list(np.flatnonzero(core))
        while frontier:
            k = frontier.pop()
            nxt = np.flatnonzero(adj[k] & ~reached)
            reached[nxt] = True
            frontier.extend(nxt.tolist())

    unsupported = sel & ~reached
    return {
        "n_selected": int(sel.sum()),
        "n_core": int(core.sum()),
        "n_supported": int((sel & reached).sum()),
        "n_unsupported": int(unsupported.sum()),
        "unsupported_indices": np.flatnonzero(unsupported).tolist(),
        **cfg.as_dict(),
    }


def pairwise_phase_penalty(xy: np.ndarray, config: PhaseConfig | None = None) -> np.ndarray:
    """Quadratic surrogate: ``(n, n)`` penalty matrix for long uncalibratable pairs.

    ``P[i, j] = 1`` when ``d_ij > coherence_length_m`` and no third pad ``k``
    lies within ``coherence_length_m`` of both ``i`` and ``j`` -- i.e. the pair
    has no available single stepping stone at all.  It is an upper-bound-style
    relaxation of the path condition: it ignores whether such a ``k`` is itself
    selected, and it ignores chains longer than one hop.
    """
    cfg = config or PhaseConfig()
    adj = connectivity_graph(xy, cfg)
    d = baseline_matrix(xy)
    two_hop = (adj.astype(int) @ adj.astype(int)) > 0
    p = ((d > cfg.coherence_length_m) & ~two_hop).astype(float)
    np.fill_diagonal(p, 0.0)
    return p
