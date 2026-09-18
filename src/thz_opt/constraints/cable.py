"""Cable and trenching cost -- the linear part, and the honest caveat.

A real array is not chosen on UV coverage alone: every occupied pad needs power,
fibre and a trench to reach it, and that cost is a large fraction of the
infrastructure budget. Cohanim et al. (2004) treated imaging performance against
cable length as an explicit two-objective problem for exactly this reason.

What fits in a QUBO
-------------------
A **star** (hub-and-spoke) layout, where each pad is trenched independently back
to a central correlator, costs

    C_star(x) = sum_i c_i x_i,   c_i = |r_i - r_hub|

which is linear -- one diagonal entry per pad, no approximation.

What does not
-------------
The true minimum trench length for a set of pads is a **minimum spanning tree**
over the selected pads. That is a global connectivity property of the selected
set, not a sum over pads or pairs, and it is not expressible as a quadratic
form -- the same failure mode as the "stepping-stone" phase condition. Attempts
to write it quadratically end up with either a relaxation that ignores
connectivity or an auxiliary-variable count that dwarfs the problem.

So: optimise with the star cost, and evaluate the real MST afterwards with
:func:`mst_length`. The star cost is an upper bound on the MST (a star *is* a
spanning tree), so it never flatters a layout -- it over-charges it, and by a
measurable amount that :func:`cable_summary` reports.

``ponytail: star cost in the objective, MST measured post hoc. A Steiner-tree
model would be more realistic still, and is not quadratic either.``
"""

from __future__ import annotations

import numpy as np

from ..arrays.validation import validate_layout

__all__ = ["star_cable_cost", "mst_length", "cable_summary"]


def star_cable_cost(pads: np.ndarray, hub: np.ndarray | None = None) -> np.ndarray:
    """Per-pad trench length to a hub, in metres. Linear QUBO coefficients.

    ``hub`` defaults to the centroid of the candidate pads, which is the
    cheapest hub for a star and keeps the cost independent of which subset is
    eventually selected -- required, since the coefficients are computed offline.
    """
    arr = validate_layout(pads)
    centre = arr.mean(axis=0) if hub is None else np.asarray(hub, dtype=float)
    return np.linalg.norm(arr - centre, axis=1)


def mst_length(pads: np.ndarray, selection=None, hub: np.ndarray | None = None) -> float:
    """Minimum spanning tree length over the selected pads plus the hub, in metres.

    This is the quantity a costing exercise actually wants, and the reason it is
    here rather than in the objective: it cannot be written quadratically.
    """
    from scipy.sparse.csgraph import minimum_spanning_tree

    arr = validate_layout(pads)
    sel = arr if selection is None else arr[np.asarray(selection, dtype=int)]
    centre = arr.mean(axis=0) if hub is None else np.asarray(hub, dtype=float)
    nodes = np.vstack((centre[None, :], sel))
    d = np.linalg.norm(nodes[:, None, :] - nodes[None, :, :], axis=2)
    return float(minimum_spanning_tree(d).sum())


def cable_summary(pads: np.ndarray, selection=None, hub: np.ndarray | None = None) -> dict:
    """Star cost, true MST cost, and how much the star over-charges.

    A star is one particular spanning tree, so ``star >= mst`` always; the ratio
    says how much realism the linear term gives up.
    """
    arr = validate_layout(pads)
    sel = np.arange(arr.shape[0]) if selection is None else np.asarray(selection, dtype=int)
    c = star_cable_cost(arr, hub)
    star = float(c[sel].sum())
    mst = mst_length(arr, sel, hub)
    return {
        "n_selected": int(sel.size),
        "star_length_m": star,
        "mst_length_m": mst,
        "star_over_mst": star / mst if mst > 0 else float("nan"),
        "note": "star cost is linear and QUBO-compatible; MST is not",
    }
