"""Science cases: letting the observation choose the objective.

Every "best layout" result so far carried the caveat that the objectives
disagree and nothing fixed which one to use. That caveat becomes a *result* as
soon as two concrete observations are specified and shown to prefer different
arrays.

Two deliberately opposed cases:

``compact_source`` -- e.g. an EHT-style black-hole shadow
    A small, bright, high-contrast object. Wants the **longest** baselines it
    can coherently use (resolution), a clean beam (low sidelobes, because the
    dynamic range is what limits the measurement), and it is the case where
    atmospheric decorrelation bites hardest. Short baselines measure almost
    nothing about a source smaller than their fringe spacing.

``extended_emission`` -- e.g. a molecular cloud or protoplanetary disc
    A large, low-contrast structure. Wants **short**-baseline completeness --
    missing short spacings produce the classic "negative bowl" and resolve out
    large-scale flux. Long baselines mostly add noise.

Both are expressed the same way: a **radial weight on UV cells**, `w(r)`, plus
a choice of what to do with the weighted occupancy. That keeps everything
compatible with the incremental state and with the quadratic QUBO form -- a
per-cell weight multiplies `a_ck`, and `n_c = sum_k a_ck y_k` stays linear.

The weight profiles below are *design choices*, not derived optima. They encode
"this observation cares more about this part of the UV plane", which is the
honest content of a science case at this stage. A real case would derive the
weights from the expected source structure (the Fourier transform of a model
brightness distribution is exactly the right weight).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "cell_radii",
    "radial_weights",
    "compact_source_weights",
    "extended_emission_weights",
    "science_case_objective",
    "SCIENCE_CASES",
]


def cell_radii(grid) -> np.ndarray:
    """UV radius (wavelengths) of every flat cell id, in grid order."""
    k = np.arange(grid.n_cells) - grid.n_cells / 2.0
    c = k * grid.cell_size
    ru, rv = np.meshgrid(c, c, indexing="ij")
    return np.hypot(ru, rv).ravel()


def radial_weights(grid, profile: str = "flat", r_break: float | None = None,
                   power: float = 1.0) -> np.ndarray:
    """Per-cell weight as a function of UV radius.

    ``profile``:

    ``"flat"``        every cell equal -- the implicit assumption so far
    ``"outer"``       ``(r / r_break) ** power``, a smooth tilt to long baselines
    ``"inner"``       ``exp(-r / r_break)``, a smooth tilt to short baselines
    ``"annulus"``     a soft band centred on ``r_break``
    ``"outer_cut"``   1 beyond ``r_break``, 0 inside it
    ``"inner_cut"``   1 inside ``r_break``, 0 beyond it

    ``r_break`` defaults to a quarter of the grid extent. Weights are
    normalised to a mean of 1 so that different profiles give comparable
    magnitudes.

    **Measured behaviour worth knowing before choosing one.** The smooth
    profiles barely move the optimum: at M=60, N=15 the layouts they select
    share 9-13 of 15 pads with the unweighted optimum, and the median baseline
    changes by under 13 %. The *cut* profiles do move it -- ``"inner_cut"``
    drops the median baseline by 27 %. The reason is that the dominant lever is
    still "touch as many cells as possible", and a smooth tilt does not
    overcome it; only removing a region from consideration does. So a science
    case expressed as a mild preference will not change the design, and one
    expressed as "these spatial frequencies are the measurement" will.
    """
    r = cell_radii(grid)
    rb = grid.uv_max / 4.0 if r_break is None else float(r_break)

    if profile == "flat":
        w = np.ones_like(r)
    elif profile == "outer":
        w = (r / rb) ** power
    elif profile == "inner":
        w = np.exp(-r / rb)
    elif profile == "annulus":
        w = np.exp(-0.5 * ((r - rb) / (0.35 * rb)) ** 2)
    elif profile == "outer_cut":
        w = (r >= rb).astype(float)
    elif profile == "inner_cut":
        w = (r <= rb).astype(float)
    else:
        raise ValueError(f"unknown profile {profile!r}")

    m = float(w.mean())
    return w / m if m > 0 else w


def compact_source_weights(grid, r_break: float | None = None,
                           hard: bool = True) -> np.ndarray:
    """Only long baselines count: a source smaller than the fringe spacing of a
    short baseline is unresolved by it, so that baseline measures total flux and
    nothing about structure.

    ``hard=True`` (default) restricts the measurement to ``r >= r_break``, which
    is what actually changes the optimum; ``hard=False`` gives the smooth tilt,
    kept so the comparison in the docstring of :func:`radial_weights` can be
    reproduced.
    """
    rb = grid.uv_max / 2.0 if r_break is None else r_break
    return (radial_weights(grid, "outer_cut", rb) if hard
            else radial_weights(grid, "outer", rb, power=1.0))


def extended_emission_weights(grid, r_break: float | None = None,
                              hard: bool = True) -> np.ndarray:
    """Only short baselines count: large-scale flux lives at small UV radius, and
    long baselines resolve it out entirely.
    """
    rb = grid.uv_max / 4.0 if r_break is None else r_break
    return (radial_weights(grid, "inner_cut", rb) if hard
            else radial_weights(grid, "inner", rb))


def science_case_objective(case: str, w_sidelobe: float = 0.0,
                           sidelobe_scale: float = 1.0,
                           cell_scale: float = 1.0):
    """Objective callable for a named science case.

    Both cases maximise *weighted* covered cells; the compact case additionally
    penalises weighted sidelobe energy, because its limit is dynamic range
    rather than sensitivity. ``w_sidelobe`` is that trade-off, and the scales
    put the two terms on comparable footing (pass values measured on a
    reference configuration).

    The state must have been constructed with the matching ``cell_weights``.
    """
    if case not in SCIENCE_CASES:
        raise ValueError(f"unknown science case {case!r}; choose from {tuple(SCIENCE_CASES)}")

    def score(state) -> float:
        val = -float(state.weighted_cells) / cell_scale
        if w_sidelobe:
            val += w_sidelobe * float(state.weighted_sum_sq) / sidelobe_scale
        return val

    score.name = f"science:{case}(w_sidelobe={w_sidelobe})"
    return score


SCIENCE_CASES = {
    "compact_source": {
        "weights": compact_source_weights,
        "description": "small bright high-contrast source; resolution and "
                       "dynamic range dominate; decorrelation is the limit",
        "default_w_sidelobe": 1.0,
    },
    "extended_emission": {
        "weights": extended_emission_weights,
        "description": "large low-contrast structure; short-spacing "
                       "completeness dominates; long baselines add little",
        "default_w_sidelobe": 0.0,
    },
}
