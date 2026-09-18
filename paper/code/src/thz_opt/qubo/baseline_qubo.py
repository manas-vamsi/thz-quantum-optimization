"""QUBO in baseline-activation space -- where the good objectives are exact.

Why this module exists
----------------------
In pad-selection space ``x`` the only objectives that are exactly quadratic are
those that decompose over *single* baselines, ``sum_{i<j} w_ij x_i x_j``.  Every
objective the array-design literature actually uses -- unique UV-cell coverage,
Cornwell's UV-point repulsion, matching a target UV density -- couples
*different baselines to each other*, and is therefore quartic in ``x``
(see ``docs/qubo_mapping.md``).

But all three are exactly **quadratic in the baseline-activation variables**

    y_k = x_i x_j,   k indexing the M(M-1)/2 candidate pairs

because each is a sum over pairs of baselines.  So the honest formulation is:
work in ``y``, where the physics is exact, and pay for the substitution
``y_k = x_i x_j`` with the standard Rosenberg penalty, which is itself exactly
quadratic:

    P_k(x_i, x_j, y_k) = x_i x_j - 2 x_i y_k - 2 x_j y_k + 3 y_k

``P_k = 0`` when ``y_k = x_i x_j`` and ``P_k >= 1`` otherwise, for every binary
assignment.  The whole model is then a genuine QUBO over ``M + M(M-1)/2``
variables, in exchange for a larger variable count than the pairwise surrogate
of :mod:`thz_opt.qubo.coefficients` and a much smaller one than the
cell-variable construction in :mod:`thz_opt.qubo.objective`.

**Be precise about what "exact" means here**, because the objectives differ:

* PSF sidelobe energy, Cornwell repulsion and target-density matching are
  represented **exactly** -- the quadratic form equals the quantity itself.
* Unique-cell coverage is **not**.  What is represented exactly is the
  second-order Bonferroni *lower bound* on it, which equals the coverage only
  when no UV cell is touched by three or more active baselines.  Exact union
  coverage needs either cell variables (see
  :func:`thz_opt.qubo.objective.build_qubo_with_cell_variables`) or a
  higher-order inclusion-exclusion construction.

So: exact representation of a bound is not the same as exact representation of
coverage, and the distinction is measured rather than glossed -- experiment 07
reports the bound gap as a function of UV resolution (25 % at the coarsest grid
tested, 0.01 % at the finest).

Three objectives are provided, each a pure ``(linear, quadratic)`` pair over
``y``:

``bonferroni_coverage_terms``
    Second-order inclusion-exclusion (Bonferroni) **lower bound** on unique
    UV-cell coverage -- the bound, not the coverage:

        |Cov| >= sum_k |C_k| y_k - sum_{k<l} |C_k ^ C_l| y_k y_l

    This is a *provable lower bound*, exact whenever no cell is covered three
    or more times, and it is the principled version of the "discount contested
    cells" heuristic.  Its coupling matrix is naturally sparse: two baselines
    interact only if their UV tracks actually overlap.

``cornwell_terms``
    Cornwell (1988) / Karastergiou et al. (2006) UV-point repulsion energy,
    ``F = R^2 sum_{i<j} |V_i - V_j|^-2``, split into the part internal to each
    baseline's own track (linear in ``y``) and the part between tracks
    (quadratic).  Dense, but the couplings fall off as ``1/d^2`` and can be
    truncated.

``density_match_terms``
    Boone (2002) style target-density matching: penalise the squared
    difference between the achieved radial UV histogram and a target
    (e.g. Gaussian).  Expanding the square gives exactly quadratic ``y`` terms.

``sidelobe_energy_terms`` / ``coherence_weighted_sidelobe_terms``
    ``sum_c n_c^2``, which by Parseval is the dirty beam's total energy, with
    each baseline optionally weighted by the fraction of its coherence that
    survives the atmosphere.  Exact.

All return coefficients for a **minimisation** problem.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..interferometry.baselines import pair_indices

__all__ = [
    "BaselineTerms",
    "build_baseline_qubo",
    "rosenberg_penalty_floor",
    "selection_penalty_floor",
    "complete_baseline_auxiliary",
    "baseline_cell_overlaps",
    "bonferroni_coverage_terms",
    "bonferroni_bound",
    "cornwell_terms",
    "track_cell_multiplicities",
    "sidelobe_energy_terms",
    "coherence_weighted_sidelobe_terms",
    "density_match_terms",
    "objective_value",
]


# --------------------------------------------------------------------------
# container
# --------------------------------------------------------------------------

@dataclass
class BaselineTerms:
    """An objective expressed over baseline-activation variables.

    ``linear[k]`` multiplies ``y_k``; ``quadratic[(k, l)]`` with ``k < l``
    multiplies ``y_k y_l``.  ``offset`` is an additive constant.  The
    dictionary form is used because the coupling matrices of interest are
    sparse (overlap) or truncatable (repulsion).
    """

    n_pads: int
    linear: np.ndarray
    quadratic: dict = field(default_factory=dict)
    offset: float = 0.0
    name: str = "objective"

    @property
    def n_pairs(self) -> int:
        return int(self.linear.shape[0])

    def sparsity(self) -> dict:
        p = self.n_pairs
        dense = p * (p - 1) // 2
        nz = len(self.quadratic)
        return {
            "n_pairs": p,
            "n_couplings": nz,
            "n_couplings_dense": dense,
            "fill_fraction": nz / dense if dense else 0.0,
        }


def objective_value(y: np.ndarray, terms: BaselineTerms) -> float:
    """Evaluate the objective directly from its definition, for validation."""
    v = np.asarray(y, dtype=float)
    val = float(terms.linear @ v) + terms.offset
    for (k, l), c in terms.quadratic.items():
        val += c * v[k] * v[l]
    return val


def selection_to_y(x: np.ndarray, n_pads: int) -> np.ndarray:
    """``y_k = x_i x_j`` for every candidate pair, in ``i < j`` order."""
    v = np.asarray(x, dtype=int)
    i, j = pair_indices(n_pads)
    return (v[i] * v[j]).astype(int)


# --------------------------------------------------------------------------
# objective 1 -- Bonferroni coverage bound
# --------------------------------------------------------------------------

def baseline_cell_overlaps(cell_sets: list) -> dict:
    """``|C_k ^ C_l|`` for every pair of baselines that share at least one cell.

    Built by inverting the baseline -> cells map, so the cost is
    ``sum_c m_c^2`` rather than ``P^2`` set intersections.
    """
    claimants: dict = {}
    for k, cells in enumerate(cell_sets):
        for c in np.asarray(cells).tolist():
            claimants.setdefault(int(c), []).append(k)

    overlaps: dict = {}
    for ks in claimants.values():
        if len(ks) < 2:
            continue
        for a in range(len(ks)):
            for b in range(a + 1, len(ks)):
                key = (ks[a], ks[b]) if ks[a] < ks[b] else (ks[b], ks[a])
                overlaps[key] = overlaps.get(key, 0) + 1
    return overlaps


def bonferroni_coverage_terms(cell_sets: list, n_pads: int) -> BaselineTerms:
    """Minimising this maximises the Bonferroni lower bound on coverage.

    Objective ``-(sum_k |C_k| y_k - sum_{k<l} |C_k ^ C_l| y_k y_l)``, so the
    returned value is minus a lower bound on ``|Cov(x)|``.
    """
    lin = -np.array([float(len(c)) for c in cell_sets])
    quad = {kl: float(v) for kl, v in baseline_cell_overlaps(cell_sets).items()}
    return BaselineTerms(n_pads=n_pads, linear=lin, quadratic=quad, name="bonferroni_coverage")


def bonferroni_bound(y: np.ndarray, cell_sets: list, overlaps: dict | None = None) -> float:
    """``sum_k |C_k| y_k - sum_{k<l} |C_k ^ C_l| y_k y_l``.

    By the second Bonferroni inequality this never exceeds the true number of
    distinct cells covered, and equals it when no cell is covered three or more
    times.  :func:`thz_opt.qubo.coefficients.exact_cells_covered` gives the
    truth to compare against.
    """
    v = np.asarray(y, dtype=float)
    ov = baseline_cell_overlaps(cell_sets) if overlaps is None else overlaps
    total = float(sum(len(c) * v[k] for k, c in enumerate(cell_sets)))
    for (k, l), c in ov.items():
        total -= c * v[k] * v[l]
    return total


# --------------------------------------------------------------------------
# objective 2 -- Cornwell / Karastergiou UV-point repulsion
# --------------------------------------------------------------------------

def cornwell_terms(
    uv_per_baseline: list,
    n_pads: int,
    power: float = 2.0,
    floor_fraction: float = 1e-3,
    coupling_cutoff: float | None = None,
) -> BaselineTerms:
    """Repulsion energy between UV samples, split into ``y`` linear/quadratic.

    ``uv_per_baseline[k]`` is the ``(n_k, 2)`` array of UV samples contributed
    by baseline ``k`` (its snapshot points, or its whole Earth-rotation track,
    conjugates included).  The energy of a configuration is the repulsion among
    all samples of all *active* baselines, which splits exactly into

        sum_k E_kk y_k  +  sum_{k<l} E_kl y_k y_l

    with ``E_kk`` the internal energy of one track and ``E_kl`` the energy
    between two tracks.  ``coupling_cutoff`` (in wavelengths) drops pairs of
    baselines whose closest samples are farther apart than the cutoff, which
    sparsifies the model at a stated, measurable cost.
    """
    tracks = [np.asarray(t, dtype=float).reshape(-1, 2) for t in uv_per_baseline]
    all_uv = np.vstack(tracks)
    R = float(np.hypot(all_uv[:, 0], all_uv[:, 1]).max())
    floor = max(floor_fraction * R, 1e-12)
    scale = R ** power

    def energy(a: np.ndarray, b: np.ndarray, same: bool) -> float:
        d = np.hypot(a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1])
        if same:
            iu = np.triu_indices(d.shape[0], k=1)
            d = d[iu]
        return float(np.sum(1.0 / np.maximum(d, floor) ** power))

    p = len(tracks)
    lin = np.array([scale * energy(t, t, True) for t in tracks])
    quad: dict = {}
    for k in range(p):
        for l in range(k + 1, p):
            if coupling_cutoff is not None:
                d = np.hypot(
                    tracks[k][:, None, 0] - tracks[l][None, :, 0],
                    tracks[k][:, None, 1] - tracks[l][None, :, 1],
                ).min()
                if d > coupling_cutoff:
                    continue
            e = scale * energy(tracks[k], tracks[l], False)
            if e != 0.0:
                quad[(k, l)] = e
    return BaselineTerms(n_pads=n_pads, linear=lin, quadratic=quad, name="cornwell_repulsion")


# --------------------------------------------------------------------------
# objective 3 -- PSF sidelobe energy (the onboarding guide's H_uv, exactly)
# --------------------------------------------------------------------------

def track_cell_multiplicities(
    pads: np.ndarray,
    wavelength: float,
    grid,
    observation=None,
    include_conjugate: bool = True,
) -> list:
    """For each candidate pair, ``{cell id: number of samples in that cell}``.

    Unlike :func:`thz_opt.qubo.coefficients.baseline_cell_sets` this keeps the
    *multiplicity*: an Earth-rotation track can revisit the same UV cell several
    times, and the sidelobe-energy objective below depends on that count, not
    just on which cells are touched.
    """
    from ..interferometry.earth_rotation import layout_to_uv_tracks
    from ..interferometry.uv import cell_ids

    arr = np.asarray(pads, dtype=float)
    i_idx, j_idx = pair_indices(arr.shape[0])

    out = []
    for a, b in zip(i_idx, j_idx):
        pair = arr[[int(a), int(b)]]
        if observation is None:
            uv = (pair[1:2] - pair[0:1]) / wavelength
            if include_conjugate:
                uv = np.vstack((uv, -uv))
        else:
            uv = layout_to_uv_tracks(pair, observation, include_conjugate)
        counts: dict = {}
        for c in cell_ids(uv, grid).tolist():
            counts[int(c)] = counts.get(int(c), 0) + 1
        out.append(counts)
    return out


def sidelobe_energy_terms(multiplicities: list, n_pads: int) -> BaselineTerms:
    """PSF sidelobe energy as an exact quadratic in the baseline variables.

    Let ``a_ck`` be the number of samples baseline ``k`` places in UV cell ``c``,
    so the gridded sampling function of a selection is ``n_c = sum_k a_ck y_k``.
    Using ``y^2 = y``,

        sum_c n_c^2 = sum_k (sum_c a_ck^2) y_k
                    + 2 sum_{k<l} (sum_c a_ck a_cl) y_k y_l

    which is exact -- verified to machine precision against direct gridding in
    ``tests/test_baseline_qubo.py``.

    Why this is the objective the onboarding guide asks for: by Parseval's
    theorem the dirty beam ``PSF = F^-1[n]`` satisfies

        sum_{l,m} |PSF(l,m)|^2 = sum_c n_c^2 / N^2

    and the beam peak is ``sum_c n_c / N^2``.  Normalising the peak to 1, the
    total energy in the beam is ``N^2 sum_c n_c^2 / (sum_c n_c)^2``.  So
    minimising ``sum_c n_c^2`` at a fixed number of samples minimises the energy
    outside the main lobe -- exactly the guide's statement that minimising the
    auto-correlation of the UV density minimises PSF sidelobes.

    Note this is *not* the same objective as maximising unique-cell coverage
    (:func:`bonferroni_coverage_terms`): coverage wants more occupied cells,
    sidelobe energy wants the occupancy spread evenly.  They agree on hating
    redundancy and disagree on everything else, so pick one deliberately or
    combine them with a stated weight.
    """
    lin = np.array([float(sum(v * v for v in d.values())) for d in multiplicities])

    quad: dict = {}
    by_cell: dict = {}
    for k, d in enumerate(multiplicities):
        for c, v in d.items():
            by_cell.setdefault(c, []).append((k, v))
    for entries in by_cell.values():
        if len(entries) < 2:
            continue
        for a in range(len(entries)):
            ka, va = entries[a]
            for b in range(a + 1, len(entries)):
                kb, vb = entries[b]
                key = (ka, kb) if ka < kb else (kb, ka)
                quad[key] = quad.get(key, 0.0) + 2.0 * va * vb

    return BaselineTerms(n_pads=n_pads, linear=lin, quadratic=quad,
                         name="psf_sidelobe_energy")


def coherence_weighted_sidelobe_terms(
    multiplicities: list,
    gamma: np.ndarray,
    n_pads: int,
    t: float = 0.0,
    gamma_power: float = 2.0,
) -> BaselineTerms:
    """Sidelobe energy with each baseline weighted by its surviving coherence.

    This is the THz-specific objective. A baseline that the atmosphere
    decorrelates contributes sampling but not signal, so the effective gridded
    weight is

        n_c = sum_k g_k a_ck y_k,      g_k = gamma_k ** gamma_power

    and the sidelobe energy is again exactly quadratic:

        sum_c n_c^2 = sum_k g_k^2 A_k y_k + 2 sum_{k<l} g_k g_l B_kl y_k y_l

    ``gamma_power``: use **2** if ``n_c`` is meant to be an inverse-variance
    (sensitivity) weight -- decorrelation scales signal by ``gamma`` at fixed
    noise, so SNR scales as ``gamma`` and the optimal weight as ``gamma^2``.
    Use **1** if ``n_c`` is meant to be the visibility amplitude response. The
    default is 2; which one is right depends on the assumed noise model and is
    a question for the collaboration, so it is a parameter.

    The normalisation problem, and its exact fix
    --------------------------------------------
    The quantity actually proportional to unit-peak beam energy is the *ratio*
    ``sum_c n_c^2 / (sum_c n_c)^2``.  With ``gamma == 1`` and equal track
    lengths the denominator is constant at fixed cardinality, so minimising the
    numerator alone is equivalent.  With coherence weighting it is **not** --
    different baselines carry different ``g_k`` -- and a ratio is not a
    quadratic form.

    The fix is Dinkelbach's parametric method for linear-fractional
    programming: minimising ``F(y)/G(y)`` is solved by repeatedly minimising

        F(y) - t * G(y)

    and updating ``t <- F(y*)/G(y*)`` until it stops moving.  Here
    ``G(y) = (sum_k g_k T_k y_k)^2`` with ``T_k`` the number of samples in
    track ``k``, which is **also** exactly quadratic in ``y``.  So each
    Dinkelbach iteration is an ordinary QUBO, and the normalised objective is
    reachable exactly rather than approximated.  Pass the current ``t`` here;
    ``t = 0`` gives the unnormalised numerator.
    """
    g = np.asarray(gamma, dtype=float) ** float(gamma_power)
    if g.shape[0] != len(multiplicities):
        raise ValueError("gamma must have one entry per candidate pair")

    lin = np.array([g[k] ** 2 * float(sum(v * v for v in d.values()))
                    for k, d in enumerate(multiplicities)])

    quad: dict = {}
    by_cell: dict = {}
    for k, d in enumerate(multiplicities):
        for c, v in d.items():
            by_cell.setdefault(c, []).append((k, v))
    for entries in by_cell.values():
        if len(entries) < 2:
            continue
        for a in range(len(entries)):
            ka, va = entries[a]
            for b in range(a + 1, len(entries)):
                kb, vb = entries[b]
                key = (ka, kb) if ka < kb else (kb, ka)
                quad[key] = quad.get(key, 0.0) + 2.0 * g[ka] * va * g[kb] * vb

    if t != 0.0:
        # subtract t * (sum_c n_c)^2 = t * (sum_k g_k T_k y_k)^2
        tot = np.array([g[k] * float(sum(d.values())) for k, d in enumerate(multiplicities)])
        lin = lin - t * tot ** 2
        for k in range(len(tot)):
            for l in range(k + 1, len(tot)):
                c = -2.0 * t * tot[k] * tot[l]
                if c != 0.0:
                    quad[(k, l)] = quad.get((k, l), 0.0) + c

    return BaselineTerms(n_pads=n_pads, linear=lin, quadratic=quad,
                         name=f"coherence_weighted_sidelobe(t={t})")


# --------------------------------------------------------------------------
# objective 4 -- target UV density matching
# --------------------------------------------------------------------------

def density_match_terms(
    uv_per_baseline: list,
    n_pads: int,
    target_counts: np.ndarray,
    bin_edges: np.ndarray,
    bin_weights: np.ndarray | None = None,
) -> BaselineTerms:
    """Squared mismatch to a target radial UV histogram, expanded over ``y``.

    With ``a_bk`` the number of samples baseline ``k`` puts in radial bin ``b``
    and ``t_b`` the target count,

        sum_b w_b ( sum_k a_bk y_k - t_b )^2
          = sum_k [ sum_b w_b (a_bk^2 - 2 t_b a_bk) ] y_k
          + sum_{k<l} [ 2 sum_b w_b a_bk a_bl ] y_k y_l
          + sum_b w_b t_b^2

    using ``y^2 = y``.  The target is supplied as counts so the caller decides
    the shape (Boone's Gaussian, a uniform disc, anything measured).
    """
    edges = np.asarray(bin_edges, dtype=float)
    t = np.asarray(target_counts, dtype=float)
    if t.shape[0] != edges.shape[0] - 1:
        raise ValueError("target_counts must have one entry per bin")
    w = np.ones_like(t) if bin_weights is None else np.asarray(bin_weights, dtype=float)

    a = np.zeros((t.shape[0], len(uv_per_baseline)))
    for k, track in enumerate(uv_per_baseline):
        arr = np.asarray(track, dtype=float).reshape(-1, 2)
        r = np.hypot(arr[:, 0], arr[:, 1])
        a[:, k], _ = np.histogram(r, bins=edges)

    lin = (w[:, None] * (a ** 2 - 2.0 * t[:, None] * a)).sum(axis=0)

    p = a.shape[1]
    quad: dict = {}
    gram = a.T @ (w[:, None] * a)  # gram[k, l] = sum_b w_b a_bk a_bl
    for k in range(p):
        for l in range(k + 1, p):
            c = 2.0 * float(gram[k, l])
            if c != 0.0:
                quad[(k, l)] = c
    return BaselineTerms(
        n_pads=n_pads,
        linear=lin,
        quadratic=quad,
        offset=float(np.sum(w * t ** 2)),
        name="density_match",
    )


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

def _max_row_weight(terms: "BaselineTerms") -> float:
    """``max_k (|linear_k| + sum_l |quadratic_kl|)`` -- the most the objective
    can change when a single baseline variable flips."""
    row = np.abs(terms.linear).astype(float).copy()
    for (k, l), c in terms.quadratic.items():
        row[k] += abs(c)
        row[l] += abs(c)
    return float(row.max() if row.size else 1.0)


def selection_penalty_floor(terms: "BaselineTerms", n_select: int, safety: float = 2.0) -> float:
    """A sufficient weight for the cardinality constraint.

    Adding one antenna beyond ``N`` switches on at most ``N`` new baselines, so
    the objective can improve by at most ``N * max_row_weight``.  The selection
    penalty for being one antenna over is ``lambda_select``, so

        lambda_select = safety * n_select * max_row_weight

    makes an over-full selection strictly worse.  Tying ``lambda_select`` to
    ``lambda_rosenberg`` instead -- the obvious shortcut -- is *not* enough: a
    fifth antenna buys four new baselines, and four small rewards can cancel
    one penalty exactly, which produced a tie in an earlier version of this
    code.
    """
    return float(safety * max(n_select, 1) * _max_row_weight(terms))


def rosenberg_penalty_floor(terms: BaselineTerms, safety: float = 2.0) -> float:
    """A sufficient Rosenberg weight.

    Breaking ``y_k = x_i x_j`` costs at least ``lambda_R``.  The most the
    objective can gain from one broken variable is
    ``|linear_k| + sum_l |quadratic_kl|``, so any
    ``lambda_R > max_k (that quantity)`` makes cheating strictly unprofitable.
    Conservative by construction; the factor is returned so it can be tightened
    deliberately.
    """
    return float(safety * _max_row_weight(terms))


def build_baseline_qubo(
    terms: BaselineTerms,
    n_select: int,
    lambda_select: float | None = None,
    lambda_rosenberg: float | None = None,
    safety: float = 2.0,
):
    """Assemble ``(Q, offset, meta)`` over ``[x (M) | y (M(M-1)/2)]``.

    ``Q`` is upper triangular in the convention of
    :mod:`thz_opt.qubo.objective` (one entry per unordered pair).  Variable
    ``M + k`` is the activation of candidate pair ``k`` in ``i < j`` order.
    """
    M = terms.n_pads
    i_idx, j_idx = pair_indices(M)
    P = len(i_idx)
    if terms.n_pairs != P:
        raise ValueError(f"terms cover {terms.n_pairs} pairs, expected {P} for M={M}")

    lam_r = rosenberg_penalty_floor(terms, safety) if lambda_rosenberg is None else float(lambda_rosenberg)
    lam_s = (selection_penalty_floor(terms, n_select, safety)
             if lambda_select is None else float(lambda_select))

    V = M + P
    Q = np.zeros((V, V))
    offset = terms.offset

    def add(u: int, v: int, val: float) -> None:
        a, b = (u, v) if u <= v else (v, u)
        Q[a, b] += val

    # selection constraint on x
    for i in range(M):
        add(i, i, lam_s * (1.0 - 2.0 * n_select))
        for j in range(i + 1, M):
            add(i, j, 2.0 * lam_s)
    offset += lam_s * n_select * n_select

    # Rosenberg:  y_k = x_i x_j
    for k in range(P):
        i, j, y = int(i_idx[k]), int(j_idx[k]), M + k
        add(i, j, lam_r)
        add(i, y, -2.0 * lam_r)
        add(j, y, -2.0 * lam_r)
        add(y, y, 3.0 * lam_r)

    # objective, purely in y
    for k in range(P):
        add(M + k, M + k, float(terms.linear[k]))
    for (k, l), c in terms.quadratic.items():
        add(M + k, M + l, float(c))

    meta = {
        "objective": terms.name,
        "n_variables": V,
        "n_x": M,
        "n_y": P,
        "n_select": n_select,
        "lambda_select": lam_s,
        "lambda_rosenberg": lam_r,
        "safety_factor": safety,
        "n_quadratic_terms": int(np.count_nonzero(np.triu(Q, 1))),
        **terms.sparsity(),
    }
    return Q, float(offset), meta


def complete_baseline_auxiliary(x: np.ndarray, n_pads: int) -> np.ndarray:
    """Full variable vector ``[x | y]`` with ``y_k = x_i x_j`` (zero penalty)."""
    v = np.asarray(x, dtype=int)
    return np.concatenate([v, selection_to_y(v, n_pads)])
