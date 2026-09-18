"""Experiment 09 -- the science case chooses the layout.

Every earlier result carried the same caveat: the objectives disagree, and
nothing said which one to use. This experiment removes the caveat by fixing two
concrete observations and optimising for each:

**compact source** -- an EHT-style shadow measurement. Small, bright, high
contrast. Resolution and dynamic range dominate, so UV cells are weighted
*upwards* with radius and sidelobe energy is penalised.

**extended emission** -- a molecular cloud or disc. Large, low contrast.
Short-spacing completeness dominates, so UV cells are weighted *downwards* with
radius.

Same pads, same N, same grid, same observing geometry. Only the weighting
changes. If the two optima differ, "the objectives disagree" stops being a
caveat and becomes a result.

Also measured here, because both change every number:

* **multi-frequency synthesis** -- how much coverage the receiver bandwidth
  adds for free, at fixed N and fixed pads
* **grid robustness** -- how much of a coverage number survives shifting the
  UV grid by fractions of a cell
* **the Pareto front** between weighted coverage and sidelobe energy, so no
  weight has to be guessed

Outputs
-------
figures/fig19_science_case_layouts.png
figures/fig20_pareto_front.png
figures/fig21_multifrequency_gain.png
data/results/exp09_science_cases.csv
data/results/exp09_pareto.csv
data/results/exp09_summary.json
"""

from __future__ import annotations

import time

import numpy as np

from common import (  # noqa: E402
    RESULTS,
    banner,
    ensure_dirs,
    load_config,
    new_figure,
    observation_config,
    save_csv,
    save_figure,
    save_json,
    stamp,
    wavelength,
)

from thz_opt.arrays.random_array import random_layout
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.gridding import grid_sensitivity
from thz_opt.interferometry.multifrequency import FrequencyBand, multi_frequency_uv
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, cell_ids, uv_occupancy
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_removal,
    local_search,
    pareto_indices,
    simulated_annealing,
)
from thz_opt.optimize.science_cases import (
    compact_source_weights,
    extended_emission_weights,
    science_case_objective,
)
from thz_opt.qubo.coefficients import exact_cells_covered

M_PADS, N_SELECT = 60, 15
R_MIN, R_MAX = 30.0, 1000.0
GRID_CELLS = 32
SEED = 20260916
BAND = FrequencyBand(3.0e11, 0.25, 5)


def multifreq_multiplicities(pads, obs, grid, band):
    """Per-pair cell multiplicities pooled over every channel of the band.

    This is the whole multi-frequency change: ``a_ck = sum_f a_ckf``. No new
    decision variables, so every objective stays exactly as quadratic as before.
    """
    from thz_opt.interferometry.baselines import pair_indices

    i_idx, j_idx = pair_indices(pads.shape[0])
    out = []
    for a, b in zip(i_idx, j_idx):
        uv = multi_frequency_uv(pads[[int(a), int(b)]], obs, band)
        counts: dict = {}
        for c in cell_ids(uv, grid).tolist():
            counts[int(c)] = counts.get(int(c), 0) + 1
        out.append(counts)
    return out


def cell_sets_from(multiplicities):
    return [np.fromiter(d.keys(), dtype=np.int64, count=len(d)) for d in multiplicities]


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    obs = observation_config(cfg)
    lam = wavelength(cfg)
    pads = random_layout(M_PADS, R_MIN, R_MAX, seed=SEED)
    sep = SeparationConfig(cfg["constraints"]["dish_diameter_m"],
                           cfg["constraints"]["minimum_separation_factor"])
    forbidden = shadow_pairs(pads, sep)

    uv_all = multi_frequency_uv(pads, obs, BAND)
    grid = UVGrid(uv_max=cfg["uv_grid"]["pad"] * float(np.abs(uv_all).max()),
                  n_cells=GRID_CELLS)
    ps = (f"M={M_PADS} pads, N={N_SELECT}, r=[{R_MIN:.0f}, {R_MAX:.0f}] m, "
          f"nu={BAND.centre_hz/1e9:.0f} GHz +/-{100*BAND.fractional_bandwidth/2:.0f}%, "
          f"{BAND.n_channels} channels, UV grid {GRID_CELLS}^2")
    banner(f"Experiment 09 -- science cases -- {ps}")

    # ---- multi-frequency gain -------------------------------------------
    print("\n  [1] what receiver bandwidth buys, at fixed pads and fixed N")
    mf_rows = []
    sel_ref = np.arange(N_SELECT)
    for fb, nch in ((0.0, 1), (0.05, 3), (0.10, 5), (0.20, 5), (0.30, 7)):
        band = FrequencyBand(BAND.centre_hz, fb, nch)
        uv = multi_frequency_uv(pads[sel_ref], obs, band)
        occ, _ = uv_occupancy(uv, grid)
        psf, info = dirty_beam(uv, grid, cfg["psf"]["weighting"])
        m = psf_metrics(psf, info["pixel_scale_arcsec"])
        mf_rows.append({"fractional_bandwidth": fb, "n_channels": nch,
                        "n_samples": int(uv.shape[0]),
                        "unique_cells": int((occ > 0).sum()),
                        "uv_radius_ratio": band.uv_radius_ratio,
                        "psf_peak_sidelobe": m["peak_sidelobe_level"]})
        print(f"      bandwidth {fb:5.0%} ({nch} ch): {int((occ>0).sum()):6d} cells, "
              f"PSL {m['peak_sidelobe_level']:.3f}")
    gain = mf_rows[-1]["unique_cells"] / mf_rows[0]["unique_cells"]
    print(f"      -> {gain:.2f}x more unique cells from bandwidth alone, "
          f"no extra antennas and no extra variables")

    # ---- grid robustness -------------------------------------------------
    print("\n  [2] how much of a coverage number is grid artifact")
    gs = grid_sensitivity(multi_frequency_uv(pads[sel_ref], obs, BAND), grid, n_shifts=3)
    print(f"      unique cells   : {gs['unique_cells']['mean']:.0f} "
          f"+/- {100*gs['unique_cells']['relative_spread']:.2f} % over shifted grids")
    print(f"      sidelobe energy: {100*gs['sidelobe_energy']['relative_spread']:.2f} % spread")

    # ---- the two science cases ------------------------------------------
    print("\n  [3] optimising each science case on identical inputs")
    t0 = time.time()
    mult = multifreq_multiplicities(pads, obs, grid, BAND)
    cells, counts = build_pair_tables(mult)
    cs = cell_sets_from(mult)
    print(f"      precompute: {time.time()-t0:.1f}s, {len(mult)} pairs")

    cases = {
        "compact_source": compact_source_weights(grid),
        "extended_emission": extended_emission_weights(grid),
    }
    results, rows = {}, []
    for case, w in cases.items():
        state = UVState(cells, counts, M_PADS, grid.n_cells ** 2, cell_weights=w)
        w_side = 1.0 if case == "compact_source" else 0.0
        # scales from a reference selection, so the two terms are comparable
        state.set_selection(sel_ref)
        score = science_case_objective(case, w_sidelobe=w_side,
                                       cell_scale=max(state.weighted_cells, 1.0),
                                       sidelobe_scale=max(state.weighted_sum_sq, 1.0))
        res = greedy_removal(state, N_SELECT, score, forbidden)
        state.set_selection(res.indices())
        res = local_search(state, score, forbidden)
        sa = simulated_annealing(state, N_SELECT, score, forbidden,
                                 n_iterations=6000, seed=SEED)
        best = sa if score_of(state, sa, score) <= res.score else res
        state.set_selection(best.indices())
        results[case] = best.indices()

        mask = np.zeros(M_PADS, dtype=int)
        mask[best.indices()] = 1
        uv = multi_frequency_uv(pads[best.indices()], obs, BAND)
        psf, info = dirty_beam(uv, grid, cfg["psf"]["weighting"])
        m = psf_metrics(psf, info["pixel_scale_arcsec"])
        d = np.hypot(*(pads[best.indices()][:, None, :] - pads[best.indices()][None, :, :]).T).T
        rows.append({
            "case": case, "selection": best.indices().tolist(),
            "unique_cells": exact_cells_covered(mask, cs, M_PADS),
            "weighted_cells": state.weighted_cells,
            "sum_sq": state.sum_sq,
            "median_baseline_m": float(np.median(d[np.triu_indices(N_SELECT, 1)])),
            "max_baseline_m": float(d.max()),
            "psf_peak_sidelobe": m["peak_sidelobe_level"],
            "psf_fwhm_arcsec": m["fwhm_arcsec"],
        })
        print(f"      {case:<19} cells {rows[-1]['unique_cells']:6d} | "
              f"median baseline {rows[-1]['median_baseline_m']:6.0f} m | "
              f"FWHM {m['fwhm_arcsec']:.4f}\" | PSL {m['peak_sidelobe_level']:.3f}")

    overlap = len(set(results["compact_source"].tolist())
                  & set(results["extended_emission"].tolist()))
    print(f"\n      the two optima share {overlap} of {N_SELECT} pads "
          f"({100*overlap/N_SELECT:.0f} %)")
    print(f"      median baseline differs by "
          f"{rows[0]['median_baseline_m'] / rows[1]['median_baseline_m']:.2f}x")

    # ---- how strongly must a science case be stated to matter? -----------
    print("\n  [3b] does a *preference* change the design, or only a *restriction*?")
    soft_hard = compare_soft_and_hard(pads, cells, counts, grid, forbidden)
    print(f"      {'weighting':<34}{'shared with unweighted':>24}{'median baseline':>18}")
    for r in soft_hard:
        print(f"      {r['weighting']:<34}{r['shared_with_flat']:>18} / {N_SELECT}"
              f"{r['median_baseline_m']:>18.0f} m")
    print("      -> a smooth radial preference barely moves the optimum; only")
    print("         excluding part of the UV plane does.")

    # ---- Pareto front ----------------------------------------------------
    print("\n  [4] Pareto front: weighted coverage against sidelobe energy")
    w = compact_source_weights(grid)
    state = UVState(cells, counts, M_PADS, grid.n_cells ** 2, cell_weights=w)
    state.set_selection(sel_ref)
    c_scale, s_scale = max(state.weighted_cells, 1.0), max(state.weighted_sum_sq, 1.0)

    pareto_rows = []
    for beta in np.linspace(0.0, 3.0, 13):
        sc = science_case_objective("compact_source", w_sidelobe=float(beta),
                                    cell_scale=c_scale, sidelobe_scale=s_scale)
        res = greedy_removal(state, N_SELECT, sc, forbidden)
        state.set_selection(res.indices())
        res = local_search(state, sc, forbidden)
        pareto_rows.append({"beta": float(beta),
                            "weighted_cells": state.weighted_cells,
                            "weighted_sum_sq": state.weighted_sum_sq,
                            "unique_cells": state.unique_cells,
                            "selection": res.indices().tolist()})
    obj = np.array([[-r["weighted_cells"] / c_scale, r["weighted_sum_sq"] / s_scale]
                    for r in pareto_rows])
    front = pareto_indices(obj)
    for i in front:
        pareto_rows[int(i)]["on_front"] = True
    for r in pareto_rows:
        r.setdefault("on_front", False)
    print(f"      {len(front)} of {len(pareto_rows)} weightings are non-dominated")

    save_csv(rows, RESULTS / "exp09_science_cases.csv")
    save_csv(soft_hard, RESULTS / "exp09_soft_vs_hard.csv")
    save_csv([{k: v for k, v in r.items() if k != "selection"} for r in pareto_rows],
             RESULTS / "exp09_pareto.csv")
    save_csv(mf_rows, RESULTS / "exp09_multifrequency.csv")
    save_json({"parameters": {"M": M_PADS, "N": N_SELECT, "grid_cells": GRID_CELLS,
                              "band": BAND.as_dict(), "seed": SEED},
               "science_cases": rows, "shared_pads": overlap,
               "soft_vs_hard": soft_hard,
               "multifrequency": mf_rows, "multifrequency_gain": gain,
               "grid_sensitivity": gs,
               "pareto": {"n_runs": len(pareto_rows), "n_front": int(len(front))}},
              RESULTS / "exp09_summary.json")

    _figures(pads, results, rows, pareto_rows, front, mf_rows, grid, obs, ps)
    print(f"\nWrote {RESULTS / 'exp09_science_cases.csv'} and exp09_summary.json")


def compare_soft_and_hard(pads, cells, counts, grid, forbidden) -> list:
    """How much does the optimum move as the science case is stated more strongly?

    Each profile re-weights UV cells and is optimised identically; the reported
    numbers are how many pads it shares with the unweighted optimum, and the
    median baseline of the layout it picks.
    """
    from thz_opt.optimize.science_cases import cell_radii, radial_weights

    r = cell_radii(grid)
    rb = grid.uv_max / 4.0
    profiles = [
        ("flat (no science case)", np.ones_like(r)),
        ("smooth tilt to long (r^1)", radial_weights(grid, "outer", rb, power=1.0)),
        ("smooth tilt to long (r^4)", radial_weights(grid, "outer", rb, power=4.0)),
        ("smooth tilt to short (exp)", radial_weights(grid, "inner", rb)),
        ("restricted to long (r > 2rb)", radial_weights(grid, "outer_cut", 2 * rb)),
        ("restricted to short (r < rb)", radial_weights(grid, "inner_cut", rb)),
    ]

    out, baseline_set = [], None
    for name, w in profiles:
        state = UVState(cells, counts, pads.shape[0], grid.n_cells ** 2, cell_weights=w)
        score = lambda st: -float(st.weighted_cells)  # noqa: E731
        res = greedy_removal(state, N_SELECT, score, forbidden)
        state.set_selection(res.indices())
        res = local_search(state, score, forbidden)
        sel = set(res.indices().tolist())
        if baseline_set is None:
            baseline_set = sel
        d = np.hypot(*(pads[res.indices()][:, None, :]
                       - pads[res.indices()][None, :, :]).T).T
        out.append({"weighting": name,
                    "shared_with_flat": len(sel & baseline_set),
                    "median_baseline_m": float(np.median(d[np.triu_indices(N_SELECT, 1)])),
                    "max_baseline_m": float(d.max()),
                    "unique_cells": int(state.unique_cells)})
    return out


def score_of(state, result, score):
    """Score a stored result without disturbing the caller's state."""
    keep = state.copy_selection()
    state.set_selection(result.indices())
    val = score(state)
    state.set_selection(np.flatnonzero(keep))
    return val


def _figures(pads, results, rows, pareto_rows, front, mf_rows, grid, obs, ps) -> None:
    # Figure 19 -- the two optimised layouts
    fig, axes = new_figure(1, 2, figsize=(10.4, 5.0))
    for ax, (case, sel) in zip(axes, results.items()):
        ax.scatter(pads[:, 0], pads[:, 1], s=10, c="0.8", label="candidate pads")
        ax.scatter(pads[sel][:, 0], pads[sel][:, 1], s=36, c="#4C72B0",
                   edgecolor="k", linewidth=0.3, label="selected")
        ax.set_aspect("equal")
        ax.set_title(case.replace("_", " "))
        ax.set_xlabel("East x (m)")
        ax.set_ylabel("North y (m)")
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("Figure 19. Same pads, same N -- different science case, different optimum",
                 y=1.01)
    stamp(fig, ps)
    save_figure(fig, "fig19_science_case_layouts.png")

    # Figure 20 -- Pareto front
    fig, ax = new_figure(figsize=(7.2, 5.0))
    x = [r["weighted_sum_sq"] for r in pareto_rows]
    y = [r["weighted_cells"] for r in pareto_rows]
    ax.scatter(x, y, s=26, c="0.6", label="all weightings")
    ax.scatter([x[int(i)] for i in front], [y[int(i)] for i in front],
               s=52, c="#C44E52", label="non-dominated")
    for r in pareto_rows:
        ax.annotate(f"{r['beta']:.2f}", (r["weighted_sum_sq"], r["weighted_cells"]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("weighted sidelobe energy (lower better)")
    ax.set_ylabel("weighted UV cells covered (higher better)")
    ax.set_title("Figure 20. Coverage against sidelobes\n"
                 "labels are the sidelobe weight; a scalar sweep finds only the convex front")
    ax.legend(fontsize=8)
    stamp(fig, ps)
    save_figure(fig, "fig20_pareto_front.png")

    # Figure 21 -- multi-frequency gain
    fig, axes = new_figure(1, 2, figsize=(11.0, 4.2))
    fb = [100 * r["fractional_bandwidth"] for r in mf_rows]
    axes[0].plot(fb, [r["unique_cells"] for r in mf_rows], "o-", color="#55A868")
    axes[0].set_xlabel("fractional bandwidth (%)")
    axes[0].set_ylabel("unique UV cells")
    axes[0].set_title("Figure 21a. Coverage from bandwidth alone\n(fixed pads, fixed N)")
    axes[1].plot(fb, [r["psf_peak_sidelobe"] for r in mf_rows], "s-", color="#C44E52")
    axes[1].set_xlabel("fractional bandwidth (%)")
    axes[1].set_ylabel("PSF peak sidelobe level")
    axes[1].set_title("Figure 21b. Effect on the beam")
    stamp(fig, ps)
    save_figure(fig, "fig21_multifrequency_gain.png")


if __name__ == "__main__":
    main()
