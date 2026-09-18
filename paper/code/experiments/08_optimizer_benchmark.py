"""Experiment 08 -- searching for a layout, instead of guessing one.

Until now every "best layout" in this repository was the best of a hand-made
shortlist. This experiment actually optimises.

The comparison is set up to be fair to the hand-made designs, not to the
optimiser. The candidate pad set is built as the **union of the designed
layouts themselves** (golden spiral, Reuleaux, Vogel golden angle, random) plus
extra filler pads. So every designed layout is literally a feasible selection
the optimiser could return. Whatever it finds is therefore at least as good as
the best design, and the interesting number is by how much.

Methods compared, all on the identical objective, pad set, UV grid and
observing geometry:

* the designed layouts, scored as-is
* greedy removal (Panduranga Rao et al. 2009 style)
* greedy addition
* 1-swap local search, seeded from greedy
* simulated annealing, multi-restart
* **exact MILP via HiGHS** -- a provably optimal ceiling, which the
  array-configuration literature essentially never reports

Outputs
-------
figures/fig17_optimizer_comparison.png
figures/fig18_optimized_layout.png
data/results/exp08_optimizer_benchmark.csv
data/results/exp08_summary.json
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

from thz_opt.arrays.geometries import phyllotaxis_layout, reuleaux_layout
from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.arrays.random_array import random_layout
from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation, shadow_pairs
from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_addition,
    greedy_removal,
    local_search,
    max_coverage_milp,
    max_unique_cells,
    simulated_annealing,
)
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities
from thz_opt.qubo.coefficients import baseline_cell_sets, exact_cells_covered

N_SELECT = 20
N_FILLER = 40
R_MIN, R_MAX = 30.0, 1000.0
GRID_CELLS = 32
SEED = 20260916
MILP_TIME_LIMIT = 300.0


def build_candidate_pads():
    """Pad set containing every designed layout, plus filler.

    Returns ``(pads, designs)`` where ``designs`` maps a layout name to the pad
    indices that reproduce it exactly.
    """
    designed = {
        "golden spiral": golden_spiral_layout(N_SELECT, r_min=R_MIN, r_max=R_MAX, n_turns=3.0),
        "Reuleaux": reuleaux_layout(N_SELECT, R_MIN, R_MAX),
        "Vogel golden angle": phyllotaxis_layout(N_SELECT, R_MIN, R_MAX),
        "random": random_layout(N_SELECT, R_MIN, R_MAX, seed=SEED),
    }
    filler = random_layout(N_FILLER, R_MIN, R_MAX, seed=SEED + 1)

    pads = np.vstack(list(designed.values()) + [filler])
    designs, cursor = {}, 0
    for name, xy in designed.items():
        designs[name] = np.arange(cursor, cursor + xy.shape[0])
        cursor += xy.shape[0]
    return pads, designs


def score_selection(pads, indices, cfg, grid, obs, cell_sets):
    """Full physical scoring of one selection, independent of the optimiser."""
    sel = np.asarray(indices, dtype=int)
    xy = pads[sel]
    uv = layout_to_uv_tracks(xy, obs)
    occ, _ = uv_occupancy(uv, grid)
    psf, info = dirty_beam(uv, grid, cfg["psf"]["weighting"])
    m = psf_metrics(psf, info["pixel_scale_arcsec"])
    mask = np.zeros(pads.shape[0], dtype=int)
    mask[sel] = 1
    sep = check_minimum_separation(
        xy, SeparationConfig(cfg["constraints"]["dish_diameter_m"],
                             cfg["constraints"]["minimum_separation_factor"])
    )
    return {
        "unique_cells_direct": exact_cells_covered(mask, cell_sets, pads.shape[0]),
        "unique_cells_grid": int((occ > 0).sum()),
        "sum_sq": int((occ.astype(np.int64) ** 2).sum()),
        "psf_peak_sidelobe": m["peak_sidelobe_level"],
        "psf_sidelobe_rms": m["sidelobe_rms"],
        "psf_fwhm_arcsec": m["fwhm_arcsec"],
        "separation_violations": sep["n_violations"],
        "d_min_m": sep["min_separation_m"],
    }


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    obs = observation_config(cfg)
    lam = wavelength(cfg)

    pads, designs = build_candidate_pads()
    M = pads.shape[0]
    sep_cfg = SeparationConfig(cfg["constraints"]["dish_diameter_m"],
                               cfg["constraints"]["minimum_separation_factor"])
    forbidden = shadow_pairs(pads, sep_cfg)

    uv_all = layout_to_uv_tracks(pads, obs)
    grid = UVGrid(uv_max=cfg["uv_grid"]["pad"] * float(np.abs(uv_all).max()),
                  n_cells=GRID_CELLS)

    ps = (f"M={M} candidate pads, N={N_SELECT} selected, r=[{R_MIN:.0f}, {R_MAX:.0f}] m, "
          f"lambda={lam*1e3:.3f} mm, UV grid {GRID_CELLS}^2, "
          f"{obs.n_times} time samples over {obs.hour_angle_end_h - obs.hour_angle_start_h:.0f} h")
    banner(f"Experiment 08 -- optimiser benchmark -- {ps}")
    print(f"  shadowing: d_min = {sep_cfg.d_min:.1f} m, "
          f"{int(np.triu(forbidden, 1).sum())} forbidden pad pairs of {M*(M-1)//2}")

    print("  precomputing UV tracks and cell occupancy for all candidate pairs...")
    t0 = time.time()
    cell_sets = baseline_cell_sets(pads, lam, grid, observation=obs)
    cells, counts = build_pair_tables(track_cell_multiplicities(pads, lam, grid, obs))
    print(f"  done in {time.time() - t0:.1f}s "
          f"({len(cell_sets)} pairs, {sum(len(c) for c in cells)} cell entries)")

    state = UVState(cells, counts, M, grid.n_cells ** 2)
    score = max_unique_cells()
    rows = []

    # -- the designed layouts, scored as they are --------------------------
    print(f"\n  {'method':<26}{'cells':>8}{'PSL':>9}{'evals':>9}{'time (s)':>10}")
    print("  " + "-" * 62)
    best_design = 0
    for name, idx in designs.items():
        s = score_selection(pads, idx, cfg, grid, obs, cell_sets)
        best_design = max(best_design, s["unique_cells_direct"])
        rows.append({"method": f"design: {name}", "kind": "design",
                     "n_evaluations": 0, "seconds": 0.0,
                     "selection": idx.tolist(), **s})
        print(f"  {'design: ' + name:<26}{s['unique_cells_direct']:>8}"
              f"{s['psf_peak_sidelobe']:>9.3f}{0:>9}{0.0:>10.2f}")

    # -- the optimisers ----------------------------------------------------
    def run(label, fn):
        t = time.time()
        res = fn()
        dt = time.time() - t
        s = score_selection(pads, res.indices(), cfg, grid, obs, cell_sets)
        rows.append({"method": label, "kind": "optimiser",
                     "n_evaluations": res.n_evaluations, "seconds": dt,
                     "selection": res.indices().tolist(), **s, **res.meta})
        print(f"  {label:<26}{s['unique_cells_direct']:>8}"
              f"{s['psf_peak_sidelobe']:>9.3f}{res.n_evaluations:>9}{dt:>10.2f}")
        return res

    gr = run("greedy removal", lambda: greedy_removal(state, N_SELECT, score, forbidden))
    run("greedy addition", lambda: greedy_addition(state, N_SELECT, score, forbidden))

    def _ls():
        state.set_selection(gr.indices())
        return local_search(state, score, forbidden)
    run("local search (from greedy)", _ls)

    sa_best, sa_hist = None, []
    for k in range(3):
        res = run(f"annealing (restart {k + 1})",
                  lambda k=k: simulated_annealing(state, N_SELECT, score, forbidden,
                                                  n_iterations=8000, seed=SEED + k,
                                                  track_history=(k == 0)))
        if k == 0:
            sa_hist = res.history
        if sa_best is None or res.unique_cells > sa_best.unique_cells:
            sa_best = res

    # -- the provably optimal ceiling --------------------------------------
    print("\n  running exact MILP (HiGHS)...")
    t = time.time()
    ex = max_coverage_milp(cell_sets, M, N_SELECT, forbidden_pairs=forbidden,
                           time_limit=MILP_TIME_LIMIT)
    dt = time.time() - t
    if ex["selection"] is not None:
        s = score_selection(pads, ex["selection"], cfg, grid, obs, cell_sets)
        rows.append({"method": "MILP (exact)", "kind": "exact",
                     "n_evaluations": ex["n_variables"], "seconds": dt,
                     "selection": ex["selection"], "proven_optimal": ex["proven_optimal"],
                     **s})
        print(f"  {'MILP (exact)':<26}{s['unique_cells_direct']:>8}"
              f"{s['psf_peak_sidelobe']:>9.3f}{ex['n_variables']:>9}{dt:>10.2f}")
    print(f"  status: {ex['message']} | proven optimal: {ex['proven_optimal']} | "
          f"{ex['n_variables']} vars, {ex['n_constraints']} constraints")

    # -- the headline numbers ---------------------------------------------
    best_opt = max(r["unique_cells_direct"] for r in rows if r["kind"] != "design")
    ceiling = ex["coverage"] if ex["coverage"] is not None else best_opt
    print("\n  " + "=" * 62)
    print(f"  best hand-designed layout : {best_design:5d} cells")
    print(f"  best searched layout      : {best_opt:5d} cells  "
          f"(+{100 * (best_opt - best_design) / best_design:.1f} %)")
    if ex["proven_optimal"]:
        print(f"  proven optimum            : {ceiling:5d} cells")
        print(f"  best design is            : {100 * best_design / ceiling:.1f} % of optimal")
        for r in rows:
            if r["kind"] == "optimiser":
                r["gap_to_optimum_percent"] = 100 * (ceiling - r["unique_cells_direct"]) / ceiling
    print("  " + "=" * 62)

    # -- a smaller instance where optimality can actually be certified ------
    certified = certified_instance(cfg, obs, lam)

    save_csv([{k: v for k, v in r.items() if k != "selection"} for r in rows],
             RESULTS / "exp08_optimizer_benchmark.csv")
    save_csv(certified["rows"], RESULTS / "exp08_certified_instance.csv")
    save_json({"parameters": {"M": M, "N": N_SELECT, "grid_cells": GRID_CELLS,
                              "r_min_m": R_MIN, "r_max_m": R_MAX, "seed": SEED,
                              "observation": obs.as_dict()},
               "rows": rows, "milp": ex,
               "best_design": best_design, "best_searched": best_opt,
               "proven_optimum": ceiling if ex["proven_optimal"] else None,
               "certified_instance": certified},
              RESULTS / "exp08_summary.json")

    _figures(pads, designs, rows, sa_hist, ex, grid, obs, cfg, ps)
    print(f"\nWrote {RESULTS / 'exp08_optimizer_benchmark.csv'} and exp08_summary.json")


def certified_instance(cfg, obs, lam, m_pads: int = 24, n_select: int = 8,
                       grid_cells: int = 16, time_limit: float = 600.0) -> dict:
    """A small problem where the MILP can certify how good the heuristics are.

    Two things are reported, and the distinction matters:

    * **proven optimum** -- the solver closed the gap, so the value *is* the
      maximum. Only happens at small sizes.
    * **dual bound** -- an upper limit on what *any* selection could achieve,
      available even when the solver times out. A heuristic within 1 % of the
      dual bound is within 1 % of optimal, proven, without the gap ever closing.

    Measured scaling: at M = 40, N = 10 on a 24-cell grid, HiGHS does not close
    within 600 s and simulated annealing beats its incumbent. So exact
    certification is a small-instance tool, and the dual bound is what carries
    the argument at realistic sizes.
    """
    banner(f"Certified instance -- M={m_pads}, N={n_select}, grid {grid_cells}^2")
    pads = random_layout(m_pads, R_MIN, R_MAX, seed=SEED + 99)
    sep_cfg = SeparationConfig(cfg["constraints"]["dish_diameter_m"],
                               cfg["constraints"]["minimum_separation_factor"])
    forbidden = shadow_pairs(pads, sep_cfg)
    grid = UVGrid(uv_max=cfg["uv_grid"]["pad"] * float(np.abs(layout_to_uv_tracks(pads, obs)).max()),
                  n_cells=grid_cells)
    cell_sets = baseline_cell_sets(pads, lam, grid, observation=obs)
    cells, counts = build_pair_tables(track_cell_multiplicities(pads, lam, grid, obs))
    state = UVState(cells, counts, m_pads, grid.n_cells ** 2)
    score = max_unique_cells()

    t = time.time()
    ex = max_coverage_milp(cell_sets, m_pads, n_select, forbidden_pairs=forbidden,
                           time_limit=time_limit)
    milp_time = time.time() - t
    ceiling = ex["coverage"]
    bound = ex.get("upper_bound")
    print(f"  MILP incumbent: {ceiling} cells | proven optimal = {ex['proven_optimal']} | "
          f"{milp_time:.1f}s | {ex['n_variables']} variables")
    if bound is not None:
        print(f"  MILP dual bound: no selection can exceed {bound:.1f} cells")

    def _gap(value):
        """Gap to the certified ceiling (dual bound), not to the incumbent."""
        ref = bound if bound is not None else ceiling
        return 100.0 * (ref - value) / ref if ref else float("nan")

    rows = [{"method": "MILP (incumbent)", "unique_cells": ceiling,
             "proven_optimal": ex["proven_optimal"], "seconds": milp_time,
             "gap_percent": _gap(ceiling), "upper_bound": bound}]

    runs = {
        "greedy removal": lambda: greedy_removal(state, n_select, score, forbidden),
        "greedy addition": lambda: greedy_addition(state, n_select, score, forbidden),
        "annealing": lambda: simulated_annealing(state, n_select, score, forbidden,
                                                 n_iterations=8000, seed=SEED),
    }
    for label, fn in runs.items():
        t = time.time()
        res = fn()
        dt = time.time() - t
        gap = _gap(res.unique_cells)
        rows.append({"method": label, "unique_cells": res.unique_cells,
                     "proven_optimal": False, "seconds": dt, "gap_percent": gap,
                     "upper_bound": bound})
        print(f"  {label:<22} {res.unique_cells:5d} cells | gap to certified ceiling "
              f"{gap:5.2f} % | {dt:.1f}s")

    # local search seeded from greedy, which is the practical recipe
    state.set_selection(greedy_removal(state, n_select, score, forbidden).indices())
    t = time.time()
    res = local_search(state, score, forbidden)
    dt = time.time() - t
    gap = _gap(res.unique_cells)
    rows.append({"method": "greedy + local search", "unique_cells": res.unique_cells,
                 "proven_optimal": False, "seconds": dt, "gap_percent": gap,
                 "upper_bound": bound})
    print(f"  {'greedy + local search':<22} {res.unique_cells:5d} cells | "
          f"gap to certified ceiling {gap:5.2f} % | {dt:.1f}s")

    best_heuristic = max(r["unique_cells"] for r in rows[1:])
    print(f"\n  best heuristic {best_heuristic} cells is within "
          f"{_gap(best_heuristic):.2f} % of a proven ceiling"
          f"{' (optimum proven)' if ex['proven_optimal'] else ''}")

    return {"m_pads": m_pads, "n_select": n_select, "grid_cells": grid_cells,
            "milp": {k: v for k, v in ex.items() if k != "selection"},
            "incumbent": ceiling, "upper_bound": bound,
            "proven_optimal": ex["proven_optimal"],
            "best_heuristic": best_heuristic, "rows": rows}


def _figures(pads, designs, rows, sa_hist, ex, grid, obs, cfg, ps) -> None:
    # Figure 17 -- method comparison + annealing convergence
    fig, axes = new_figure(1, 2, figsize=(13.0, 5.0))

    labels = [r["method"] for r in rows]
    vals = [r["unique_cells_direct"] for r in rows]
    colors = {"design": "#C44E52", "optimiser": "#4C72B0", "exact": "#55A868"}
    axes[0].barh(range(len(labels)), vals,
                 color=[colors[r["kind"]] for r in rows])
    axes[0].set_yticks(range(len(labels)))
    axes[0].set_yticklabels(labels, fontsize=7)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("unique UV cells covered (higher is better)")
    axes[0].set_title("Figure 17a. Hand-designed vs searched vs provably optimal\n"
                      "red = design, blue = heuristic search, green = exact MILP")
    if ex["coverage"] is not None and ex["proven_optimal"]:
        axes[0].axvline(ex["coverage"], color="k", ls="--", lw=0.9,
                        label=f"proven optimum = {ex['coverage']}")
        axes[0].legend(loc="lower right", fontsize=7)

    if sa_hist:
        it = [h["iteration"] for h in sa_hist]
        axes[1].plot(it, [-h["current"] for h in sa_hist], lw=0.7, color="0.6",
                     label="current")
        axes[1].plot(it, [-h["best"] for h in sa_hist], lw=1.4, color="#4C72B0",
                     label="best so far")
        if ex["coverage"] is not None and ex["proven_optimal"]:
            axes[1].axhline(ex["coverage"], color="k", ls="--", lw=0.9,
                            label="proven optimum")
        axes[1].set_xlabel("annealing iteration")
        axes[1].set_ylabel("unique UV cells")
        axes[1].set_title("Figure 17b. Simulated-annealing convergence")
        axes[1].legend(loc="lower right", fontsize=7)
    stamp(fig, ps)
    save_figure(fig, "fig17_optimizer_comparison.png")

    # Figure 18 -- the optimised layout and its UV coverage
    best = max((r for r in rows if r["kind"] != "design"),
               key=lambda r: r["unique_cells_direct"])
    sel = np.asarray(best["selection"], dtype=int)
    xy = pads[sel]
    uv = layout_to_uv_tracks(xy, obs)
    occ, _ = uv_occupancy(uv, grid)
    gold = pads[designs["golden spiral"]]
    uv_gold = layout_to_uv_tracks(gold, obs)
    occ_gold, _ = uv_occupancy(uv_gold, grid)

    fig, axes = new_figure(1, 3, figsize=(14.0, 4.6))
    axes[0].scatter(pads[:, 0], pads[:, 1], s=12, c="0.75", label="candidate pads")
    axes[0].scatter(xy[:, 0], xy[:, 1], s=34, c="#4C72B0", edgecolor="k",
                    linewidth=0.3, label=f"selected ({best['method']})")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("East x (m)")
    axes[0].set_ylabel("North y (m)")
    axes[0].set_title("Figure 18a. Optimised selection")
    axes[0].legend(fontsize=7, loc="upper right")

    for ax, o, ttl in ((axes[1], occ, f"Figure 18b. Optimised UV ({int((occ>0).sum())} cells)"),
                       (axes[2], occ_gold,
                        f"Figure 18c. Golden spiral UV ({int((occ_gold>0).sum())} cells)")):
        ax.imshow((o.T > 0).astype(float), origin="lower", cmap="Greys",
                  extent=[-grid.uv_max / 1e3, grid.uv_max / 1e3,
                          -grid.uv_max / 1e3, grid.uv_max / 1e3])
        ax.set_aspect("equal")
        ax.set_xlabel("u (kilo-wavelengths)")
        ax.set_ylabel("v (kilo-wavelengths)")
        ax.set_title(ttl)
    stamp(fig, ps)
    save_figure(fig, "fig18_optimized_layout.png")


if __name__ == "__main__":
    main()
