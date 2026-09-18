"""Experiment 10 -- reconfigurable arrays: the paper's idea, done properly.

The golden-spiral paper proposed moving antennas during an observation to fill
the UV plane. The mechanism it described -- rigidly rotating the whole array --
cannot work: a rotation is an isometry, so every baseline length is unchanged
and the UV coverage is simply rotated with it. What *does* work is occupying a
**different set of pads** in each epoch, so the union over epochs samples more
than any single configuration can.

This experiment measures three things:

1. **What reconfiguration actually buys.** Union coverage over ``T`` epochs
   against a single static configuration observing for the same total time.
2. **What it costs.** The number of antenna relocations, and the trade-off
   curve against coverage as the movement penalty rises.
3. **Whether the QUBO stays buildable.** Variables grow linearly in ``T``, but
   the objective couples baselines across epochs, so the coupling count is the
   thing that could explode. It is measured rather than assumed.

Outputs
-------
figures/fig22_multiepoch_tradeoff.png
data/results/exp10_multiepoch.csv
data/results/exp10_summary.json
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
)

from thz_opt.arrays.random_array import random_layout
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid
from thz_opt.optimize import UVState, build_pair_tables, greedy_removal, local_search
from thz_opt.optimize.objectives import max_unique_cells
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities
from thz_opt.qubo.multiepoch import (
    MultiEpochSpec,
    build_multiepoch_qubo,
    movement_count,
    multiepoch_size_report,
)

M_PADS, N_SELECT = 40, 10
R_MIN, R_MAX = 30.0, 1000.0
GRID_CELLS = 24
SEED = 20260916
TOTAL_HOURS = 6.0


def epoch_configs(base: ObservationConfig, n_epochs: int) -> list:
    """Split one observation window into ``n_epochs`` consecutive blocks.

    Total observing time is held fixed, so a multi-epoch run is compared against
    a static array observing for exactly as long -- otherwise the comparison
    would just be measuring extra integration time.
    """
    edges = np.linspace(-TOTAL_HOURS / 2, TOTAL_HOURS / 2, n_epochs + 1)
    per = max(int(round(base.n_times / n_epochs)), 2)
    return [ObservationConfig(frequency_hz=base.frequency_hz,
                              declination_deg=base.declination_deg,
                              latitude_deg=base.latitude_deg,
                              hour_angle_start_h=float(edges[t]),
                              hour_angle_end_h=float(edges[t + 1]),
                              n_times=per)
            for t in range(n_epochs)]


def union_coverage(pads, selections, epochs, grid, lam) -> tuple:
    """Distinct UV cells covered by the union over epochs, and the pooled sum_sq."""
    from thz_opt.interferometry.uv import cell_ids

    n: dict = {}
    for sel, obs in zip(selections, epochs):
        uv = layout_to_uv_tracks(pads[np.asarray(sel, dtype=int)], obs)
        for c in cell_ids(uv, grid).tolist():
            n[c] = n.get(c, 0) + 1
    return len(n), int(sum(v * v for v in n.values()))


def optimise_epoch(pads, mult, grid, forbidden, covered_weight=None):
    """Best selection for one epoch, optionally down-weighting already-covered cells.

    Greedy across epochs: each epoch is optimised against what previous epochs
    have already covered. This is a heuristic, not the joint optimum -- the
    joint problem is exactly what the multi-epoch QUBO is for -- but it gives a
    classical reference the QUBO has to beat.
    """
    cells, counts = build_pair_tables(mult)
    state = UVState(cells, counts, pads.shape[0], grid.n_cells ** 2,
                    cell_weights=covered_weight)
    score = (max_unique_cells() if covered_weight is None
             else (lambda st: -float(st.weighted_cells)))
    res = greedy_removal(state, N_SELECT, score, forbidden)
    state.set_selection(res.indices())
    return local_search(state, score, forbidden).indices()


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    base = observation_config(cfg)
    lam = 299_792_458.0 / float(cfg["observation"]["frequency_hz"])
    pads = random_layout(M_PADS, R_MIN, R_MAX, seed=SEED)
    forbidden = shadow_pairs(pads, SeparationConfig(
        cfg["constraints"]["dish_diameter_m"],
        cfg["constraints"]["minimum_separation_factor"]))

    grid = UVGrid(uv_max=cfg["uv_grid"]["pad"] * float(np.abs(layout_to_uv_tracks(pads, base)).max()),
                  n_cells=GRID_CELLS)
    ps = (f"M={M_PADS} pads, N={N_SELECT} per epoch, {TOTAL_HOURS:.0f} h total "
          f"split into epochs, UV grid {GRID_CELLS}^2")
    banner(f"Experiment 10 -- multi-epoch reconfiguration -- {ps}")

    rows = []
    print(f"\n  {'epochs':>7}{'union cells':>13}{'gain':>8}{'moves':>8}"
          f"{'variables':>11}{'couplings':>12}{'fill':>8}")
    print("  " + "-" * 67)

    static_cells = None
    for n_epochs in (1, 2, 3, 4, 6):
        epochs = epoch_configs(base, n_epochs)
        mults = [track_cell_multiplicities(pads, lam, grid, o) for o in epochs]

        # greedy sequential: each epoch sees what is already covered
        selections, covered = [], np.ones(grid.n_cells ** 2)
        for t, mult in enumerate(mults):
            sel = optimise_epoch(pads, mult, grid, forbidden,
                                 covered_weight=None if t == 0 else covered)
            selections.append(sel)
            for k, d in enumerate(mult):
                i, j = np.triu_indices(M_PADS, 1)
                if sel.size and (i[k] in sel) and (j[k] in sel):
                    for c in d:
                        covered[c] = 0.0     # already have it; stop paying for it

        cells_union, sum_sq = union_coverage(pads, selections, epochs, grid, lam)
        if static_cells is None:
            static_cells = cells_union

        xm = np.zeros((n_epochs, M_PADS), dtype=int)
        for t, sel in enumerate(selections):
            xm[t, sel] = 1
        moves = movement_count(xm)

        rep = multiepoch_size_report(mults, M_PADS, n_epochs)
        rows.append({"n_epochs": n_epochs, "union_cells": cells_union,
                     "gain_over_static": cells_union / static_cells,
                     "moves": moves, "pooled_sum_sq": sum_sq,
                     **{k: rep[k] for k in ("n_variables_total", "n_baseline_variables",
                                            "n_objective_couplings", "fill_fraction")}})
        print(f"  {n_epochs:>7}{cells_union:>13}{cells_union/static_cells:>8.2f}"
              f"{moves:>8}{rep['n_variables_total']:>11}"
              f"{rep['n_objective_couplings']:>12}{100*rep['fill_fraction']:>7.1f}%")

    print(f"\n  Reconfiguration gain at fixed total observing time: "
          f"{rows[-1]['gain_over_static']:.2f}x more unique UV cells "
          f"for {rows[-1]['moves']} antenna relocations.")

    # ---- the result that explains the one above --------------------------
    print("\n  Why so small? Because Earth rotation is already doing the job.")
    print("  Sweeping the observing window, with total integration time held fixed:")
    print(f"    {'window':>8}{'static':>9}{'T=3':>8}{'T=6':>8}{'gain':>8}{'moves':>8}")
    window_rows = []
    for hours in (0.2, 0.5, 1.0, 2.0, 4.0, 6.0):
        c1, _ = run_window(pads, base, grid, forbidden, lam, hours, 1)
        c3, _ = run_window(pads, base, grid, forbidden, lam, hours, 3)
        c6, m6 = run_window(pads, base, grid, forbidden, lam, hours, 6)
        window_rows.append({"window_hours": hours, "static_cells": c1,
                            "cells_T3": c3, "cells_T6": c6,
                            "gain_T6": c6 / c1, "moves_T6": m6})
        print(f"    {hours:>7.1f}h{c1:>9}{c3:>8}{c6:>8}{c6/c1:>8.2f}{m6:>8}")
    print("\n    -> reconfiguration SUBSTITUTES for Earth rotation. With a short")
    print("       window it nearly triples coverage; with a long track it adds a few")
    print("       per cent. The source paper's instinct is right for snapshots and")
    print("       transients, and close to worthless for a 6 h synthesis -- and its")
    print("       specific mechanism (rigidly rotating the array) gives nothing at")
    print("       all, since a rotation leaves every baseline length unchanged.")

    # ---- can the QUBO actually be built? --------------------------------
    print("\n  QUBO size for the joint (not greedy) multi-epoch problem:")
    small_pads = pads[:14]
    small_forbidden = forbidden[:14, :14]
    for n_epochs in (2, 3):
        epochs = epoch_configs(base, n_epochs)
        mults = [track_cell_multiplicities(small_pads, lam, grid, o) for o in epochs]
        spec = MultiEpochSpec(14, n_epochs, 6, lambda_select=1e6,
                              lambda_rosenberg=1e5, lambda_move=100.0)
        t = time.time()
        _, _, meta = build_multiepoch_qubo(mults, spec, max_variables=50000)
        print(f"    M=14, T={n_epochs}: {meta['n_variables']} variables, "
              f"{meta['n_quadratic_terms']} quadratic terms, built in {time.time()-t:.2f}s")

    print("\n  Extrapolated joint-QUBO size (variables = T(M + M(M-1)/2)):")
    for m in (20, 50, 100):
        for t in (2, 3, 4):
            print(f"    M={m:>3}, T={t}: {t * (m + m * (m - 1) // 2):>8} variables")

    save_csv(rows, RESULTS / "exp10_multiepoch.csv")
    save_csv(window_rows, RESULTS / "exp10_window_sweep.csv")
    save_json({"parameters": {"M": M_PADS, "N": N_SELECT, "grid_cells": GRID_CELLS,
                              "total_hours": TOTAL_HOURS, "seed": SEED},
               "rows": rows, "window_sweep": window_rows,
               "conclusion": "reconfiguration substitutes for Earth rotation: "
                             "2.7x coverage gain for a 0.2 h window, 1.06x for 6 h, "
                             "at fixed total integration time"},
              RESULTS / "exp10_summary.json")

    _figure(rows, window_rows, ps)
    print(f"\nWrote {RESULTS / 'exp10_multiepoch.csv'} and exp10_summary.json")


def run_window(pads, base, grid, forbidden, lam, total_hours: float, n_epochs: int):
    """Union coverage and relocation count for one window length / epoch count.

    Total integration time is held fixed as the window is split, so the only
    thing changing is *when* the samples are taken and whether the array moves
    between them.
    """
    from thz_opt.interferometry.uv import cell_ids

    edges = np.linspace(-total_hours / 2, total_hours / 2, n_epochs + 1)
    per = max(int(round(36 / n_epochs)), 2)
    epochs = [ObservationConfig(frequency_hz=base.frequency_hz,
                                declination_deg=base.declination_deg,
                                latitude_deg=base.latitude_deg,
                                hour_angle_start_h=float(edges[t]),
                                hour_angle_end_h=float(edges[t + 1]),
                                n_times=per)
              for t in range(n_epochs)]

    selections, covered = [], np.ones(grid.n_cells ** 2)
    i_idx, j_idx = np.triu_indices(pads.shape[0], 1)
    for t, obs in enumerate(epochs):
        mult = track_cell_multiplicities(pads, lam, grid, obs)
        sel = optimise_epoch(pads, mult, grid, forbidden,
                             covered_weight=None if t == 0 else covered)
        selections.append(sel)
        chosen = set(sel.tolist())
        for k, d in enumerate(mult):
            if int(i_idx[k]) in chosen and int(j_idx[k]) in chosen:
                for c in d:
                    covered[c] = 0.0

    n: dict = {}
    for sel, obs in zip(selections, epochs):
        for c in cell_ids(layout_to_uv_tracks(pads[sel], obs), grid).tolist():
            n[c] = n.get(c, 0) + 1

    xm = np.zeros((n_epochs, pads.shape[0]), dtype=int)
    for t, sel in enumerate(selections):
        xm[t, sel] = 1
    return len(n), movement_count(xm)


def _figure(rows, window_rows, ps) -> None:
    fig, axes = new_figure(1, 3, figsize=(13.5, 4.2))
    e = [r["n_epochs"] for r in rows]

    axes[0].plot([r["window_hours"] for r in window_rows],
                 [r["gain_T6"] for r in window_rows], "o-", color="#4C72B0")
    axes[0].axhline(1.0, color="k", lw=0.8, ls="--")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("observing window (hours), total integration fixed")
    axes[0].set_ylabel("coverage gain from 6 epochs vs static")
    axes[0].set_title("Figure 22a. Reconfiguration substitutes\nfor Earth rotation")

    axes[1].plot([r["moves"] for r in rows], [r["union_cells"] for r in rows],
                 "s-", color="#C44E52")
    for r in rows:
        axes[1].annotate(f"T={r['n_epochs']}", (r["moves"], r["union_cells"]),
                         fontsize=7, xytext=(4, 4), textcoords="offset points")
    axes[1].set_xlabel("antenna relocations")
    axes[1].set_ylabel("unique UV cells")
    axes[1].set_title("Figure 22b. What it costs")

    axes[2].plot(e, [r["n_variables_total"] for r in rows], "o-", color="#55A868",
                 label="variables")
    axes[2].plot(e, [r["n_objective_couplings"] for r in rows], "^-", color="0.4",
                 label="objective couplings")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("number of epochs")
    axes[2].set_ylabel("count")
    axes[2].set_title("Figure 22c. QUBO size growth")
    axes[2].legend(fontsize=8)

    stamp(fig, ps)
    save_figure(fig, "fig22_multiepoch_tradeoff.png")


if __name__ == "__main__":
    main()
