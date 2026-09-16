"""Experiment 07 -- three ways to put UV coverage into a QUBO, measured.

The objective a solver sees is a modelling choice, and the three available
choices trade accuracy against size very differently:

A. **Pairwise surrogate in pad space** (``thz_opt.qubo.coefficients``)
   ``M`` variables.  Cheapest possible.  Approximates unique coverage by a sum
   over single baselines, which cannot represent overlap between baselines.

B. **Bonferroni bound in baseline space** (``thz_opt.qubo.baseline_qubo``)
   ``M + M(M-1)/2`` variables.  Second-order inclusion-exclusion:

       |Cov| >= sum_k |C_k| y_k - sum_{k<l} |C_k ^ C_l| y_k y_l

   A *provable lower bound*, exactly quadratic in ``y``, exact whenever no UV
   cell is covered three or more times, and sparse because most baseline pairs
   never touch the same cell.

C. **Exact coverage with cell variables** (``thz_opt.qubo.objective``)
   ``M + pairs + cells + slack`` variables.  No approximation at all, and far
   too large to be practical beyond toy sizes.

This script measures all three on the same toy problem over a sweep of UV grid
resolutions, under Earth-rotation sampling: rank correlation against exact
coverage, whether the formulation's own optimum is the true optimum, how loose
the bound is, and what each costs in variables and couplings.

Outputs
-------
figures/fig16_formulation_comparison.png
data/results/exp07_formulations.csv
data/results/exp07_summary.json
"""

from __future__ import annotations

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

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid
from thz_opt.qubo.baseline_qubo import (
    baseline_cell_overlaps,
    bonferroni_bound,
    bonferroni_coverage_terms,
    build_baseline_qubo,
    complete_baseline_auxiliary,
    selection_to_y,
)
from thz_opt.qubo.coefficients import (
    baseline_cell_sets,
    cell_claim_counts,
    exact_cells_covered,
    pairwise_uv_weights,
    surrogate_value,
)
from thz_opt.qubo.exhaustive import feasible_configurations
from thz_opt.qubo.objective import build_qubo_with_cell_variables, qubo_energy

M, N_SEL = 8, 4
GRIDS = [8, 12, 16, 24, 32, 48, 64, 128]


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    obs = observation_config(cfg)
    lam = 299_792_458.0 / float(cfg["observation"]["frequency_hz"])
    pads = golden_spiral_layout(M, r_min=20.0, r_max=400.0)
    uv = layout_to_uv_tracks(pads, obs)
    configs = feasible_configurations(M, N_SEL)

    ps = (f"M={M} pads, N={N_SEL} selected, golden-spiral candidates, "
          f"lambda={lam*1e3:.3f} mm, Earth-rotation sampling, {obs.n_times} time samples")
    banner(f"Experiment 07 -- formulation comparison -- {ps}")
    print(f"  {len(configs)} feasible selections enumerated exhaustively\n")
    print(f"  {'cells/axis':>10}{'contested':>10}{'A rho':>9}{'B rho':>9}"
          f"{'A opt':>7}{'B opt':>7}{'bound gap %':>12}{'B fill %':>10}")
    print("  " + "-" * 76)

    rows = []
    for n_cells in GRIDS:
        grid = UVGrid.from_uv(uv, n_cells)
        cell_sets = baseline_cell_sets(pads, lam, grid, observation=obs)
        claims = cell_claim_counts(cell_sets)
        contested = sum(1 for v in claims.values() if v > 1)
        triple = sum(1 for v in claims.values() if v > 2)

        exact = np.array([exact_cells_covered(x, cell_sets, M) for x in configs], dtype=float)

        w = pairwise_uv_weights(cell_sets, grid, "shared")
        surr = np.array([surrogate_value(x, w, M) for x in configs])

        terms = bonferroni_coverage_terms(cell_sets, M)
        overlaps = baseline_cell_overlaps(cell_sets)
        bonf = np.array([bonferroni_bound(selection_to_y(x, M), cell_sets, overlaps)
                         for x in configs])

        Qb, offb, metab = build_baseline_qubo(terms, N_SEL)
        # the QUBO must reproduce the bound exactly at the intended completion
        resid = max(abs(qubo_energy(complete_baseline_auxiliary(x, M), Qb, offb) + b)
                    for x, b in zip(configs, bonf))

        rho_a, rho_b = _rho(surr, exact), _rho(bonf, exact)
        opt_a = bool(exact[int(np.argmax(surr))] >= exact.max() - 1e-9)
        opt_b = bool(exact[int(np.argmax(bonf))] >= exact.max() - 1e-9)
        gap = float(np.mean((exact - bonf) / np.maximum(exact, 1.0)))
        fill = metab["fill_fraction"]

        rows.append({
            "n_cells": n_cells, "cell_size_lambda": grid.cell_size,
            "contested_cells": contested, "triple_covered_cells": triple,
            "A_variables": M, "A_spearman": rho_a, "A_finds_optimum": opt_a,
            "B_variables": metab["n_variables"], "B_couplings": metab["n_couplings"],
            "B_fill_fraction": fill, "B_spearman": rho_b, "B_finds_optimum": opt_b,
            "B_bound_gap_fraction": gap, "B_qubo_residual": resid,
            "B_exact_here": bool(triple == 0),
        })
        print(f"  {n_cells:>10}{contested:>10}{rho_a:>9.4f}{rho_b:>9.4f}"
              f"{str(opt_a):>7}{str(opt_b):>7}{100*gap:>12.2f}{100*fill:>10.1f}")

    # cost of the exact cell-variable construction, for one representative grid
    grid = UVGrid.from_uv(uv, 16)
    cell_sets = baseline_cell_sets(pads, lam, grid, observation=obs)
    _, _, metac = build_qubo_with_cell_variables(
        cell_sets, M, N_SEL, lambda_select=100.0, lambda_aux=100.0, max_variables=200000
    )
    terms = bonferroni_coverage_terms(cell_sets, M)
    _, _, metab = build_baseline_qubo(terms, N_SEL)

    print("\n  Cost at 16 cells/axis, Earth-rotation sampling:")
    print(f"    A  pairwise surrogate (pad space)     {M:>7} variables")
    print(f"    B  Bonferroni bound (baseline space)  {metab['n_variables']:>7} variables, "
          f"{metab['n_couplings']:>6} couplings ({100*metab['fill_fraction']:.1f}% fill)")
    print(f"    C  exact coverage (cell variables)    {metac['n_variables']:>7} variables "
          f"(x={metac['n_x']}, y={metac['n_y']}, z={metac['n_z']}, slack={metac['n_slack']})")

    print("\n  Scaling of B, the formulation worth using, with pad count:")
    for m in (8, 20, 50, 100, 200):
        print(f"    M={m:>4}  ->  {m + m*(m-1)//2:>7} binary variables")

    save_csv(rows, RESULTS / "exp07_formulations.csv")
    save_json({"parameters": {"M": M, "N": N_SEL, "wavelength_m": lam,
                              "observation": obs.as_dict(), "grids": GRIDS},
               "rows": rows,
               "cost_at_16_cells": {"A_variables": M, "B": metab, "C": metac}},
              RESULTS / "exp07_summary.json")
    _figure(rows, ps)
    print(f"\nWrote {RESULTS / 'exp07_formulations.csv'} and exp07_summary.json")


def _figure(rows: list, ps: str) -> None:
    n = [r["n_cells"] for r in rows]
    fig, axes = new_figure(1, 3, figsize=(13.5, 4.2))

    axes[0].plot(n, [r["A_spearman"] for r in rows], "o-", ms=4,
                 label="A: pairwise surrogate (M vars)")
    axes[0].plot(n, [r["B_spearman"] for r in rows], "s-", ms=4,
                 label="B: Bonferroni bound (M + pairs)")
    axes[0].set_xscale("log", base=2)
    axes[0].set_xlabel("UV grid cells per axis")
    axes[0].set_ylabel("Spearman rho vs exact unique-cell coverage")
    axes[0].set_title("Figure 16a. How well each formulation ranks selections")
    axes[0].legend(loc="lower right")

    axes[1].plot(n, [100 * r["B_bound_gap_fraction"] for r in rows], "s-", ms=4, color="#C44E52")
    axes[1].set_xscale("log", base=2)
    axes[1].set_xlabel("UV grid cells per axis")
    axes[1].set_ylabel("mean (exact - bound) / exact  [%]")
    axes[1].set_title("Figure 16b. Tightness of the Bonferroni bound\n"
                      "(zero once no cell is triple-covered)")

    axes[2].plot(n, [100 * r["B_fill_fraction"] for r in rows], "s-", ms=4, color="#55A868")
    axes[2].set_xscale("log", base=2)
    axes[2].set_xlabel("UV grid cells per axis")
    axes[2].set_ylabel("coupling matrix fill [%]")
    axes[2].set_title("Figure 16c. Sparsity of the baseline-space QUBO")

    stamp(fig, ps)
    save_figure(fig, "fig16_formulation_comparison.png")


if __name__ == "__main__":
    main()
