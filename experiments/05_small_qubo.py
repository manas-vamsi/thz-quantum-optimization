"""Experiment 05 -- toy QUBO: build, validate, enumerate, and measure the
pairwise approximation.

Toy setup (illustrative, NOT an instrument design):

    M = 8 candidate pads on a golden spiral, r in [20, 400] m
    N = 4 antennas to select
    wavelength = 1 mm
    snapshot UV for the QUBO itself; the resolution sweep repeats under
    Earth-rotation sampling as well

Four things are done, in order:

1. **Validation (mandatory).**  The matrix QUBO is compared with a direct,
   loop-based implementation of the term definitions over all 2^8 = 256
   configurations, and the upper-triangular and symmetric conventions are
   cross-checked.  Nothing else in this script is meaningful if this fails.

2. **Exhaustive optimisation.**  All 256 configurations are ranked, the
   constraint-satisfying optimum is identified, and it is confirmed that the
   penalty rule of ``suggest_penalty_weights`` keeps the global optimum
   feasible.

3. **The central research measurement.**  The UV grid is swept from coarse to
   fine.  For each resolution the script reports whether the exact
   unique-cell coverage is representable as a quadratic form at all (no cell
   contested by two candidate pairs), and how well each pairwise weighting
   ranks the C(8,4) = 70 feasible selections against exact coverage.

4. **The exact alternative.**  The auxiliary-variable construction is built for
   the same problem and its variable count is reported, quantifying the price
   of exactness.

Outputs
-------
figures/fig10_qubo_energy_landscape.png
data/results/exp05_qubo_validation.json
data/results/exp05_surrogate_sweep.csv
data/generated/exp05_qubo_matrix.npz
"""

from __future__ import annotations

import numpy as np

from common import (  # noqa: E402
    GENERATED,
    RESULTS,
    banner,
    ensure_dirs,
    load_config,
    new_figure,
    observation_config,
    save_csv,
    save_figure,
    save_json,
    save_npz,
    stamp,
)

from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation, shadow_pairs
from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid, layout_to_uv
from thz_opt.qubo.coefficients import (
    baseline_cell_sets,
    cell_claim_counts,
    coverage_is_exactly_quadratic,
    exact_cells_covered,
    pairwise_uv_weights,
    suggest_penalty_weights,
)
from thz_opt.qubo.exhaustive import enumerate_energies, feasible_configurations, rank_table
from thz_opt.qubo.objective import QUBOTerms, build_qubo, build_qubo_with_cell_variables, direct_terms
from thz_opt.qubo.validation import check_convention_roundtrip, compare_surrogate_to_exact, validate_qubo


def main() -> None:
    ensure_dirs()
    cfg = load_config()["toy_qubo"]
    M, N = cfg["n_pads"], cfg["n_select"]
    lam = float(cfg["wavelength_m"])
    pads = golden_spiral_layout(M, r_min=cfg["r_min_m"], r_max=cfg["r_max_m"])
    ps = (f"M={M} pads, N={N} selected, golden-spiral candidates, "
          f"r=[{cfg['r_min_m']:.0f}, {cfg['r_max_m']:.0f}] m, lambda={lam*1e3:.1f} mm, snapshot UV")
    banner(f"Experiment 05 -- toy QUBO -- {ps}")

    uv = layout_to_uv(pads, lam)
    grid = UVGrid.from_uv(uv, cfg["grid_cells_report"])
    cell_sets = baseline_cell_sets(pads, lam, grid)
    w = pairwise_uv_weights(cell_sets, grid, "shared")

    sep = SeparationConfig(6.0, 1.5)
    sepr = check_minimum_separation(pads, sep)
    shadow = shadow_pairs(pads, sep).astype(float)
    print(f"  d_min = {sepr['d_min_m']:.1f} m; closest candidate pair "
          f"{sepr['min_separation_m']:.1f} m; shadowed pairs: {int(np.triu(shadow, 1).sum())}")

    weights = suggest_penalty_weights(w, N, int(np.triu(shadow, 1).sum()),
                                      cfg["penalty_safety_factor"])
    print(f"  penalty rule: {weights['rule']}")
    print(f"  U_max = {weights['u_max']:.3f} -> lambda_select = lambda_shadow = "
          f"{weights['lambda_select']:.3f}  (toy values, labelled as such)")

    terms = QUBOTerms(n_pads=M, n_select=N, w_uv=w, shadow=shadow,
                      lambda_select=weights["lambda_select"],
                      lambda_shadow=weights["lambda_shadow"])
    Q, offset = build_qubo(terms)

    # ---- 1. validation -------------------------------------------------
    v = validate_qubo(terms)
    c = check_convention_roundtrip(Q, offset)
    print(f"\n  [1] direct objective vs QUBO over {v['n_configurations']} configurations: "
          f"max |difference| = {v['max_abs_difference']:.3e} -> "
          f"{'PASS' if v['within_tolerance'] else 'FAIL'}")
    print(f"      upper-triangular vs symmetric convention: "
          f"max |difference| = {c['max_abs_difference']:.3e} -> "
          f"{'PASS' if c['within_tolerance'] else 'FAIL'}")
    if not (v["within_tolerance"] and c["within_tolerance"]):
        raise SystemExit("QUBO validation failed -- no further results are meaningful")

    # ---- 2. exhaustive optimisation ------------------------------------
    configs, direct, quad = enumerate_energies(terms, Q, offset)
    feasible = configs.sum(axis=1) == N
    best_idx = int(np.argmin(quad))
    best = configs[best_idx]
    print(f"\n  [2] global optimum: {''.join(map(str, best))} "
          f"(E = {quad[best_idx]:.4f}, |x| = {int(best.sum())})")
    print(f"      best feasible {quad[feasible].min():.4f} vs best infeasible "
          f"{quad[~feasible].min():.4f} -> constraint "
          f"{'holds' if quad[feasible].min() < quad[~feasible].min() else 'VIOLATED'}")
    print("      lowest-energy configurations:")
    for bits, e in rank_table(configs, quad, 5):
        print(f"        {bits}  E = {e:>10.4f}  cells covered = "
              f"{exact_cells_covered(np.array(list(bits), dtype=int), cell_sets, M)}")
    print("      term breakdown at the optimum: " +
          ", ".join(f"{k}={val:.3f}" for k, val in direct_terms(best, terms).items()))

    # ---- 3. surrogate vs exact coverage, swept over UV resolution -------
    print("\n  [3] pairwise surrogate vs exact unique-cell coverage")
    print("      'flat-exact' marks a resolution where every feasible selection covers the")
    print("      same number of cells; 'flat-surr' where the surrogate cannot tell them apart.")
    print(f"      {'sampling':>14} {'cells/axis':>10} {'contested':>10} {'exactly quad':>13} "
          f"{'weighting':>10} {'spearman':>9} {'argmax ok':>10} {'shortfall %':>12}")
    obs_cfg = observation_config(load_config())
    sweep = []
    for sampling, obs in (("snapshot", None), ("earth_rotation", obs_cfg)):
        for n_cells in cfg["grid_cells_sweep"]:
            base_uv = uv if obs is None else layout_to_uv_tracks(pads, obs)
            g = UVGrid.from_uv(base_uv, n_cells)
            cs = baseline_cell_sets(pads, lam, g, observation=obs)
            claims = cell_claim_counts(cs)
            contested = sum(1 for v_ in claims.values() if v_ > 1)
            exact_ok = coverage_is_exactly_quadratic(cs)
            for weighting in cfg["weightings"]:
                ww = pairwise_uv_weights(cs, g, weighting)
                r = compare_surrogate_to_exact(cs, ww, M, N)
                sweep.append({"sampling": sampling, "n_cells": n_cells,
                              "cell_size_lambda": g.cell_size,
                              "contested_cells": contested, "exactly_quadratic": exact_ok,
                              "weighting": weighting, **r})
                if r["degenerate_surrogate"]:
                    tag = "flat-surr"
                elif r["degenerate_exact"]:
                    tag = "flat-exact"
                else:
                    tag = f"{r['spearman_rho']:.4f}"
                print(f"      {sampling:>14} {n_cells:>10} {contested:>10} {str(exact_ok):>13} "
                      f"{weighting:>10} {tag:>9} "
                      f"{str(r['surrogate_argmax_is_exact_argmax']):>10} "
                      f"{r['coverage_shortfall_percent']:>12.2f}")

    # ---- 4. the exact auxiliary-variable construction -------------------
    Qx, offx, meta = build_qubo_with_cell_variables(
        cell_sets, M, N, lambda_select=weights["lambda_select"], lambda_aux=weights["lambda_select"]
    )
    print(f"\n  [4] exact auxiliary construction on the {grid.n_cells}^2 grid: "
          f"{meta['n_variables']} binary variables "
          f"(x={meta['n_x']}, y={meta['n_y']}, z={meta['n_z']}, slack={meta['n_slack']}) "
          f"vs {M} for the pairwise QUBO")

    # ---- outputs -------------------------------------------------------
    save_npz(GENERATED / "exp05_qubo_matrix.npz", Q=Q, offset=np.array([offset]),
             w_uv=w, shadow=shadow, pads_xy_m=pads,
             configs=configs, energy_direct=direct, energy_qubo=quad)
    save_json({"parameters": {**cfg, "note": "toy values for illustration only"},
               "penalty_weights": weights, "validation": v, "convention_check": c,
               "optimum": {"bits": "".join(map(str, best)), "energy": float(quad[best_idx]),
                           "terms": direct_terms(best, terms),
                           "cells_covered": exact_cells_covered(best, cell_sets, M)},
               "exact_construction": {k: meta[k] for k in
                                      ("n_variables", "n_x", "n_y", "n_z", "n_slack")},
               "separation": sepr},
              RESULTS / "exp05_qubo_validation.json")
    save_csv(sweep, RESULTS / "exp05_surrogate_sweep.csv")

    _figure(configs, quad, feasible, cell_sets, M, N, w, ps, sweep)
    print(f"\nWrote {RESULTS / 'exp05_qubo_validation.json'} and exp05_surrogate_sweep.csv")


def _figure(configs, quad, feasible, cell_sets, M, N, w, ps, sweep) -> None:
    fig, axes = new_figure(1, 3, figsize=(13.5, 4.2))

    order = np.argsort(quad)
    axes[0].plot(np.arange(len(quad)), quad[order], lw=0.8, color="0.5",
                 label="all 2^8 configurations")
    fs = feasible[order]
    axes[0].scatter(np.flatnonzero(fs), quad[order][fs], s=10, color="#C44E52",
                    label=f"feasible (|x| = {N})", zorder=3)
    axes[0].set_title("Figure 10a. Sorted QUBO energy landscape")
    axes[0].set_xlabel("configuration rank (sorted by energy)")
    axes[0].set_ylabel("QUBO energy E(x)")
    axes[0].legend()

    pop = configs.sum(axis=1)
    axes[1].scatter(pop + np.random.default_rng(0).uniform(-0.15, 0.15, len(pop)), quad,
                    s=8, alpha=0.5, color="#4C72B0")
    axes[1].set_title("Figure 10b. Energy vs number of selected pads")
    axes[1].set_xlabel("number of selected pads |x|")
    axes[1].set_ylabel("QUBO energy E(x)")

    feas = feasible_configurations(M, N)
    exact = np.array([exact_cells_covered(c, cell_sets, M) for c in feas], dtype=float)
    surr = np.array([float(np.sum(w * c[np.triu_indices(M, 1)[0]] * c[np.triu_indices(M, 1)[1]]))
                     for c in feas])
    axes[2].scatter(surr, exact, s=14, color="#55A868")
    axes[2].set_title("Figure 10c. Pairwise surrogate vs exact coverage\n"
                      f"({len(feas)} feasible selections, shared weighting)")
    axes[2].set_xlabel("sum w_ij x_i x_j  (quadratic surrogate)")
    axes[2].set_ylabel("|Cov(x)|  (exact unique UV cells)")

    stamp(fig, ps)
    save_figure(fig, "fig10_qubo_energy_landscape.png")


if __name__ == "__main__":
    main()
