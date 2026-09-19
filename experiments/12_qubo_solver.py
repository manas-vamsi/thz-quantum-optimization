"""Experiment 12 -- solving the QUBO, and separating the two gaps that matter.

Every earlier experiment either *formulated* a QUBO or searched the selection
problem directly with a classical heuristic. Neither tells you whether the QUBO
route actually works. This one solves the QUBO with a real annealing sampler
(``dwave-samplers``, the same software interface used to drive D-Wave hardware)
and measures how far its answer lands from the proven optimum.

The distinction this experiment exists to make
----------------------------------------------
"The QUBO gave a worse answer" can mean two completely different things, and
conflating them is the standard way this kind of study misleads:

**Formulation gap.** The QUBO does not optimise unique-cell coverage; it
optimises the second-order Bonferroni *lower bound* on coverage, because that
is the part of coverage which is quadratic. Even a perfect solver returns the
minimiser of the bound, not of coverage. This gap is a property of the
mathematics and no solver can close it.

**Solver gap.** Given the QUBO, how close does the sampler get to that QUBO's
own optimum? This is a property of the solver and closes with more effort.

Part A measures the first by brute force: on a 16-pad instance every feasible
selection is enumerated, so both the true coverage optimum and the QUBO's own
optimum are known exactly. Part B measures the second against the MILP
certificates on larger real instances.

A third number decides whether any of this is worth doing: the same instances
are also solved by the direct classical heuristic, which optimises true
coverage and ignores the QUBO entirely.

What this is not
----------------
No quantum hardware was used and no quantum advantage is claimed or implied.
Simulated annealing on a QUBO is a classical algorithm. What it establishes is
that the formulation is *solvable by a sampler of the kind quantum annealers
expose*, and what accuracy that costs -- which is the prerequisite for a
hardware run, not a substitute for one.

Outputs
-------
figures/fig25_qubo_gap_decomposition.png
figures/fig26_qubo_solver_convergence.png
data/results/exp12_gap_decomposition.csv
data/results/exp12_solver_benchmark.csv
data/results/exp12_summary.json
"""

from __future__ import annotations

import itertools
import time

import numpy as np

from common import (  # noqa: E402
    RESULTS,
    ROOT,
    banner,
    ensure_dirs,
    load_config,
    new_figure,
    save_csv,
    save_figure,
    save_json,
    stamp,
)

import dimod
from dwave.samplers import SimulatedAnnealingSampler

from thz_opt.arrays.real_arrays import REAL_SITES, load_cfg
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_removal,
    local_search,
    max_coverage_milp,
    max_unique_cells,
)
from thz_opt.qubo.baseline_qubo import (
    bonferroni_coverage_terms,
    track_cell_multiplicities,
    build_baseline_qubo,
    complete_baseline_auxiliary,
    objective_value,
    selection_to_y,
)
from thz_opt.qubo.coefficients import baseline_cell_sets, exact_cells_covered
from thz_opt.qubo.objective import qubo_energy
from thz_opt.interferometry.baselines import pair_indices

FREQ_HZ = 3.0e11
SEED = 20260919
CFG_FILE = ROOT / "data" / "external" / "alma.all.cfg"

#: (M, N, uv grid cells) for each instance. The first is small enough to
#: enumerate exhaustively; the other two match the MILP certificates of
#: experiment 11, so the proven optimum is the yardstick.
INSTANCES = [(16, 6, 12), (24, 8, 16), (40, 10, 16)]
READ_SWEEP = (1, 5, 20, 100, 500)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def build_instance(pads, lam, obs, sep_cfg, m_sub, n_cells, order):
    """A real contiguous sub-array, its UV cell sets and its shadow rule."""
    sub = order[:m_sub]
    sub_pads = pads[sub]
    grid = UVGrid(
        uv_max=1.02 * float(np.abs(layout_to_uv_tracks(sub_pads, obs)).max()),
        n_cells=n_cells,
    )
    cell_sets = baseline_cell_sets(sub_pads, lam, grid, observation=obs)
    mult = track_cell_multiplicities(sub_pads, lam, grid, obs)
    return sub_pads, grid, cell_sets, shadow_pairs(sub_pads, sep_cfg), mult


def decode(sample, m_pads, n_select):
    """Split a sampler assignment into ``x``, and audit the penalty terms.

    A QUBO enforces its constraints by penalty, so a returned sample is *not*
    guaranteed feasible. Reporting the energy without checking this is how a
    QUBO study quietly reports infeasible solutions as wins.
    """
    v = np.array([sample[i] for i in range(len(sample))], dtype=int)
    x, y = v[:m_pads], v[m_pads:]
    i_idx, j_idx = pair_indices(m_pads)
    return {
        "x": x,
        "n_selected": int(x.sum()),
        "cardinality_ok": int(x.sum()) == n_select,
        "rosenberg_violations": int(np.sum(y != (x[i_idx] * x[j_idx]))),
    }


def solve_qubo(Q, offset, m_pads, n_select, num_reads, seed):
    bqm = dimod.BinaryQuadraticModel(Q, "BINARY")
    bqm.offset = offset
    t0 = time.time()
    sampleset = SimulatedAnnealingSampler().sample(
        bqm, num_reads=num_reads, seed=seed, num_sweeps=1000
    )
    elapsed = time.time() - t0
    best = sampleset.first
    info = decode(best.sample, m_pads, n_select)
    info.update(energy=float(best.energy), seconds=elapsed, num_reads=num_reads)
    # how often the sampler lands on a feasible assignment at all
    feasible = sum(
        1
        for s in sampleset.samples()
        if decode(s, m_pads, n_select)["cardinality_ok"]
        and decode(s, m_pads, n_select)["rosenberg_violations"] == 0
    )
    info["feasible_fraction"] = feasible / len(sampleset)
    return info, sampleset


# --------------------------------------------------------------------------
# Part A -- exhaustive gap decomposition
# --------------------------------------------------------------------------

def gap_decomposition(cell_sets, m_pads, n_select, forbidden):
    """Enumerate every feasible selection; report both optima exactly."""
    terms = bonferroni_coverage_terms(cell_sets, m_pads)
    rows = []
    best_true = (-1, None)
    best_qubo = (np.inf, None)
    n_feasible = 0

    for combo in itertools.combinations(range(m_pads), n_select):
        x = np.zeros(m_pads, dtype=int)
        x[list(combo)] = 1
        if any(forbidden[a, b] for a, b in itertools.combinations(combo, 2)):
            continue
        n_feasible += 1
        true_cov = exact_cells_covered(x, cell_sets, m_pads)
        q = objective_value(selection_to_y(x, m_pads), terms)
        rows.append((true_cov, q))
        if true_cov > best_true[0]:
            best_true = (true_cov, x.copy())
        if q < best_qubo[0]:
            best_qubo = (q, x.copy())

    arr = np.array(rows, dtype=float)
    cov_of_qubo_opt = exact_cells_covered(best_qubo[1], cell_sets, m_pads)
    return {
        "n_feasible": n_feasible,
        "true_optimum": int(best_true[0]),
        "qubo_optimum_energy": float(best_qubo[0]),
        "coverage_at_qubo_optimum": int(cov_of_qubo_opt),
        "formulation_gap_cells": int(best_true[0] - cov_of_qubo_opt),
        "formulation_gap_percent": 100.0 * (best_true[0] - cov_of_qubo_opt) / best_true[0],
        "rank_correlation": float(
            np.corrcoef(arr[:, 0], -arr[:, 1])[0, 1]
        ),
        "coverage_mean": float(arr[:, 0].mean()),
        "coverage_worst": int(arr[:, 0].min()),
        "scatter": arr,
        "x_true": best_true[1],
        "x_qubo": best_qubo[1],
    }


# --------------------------------------------------------------------------

def main() -> None:
    ensure_dirs()
    cfg = load_config()
    site = REAL_SITES["alma"]
    if not CFG_FILE.exists():
        raise SystemExit(f"missing real pad file: {CFG_FILE}")

    pads_cfg = load_cfg(CFG_FILE, dish_diameter_m=site["dish_diameter_m"])
    pads = pads_cfg.xy
    lam = 299_792_458.0 / FREQ_HZ
    obs = ObservationConfig(
        frequency_hz=FREQ_HZ,
        declination_deg=-30.0,
        latitude_deg=site["latitude_deg"],
        hour_angle_start_h=-2.0, hour_angle_end_h=2.0, n_times=41,
    )
    sep_cfg = SeparationConfig(site["dish_diameter_m"],
                               cfg["constraints"]["minimum_separation_factor"])
    order = np.argsort(np.hypot(pads[:, 0], pads[:, 1]))

    ps = (f"real ALMA pads, nu={FREQ_HZ/1e9:.0f} GHz, 4 h rotation, "
          f"Bonferroni coverage QUBO, simulated annealing")
    banner(f"Experiment 12 -- solving the QUBO -- {ps}")

    # ---- Part A ---------------------------------------------------------
    m_a, n_a, g_a = INSTANCES[0]
    sub_pads, grid, cell_sets, forbidden, _ = build_instance(
        pads, lam, obs, sep_cfg, m_a, g_a, order)
    print(f"\n  [A] formulation gap by exhaustive enumeration, M={m_a}, N={n_a}")
    t0 = time.time()
    gap = gap_decomposition(cell_sets, m_a, n_a, forbidden)
    print(f"      enumerated {gap['n_feasible']} feasible selections in "
          f"{time.time()-t0:.1f}s")
    print(f"      true coverage optimum              : {gap['true_optimum']} cells")
    print(f"      coverage at the QUBO's own optimum : "
          f"{gap['coverage_at_qubo_optimum']} cells")
    print(f"      formulation gap                    : "
          f"{gap['formulation_gap_cells']} cells "
          f"({gap['formulation_gap_percent']:.1f} %)")
    print(f"      correlation, QUBO score vs true coverage: "
          f"{gap['rank_correlation']:.4f}")
    print(f"      (mean over all feasible selections: {gap['coverage_mean']:.1f}, "
          f"worst: {gap['coverage_worst']})")

    # ---- Part B ---------------------------------------------------------
    print("\n  [B] solver gap against the proven optimum")
    bench = []
    convergence = {}
    for m_sub, n_sub, g_cells in INSTANCES:
        sub_pads, grid, cs, forb, mult = build_instance(
            pads, lam, obs, sep_cfg, m_sub, g_cells, order)
        terms = bonferroni_coverage_terms(cs, m_sub)
        Q, off, meta = build_baseline_qubo(terms, n_sub)

        t0 = time.time()
        ex = max_coverage_milp(cs, m_sub, n_sub, forbidden_pairs=forb, time_limit=300.0)
        milp_s = time.time() - t0
        opt, proven = ex["coverage"], ex["proven_optimal"]

        print(f"\n      M={m_sub}, N={n_sub}: {meta['n_variables']} QUBO variables "
              f"({meta['n_x']} x + {meta['n_y']} y), "
              f"{meta['n_quadratic_terms']} couplings")
        print(f"      MILP optimum {opt} cells (proven={proven}, {milp_s:.0f}s), "
              f"lambda_R={meta['lambda_rosenberg']:.3g}, "
              f"lambda_sel={meta['lambda_select']:.3g}")
        print(f"        {'reads':>6}{'energy':>14}{'cells':>8}{'vs opt':>9}"
              f"{'feas %':>9}{'time s':>9}")

        curve = []
        for reads in READ_SWEEP:
            info, _ = solve_qubo(Q, off, m_sub, n_sub, reads, SEED)
            cov = (exact_cells_covered(info["x"], cs, m_sub)
                   if info["cardinality_ok"] else 0)
            rel = 100.0 * cov / opt if opt else float("nan")
            curve.append((reads, cov, rel, info["seconds"]))
            bench.append({
                "m_pads": m_sub, "n_select": n_sub, "num_reads": reads,
                "qubo_variables": meta["n_variables"],
                "energy": info["energy"], "cells": cov,
                "milp_optimum": opt, "milp_proven": proven,
                "percent_of_optimum": rel,
                "cardinality_ok": info["cardinality_ok"],
                "rosenberg_violations": info["rosenberg_violations"],
                "feasible_fraction": info["feasible_fraction"],
                "seconds": info["seconds"],
            })
            flag = "" if info["cardinality_ok"] and not info["rosenberg_violations"] else "  INFEASIBLE"
            print(f"        {reads:>6}{info['energy']:>14.1f}{cov:>8}"
                  f"{rel:>8.1f}%{100*info['feasible_fraction']:>8.0f}%"
                  f"{info['seconds']:>9.2f}{flag}")
        convergence[f"M{m_sub}_N{n_sub}"] = {"curve": curve, "optimum": opt}

        # the classical heuristic, for scale
        cells_t, counts_t = build_pair_tables(mult)
        state = UVState(cells_t, counts_t, m_sub, grid.n_cells ** 2)
        score = max_unique_cells()
        t0 = time.time()
        res = greedy_removal(state, n_sub, score, forb)
        state.set_selection(res.selection)
        res = local_search(state, score, forb)
        h_s = time.time() - t0
        x_h = np.asarray(res.selection, dtype=int)   # boolean mask over pads
        cov_h = exact_cells_covered(x_h, cs, m_sub)
        print(f"        direct classical heuristic: {cov_h} cells "
              f"({100*cov_h/opt:.1f} % of optimum) in {h_s:.2f}s")
        bench.append({
            "m_pads": m_sub, "n_select": n_sub, "num_reads": None,
            "qubo_variables": m_sub, "energy": None, "cells": cov_h,
            "milp_optimum": opt, "milp_proven": proven,
            "percent_of_optimum": 100.0 * cov_h / opt,
            "cardinality_ok": True, "rosenberg_violations": 0,
            "feasible_fraction": 1.0, "seconds": h_s,
            "method": "direct_heuristic",
        })
        convergence[f"M{m_sub}_N{n_sub}"]["heuristic"] = cov_h

    # ---- outputs --------------------------------------------------------
    save_csv([{k: v for k, v in r.items()} for r in bench],
             RESULTS / "exp12_solver_benchmark.csv")
    save_csv([{"true_cells": int(a), "qubo_energy": float(b)}
              for a, b in gap["scatter"]],
             RESULTS / "exp12_gap_decomposition.csv")
    summary = {
        "instances": INSTANCES,
        "read_sweep": list(READ_SWEEP),
        "sampler": "dwave-samplers SimulatedAnnealingSampler, 1000 sweeps",
        "quantum_hardware_used": False,
        "formulation_gap": {k: v for k, v in gap.items()
                            if k not in ("scatter", "x_true", "x_qubo")},
        "benchmark": bench,
    }
    save_json(summary, RESULTS / "exp12_summary.json")
    _figures(gap, convergence, ps)
    print(f"\nWrote {RESULTS/'exp12_solver_benchmark.csv'} and exp12_summary.json")


def _figures(gap, convergence, ps) -> None:
    fig, ax = new_figure(1, 2, figsize=(11, 4.4))
    s = gap["scatter"]
    ax[0].scatter(-s[:, 1], s[:, 0], s=4, alpha=0.25, color="#2b6cb0",
                  edgecolors="none")
    ax[0].axhline(gap["true_optimum"], color="#c53030", lw=1.2, ls="--",
                  label=f"true optimum ({gap['true_optimum']})")
    ax[0].axhline(gap["coverage_at_qubo_optimum"], color="#2f855a", lw=1.2,
                  ls=":", label=f"QUBO's choice ({gap['coverage_at_qubo_optimum']})")
    ax[0].set_xlabel("Bonferroni QUBO score (higher = better)")
    ax[0].set_ylabel("true unique UV cells")
    ax[0].set_title(f"Formulation gap, every feasible selection\n"
                    f"r = {gap['rank_correlation']:.3f}, "
                    f"gap = {gap['formulation_gap_percent']:.1f} %")
    ax[0].legend(fontsize=7, loc="lower right")

    for name, d in convergence.items():
        curve = np.array([(c[0], c[2]) for c in d["curve"]], dtype=float)
        ax[1].plot(curve[:, 0], curve[:, 1], "o-", lw=1.3, ms=4, label=name)
        if "heuristic" in d:
            ax[1].axhline(100.0 * d["heuristic"] / d["optimum"], lw=0.8,
                          ls=":", color="grey")
    ax[1].axhline(100.0, color="#c53030", lw=1.2, ls="--",
                  label="proven optimum")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("annealing reads")
    ax[1].set_ylabel("% of proven optimum")
    ax[1].set_title("Solver gap: QUBO annealing vs certified optimum\n"
                    "(dotted grey: direct classical heuristic)")
    ax[1].legend(fontsize=7, loc="lower right")
    stamp(fig, ps)
    save_figure(fig, "fig25_qubo_gap_decomposition")

    fig2, ax2 = new_figure(figsize=(6.2, 4.2))
    for name, d in convergence.items():
        curve = np.array([(c[3], c[2]) for c in d["curve"]], dtype=float)
        ax2.plot(curve[:, 0], curve[:, 1], "o-", lw=1.3, ms=4, label=name)
    ax2.axhline(100.0, color="#c53030", lw=1.2, ls="--", label="proven optimum")
    ax2.set_xscale("log")
    ax2.set_xlabel("sampler wall-clock time (s)")
    ax2.set_ylabel("% of proven optimum")
    ax2.set_title("Accuracy bought per second of annealing")
    ax2.legend(fontsize=7, loc="lower right")
    stamp(fig2, ps)
    save_figure(fig2, "fig26_qubo_solver_convergence")


if __name__ == "__main__":
    main()
