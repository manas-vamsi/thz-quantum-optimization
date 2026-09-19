"""Experiment 13 -- an array optimised for a measurement, not for a metric.

Every other experiment in this project maximises occupied UV cells, or occupied
cells times a radial weight that was chosen by hand. Nobody builds a telescope
to occupy UV cells. This experiment replaces the invented figure of merit with
one derived from a science case, and reports the result in a physical unit: the
smallest error bar the array can place on a property of a real source.

The chain, with nothing chosen in the middle
--------------------------------------------
source model (HL Tau, measured)
    -> visibility V(q), by Hankel transform
    -> dV/dtheta for the parameter being measured
    -> Fisher information F = sum_sampled |dV/dtheta|^2      [the UV weight]
    -> Cramer-Rao bound sigma(theta) >= 1/sqrt(F)            [the answer]

The per-cell weight is ``|dV/dtheta|^2``, which is exactly the shape the
existing optimiser consumes, so no solver changes are needed. See
``thz_opt.science.disk_gap`` for the derivation.

The question this settles
-------------------------
Does it matter? If the cell-count-optimal array happens to also be the
Fisher-optimal array, then the hand-chosen metric was fine all along and this
module is decoration. The experiment is built to be able to return that answer.

Three parameters are compared -- gap depth, radius and width -- because they
are informed by different spatial frequencies, so "optimise for HL Tau" is not
yet a well-posed instruction. Quantifying that is part of the result.

Outputs
-------
figures/fig27_science_weight_derivation.png
figures/fig28_science_vs_cellcount.png
data/results/exp13_science_case.csv
data/results/exp13_summary.json
"""

from __future__ import annotations

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

from thz_opt.arrays.real_arrays import REAL_SITES, load_cfg
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_removal,
    local_search,
    max_unique_cells,
)
from thz_opt.interferometry.baselines import pair_indices
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities
from thz_opt.science.disk_gap import (
    HL_TAU,
    PARAMETERS,
    brightness_profile,
    depth_error_bar,
    fisher_weight_profile,
    gap_derivative_profile,
    resolution_requirement,
    uv_cell_weights,
    visibility_profile,
)

FREQ_HZ = 3.0e11
N_SELECT = 20
M_PADS = 80
GRID_CELLS = 48
SEED = 20260919
CFG_FILE = ROOT / "data" / "external" / "alma.all.cfg"


def weighted_cells_objective():
    """Maximise Fisher weight over *distinct* occupied cells.

    This is the imaging-flavoured reading of the science case: reward reaching
    informative parts of the UV plane, but count each cell once.
    """
    def score(state) -> float:
        return -float(state.weighted_cells)
    score.name = "fisher_weighted_cells"
    return score


def pair_fisher_weights(mult, cell_w):
    """``W_k = sum_c w_c * a_ck`` -- the Fisher information one baseline carries.

    Equation (1) sums ``|dV/dtheta|^2`` over *samples*, not over distinct
    cells, because with independent noise a repeated visibility measurement is
    a genuine second measurement. That makes the Fisher objective

        F(x) = sum_{i<j} W_ij x_i x_j

    exactly linear in the baseline-activation variables, so unlike coverage it
    needs no Rosenberg quadratization at all: it is already a plain QUBO over
    the pad variables. That is a real structural advantage of scoring an array
    by what it measures rather than by how many cells it touches.
    """
    return np.array([sum(cell_w[c] * a for c, a in d.items()) if d else 0.0
                     for d in mult], dtype=float)


def true_fisher_objective(pair_w, n_pads):
    """Maximise ``sum_{i<j} W_ij x_i x_j``, redundancy included."""
    i_idx, j_idx = pair_indices(n_pads)

    def score(state) -> float:
        sel = state.selected
        return -float(pair_w[sel[i_idx] & sel[j_idx]].sum())
    score.name = "fisher_information"
    return score


def redundancy(state) -> float:
    """Mean samples per occupied cell -- how much the array repeats itself."""
    n = state._n
    occ = n[n != 0]
    return float(occ.mean()) if occ.size else 0.0


def sampled_q(pads, selection, obs, grid):
    """Spatial frequencies actually sampled by a selection, one per UV point."""
    sel = np.asarray(selection, dtype=bool)
    uv = layout_to_uv_tracks(pads[sel], obs)
    return np.hypot(uv[..., 0], uv[..., 1]).ravel()


def optimise(state, n_select, forbidden, score, seed):
    state.set_selection([])
    t0 = time.time()
    res = greedy_removal(state, n_select, score, forbidden)
    state.set_selection(res.selection)
    res = local_search(state, score, forbidden)
    return res, time.time() - t0


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    site = REAL_SITES["alma"]
    if not CFG_FILE.exists():
        raise SystemExit(f"missing real pad file: {CFG_FILE}")

    pads_cfg = load_cfg(CFG_FILE, dish_diameter_m=site["dish_diameter_m"])
    all_pads = pads_cfg.xy
    order = np.argsort(np.hypot(all_pads[:, 0], all_pads[:, 1]))
    pads = all_pads[order[:M_PADS]]
    lam = 299_792_458.0 / FREQ_HZ
    obs = ObservationConfig(
        frequency_hz=FREQ_HZ,
        declination_deg=-30.0,
        latitude_deg=site["latitude_deg"],
        hour_angle_start_h=-2.0, hour_angle_end_h=2.0, n_times=41,
    )
    sep_cfg = SeparationConfig(site["dish_diameter_m"],
                               cfg["constraints"]["minimum_separation_factor"])
    forbidden = shadow_pairs(pads, sep_cfg)

    ps = (f"{M_PADS} real ALMA pads, N={N_SELECT}, nu={FREQ_HZ/1e9:.0f} GHz, "
          f"HL Tau gap model, Fisher-derived UV weight")
    banner(f"Experiment 13 -- science-derived objective -- {ps}")

    # ---- 1. what the science demands, before any optimisation -----------
    req = resolution_requirement(HL_TAU, lam)
    print(f"\n  [1] the science case: {HL_TAU.name}")
    print(f"      distance {HL_TAU.distance_pc:.0f} pc, disc radius "
          f"{HL_TAU.r_out_au:.0f} au, gap at {HL_TAU.gap_radius_au:.1f} au "
          f"(width {HL_TAU.gap_width_au:.1f} au)")
    print(f"      gap subtends {req['theta_gap_arcsec']*1000:.1f} mas, disc "
          f"{req['theta_disc_arcsec']:.2f} arcsec")
    print(f"      to resolve the gap:   b_max >= {req['b_max_required_m']:.0f} m")
    print(f"      to keep the disc:     b_min <= {req['b_min_allowed_m']:.0f} m")

    # ---- 2. the derived weight, per parameter ---------------------------
    q_scan = np.logspace(3.5, 7.0, 900)
    print("\n  [2] Fisher weight derived per gap parameter (flux-conserving)")
    print(f"      {'parameter':<10}{'peak baseline':>15}{'half-max band':>22}")
    bands = {}
    for par in PARAMETERS:
        w = fisher_weight_profile(q_scan, parameter=par)
        hi = q_scan[w >= 0.5]
        bands[par] = {
            "peak_baseline_m": float(q_scan[int(np.argmax(w))] * lam),
            "band_lo_m": float(hi.min() * lam), "band_hi_m": float(hi.max() * lam),
        }
        print(f"      {par:<10}{bands[par]['peak_baseline_m']:>12.0f} m"
              f"{bands[par]['band_lo_m']:>14.0f} -{bands[par]['band_hi_m']:>6.0f} m")
    print("      The three differ, so 'optimise for HL Tau' is under-specified")
    print("      until the parameter of interest is named.")

    # ---- 3. optimise for cells, then for each parameter -----------------
    uv_all = layout_to_uv_tracks(pads, obs)
    grid = UVGrid(uv_max=1.02 * float(np.abs(uv_all).max()), n_cells=GRID_CELLS)
    mult = track_cell_multiplicities(pads, lam, grid, obs)
    cells_t, counts_t = build_pair_tables(mult)

    print(f"\n  [3] optimising {N_SELECT} of {M_PADS} real pads")
    print(f"      {'objective':<20}{'cells':>8}{'redun':>7}{'sig(depth)':>13}"
          f"{'sig(radius)':>13}{'sig(width)':>13}{'s':>7}")

    runs = {}
    plain = UVState(cells_t, counts_t, len(pads), grid.n_cells ** 2)
    res, dt = optimise(plain, N_SELECT, forbidden, max_unique_cells(), SEED)
    runs["max cells"] = (res.selection, dt)

    for par in PARAMETERS:
        w = uv_cell_weights(grid, HL_TAU, parameter=par)
        st = UVState(cells_t, counts_t, len(pads), grid.n_cells ** 2,
                     cell_weights=w)
        res, dt = optimise(st, N_SELECT, forbidden, weighted_cells_objective(), SEED)
        runs[f"w-cells: {par}"] = (res.selection, dt)

        pw = pair_fisher_weights(mult, w)
        st2 = UVState(cells_t, counts_t, len(pads), grid.n_cells ** 2,
                      cell_weights=w)
        res, dt = optimise(st2, N_SELECT, forbidden,
                           true_fisher_objective(pw, len(pads)), SEED)
        runs[f"Fisher: {par}"] = (res.selection, dt)

    rows = []
    for name, (sel, dt) in runs.items():
        q = sampled_q(pads, sel, obs, grid)
        plain.set_selection(np.flatnonzero(sel))
        errs = {p: depth_error_bar(q, HL_TAU, 1.0, p) for p in PARAMETERS}
        rows.append({"objective": name, "cells": plain.unique_cells,
                     "redundancy": redundancy(plain),
                     "sigma_depth": errs["depth"], "sigma_radius": errs["radius"],
                     "sigma_width": errs["width"], "seconds": dt,
                     "b_max_m": float(np.hypot(*(pads[sel] - pads[sel].mean(0)).T).max() * 2)})
        print(f"      {name:<20}{plain.unique_cells:>8}{redundancy(plain):>7.1f}"
              f"{errs['depth']:>13.3e}{errs['radius']:>13.3e}"
              f"{errs['width']:>13.3e}{dt:>7.1f}")

    # ---- 4. does the science-derived objective actually buy anything? ---
    print("\n  [4] does deriving the weight change the answer?")
    base = next(r for r in rows if r["objective"] == "max cells")
    verdict = {}
    for par in PARAMETERS:
        tuned = next(r for r in rows if r["objective"] == f"Fisher: {par}")
        key = f"sigma_{par}"
        gain = 100.0 * (base[key] - tuned[key]) / base[key]
        cell_cost = 100.0 * (base["cells"] - tuned["cells"]) / base["cells"]
        verdict[par] = {"error_bar_improvement_percent": gain,
                        "cell_count_cost_percent": cell_cost,
                        "redundancy_cellcount": base["redundancy"],
                        "redundancy_fisher": tuned["redundancy"],
                        "sigma_cellcount_array": base[key],
                        "sigma_fisher_array": tuned[key]}
        print(f"      {par:<8}: error bar {gain:+6.1f} % vs the cell-count array, "
              f"costing {cell_cost:+5.1f} % of cells, "
              f"redundancy {base['redundancy']:.1f} -> {tuned['redundancy']:.1f}")

    save_csv(rows, RESULTS / "exp13_science_case.csv")
    save_json({"source": HL_TAU.name, "model": HL_TAU.__dict__,
               "requirement": req, "weight_bands_m": bands,
               "runs": rows, "verdict": verdict,
               "wavelength_m": lam, "n_select": N_SELECT, "m_pads": M_PADS},
              RESULTS / "exp13_summary.json")
    _figures(q_scan, lam, rows, verdict, pads, runs, ps)
    print(f"\nWrote {RESULTS/'exp13_science_case.csv'} and exp13_summary.json")


def _figures(q_scan, lam, rows, verdict, pads, runs, ps) -> None:
    fig, ax = new_figure(1, 3, figsize=(13.5, 4.0))
    r = np.linspace(HL_TAU.r_in_au, HL_TAU.r_out_au, 2000)
    ax[0].plot(r, brightness_profile(r), color="#2b6cb0", lw=1.5, label="I(r)")
    ax[0].plot(r, gap_derivative_profile(r, HL_TAU, True, "depth"),
               color="#c53030", lw=1.2, label="dI/d(depth), flux conserving")
    ax[0].axvline(HL_TAU.gap_radius_au, color="grey", ls=":", lw=0.9)
    ax[0].set_xlabel("radius (au)")
    ax[0].set_ylabel("normalised brightness")
    ax[0].set_title(f"Source model\n{HL_TAU.name}")
    ax[0].legend(fontsize=7)

    for par, c in zip(PARAMETERS, ("#2b6cb0", "#c53030", "#2f855a")):
        ax[1].semilogx(q_scan * lam, fisher_weight_profile(q_scan, parameter=par),
                       lw=1.4, color=c, label=f"d/d({par})")
    ax[1].set_xlabel("baseline length (m)")
    ax[1].set_ylabel(r"$|dV/d\theta|^2$  (normalised)")
    ax[1].set_title("Derived UV weight\n(not chosen: Fisher information)")
    ax[1].legend(fontsize=7)

    names = [r["objective"] for r in rows]
    x = np.arange(len(names))
    width = 0.27
    for k, (par, c) in enumerate(zip(PARAMETERS, ("#2b6cb0", "#c53030", "#2f855a"))):
        vals = np.array([r[f"sigma_{par}"] for r in rows], dtype=float)
        ax[2].bar(x + (k - 1) * width, vals / vals.max(), width, color=c,
                  label=f"sigma({par})")
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(names, rotation=20, ha="right", fontsize=7)
    ax[2].set_ylabel("error bar, relative to worst")
    ax[2].set_title("Cramer-Rao bound by objective\n(lower is better)")
    ax[2].legend(fontsize=7)
    stamp(fig, ps)
    save_figure(fig, "fig27_science_weight_derivation")

    fig2, ax2 = new_figure(1, len(runs), figsize=(3.2 * len(runs), 3.4))
    ax2 = np.atleast_1d(ax2)
    for a, (name, (sel, _)) in zip(ax2, runs.items()):
        a.scatter(pads[:, 0] / 1e3, pads[:, 1] / 1e3, s=5, color="#cbd5e0")
        a.scatter(pads[sel, 0] / 1e3, pads[sel, 1] / 1e3, s=18, color="#c53030")
        a.set_aspect("equal")
        a.set_title(name, fontsize=8)
        a.set_xlabel("east (km)")
    ax2[0].set_ylabel("north (km)")
    stamp(fig2, ps)
    save_figure(fig2, "fig28_science_vs_cellcount")


if __name__ == "__main__":
    main()
