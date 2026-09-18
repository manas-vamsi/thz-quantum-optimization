"""Experiment 11 -- the selection problem on real ALMA pads.

Every other experiment here places candidate pads with a random generator. This
one uses the **actual 174 twelve-metre pads on the Chajnantor plateau**, as
distributed in the CASA observatory configuration files, at ALMA's real
latitude, with the measured phase-stability constants for that site. Geometry,
dish diameter, site coordinates, precipitable water vapour and phase structure
function are all measured quantities here.

That makes the comparison concrete in a way a synthetic pad field cannot be:
ALMA itself selects roughly 40-50 antennas from this pad list each cycle, which
is precisely the problem this project formulates.

Reported:

1. What the real pad field looks like, and whether the shadowing rule taken
   from the project note is consistent with it.
2. Searched selections against the analytic reference layouts, on real pads.
3. Certification: how large an instance can be proven optimal on real geometry.
4. Coherence at ALMA's own site latitude and a terahertz observing frequency, using the measured
   phase-stability constants for the site.

Outputs
-------
figures/fig23_real_alma_pads.png
figures/fig24_real_alma_selection.png
data/results/exp11_real_alma.csv
data/results/exp11_summary.json
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

from thz_opt.arrays.geometries import reuleaux_layout
from thz_opt.arrays.golden_spiral import golden_spiral_layout
from thz_opt.arrays.real_arrays import REAL_SITES, load_cfg
from thz_opt.constraints.atmosphere_data import (
    ALMA_PHASE,
    ALMA_PWV_SUMMARY,
    monthly_pwv,
)
from thz_opt.constraints.coherence import CoherenceConfig, coherence_summary
from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation, shadow_pairs
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_removal,
    local_search,
    max_coverage_milp,
    max_unique_cells,
    simulated_annealing,
)
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities
from thz_opt.qubo.coefficients import baseline_cell_sets, exact_cells_covered

N_SELECT = 20
GRID_CELLS = 32
FREQ_HZ = 3.0e11          # 300 GHz, ALMA band 7 territory
SEED = 20260918
CFG_FILE = ROOT / "data" / "external" / "alma.all.cfg"


def snap_to_pads(target, pads, forbidden=None):
    """Nearest available real pad for each point of an idealised layout.

    An analytic curve is not a feasible answer on a real site: the antennas have
    to stand on pads that exist. Snapping is how a designed layout is actually
    deployed, and it is the only way to compare an analytic construction against
    a searched selection on equal terms -- both must then be subsets of the same
    real pad list. Each pad is used at most once, and pads that would break the
    shadowing rule against an already-chosen pad are skipped.
    """
    chosen: list = []
    used = np.zeros(pads.shape[0], dtype=bool)
    for t in np.asarray(target, dtype=float):
        order = np.argsort(np.hypot(pads[:, 0] - t[0], pads[:, 1] - t[1]))
        for k in order:
            if used[k]:
                continue
            if forbidden is not None and chosen and np.any(forbidden[k, chosen]):
                continue
            used[k] = True
            chosen.append(int(k))
            break
    return np.asarray(chosen, dtype=int)


def score_layout(xy, cfg, grid, obs, sep_cfg):
    uv = layout_to_uv_tracks(xy, obs)
    occ, _ = uv_occupancy(uv, grid)
    psf, info = dirty_beam(uv, grid, cfg["psf"]["weighting"])
    m = psf_metrics(psf, info["pixel_scale_arcsec"])
    sep = check_minimum_separation(xy, sep_cfg)
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    off = d[~np.eye(xy.shape[0], dtype=bool)]
    return {
        "unique_cells": int((occ > 0).sum()),
        "psf_peak_sidelobe": m["peak_sidelobe_level"],
        "psf_fwhm_arcsec": m["fwhm_arcsec"],
        "max_baseline_m": float(off.max()),
        "min_separation_m": float(off.min()),
        "separation_violations": sep["n_violations"],
    }


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    site = REAL_SITES["alma"]

    if not CFG_FILE.exists():
        raise SystemExit(f"missing real pad file: {CFG_FILE}\n"
                         f"see data/external/SOURCES.md for provenance")

    pads_cfg = load_cfg(CFG_FILE, dish_diameter_m=site["dish_diameter_m"])
    pads = pads_cfg.xy
    M = len(pads_cfg)
    lam = 299_792_458.0 / FREQ_HZ

    obs = ObservationConfig(
        frequency_hz=FREQ_HZ,
        declination_deg=-30.0,
        latitude_deg=site["latitude_deg"],     # the real site latitude
        hour_angle_start_h=-2.0, hour_angle_end_h=2.0, n_times=41,
    )
    sep_cfg = SeparationConfig(site["dish_diameter_m"],
                               cfg["constraints"]["minimum_separation_factor"])

    ps = (f"{M} real ALMA 12 m pads, N={N_SELECT}, lat={site['latitude_deg']:.4f} deg, "
          f"nu={FREQ_HZ/1e9:.0f} GHz, UV grid {GRID_CELLS}^2, 4 h Earth rotation")
    banner(f"Experiment 11 -- real ALMA pads -- {ps}")

    # ---- 1. the real pad field -----------------------------------------
    s = pads_cfg.summary()
    print(f"\n  [1] the measured pad field ({pads_cfg.source})")
    print(f"      {M} pads, radius {s['radius_min_m']:.0f} to {s['radius_max_m']:.0f} m, "
          f"longest baseline {s['baseline_max_m']/1e3:.1f} km")
    print(f"      closest pad pair: {s['separation_min_m']:.1f} m")
    print(f"      shadowing rule d_min = {cfg['constraints']['minimum_separation_factor']} x "
          f"{site['dish_diameter_m']:.0f} m = {sep_cfg.d_min:.0f} m")

    forbidden = shadow_pairs(pads, sep_cfg)
    n_forbidden = int(np.triu(forbidden, 1).sum())
    print(f"      pad pairs the rule forbids: {n_forbidden} of {M*(M-1)//2} "
          f"({100*n_forbidden/(M*(M-1)//2):.2f} %)")
    if s["separation_min_m"] < sep_cfg.d_min:
        print(f"      note: the real pad field contains pairs closer than the rule allows,")
        print(f"      so the 1.5 x D figure is stricter than ALMA's own pad spacing.")

    # ---- 2. search against analytic references on real pads -------------
    print(f"\n  [2] selecting {N_SELECT} of {M} real pads")
    t0 = time.time()
    uv_all = layout_to_uv_tracks(pads, obs)
    grid = UVGrid(uv_max=cfg["uv_grid"]["pad"] * float(np.abs(uv_all).max()),
                  n_cells=GRID_CELLS)
    cell_sets = baseline_cell_sets(pads, lam, grid, observation=obs)
    cells, counts = build_pair_tables(track_cell_multiplicities(pads, lam, grid, obs))
    print(f"      precomputed {len(cell_sets)} candidate baselines in {time.time()-t0:.1f}s")

    state = UVState(cells, counts, M, grid.n_cells ** 2)
    score = max_unique_cells()
    rows = []

    # Analytic curves at the same radial extent. These are NOT feasible answers:
    # they place antennas wherever the formula says, not on pads that exist. They
    # are reported only to show the gap between an idealised curve and what the
    # site can actually support.
    r_min, r_max = s["radius_min_m"], s["radius_max_m"]
    ideal = {
        "golden spiral (ideal, infeasible)": golden_spiral_layout(
            N_SELECT, r_min=r_min, r_max=r_max, n_turns=3.0),
        "Reuleaux (ideal, infeasible)": reuleaux_layout(N_SELECT, r_min, r_max),
    }
    for name, xy in ideal.items():
        rows.append({"method": name, "kind": "ideal", "seconds": 0.0,
                     **score_layout(xy, cfg, grid, obs, sep_cfg)})

    # The same curves snapped to real pads: this is how a designed layout is
    # actually deployed, and it is the fair comparison for a searched selection.
    for name, xy in ideal.items():
        idx = snap_to_pads(xy, pads, forbidden)
        rows.append({"method": name.replace(" (ideal, infeasible)", " (snapped to pads)"),
                     "kind": "analytic", "seconds": 0.0,
                     "selection": [pads_cfg.names[i] for i in idx],
                     **score_layout(pads[idx], cfg, grid, obs, sep_cfg)})

    # Random feasible selections, as a null model on the same pad list.
    rng = np.random.default_rng(SEED)
    best_rand, rand_cells = None, []
    for _ in range(25):
        pick = rng.choice(M, N_SELECT, replace=False)
        if np.any(np.triu(forbidden[np.ix_(pick, pick)], 1)):
            continue
        r_ = score_layout(pads[pick], cfg, grid, obs, sep_cfg)
        rand_cells.append(r_["unique_cells"])
        if best_rand is None or r_["unique_cells"] > best_rand["unique_cells"]:
            best_rand = r_
    if best_rand is not None:
        rows.append({"method": f"best of {len(rand_cells)} random feasible",
                     "kind": "random", "seconds": 0.0, **best_rand})

    # The innermost pads: a real compact configuration, scored on the grid sized
    # for the extended array, which is why its coverage is so small.
    order = np.argsort(np.hypot(pads[:, 0], pads[:, 1]))
    rows.append({"method": "innermost pads (compact config)", "kind": "real subset",
                 "seconds": 0.0,
                 **score_layout(pads[order[:N_SELECT]], cfg, grid, obs, sep_cfg)})

    for label, fn in (
        ("greedy removal", lambda: greedy_removal(state, N_SELECT, score, forbidden)),
        ("annealing", lambda: simulated_annealing(state, N_SELECT, score, forbidden,
                                                  n_iterations=8000, seed=SEED)),
    ):
        t = time.time()
        res = fn()
        dt = time.time() - t
        rows.append({"method": label, "kind": "searched", "seconds": dt,
                     "selection": [pads_cfg.names[i] for i in res.indices()],
                     **score_layout(pads[res.indices()], cfg, grid, obs, sep_cfg)})

    state.set_selection(rows[-2].get("selection_idx", greedy_removal(
        state, N_SELECT, score, forbidden).indices()))
    t = time.time()
    res = local_search(state, score, forbidden)
    rows.append({"method": "greedy + local search", "kind": "searched",
                 "seconds": time.time() - t,
                 "selection": [pads_cfg.names[i] for i in res.indices()],
                 **score_layout(pads[res.indices()], cfg, grid, obs, sep_cfg)})
    best_real = res

    print(f"\n      {'method':<30}{'cells':>8}{'PSL':>8}{'b_max/km':>10}{'time s':>9}")
    print("      " + "-" * 65)
    for r in rows:
        print(f"      {r['method']:<30}{r['unique_cells']:>8}"
              f"{r['psf_peak_sidelobe']:>8.3f}{r['max_baseline_m']/1e3:>10.1f}"
              f"{r['seconds']:>9.1f}")

    analytic = max(r["unique_cells"] for r in rows if r["kind"] == "analytic")
    searched = max(r["unique_cells"] for r in rows if r["kind"] == "searched")
    ideal_best = max(r["unique_cells"] for r in rows if r["kind"] == "ideal")
    rand_best = max([r["unique_cells"] for r in rows if r["kind"] == "random"],
                    default=0)

    print(f"\n      Among feasible selections of real pads:")
    print(f"        best analytic layout snapped to pads : {analytic}")
    print(f"        best random feasible selection       : {rand_best}")
    print(f"        best searched selection              : {searched}"
          f"   ({100*(searched-analytic)/analytic:+.1f} % vs snapped analytic)")
    print(f"      An idealised curve ignoring the pad list would score {ideal_best},")
    print(f"      but it is not a layout this site can build: it places antennas")
    print(f"      where no pad exists.")

    # ---- 3. certification on real geometry ------------------------------
    print("\n  [3] how much of the real problem can be proven optimal")
    cert = []
    for m_sub, n_sub in ((24, 8), (40, 10)):
        sub = order[:m_sub]                      # a contiguous real sub-array
        sub_pads = pads[sub]
        g = UVGrid(uv_max=1.02 * float(np.abs(layout_to_uv_tracks(sub_pads, obs)).max()),
                   n_cells=16)
        cs = baseline_cell_sets(sub_pads, lam, g, observation=obs)
        t = time.time()
        ex = max_coverage_milp(cs, m_sub, n_sub,
                               forbidden_pairs=shadow_pairs(sub_pads, sep_cfg),
                               time_limit=180.0)
        dt = time.time() - t
        cert.append({"m_pads": m_sub, "n_select": n_sub, "incumbent": ex["coverage"],
                     "upper_bound": ex.get("upper_bound"),
                     "proven_optimal": ex["proven_optimal"], "seconds": dt})
        print(f"      M={m_sub}, N={n_sub}: incumbent {ex['coverage']}, "
              f"bound {ex.get('upper_bound')}, proven={ex['proven_optimal']}, {dt:.0f}s")

    # ---- 4. coherence at the real site ----------------------------------
    print(f"\n  [4] atmospheric coherence at {FREQ_HZ/1e9:.0f} GHz, measured constants")
    print(f"      {'observing condition':<22}{'sigma(1km)':>12}{'mean gamma':>12}"
          f"{'below 0.5':>11}")
    coh = {}
    for cond in ("uncorrected", "uncorrected_low_pwv", "wvr_corrected",
                 "wvr_corrected_windy"):
        cc = CoherenceConfig.from_measured(cond)
        c = coherence_summary(pads[best_real.indices()], lam, cc)
        coh[cond] = c
        print(f"      {cond:<22}{cc.sigma_1km_m*1e6:>9.0f} um{c['coherence_mean']:>12.3f}"
              f"{100*c['fraction_below_0p5']:>10.0f} %")
    print(f"      source: {ALMA_PHASE['wvr_corrected']['source']}")

    pwv = {m: monthly_pwv(m) for m in ("Jan", "Jun", "Aug", "Dec")}
    print("\n      measured PWV medians (mm): " +
          ", ".join(f"{m} {v}" for m, v in pwv.items()) +
          f"  [{ALMA_PWV_SUMMARY['span_years']}-year year-round median "
          f"{ALMA_PWV_SUMMARY['median_year_round_mm']} mm]")

    save_csv([{k: v for k, v in r.items() if k != "selection"} for r in rows],
             RESULTS / "exp11_real_alma.csv")
    save_json({"site": site, "pad_file": pads_cfg.source, "pad_summary": s,
               "n_pads": M, "n_select": N_SELECT, "frequency_hz": FREQ_HZ,
               "observation": obs.as_dict(), "forbidden_pairs": n_forbidden,
               "rows": rows, "certification": cert,
               "coherence": coh,
               "pwv_measured_mm": {m: monthly_pwv(m) for m in
                                   ("Jan", "Jun", "Aug", "Dec")},
               "selected_pads": [pads_cfg.names[i] for i in best_real.indices()]},
              RESULTS / "exp11_summary.json")

    _figures(pads, pads_cfg, best_real.indices(), grid, obs, cfg, ps, sep_cfg)
    print(f"\nWrote {RESULTS / 'exp11_real_alma.csv'} and exp11_summary.json")


def _figures(pads, pads_cfg, selection, grid, obs, cfg, ps, sep_cfg) -> None:
    fig, axes = new_figure(1, 2, figsize=(11.5, 5.2))
    r = np.hypot(pads[:, 0], pads[:, 1])
    axes[0].scatter(pads[:, 0] / 1e3, pads[:, 1] / 1e3, s=9, c="#1f4e79")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("East (km)")
    axes[0].set_ylabel("North (km)")
    axes[0].set_title(f"Figure 23a. All {len(pads)} real ALMA 12 m pads")

    inner = r < 500
    axes[1].scatter(pads[inner, 0], pads[inner, 1], s=16, c="#1f4e79")
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("East (m)")
    axes[1].set_ylabel("North (m)")
    axes[1].set_title("Figure 23b. Inner 500 m of the pad field")
    stamp(fig, ps)
    save_figure(fig, "fig23_real_alma_pads.png")

    sel = np.asarray(selection, dtype=int)
    uv = layout_to_uv_tracks(pads[sel], obs)
    occ, _ = uv_occupancy(uv, grid)
    fig, axes = new_figure(1, 2, figsize=(11.5, 5.0))
    axes[0].scatter(pads[:, 0] / 1e3, pads[:, 1] / 1e3, s=8, c="0.78",
                    label="available pads")
    axes[0].scatter(pads[sel, 0] / 1e3, pads[sel, 1] / 1e3, s=34, c="#C44E52",
                    edgecolor="k", linewidth=0.3, label="selected")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("East (km)")
    axes[0].set_ylabel("North (km)")
    axes[0].set_title("Figure 24a. Selection on real ALMA pads")
    axes[0].legend(fontsize=7)

    axes[1].imshow((occ.T > 0).astype(float), origin="lower", cmap="Greys",
                   extent=[-grid.uv_max / 1e6, grid.uv_max / 1e6,
                           -grid.uv_max / 1e6, grid.uv_max / 1e6])
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("u (mega-wavelengths)")
    axes[1].set_ylabel("v (mega-wavelengths)")
    axes[1].set_title(f"Figure 24b. Resulting UV coverage ({int((occ>0).sum())} cells)")
    stamp(fig, ps)
    save_figure(fig, "fig24_real_alma_selection.png")


if __name__ == "__main__":
    main()
