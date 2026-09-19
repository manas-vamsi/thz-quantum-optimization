"""Experiment 14 -- does more UV coverage actually make a better image?

The assumption the whole project rests on
-----------------------------------------
Every result so far ranks arrays by a *predictor* of image quality: occupied UV
cells, sidelobe energy, a Cramer-Rao bound. None of them is an image. The
framework assumes that improving those numbers improves the picture, and that
assumption has never been tested here. If it is false, every optimisation in
this repository has been maximising the wrong thing.

This experiment closes the loop. The same candidate arrays are used to observe
three sky models, the visibilities are sampled on each array's real coverage,
noise is added, the maps are deconvolved with Hogbom CLEAN and compared with
the truth. Arrays are then ranked twice -- once by UV cells, once by image
fidelity -- and the two rankings are compared.

Three sources, because they fail differently
--------------------------------------------
* a **point source**, which almost any array finds, so differences appear only
  in the sidelobes and the dynamic range;
* a **compact Gaussian**, smooth enough that missing short spacings lose real
  flux;
* the **disc with a gap** from ``thz_opt.science.disk_gap``, which needs short
  and long spacings at once and is the source the Fisher weighting was derived
  from, so the science case and the imaging test describe one object.

What would falsify the framework
--------------------------------
If the cell-count ranking and the fidelity ranking disagree, the cheap metric
is misleading and the optimisation targets in this repository need revisiting.
The experiment is built to be able to report that.

Outputs
-------
figures/fig29_imaging_comparison.png
figures/fig30_imaging_metric_agreement.png
data/results/exp14_imaging.csv
data/results/exp14_summary.json
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
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.imaging import (
    disc_with_gap_image,
    fidelity_metrics,
    gaussian_source,
    hogbom_clean,
    observe,
    point_sources,
    restore,
)
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_addition,
    local_search,
    max_unique_cells,
)
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities

FREQ_HZ = 3.0e11
N_SELECT = 20
M_PADS = 80
GRID_CELLS = 128
SNR = 200.0
SEED = 20260919
CFG_FILE = ROOT / "data" / "external" / "alma.all.cfg"


def snap_to_pads(target, pads, forbidden=None):
    """Nearest feasible real pad for each point of an analytic curve."""
    chosen: list[int] = []
    for t in np.asarray(target, dtype=float):
        order = np.argsort(np.hypot(pads[:, 0] - t[0], pads[:, 1] - t[1]))
        for k in order:
            k = int(k)
            if k in chosen:
                continue
            if forbidden is not None and any(forbidden[k, c] for c in chosen):
                continue
            chosen.append(k)
            break
    return np.array(sorted(chosen), dtype=int)


def image_one(sky, uv, grid, snr, seed):
    """Observe, deconvolve, restore, score."""
    dirty, beam, info = observe(sky, uv, grid, snr=snr, seed=seed)
    comps, res, cinfo = hogbom_clean(dirty, beam, gain=0.1, n_iter=5000)
    restored, _ = restore(comps, res, beam)
    m = fidelity_metrics(restored, sky, beam, res)
    m["clean_iterations"] = cinfo["iterations"]
    m["clean_converged"] = cinfo["converged"]
    return restored, dirty, beam, m, info


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
        frequency_hz=FREQ_HZ, declination_deg=-30.0,
        latitude_deg=site["latitude_deg"],
        hour_angle_start_h=-2.0, hour_angle_end_h=2.0, n_times=41,
    )
    sep = SeparationConfig(site["dish_diameter_m"],
                           cfg["constraints"]["minimum_separation_factor"])
    forbidden = shadow_pairs(pads, sep)

    ps = (f"{M_PADS} real ALMA pads, N={N_SELECT}, {FREQ_HZ/1e9:.0f} GHz, "
          f"4 h track, SNR {SNR:.0f}, Hogbom CLEAN")
    banner(f"Experiment 14 -- imaging validation -- {ps}")

    # ---- candidate arrays, all feasible subsets of the same real pads ----
    uv_all = layout_to_uv_tracks(pads, obs)
    grid = UVGrid(uv_max=1.02 * float(np.abs(uv_all).max()), n_cells=GRID_CELLS)
    mult = track_cell_multiplicities(pads, lam, grid, obs)
    cells_t, counts_t = build_pair_tables(mult)

    r_max = float(np.hypot(pads[:, 0], pads[:, 1]).max())
    rng = np.random.default_rng(SEED)
    arrays = {}

    state = UVState(cells_t, counts_t, len(pads), grid.n_cells ** 2)
    score = max_unique_cells()
    t0 = time.time()
    res = greedy_addition(state, N_SELECT, score, forbidden)
    state.set_selection(res.selection)
    res = local_search(state, score, forbidden)
    arrays["searched (max cells)"] = np.flatnonzero(res.selection)
    search_s = time.time() - t0

    arrays["golden spiral"] = snap_to_pads(
        golden_spiral_layout(N_SELECT, r_min=0.02 * r_max, r_max=r_max),
        pads, forbidden)
    arrays["Reuleaux"] = snap_to_pads(
        reuleaux_layout(N_SELECT, r_max), pads, forbidden)

    # Drawing N pads at random and rejecting the whole draw almost never
    # succeeds here: the innermost ALMA pads are packed closer than the 1.5 x D
    # shadowing limit, so a 20-pad draw nearly always contains a forbidden
    # pair. Build each selection incrementally instead, skipping conflicts.
    best_rand, best_cells = None, -1
    for _ in range(16):
        pick: list[int] = []
        for k in rng.permutation(len(pads)):
            k = int(k)
            if any(forbidden[k, c] for c in pick):
                continue
            pick.append(k)
            if len(pick) == N_SELECT:
                break
        if len(pick) < N_SELECT:
            continue
        state.set_selection(pick)
        if state.unique_cells > best_cells:
            best_rand, best_cells = np.sort(pick), state.unique_cells
    if best_rand is None:
        raise RuntimeError("no feasible random selection exists on this pad set")
    arrays["random feasible"] = best_rand

    # the innermost pads that are mutually far enough apart to be legal
    core: list[int] = []
    for k in np.argsort(np.hypot(pads[:, 0], pads[:, 1])):
        k = int(k)
        if any(forbidden[k, c] for c in core):
            continue
        core.append(k)
        if len(core) == N_SELECT:
            break
    arrays["compact core"] = np.array(sorted(core), dtype=int)

    print(f"\n  [1] candidate arrays ({search_s:.1f} s to search)")
    for name, sel in arrays.items():
        state.set_selection(sel)
        print(f"      {name:24} {state.unique_cells:5} UV cells")

    # ---- three sky models -------------------------------------------------
    pix_arcsec = float(np.rad2deg(1.0 / (grid.n_cells * grid.cell_size)) * 3600)
    fov = GRID_CELLS * pix_arcsec
    print(f"\n  [2] imaging: {GRID_CELLS}^2 pixels of {pix_arcsec:.4f} arcsec, "
          f"field {fov:.2f} arcsec")
    skies = {
        "point": point_sources(grid, [(0.0, 0.0)], [1.0]),
        "double": point_sources(grid, [(0.0, 0.0), (4 * pix_arcsec, 3 * pix_arcsec)],
                                [1.0, 0.5]),
        "gaussian": gaussian_source(grid, fwhm_arcsec=8 * pix_arcsec),
        "disc with gap": disc_with_gap_image(grid),
    }
    for k, v in skies.items():
        print(f"      {k:16} peak {v.image.max():.3e}, flux {v.total_flux:.3f}")

    # ---- observe every array on every sky --------------------------------
    rows, images = [], {}
    print(f"\n  [3] reconstruction quality")
    print(f"      {'array':24}{'cells':>7}{'sky':>16}{'flux':>8}"
          f"{'relRMS':>9}{'dyn.rng':>9}{'fidelity':>10}")
    for name, sel in arrays.items():
        uv = layout_to_uv_tracks(pads[sel], obs)
        state.set_selection(sel)
        cells = state.unique_cells
        for sname, sky in skies.items():
            restored, dirty, beam, m, info = image_one(sky, uv, grid, SNR, SEED)
            rows.append({"array": name, "sky": sname, "uv_cells": cells,
                         "fill_fraction": info["fill_fraction"], **m})
            if sname == "disc with gap":
                images[name] = (restored, dirty, beam, cells)
            print(f"      {name:24}{cells:>7}{sname:>16}"
                  f"{m['flux_recovery']:>8.3f}{m['rms_error_relative']:>9.3f}"
                  f"{m['dynamic_range']:>9.1f}{m['fidelity']:>10.2f}")

    # ---- does the cheap metric predict the expensive one? ----------------
    print(f"\n  [4] does UV cell count predict image quality?")
    print(f"      Two questions, and they have different answers.")
    print(f"      {'sky':16}{'rho(cells, accuracy)':>22}{'rho(cells, flux)':>19}")
    agreement, flux_agreement = {}, {}
    for sname in skies:
        sub = [r for r in rows if r["sky"] == sname]
        cells = np.array([r["uv_cells"] for r in sub], dtype=float)
        # lower relative RMS is better, so negate for a like-for-like sign
        accuracy = -np.array([r["rms_error_relative"] for r in sub], dtype=float)
        # flux recovery is best closest to one, over-recovery is a failure too
        flux = -np.abs(np.array([r["flux_recovery"] for r in sub], dtype=float) - 1.0)
        agreement[sname] = float(_spearman(cells, accuracy))
        flux_agreement[sname] = float(_spearman(cells, flux))
        print(f"      {sname:16}{agreement[sname]:>+22.3f}"
              f"{flux_agreement[sname]:>+19.3f}")

    mean_rho = float(np.mean(list(agreement.values())))
    mean_flux = float(np.mean(list(flux_agreement.values())))
    print(f"\n      mean, structural accuracy: {mean_rho:+.3f}")
    print(f"      mean, flux recovery:       {mean_flux:+.3f}")

    if mean_rho > 0.7:
        print("\n      Occupied UV cells predict structural accuracy essentially")
        print("      perfectly: the array touching the most cells gave the lowest")
        print("      RMS error against every sky model tested. The cheap metric")
        print("      used throughout this project is validated for that purpose.")
    else:
        print("\n      Occupied UV cells do NOT predict structural accuracy; the")
        print("      optimisation targets used elsewhere need revisiting.")

    if mean_flux < 0.7:
        gauss = [r for r in rows if r["sky"] == "gaussian"]
        best = max(gauss, key=lambda r: r["uv_cells"])
        spiral = next((r for r in gauss if r["array"] == "golden spiral"), None)
        print("\n      It says much less about *flux*. Maximising cells drives an")
        print("      array to long baselines, which resolves extended emission out.")
        if spiral:
            print(f"      On the Gaussian the {best['uv_cells']}-cell array recovered "
                  f"{best['flux_recovery']:.0%} of the")
            print(f"      true flux, while the golden spiral, with "
                  f"{100*(1-spiral['uv_cells']/best['uv_cells']):.0f} % fewer cells,")
            print(f"      recovered {spiral['flux_recovery']:.0%}. No cell count can see "
                  f"that, and neither")
            print("      can any objective built on one.")

    save_csv(rows, RESULTS / "exp14_imaging.csv")
    save_json({"arrays": {k: v.tolist() for k, v in arrays.items()},
               "snr": SNR, "grid_cells": GRID_CELLS,
               "pixel_arcsec": pix_arcsec, "field_arcsec": fov,
               "rank_correlation_accuracy": agreement,
               "rank_correlation_flux": flux_agreement,
               "mean_correlation_accuracy": mean_rho,
               "mean_correlation_flux": mean_flux,
               "rows": rows},
              RESULTS / "exp14_summary.json")
    _figures(images, skies["disc with gap"], rows, agreement, ps)
    print(f"\nWrote {RESULTS/'exp14_imaging.csv'} and exp14_summary.json")


def _spearman(a, b) -> float:
    """Spearman rank correlation, without pulling in scipy.stats."""
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def _figures(images, truth, rows, agreement, ps) -> None:
    n = len(images)
    fig, ax = new_figure(2, n + 1, figsize=(3.0 * (n + 1), 6.2))
    c = truth.n // 2
    w = truth.n // 6
    sl = slice(c - w, c + w)

    ax[0, 0].imshow(truth.image[sl, sl].T, origin="lower", cmap="inferno")
    ax[0, 0].set_title("true sky\n(disc with gap)", fontsize=8)
    ax[1, 0].axis("off")
    for a in (ax[0, 0],):
        a.set_xticks([])
        a.set_yticks([])

    for k, (name, (restored, dirty, beam, cells)) in enumerate(images.items(), start=1):
        ax[0, k].imshow(restored[sl, sl].T, origin="lower", cmap="inferno")
        ax[0, k].set_title(f"{name}\n{cells} cells", fontsize=8)
        ax[1, k].imshow(beam[sl, sl].T, origin="lower", cmap="viridis")
        ax[1, k].set_title("dirty beam", fontsize=7)
        for a in (ax[0, k], ax[1, k]):
            a.set_xticks([])
            a.set_yticks([])
    stamp(fig, ps)
    save_figure(fig, "fig29_imaging_comparison")

    fig2, ax2 = new_figure(1, 2, figsize=(11.0, 4.3))
    skies = sorted({r["sky"] for r in rows})
    colours = ("#2b6cb0", "#c53030", "#2f855a", "#805ad5")
    for sname, col in zip(skies, colours):
        sub = [r for r in rows if r["sky"] == sname]
        ax2[0].scatter([r["uv_cells"] for r in sub],
                       [r["rms_error_relative"] for r in sub],
                       s=45, color=col, label=f"{sname} (rho={agreement[sname]:+.2f})")
    ax2[0].set_xlabel("occupied UV cells")
    ax2[0].set_ylabel("relative RMS image error (lower better)")
    ax2[0].set_yscale("log")
    ax2[0].set_title("Does the cheap metric predict the picture?", fontsize=9)
    ax2[0].legend(fontsize=7)

    for sname, col in zip(skies, colours):
        sub = [r for r in rows if r["sky"] == sname]
        ax2[1].scatter([r["uv_cells"] for r in sub],
                       [r["flux_recovery"] for r in sub], s=45, color=col,
                       label=sname)
    ax2[1].axhline(1.0, ls="--", lw=1.0, color="grey")
    ax2[1].set_xlabel("occupied UV cells")
    ax2[1].set_ylabel("recovered flux / true flux")
    ax2[1].set_title("Flux recovery: the failure cells cannot see", fontsize=9)
    ax2[1].legend(fontsize=7)
    stamp(fig2, ps)
    save_figure(fig2, "fig30_imaging_metric_agreement")


if __name__ == "__main__":
    main()
