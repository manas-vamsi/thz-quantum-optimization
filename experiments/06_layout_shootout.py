"""Experiment 06 -- the golden spiral against the geometries the field uses.

Experiment 01 asked "golden spiral or Fibonacci?".  That is the wrong pair of
options: both are closed-form curves picked for their mathematics rather than
for an imaging objective.  This experiment puts them alongside the geometries
the radio-array literature actually converged on, and scores everything on
three objectives at once so it is visible that the objectives disagree:

* **unique UV cells** -- the cell-counting metric used in experiments 01-03
* **Cornwell F** -- UV-point repulsion energy (Cornwell 1988), lower = more
  uniformly spread UV samples
* **Gaussian density chi2** -- mismatch to a Gaussian radial UV density
  (Boone 2002), lower = beam closer to a sidelobe-free Gaussian

plus the PSF diagnostics and the minimum-separation feasibility check.

It also runs the experiment that the golden-ratio literature never runs: a
sweep over the phyllotaxis angle ``alpha``, including the golden angle, other
badly approximable irrationals, and rationals.  If the golden ratio is special
for array design, the golden angle should stand out.

Outputs
-------
figures/fig13_layout_gallery.png
figures/fig14_objective_disagreement.png
figures/fig15_alpha_sweep.png
data/results/exp06_shootout.csv
data/results/exp06_alpha_sweep.csv
"""

from __future__ import annotations

import numpy as np

from common import (  # noqa: E402
    RESULTS,
    banner,
    ensure_dirs,
    load_config,
    merge,
    new_figure,
    observation_config,
    param_string,
    save_csv,
    save_figure,
    stamp,
    wavelength,
)

from thz_opt.arrays.fibonacci import fibonacci_layout
from thz_opt.arrays.geometries import (
    gaussian_layout,
    hierarchical_layout,
    multi_arm_spiral,
    phyllotaxis_layout,
    reuleaux_base,
    reuleaux_layout,
)
from thz_opt.arrays.golden_spiral import PHI, golden_spiral_layout
from thz_opt.arrays.random_array import random_layout, uniform_ring_layout
from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation
from thz_opt.interferometry.baselines import baseline_lengths
from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.metrics.density_matching import cornwell_energy, density_match_chi2
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.metrics.redundancy import redundancy_metrics
from thz_opt.metrics.uv_coverage import unique_cells

N_ANT = 20
EXTENT = dict(r_min=20.0, r_max=1000.0)


def build_catalogue(n: int) -> dict:
    """Every layout family, all at the same N and the same radial extent."""
    base3 = reuleaux_base(1.0, 3)
    return {
        "golden spiral (1 arm)": golden_spiral_layout(n, n_turns=3.0, **EXTENT),
        "log spiral (3 arms)": multi_arm_spiral(n, n_arms=3, n_turns=1.0, **EXTENT),
        "log spiral (5 arms)": multi_arm_spiral(n, n_arms=5, n_turns=1.0, **EXTENT),
        "log spiral (3 arms, growth 4)": multi_arm_spiral(n, n_arms=3, growth=4.0, **EXTENT),
        "Fibonacci radial": fibonacci_layout(n, variant="radial", skip=2, **EXTENT),
        "Vogel golden angle": phyllotaxis_layout(n, **EXTENT),
        "Reuleaux (1 ring)": reuleaux_layout(n, **EXTENT),
        "Reuleaux (3 rings)": reuleaux_layout(n, n_rings=3, **EXTENT),
        "hierarchical (3^3)": hierarchical_layout(base3, levels=3, scale=0.35,
                                                  r_max=EXTENT["r_max"]),
        "Gaussian (Boone target)": gaussian_layout(n, seed=20250916, **EXTENT),
        "random (area-uniform)": random_layout(n, seed=20250916, **EXTENT),
        "uniform rings": uniform_ring_layout(n, n_rings=3, **EXTENT),
    }


def score(name: str, xy: np.ndarray, cfg: dict, grid: UVGrid, obs) -> dict:
    uv = layout_to_uv_tracks(xy, obs)
    occ, outside = uv_occupancy(uv, grid)
    psf, info = dirty_beam(uv, grid, cfg["psf"]["weighting"])
    lengths = baseline_lengths(xy)
    sep = check_minimum_separation(
        xy, SeparationConfig(cfg["constraints"]["dish_diameter_m"],
                             cfg["constraints"]["minimum_separation_factor"])
    )
    row = {
        "layout": name,
        "n_antennas": int(xy.shape[0]),
        "unique_cells": unique_cells(occ),
        "samples_outside_grid": outside,
        "redundancy_fraction": redundancy_metrics(occ)["redundancy_fraction"],
        "cornwell_F": cornwell_energy(uv, max_samples=3000)["cornwell_F"],
        "gaussian_chi2": density_match_chi2(uv, n_bins=30)["normalised_chi2"],
        "d_min_m": float(lengths.min()),
        "d_max_m": float(lengths.max()),
        "separation_violations": sep["n_violations"],
        "feasible": sep["valid"],
    }
    row.update({f"psf_{k}": v for k, v in psf_metrics(psf, info["pixel_scale_arcsec"]).items()})
    return row


def main() -> None:
    ensure_dirs()
    cfg_all = load_config()
    cfg = merge(cfg_all, {"array": {"n_antennas": N_ANT}})
    obs = observation_config(cfg)
    ps = param_string(cfg, "Earth-rotation sampling, shared grid, uniform weighting")
    banner(f"Experiment 06 -- layout shootout -- {ps}")

    layouts = build_catalogue(N_ANT)
    uv_max = cfg["uv_grid"]["pad"] * max(
        float(np.abs(layout_to_uv_tracks(xy, obs)).max()) for xy in layouts.values()
    )
    grid = UVGrid(uv_max=uv_max, n_cells=cfg["uv_grid"]["n_cells"])

    rows = [score(name, xy, cfg, grid, obs) for name, xy in layouts.items()]
    rows.sort(key=lambda r: -r["unique_cells"])

    print(f"\n  {'layout':<30}{'cells':>7}{'redund':>8}{'CornwellF':>11}"
          f"{'gaussX2':>9}{'PSL':>8}{'d_min':>8}{'ok':>4}")
    print("  " + "-" * 87)
    for r in rows:
        print(f"  {r['layout']:<30}{r['unique_cells']:>7}{r['redundancy_fraction']:>8.3f}"
              f"{r['cornwell_F']:>11.2f}{r['gaussian_chi2']:>9.3f}"
              f"{r['psf_peak_sidelobe_level']:>8.3f}{r['d_min_m']:>8.1f}"
              f"{'y' if r['feasible'] else 'n':>4}")

    best = {
        "unique cells (max)": max(rows, key=lambda r: r["unique_cells"])["layout"],
        "Cornwell F (min)": min(rows, key=lambda r: r["cornwell_F"])["layout"],
        "Gaussian chi2 (min)": min(rows, key=lambda r: r["gaussian_chi2"])["layout"],
        "PSF sidelobe (min)": min(rows, key=lambda r: r["psf_peak_sidelobe_level"])["layout"],
        "PSF FWHM (min)": min(rows, key=lambda r: r["psf_fwhm_arcsec"])["layout"],
    }
    print("\n  Best layout per objective -- note that they disagree:")
    for k, v in best.items():
        print(f"    {k:<22} {v}")

    save_csv(rows, RESULTS / "exp06_shootout.csv")
    alpha_rows = alpha_sweep(cfg, grid, obs)
    _figures(cfg, layouts, rows, alpha_rows, ps)
    print(f"\nWrote {RESULTS / 'exp06_shootout.csv'} and exp06_alpha_sweep.csv")


def alpha_sweep(cfg: dict, grid: UVGrid, obs) -> list:
    """Is the golden angle special?  Scan the phyllotaxis angle."""
    named = {
        "1/PHI^2 (golden angle)": 1.0 / PHI ** 2,
        "sqrt(2)-1": np.sqrt(2.0) - 1.0,
        "sqrt(3)-1": np.sqrt(3.0) - 1.0,
        "1/e": 1.0 / np.e,
        "1/pi": 1.0 / np.pi,
        "1/2 (rational)": 0.5,
        "1/3 (rational)": 1.0 / 3.0,
        "3/8 (Fibonacci ratio)": 3.0 / 8.0,
        "5/13 (Fibonacci ratio)": 5.0 / 13.0,
        "21/55 (Fibonacci ratio)": 21.0 / 55.0,
    }
    rows = []
    for label, a in named.items():
        xy = phyllotaxis_layout(N_ANT, alpha=a, **EXTENT)
        r = score(label, xy, cfg, grid, obs)
        r["alpha"] = a
        rows.append(r)

    print("\n  Phyllotaxis angle sweep (same N, extent, grid):")
    print(f"    {'alpha':<26}{'value':>10}{'cells':>8}{'CornwellF':>11}{'PSL':>8}")
    for r in sorted(rows, key=lambda q: -q["unique_cells"]):
        print(f"    {r['layout']:<26}{r['alpha']:>10.5f}{r['unique_cells']:>8}"
              f"{r['cornwell_F']:>11.2f}{r['psf_peak_sidelobe_level']:>8.3f}")

    # continuous scan for the figure
    fine = []
    for a in np.linspace(0.05, 0.95, 181):
        xy = phyllotaxis_layout(N_ANT, alpha=float(a), **EXTENT)
        uv = layout_to_uv_tracks(xy, obs)
        occ, _ = uv_occupancy(uv, grid)
        fine.append({"alpha": float(a), "unique_cells": unique_cells(occ)})
    save_csv(rows + fine, RESULTS / "exp06_alpha_sweep.csv")
    return rows + fine


def _figures(cfg: dict, layouts: dict, rows: list, alpha_rows: list, ps: str) -> None:
    # Figure 13 -- layout gallery
    names = list(layouts)
    ncol = 4
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = new_figure(nrow, ncol, figsize=(3.0 * ncol, 3.0 * nrow))
    for ax, name in zip(axes.ravel(), names):
        xy = layouts[name]
        ax.scatter(xy[:, 0], xy[:, 1], s=9, c="#1f4e79")
        ax.set_aspect("equal")
        ax.set_title(name, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.ravel()[len(names):]:
        ax.axis("off")
    fig.suptitle("Figure 13. Candidate layout families, same N and radial extent", y=1.0)
    stamp(fig, ps)
    save_figure(fig, "fig13_layout_gallery.png")

    # Figure 14 -- objectives disagree
    labels = [r["layout"] for r in rows]
    cells = np.array([r["unique_cells"] for r in rows], dtype=float)
    corn = np.array([r["cornwell_F"] for r in rows], dtype=float)
    chi2 = np.array([r["gaussian_chi2"] for r in rows], dtype=float)
    psl = np.array([r["psf_peak_sidelobe_level"] for r in rows], dtype=float)

    fig, axes = new_figure(1, 2, figsize=(12.5, 5.0))
    axes[0].scatter(cells, corn, s=30, c="#4C72B0")
    for x, y, t in zip(cells, corn, labels):
        axes[0].annotate(t, (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("unique UV cells (higher = more of the UV plane touched)")
    axes[0].set_ylabel("Cornwell F (lower = more uniformly spread)")
    axes[0].set_title("Figure 14a. Cell count vs UV-point repulsion")

    axes[1].scatter(chi2, psl, s=30, c="#C44E52")
    for x, y, t in zip(chi2, psl, labels):
        axes[1].annotate(t, (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Gaussian density mismatch chi2 (lower = closer to Boone target)")
    axes[1].set_ylabel("PSF peak sidelobe level")
    axes[1].set_title("Figure 14b. Density-shape mismatch vs measured sidelobes")
    stamp(fig, ps)
    save_figure(fig, "fig14_objective_disagreement.png")

    # Figure 15 -- alpha sweep
    fine = [r for r in alpha_rows if "layout" not in r]
    named = [r for r in alpha_rows if "layout" in r]
    fig, ax = new_figure(figsize=(9.0, 4.4))
    ax.plot([r["alpha"] for r in fine], [r["unique_cells"] for r in fine],
            lw=0.9, color="0.4", label="phyllotaxis scan")
    for r in named:
        ax.scatter([r["alpha"]], [r["unique_cells"]], s=26, zorder=3)
        ax.annotate(r["layout"], (r["alpha"], r["unique_cells"]), fontsize=6,
                    xytext=(3, 4), textcoords="offset points")
    ax.set_xlabel("phyllotaxis angle parameter alpha (turns per antenna)")
    ax.set_ylabel("unique UV cells")
    ax.set_title("Figure 15. Is the golden angle special?\n"
                 "UV cells vs angular step; dips are rational (commensurate) angles")
    ax.legend(loc="lower right")
    stamp(fig, ps)
    save_figure(fig, "fig15_alpha_sweep.png")


if __name__ == "__main__":
    main()
