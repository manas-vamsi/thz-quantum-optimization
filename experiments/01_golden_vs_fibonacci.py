"""Experiment 01 -- controlled golden-spiral vs Fibonacci comparison.

Runs the whole classical pipeline for every layout model at N = 10, 20, 50, 100:

    layout -> baselines -> snapshot UV -> Earth-rotation UV -> UV occupancy
           -> coverage / redundancy / baseline-distribution / PSF metrics

Every model in a given N shares the radial extent, the wavelength, the
observing geometry, the UV grid and the metric definitions.  The UV grid extent
is taken from the largest |uv| across all models at that N, so no model is
measured on a grid tailored to itself.

The script prints a Metric | Golden | Fibonacci | Difference table and declares
no winner: the numbers are what they are, and which of them matters depends on
a science case this repository has not yet fixed.

Outputs
-------
data/generated/exp01_<label>_<model>.npz   coordinates, baselines, uv, occupancy
data/results/exp01_metrics.csv             one row per (N, model)
data/results/exp01_metrics.json            same, plus the full configuration
figures/fig01_golden_positions.png
figures/fig02_fibonacci_positions.png
figures/fig05_baseline_distributions.png
figures/fig09_metric_comparison.png
"""

from __future__ import annotations

import numpy as np

from common import (  # noqa: E402
    GENERATED,
    LAYOUT_LABELS,
    RESULTS,
    banner,
    build_layout,
    ensure_dirs,
    load_config,
    merge,
    new_figure,
    observation_config,
    param_string,
    save_csv,
    save_figure,
    save_json,
    save_npz,
    stamp,
    wavelength,
)

from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation
from thz_opt.interferometry.baselines import baseline_lengths, compute_baselines
from thz_opt.interferometry.earth_rotation import elevation_deg, layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, layout_to_uv, uv_occupancy
from thz_opt.metrics.baseline_distribution import baseline_summary, empirical_cdf, length_histogram
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.metrics.redundancy import redundancy_metrics
from thz_opt.metrics.uv_coverage import coverage_summary


def analyse(name: str, xy: np.ndarray, cfg: dict, grid_snap: UVGrid, grid_rot: UVGrid) -> dict:
    lam = wavelength(cfg)
    obs = observation_config(cfg)
    sep = SeparationConfig(
        cfg["constraints"]["dish_diameter_m"], cfg["constraints"]["minimum_separation_factor"]
    )

    uv_snap = layout_to_uv(xy, lam)
    occ_snap, _ = uv_occupancy(uv_snap, grid_snap)
    uv_rot = layout_to_uv_tracks(xy, obs)
    occ_rot, _ = uv_occupancy(uv_rot, grid_rot)

    psf, psf_info = dirty_beam(uv_rot, grid_rot, cfg["psf"]["weighting"])
    lengths = baseline_lengths(xy)
    b, _, _ = compute_baselines(xy)

    row = {"model": name, "n_antennas": xy.shape[0]}
    row.update({f"snap_{k}": v for k, v in coverage_summary(uv_snap, occ_snap, grid_snap).items()})
    row.update({f"rot_{k}": v for k, v in coverage_summary(uv_rot, occ_rot, grid_rot).items()})
    row.update({f"rot_red_{k}": v for k, v in redundancy_metrics(occ_rot).items()})
    row.update(baseline_summary(lengths, cfg["analysis"]["n_length_bins"]))
    row.update({f"psf_{k}": v for k, v in
                psf_metrics(psf, psf_info["pixel_scale_arcsec"]).items()})
    sepr = check_minimum_separation(xy, sep)
    row.update({"sep_valid": sepr["valid"], "sep_violations": sepr["n_violations"],
                "sep_min_m": sepr["min_separation_m"], "sep_d_min_m": sepr["d_min_m"]})

    save_npz(
        GENERATED / f"exp01_{cfg['_label']}_{name}.npz",
        antennas_xy_m=xy, baselines_m=b, baseline_lengths_m=lengths,
        uv_snapshot_lambda=uv_snap, uv_rotation_lambda=uv_rot,
        occupancy_snapshot=occ_snap, occupancy_rotation=occ_rot, psf=psf,
    )
    return row


def comparison_table(rows: list, a: str, b: str) -> list:
    """Metric | A | B | Difference rows, in the order a scientist reads them."""
    ra = next(r for r in rows if r["model"] == a)
    rb = next(r for r in rows if r["model"] == b)
    keys = [
        ("snap_unique_cells", "Snapshot unique UV cells", "cells"),
        ("snap_occupied_fraction_annulus", "Snapshot filled fraction (annulus)", "-"),
        ("rot_unique_cells", "Earth-rotation unique UV cells", "cells"),
        ("rot_occupied_fraction_annulus", "Earth-rotation filled fraction (annulus)", "-"),
        ("rot_red_redundancy_fraction", "Earth-rotation redundancy fraction", "-"),
        ("rot_density_uniformity", "UV density uniformity (entropy)", "-"),
        ("rot_angular_uniformity", "UV angular uniformity (entropy)", "-"),
        ("d_min_m", "Shortest baseline", "m"),
        ("d_max_m", "Longest baseline", "m"),
        ("dynamic_range", "Baseline dynamic range", "-"),
        ("log_spacing_uniformity", "Baseline log-spacing uniformity", "-"),
        ("powerlaw_alpha", "Baseline power-law exponent", "-"),
        ("powerlaw_r_squared", "Power-law fit R^2", "-"),
        ("psf_peak_sidelobe_level", "PSF peak sidelobe level", "-"),
        ("psf_sidelobe_rms", "PSF sidelobe RMS", "-"),
        ("psf_fwhm_arcsec", "PSF FWHM", "arcsec"),
        ("psf_ellipticity", "PSF ellipticity", "-"),
    ]
    out = []
    for k, label, unit in keys:
        va, vb = float(ra[k]), float(rb[k])
        out.append({"metric": label, "unit": unit, a: va, b: vb, "difference": vb - va,
                    "relative_difference": (vb - va) / va if va else float("nan")})
    return out


def print_table(table: list, a: str, b: str) -> None:
    print(f"\n{'Metric':<42}{'unit':<8}{a:>14}{b:>14}{'difference':>14}")
    print("-" * 92)
    for r in table:
        print(f"{r['metric']:<42}{r['unit']:<8}{r[a]:>14.4g}{r[b]:>14.4g}"
              f"{r['difference']:>14.4g}")
    print("-" * 92)
    print("No winner is declared: these are descriptive statistics of one "
          "controlled comparison,\nnot a ranking of imaging performance.")


def main() -> None:
    ensure_dirs()
    cfg_all = load_config()
    models = cfg_all["layouts"]
    all_rows, all_tables = [], {}

    for key, spec in cfg_all["experiments"].items():
        cfg = merge(cfg_all, {k: v for k, v in spec.items() if k != "label"})
        cfg["_label"] = spec["label"]
        banner(f"Experiment {key} -- {spec['label']} -- {param_string(cfg)}")

        layouts = {m: build_layout(m, cfg) for m in models}
        lam = wavelength(cfg)
        obs = observation_config(cfg)
        el = elevation_deg(obs)
        print(f"  source elevation over the window: {el.min():.1f} to {el.max():.1f} deg")

        # one shared grid per experiment, sized by the most extended model
        pad, n_cells = cfg["uv_grid"]["pad"], cfg["uv_grid"]["n_cells"]
        uv_max_snap = max(np.abs(layout_to_uv(xy, lam)).max() for xy in layouts.values())
        uv_max_rot = max(np.abs(layout_to_uv_tracks(xy, obs)).max() for xy in layouts.values())
        grid_snap = UVGrid(uv_max=pad * uv_max_snap, n_cells=n_cells)
        grid_rot = UVGrid(uv_max=pad * uv_max_rot, n_cells=n_cells)
        print(f"  shared snapshot grid: cell {grid_snap.cell_size:.3g} lambda, "
              f"rotation grid: cell {grid_rot.cell_size:.3g} lambda")

        rows = []
        for m, xy in layouts.items():
            row = analyse(m, xy, cfg, grid_snap, grid_rot)
            row["experiment"] = key
            row["label"] = spec["label"]
            rows.append(row)
            print(f"  {LAYOUT_LABELS[m]:<28} unique cells (rot): {row['rot_unique_cells']:>7} "
                  f"| PSL: {row['psf_peak_sidelobe_level']:.3f} "
                  f"| d_min: {row['d_min_m']:.1f} m")

        table = comparison_table(rows, "golden", "fibonacci_radial")
        all_tables[spec["label"]] = table
        if spec["label"] == "N20":
            print_table(table, "golden", "fibonacci_radial")
            _figures(cfg, layouts, rows, table)
        all_rows.extend(rows)

    save_csv(all_rows, RESULTS / "exp01_metrics.csv")
    save_json({"config": {k: v for k, v in cfg_all.items()},
               "rows": all_rows, "tables": all_tables},
              RESULTS / "exp01_metrics.json")
    print(f"\nWrote {RESULTS / 'exp01_metrics.csv'} and exp01_metrics.json")


def _figures(cfg: dict, layouts: dict, rows: list, table: list) -> None:
    ps = param_string(cfg)

    for fname, model, fignum in (
        ("fig01_golden_positions.png", "golden", "Figure 1"),
        ("fig02_fibonacci_positions.png", "fibonacci_radial", "Figure 2"),
    ):
        xy = layouts[model]
        fig, ax = new_figure(figsize=(5.2, 5.0))
        ax.scatter(xy[:, 0], xy[:, 1], s=26, c=np.arange(len(xy)), cmap="viridis",
                   edgecolor="k", linewidth=0.3, zorder=3, label="antenna pad")
        ax.plot(xy[:, 0], xy[:, 1], lw=0.5, c="0.6", zorder=2, label="index order")
        ax.set_aspect("equal")
        ax.set_title(f"{fignum}. {LAYOUT_LABELS[model]} antenna positions")
        ax.set_xlabel("East x (m)")
        ax.set_ylabel("North y (m)")
        ax.legend(loc="upper right")
        stamp(fig, ps)
        save_figure(fig, fname)

    # Figure 5 -- baseline length distributions
    fig, axes = new_figure(1, 3, figsize=(12.5, 3.8))
    for m in ("golden", "fibonacci_radial", "fibonacci_golden_angle", "random"):
        d = baseline_lengths(layouts[m])
        c, n, _ = length_histogram(d, cfg["analysis"]["n_length_bins"])
        axes[0].step(c, n, where="mid", label=LAYOUT_LABELS[m])
        cl, nl, _ = length_histogram(d, cfg["analysis"]["n_length_bins"], log_bins=True)
        axes[1].step(cl, nl, where="mid", label=LAYOUT_LABELS[m])
        x, cdf = empirical_cdf(d)
        axes[2].plot(x, cdf, label=LAYOUT_LABELS[m])
    axes[0].set_title("Figure 5a. Baseline lengths (linear bins)")
    axes[0].set_xlabel("baseline length (m)")
    axes[0].set_ylabel("number of baselines")
    axes[1].set_title("Figure 5b. Baseline lengths (log bins)")
    axes[1].set_xlabel("baseline length (m)")
    axes[1].set_ylabel("number of baselines per log bin")
    axes[1].set_xscale("log")
    axes[2].set_title("Figure 5c. Empirical CDF")
    axes[2].set_xlabel("baseline length (m)")
    axes[2].set_ylabel("fraction of baselines <= d")
    axes[2].set_xscale("log")
    for ax in axes:
        ax.legend()
    stamp(fig, ps)
    save_figure(fig, "fig05_baseline_distributions.png")

    # Figure 9 -- metric comparison
    keys = [t["metric"] for t in table]
    rel = [t["relative_difference"] for t in table]
    fig, ax = new_figure(figsize=(8.0, 5.2))
    colors = ["#4C72B0" if v >= 0 else "#C44E52" for v in rel]
    ax.barh(range(len(keys)), rel, color=colors)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(keys, fontsize=7)
    ax.axvline(0.0, color="k", lw=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("(Fibonacci radial - Golden spiral) / Golden spiral")
    ax.set_title("Figure 9. Relative metric difference, Fibonacci vs golden spiral\n"
                 "(sign only shows direction; no metric is a quality score)")
    stamp(fig, ps)
    save_figure(fig, "fig09_metric_comparison.png")


if __name__ == "__main__":
    main()
