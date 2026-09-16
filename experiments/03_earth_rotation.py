"""Experiment 03 -- Earth-rotation aperture synthesis.

Adds the tracked-observation model of
:mod:`thz_opt.interferometry.earth_rotation` on top of the same layouts, and
measures how UV coverage grows with integration time.  The snapshot result of
experiment 02 is the ``n_times = 1`` limit of this experiment, which is checked
numerically at the start of the run.

Outputs
-------
figures/fig11_uv_tracks.png
figures/fig12_coverage_vs_time.png
data/results/exp03_rotation_metrics.csv
"""

from __future__ import annotations

import numpy as np

from common import (  # noqa: E402
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
    stamp,
    wavelength,
)

from thz_opt.interferometry.earth_rotation import ObservationConfig, elevation_deg, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid, layout_to_uv, uv_occupancy
from thz_opt.metrics.redundancy import redundancy_metrics
from thz_opt.metrics.uv_coverage import coverage_summary

COMPARE = ("golden", "fibonacci_radial", "fibonacci_golden_angle", "random")


def main() -> None:
    ensure_dirs()
    cfg_all = load_config()
    cfg = merge(cfg_all, cfg_all["experiments"]["B"])  # N = 20
    lam = wavelength(cfg)
    obs = observation_config(cfg)
    el = elevation_deg(obs)
    ps = param_string(
        cfg,
        f"dec={obs.declination_deg} deg, lat={obs.latitude_deg} deg, "
        f"HA in [{obs.hour_angle_start_h}, {obs.hour_angle_end_h}] h, {obs.n_times} samples",
    )
    banner(f"Experiment 03 -- Earth rotation -- {ps}")
    print(f"  source elevation {el.min():.1f} to {el.max():.1f} deg (no elevation cut applied)")

    layouts = {m: build_layout(m, cfg) for m in cfg_all["layouts"]}

    # consistency check against the snapshot module
    zen = ObservationConfig(frequency_hz=obs.frequency_hz, declination_deg=obs.latitude_deg,
                            latitude_deg=obs.latitude_deg, hour_angle_start_h=0.0,
                            hour_angle_end_h=0.0, n_times=1)
    err = np.abs(layout_to_uv_tracks(layouts["golden"], zen) - layout_to_uv(layouts["golden"], lam)).max()
    print(f"  zenith snapshot consistency: max |uv difference| = {err:.3e} wavelengths")

    tracks = {m: layout_to_uv_tracks(xy, obs) for m, xy in layouts.items()}
    uv_max = cfg["uv_grid"]["pad"] * max(np.abs(t).max() for t in tracks.values())
    grid = UVGrid(uv_max=uv_max, n_cells=cfg["uv_grid"]["n_cells"])

    rows = []
    for m, uv in tracks.items():
        occ, outside = uv_occupancy(uv, grid)
        row = {"model": m, "n_uv_samples": uv.shape[0], "samples_outside_grid": outside}
        row.update(coverage_summary(uv, occ, grid))
        row.update(redundancy_metrics(occ))
        rows.append(row)
        print(f"  {LAYOUT_LABELS[m]:<28} samples {uv.shape[0]:>7} | unique cells "
              f"{row['unique_cells']:>6} | redundancy {row['redundancy_fraction']:.3f}")

    # Figure 11 -- UV tracks
    fig, axes = new_figure(1, len(COMPARE), figsize=(4.0 * len(COMPARE), 4.2), sharex=True, sharey=True)
    for ax, m in zip(axes, COMPARE):
        uv = tracks[m]
        ax.scatter(uv[:, 0] / 1e3, uv[:, 1] / 1e3, s=0.6, alpha=0.35, c="#1f4e79", linewidths=0)
        ax.set_aspect("equal")
        ax.set_title(f"{LAYOUT_LABELS[m]}")
        ax.set_xlabel("u (kilo-wavelengths)")
    axes[0].set_ylabel("v (kilo-wavelengths)")
    fig.suptitle("Figure 11. Earth-rotation UV tracks (conjugate points included)", y=1.02)
    stamp(fig, ps)
    save_figure(fig, "fig11_uv_tracks.png")

    # Figure 12 -- coverage growth with integration time
    n_steps = [1, 3, 5, 9, 17, 25, 33, 41]
    fig, axes = new_figure(1, 2, figsize=(11.0, 4.2))
    growth_rows = []
    for m in COMPARE:
        uniq, redun = [], []
        for n in n_steps:
            o = ObservationConfig(frequency_hz=obs.frequency_hz,
                                  declination_deg=obs.declination_deg,
                                  latitude_deg=obs.latitude_deg,
                                  hour_angle_start_h=obs.hour_angle_start_h,
                                  hour_angle_end_h=obs.hour_angle_end_h, n_times=n)
            occ, _ = uv_occupancy(layout_to_uv_tracks(layouts[m], o), grid)
            r = redundancy_metrics(occ)
            uniq.append(r["unique_cells"])
            redun.append(r["redundancy_fraction"])
            growth_rows.append({"model": m, "n_times": n, **r})
        axes[0].plot(n_steps, uniq, marker="o", ms=3, label=LAYOUT_LABELS[m])
        axes[1].plot(n_steps, redun, marker="o", ms=3, label=LAYOUT_LABELS[m])
    axes[0].set_title("Figure 12a. Unique UV cells vs number of time samples")
    axes[0].set_ylabel("unique occupied cells")
    axes[1].set_title("Figure 12b. Redundancy fraction vs number of time samples")
    axes[1].set_ylabel("redundant samples / total samples")
    for ax in axes:
        ax.set_xlabel("time samples over the 4 h hour-angle window")
        ax.legend()
    stamp(fig, ps)
    save_figure(fig, "fig12_coverage_vs_time.png")

    save_csv(rows, RESULTS / "exp03_rotation_metrics.csv")
    save_csv(growth_rows, RESULTS / "exp03_coverage_growth.csv")
    print(f"\nWrote {RESULTS / 'exp03_rotation_metrics.csv'} and exp03_coverage_growth.csv")


if __name__ == "__main__":
    main()
