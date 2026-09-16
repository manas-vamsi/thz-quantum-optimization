"""Experiment 02 -- snapshot UV coverage.

Instantaneous (zenith, coplanar) UV sampling only: no Earth rotation, so the
sampling is exactly ``2 * N(N-1)/2`` points and every difference between models
comes from geometry alone.  This is deliberately separated from experiment 03.

Outputs
-------
figures/fig03_golden_uv.png
figures/fig04_fibonacci_uv.png
figures/fig06_uv_radial_density.png
figures/fig07_uv_angular_density.png
data/results/exp02_snapshot_metrics.csv
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
    param_string,
    save_csv,
    save_figure,
    stamp,
    wavelength,
)

from thz_opt.interferometry.uv import UVGrid, layout_to_uv, uv_occupancy
from thz_opt.metrics.redundancy import redundancy_metrics
from thz_opt.metrics.uv_coverage import angular_profile, coverage_summary, radial_profile

COMPARE = ("golden", "fibonacci_radial", "fibonacci_golden_angle", "random")


def uv_panel(ax, uv, occ, grid, title):
    ax.scatter(uv[:, 0] / 1e3, uv[:, 1] / 1e3, s=3, alpha=0.7, c="#1f4e79")
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("u (kilo-wavelengths)")
    ax.set_ylabel("v (kilo-wavelengths)")


def main() -> None:
    ensure_dirs()
    cfg_all = load_config()
    cfg = merge(cfg_all, cfg_all["experiments"]["B"])  # N = 20
    cfg["_label"] = cfg_all["experiments"]["B"]["label"]
    lam = wavelength(cfg)
    ps = param_string(cfg, "snapshot (no Earth rotation), conjugate points included")
    banner(f"Experiment 02 -- snapshot UV -- {ps}")

    layouts = {m: build_layout(m, cfg) for m in cfg_all["layouts"]}
    uvs = {m: layout_to_uv(xy, lam) for m, xy in layouts.items()}
    uv_max = cfg["uv_grid"]["pad"] * max(np.abs(u).max() for u in uvs.values())
    grid = UVGrid(uv_max=uv_max, n_cells=cfg["uv_grid"]["n_cells"])
    print(f"  shared grid: {grid.n_cells}^2 cells of {grid.cell_size:.4g} lambda")

    rows = []
    for m, uv in uvs.items():
        occ, outside = uv_occupancy(uv, grid)
        row = {"model": m, "n_uv_samples": uv.shape[0], "samples_outside_grid": outside}
        row.update(coverage_summary(uv, occ, grid))
        row.update(redundancy_metrics(occ))
        rows.append(row)
        print(f"  {LAYOUT_LABELS[m]:<28} samples {uv.shape[0]:>5} | unique cells "
              f"{row['unique_cells']:>5} | redundancy {row['redundancy_fraction']:.3f}")

    # Figures 3 and 4 -- snapshot UV coverage
    for fname, model, fignum in (
        ("fig03_golden_uv.png", "golden", "Figure 3"),
        ("fig04_fibonacci_uv.png", "fibonacci_radial", "Figure 4"),
    ):
        uv = uvs[model]
        occ, _ = uv_occupancy(uv, grid)
        fig, axes = new_figure(1, 2, figsize=(10.0, 4.6))
        uv_panel(axes[0], uv, occ, grid,
                 f"{fignum}a. {LAYOUT_LABELS[model]} snapshot (u,v) samples")
        im = axes[1].imshow(
            (occ.T > 0).astype(float), origin="lower", cmap="Greys",
            extent=[-grid.uv_max / 1e3, grid.uv_max / 1e3, -grid.uv_max / 1e3, grid.uv_max / 1e3],
        )
        axes[1].set_title(f"{fignum}b. Occupied UV cells (grid {grid.n_cells}^2)")
        axes[1].set_xlabel("u (kilo-wavelengths)")
        axes[1].set_ylabel("v (kilo-wavelengths)")
        axes[1].set_aspect("equal")
        fig.colorbar(im, ax=axes[1], label="cell occupied (1) / empty (0)", shrink=0.8)
        stamp(fig, ps)
        save_figure(fig, fname)

    # Figure 6 -- radial UV density
    fig, axes = new_figure(1, 2, figsize=(11.0, 4.2))
    for m in COMPARE:
        r, counts, density = radial_profile(uvs[m], cfg["analysis"]["n_radial_bins"],
                                            r_max=grid.uv_max)
        axes[0].step(r / 1e3, counts, where="mid", label=LAYOUT_LABELS[m])
        axes[1].step(r / 1e3, density * 1e6, where="mid", label=LAYOUT_LABELS[m])
    axes[0].set_title("Figure 6a. Radial UV sample count")
    axes[0].set_ylabel("samples per annulus")
    axes[1].set_title("Figure 6b. Radial UV density per unit area")
    axes[1].set_ylabel("samples per (kilo-wavelength)^2")
    axes[1].set_yscale("log")
    for ax in axes:
        ax.set_xlabel("UV radius |(u,v)| (kilo-wavelengths)")
        ax.legend()
    stamp(fig, ps)
    save_figure(fig, "fig06_uv_radial_density.png")

    # Figure 7 -- angular UV density
    fig, ax = new_figure(figsize=(7.2, 4.2), subplot_kw={"projection": "polar"})
    for m in COMPARE:
        ang, counts = angular_profile(uvs[m], cfg["analysis"]["n_angular_bins"])
        ang = np.concatenate([ang, ang + np.pi, ang[:1]])
        c = np.concatenate([counts, counts, counts[:1]])
        ax.plot(ang, c, lw=1.2, label=LAYOUT_LABELS[m])
    ax.set_title("Figure 7. UV position-angle distribution (folded onto [0, pi)\n"
                 "and mirrored for display; radius = samples per angular bin)")
    ax.legend(loc="lower left", bbox_to_anchor=(1.02, 0.0))
    stamp(fig, ps)
    save_figure(fig, "fig07_uv_angular_density.png")

    save_csv(rows, RESULTS / "exp02_snapshot_metrics.csv")
    print(f"\nWrote {RESULTS / 'exp02_snapshot_metrics.csv'}")


if __name__ == "__main__":
    main()
