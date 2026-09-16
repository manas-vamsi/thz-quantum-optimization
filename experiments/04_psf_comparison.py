"""Experiment 04 -- dirty beam / PSF comparison.

The PSF is the inverse FFT of the gridded sampling function, so it inherits
every choice made about the grid.  All models here share one grid, one
weighting and one observing geometry; absolute sidelobe levels are therefore
comparable between the models in this figure and nowhere else.

Outputs
-------
figures/fig08_psf_comparison.png
data/results/exp04_psf_metrics.csv
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
)

from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid

COMPARE = ("golden", "fibonacci_radial", "fibonacci_golden_angle", "random")


def main() -> None:
    ensure_dirs()
    cfg_all = load_config()
    cfg = merge(cfg_all, cfg_all["experiments"]["B"])  # N = 20
    obs = observation_config(cfg)
    weighting = cfg["psf"]["weighting"]
    ps = param_string(cfg, f"Earth-rotation sampling, {weighting} weighting, no taper, w = 0")
    banner(f"Experiment 04 -- PSF -- {ps}")

    layouts = {m: build_layout(m, cfg) for m in COMPARE}
    tracks = {m: layout_to_uv_tracks(xy, obs) for m, xy in layouts.items()}
    uv_max = cfg["uv_grid"]["pad"] * max(np.abs(t).max() for t in tracks.values())
    grid = UVGrid(uv_max=uv_max, n_cells=cfg["uv_grid"]["n_cells"])

    from thz_opt.metrics.psf_metrics import psf_metrics

    psfs, rows = {}, []
    for m, uv in tracks.items():
        psf, info = dirty_beam(uv, grid, weighting)
        psfs[m] = psf
        met = psf_metrics(psf, info["pixel_scale_arcsec"])
        rows.append({"model": m, **info, **met})
        print(f"  {LAYOUT_LABELS[m]:<28} PSL {met['peak_sidelobe_level']:.4f} | "
              f"sidelobe RMS {met['sidelobe_rms']:.4f} | FWHM {met['fwhm_arcsec']:.4f} arcsec | "
              f"imag residual {info['imag_residual_fraction']:.2e}")

    n = grid.n_cells
    half = 48  # zoom window in pixels around the beam centre
    sl = slice(n // 2 - half, n // 2 + half)
    pix = rows[0]["pixel_scale_arcsec"]
    ext = [-half * pix, half * pix, -half * pix, half * pix]

    fig, axes = new_figure(2, len(COMPARE), figsize=(4.0 * len(COMPARE), 7.6))
    for k, m in enumerate(COMPARE):
        img = psfs[m][sl, sl]
        im = axes[0, k].imshow(img, origin="lower", extent=ext, cmap="magma", vmin=-0.2, vmax=1.0)
        axes[0, k].set_title(LAYOUT_LABELS[m])
        axes[0, k].set_xlabel("l offset (arcsec)")
        if k == 0:
            axes[0, k].set_ylabel("m offset (arcsec)")
        fig.colorbar(im, ax=axes[0, k], shrink=0.8, label="normalised PSF")

        cut = psfs[m][n // 2, :]
        x = (np.arange(n) - n // 2) * pix
        axes[1, k].plot(x, cut, lw=0.9, label="m = 0 cut")
        axes[1, k].plot(x, psfs[m][:, n // 2], lw=0.9, ls="--", label="l = 0 cut")
        axes[1, k].set_xlim(-half * pix, half * pix)
        axes[1, k].set_ylim(-0.4, 1.05)
        axes[1, k].axhline(0.0, color="k", lw=0.6)
        axes[1, k].set_xlabel("offset (arcsec)")
        if k == 0:
            axes[1, k].set_ylabel("normalised PSF amplitude")
        axes[1, k].legend()
    fig.suptitle("Figure 8. Dirty beam comparison (top: central 96x96 pixels; "
                 "bottom: principal cuts)", y=0.98)
    stamp(fig, ps)
    save_figure(fig, "fig08_psf_comparison.png")

    save_csv(rows, RESULTS / "exp04_psf_metrics.csv")
    print(f"\nWrote {RESULTS / 'exp04_psf_metrics.csv'}")


if __name__ == "__main__":
    main()
