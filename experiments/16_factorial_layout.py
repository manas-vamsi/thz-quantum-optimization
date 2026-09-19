"""Experiment 16 -- which matters, the radial law or the angular law?

The gap this closes
-------------------
Experiment 6 scores twelve named layout families against each other. It is a
fair comparison of those twelve objects and it cannot say *why* any of them
wins, because they differ in two ways simultaneously: a golden spiral changes
both how radius grows with index and how bearing advances. When it beats a
random array, the experiment cannot attribute that to the radial law, the
angular law, or only to their combination.

The distinction is practical. "Space antennas uniformly in area" and "advance
by the golden angle" are separate pieces of advice, and neither can be read off
a ranking of composite objects.

Design
------
A full factorial: five radial laws crossed with five angular laws, each cell
replicated with positional jitter so there is a genuine error term. Everything
else -- N, radial extent, UV grid, declination, track length, shadowing rule --
is held fixed, so the only things varying are the two factors.

The variance in each score is then decomposed the standard way,

    SS_total = SS_radial + SS_angular + SS_interaction + SS_residual

and reported as the fraction of explained variance attributable to each, with
an F ratio against the replication noise. If the interaction term dominates,
the two laws cannot be chosen independently and the composite families were the
right unit of comparison after all -- which is itself worth knowing.

Outputs
-------
figures/fig33_factorial_heatmap.png
figures/fig34_factorial_variance.png
data/results/exp16_factorial.csv
data/results/exp16_summary.json
"""

from __future__ import annotations

import itertools

import numpy as np

from common import (  # noqa: E402
    RESULTS,
    banner,
    ensure_dirs,
    load_config,
    new_figure,
    save_csv,
    save_figure,
    save_json,
    stamp,
)

from thz_opt.arrays.factorial import (
    ANGULAR_LAWS,
    RADIAL_LAWS,
    layout_from_laws,
)
from thz_opt.constraints.separation import SeparationConfig, check_minimum_separation
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.metrics.psf_metrics import psf_metrics

N_ANT = 24
R_MIN = 30.0
R_MAX = 2000.0
FREQ_HZ = 3.0e11
GRID_CELLS = 64
N_REPLICATES = 5
JITTER_M = 12.0
SEED = 20260919


def two_way_anova(values, row_idx, col_idx, n_rows, n_cols) -> dict:
    """Decompose variance into two main effects, interaction and residual.

    A balanced design, so the textbook sums of squares apply directly. The
    residual comes from the jitter replicates, which is why replication was
    built into the design rather than added afterwards.
    """
    v = np.asarray(values, dtype=float)
    grand = v.mean()
    ss_total = float(((v - grand) ** 2).sum())

    cell_mean = np.zeros((n_rows, n_cols))
    cell_n = np.zeros((n_rows, n_cols), dtype=int)
    for k in range(len(v)):
        cell_mean[row_idx[k], col_idx[k]] += v[k]
        cell_n[row_idx[k], col_idx[k]] += 1
    cell_mean /= np.maximum(cell_n, 1)

    n_rep = int(cell_n.max())
    row_mean = cell_mean.mean(axis=1)
    col_mean = cell_mean.mean(axis=0)

    ss_row = float(n_cols * n_rep * ((row_mean - grand) ** 2).sum())
    ss_col = float(n_rows * n_rep * ((col_mean - grand) ** 2).sum())
    inter = cell_mean - row_mean[:, None] - col_mean[None, :] + grand
    ss_inter = float(n_rep * (inter ** 2).sum())
    ss_resid = float(sum((v[k] - cell_mean[row_idx[k], col_idx[k]]) ** 2
                         for k in range(len(v))))

    df_row, df_col = n_rows - 1, n_cols - 1
    df_inter = df_row * df_col
    df_resid = n_rows * n_cols * (n_rep - 1)
    ms_resid = ss_resid / df_resid if df_resid > 0 else np.nan

    def f_ratio(ss, df):
        return float((ss / df) / ms_resid) if ms_resid and ms_resid > 0 else np.nan

    return {
        "ss_total": ss_total,
        "eta2_radial": ss_row / ss_total if ss_total else np.nan,
        "eta2_angular": ss_col / ss_total if ss_total else np.nan,
        "eta2_interaction": ss_inter / ss_total if ss_total else np.nan,
        "eta2_residual": ss_resid / ss_total if ss_total else np.nan,
        "f_radial": f_ratio(ss_row, df_row),
        "f_angular": f_ratio(ss_col, df_col),
        "f_interaction": f_ratio(ss_inter, df_inter),
        "cell_mean": cell_mean,
        "n_replicates": n_rep,
    }


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    lam = 299_792_458.0 / FREQ_HZ

    obs = ObservationConfig(
        frequency_hz=FREQ_HZ, declination_deg=-30.0, latitude_deg=-23.0229,
        hour_angle_start_h=-2.0, hour_angle_end_h=2.0, n_times=41,
    )
    sep = SeparationConfig(12.0, cfg["constraints"]["minimum_separation_factor"])

    ps = (f"N={N_ANT}, r={R_MIN:.0f}-{R_MAX:.0f} m, {FREQ_HZ/1e9:.0f} GHz, "
          f"4 h track, {N_REPLICATES} replicates, {JITTER_M:.0f} m jitter")
    banner(f"Experiment 16 -- factorial layout design -- {ps}")
    print(f"\n  {len(RADIAL_LAWS)} radial laws x {len(ANGULAR_LAWS)} angular laws "
          f"x {N_REPLICATES} replicates = "
          f"{len(RADIAL_LAWS)*len(ANGULAR_LAWS)*N_REPLICATES} layouts")

    # one common UV grid so every layout is scored on the same cells
    ref = layout_from_laws("sqrt", "golden", N_ANT, R_MIN, R_MAX, seed=SEED)
    uv_ref = layout_to_uv_tracks(ref, obs)
    grid = UVGrid(uv_max=1.05 * float(np.abs(uv_ref).max()), n_cells=GRID_CELLS)

    rows = []
    for ri, radial in enumerate(RADIAL_LAWS):
        for ci, angular in enumerate(ANGULAR_LAWS):
            for rep in range(N_REPLICATES):
                xy = layout_from_laws(radial, angular, N_ANT, R_MIN, R_MAX,
                                      jitter_m=JITTER_M,
                                      seed=SEED + 1000 * rep + 7 * ri + ci)
                uv = layout_to_uv_tracks(xy, obs)
                occ, _ = uv_occupancy(uv.reshape(-1, 2), grid)
                beam, binfo = dirty_beam(uv.reshape(-1, 2), grid)
                m = psf_metrics(beam, binfo["pixel_scale_arcsec"])
                ok = check_minimum_separation(xy, sep)["valid"]
                rows.append({
                    "radial": radial, "angular": angular, "replicate": rep,
                    "radial_idx": ri, "angular_idx": ci,
                    "unique_cells": int(np.count_nonzero(occ)),
                    "peak_sidelobe": float(m.get("peak_sidelobe_level", np.nan)),
                    "fwhm_arcsec": float(m.get("fwhm_arcsec", np.nan)),
                    "feasible": bool(ok),
                })

    n_feasible = sum(r["feasible"] for r in rows)
    print(f"  {n_feasible} of {len(rows)} layouts satisfy the shadowing rule")

    # ---- the decomposition ------------------------------------------------
    print(f"\n  [1] where does the variation come from?")
    anova = {}
    for metric, better in (("unique_cells", "higher"),
                           ("peak_sidelobe", "lower")):
        vals = [r[metric] for r in rows]
        a = two_way_anova(vals, [r["radial_idx"] for r in rows],
                          [r["angular_idx"] for r in rows],
                          len(RADIAL_LAWS), len(ANGULAR_LAWS))
        anova[metric] = {k: v for k, v in a.items() if k != "cell_mean"}
        anova[metric]["cell_mean"] = a["cell_mean"].tolist()
        print(f"\n      {metric} ({better} is better)")
        print(f"        radial law      {100*a['eta2_radial']:5.1f} % of variance"
              f"   F = {a['f_radial']:8.1f}")
        print(f"        angular law     {100*a['eta2_angular']:5.1f} % of variance"
              f"   F = {a['f_angular']:8.1f}")
        print(f"        interaction     {100*a['eta2_interaction']:5.1f} % of variance"
              f"   F = {a['f_interaction']:8.1f}")
        print(f"        jitter noise    {100*a['eta2_residual']:5.1f} % of variance")

    # ---- best and worst of each factor, averaged over the other ----------
    print(f"\n  [2] marginal means (averaged over the other factor)")
    for metric in ("unique_cells", "peak_sidelobe"):
        cm = np.array(anova[metric]["cell_mean"])
        print(f"\n      {metric}")
        rmeans = cm.mean(axis=1)
        cmeans = cm.mean(axis=0)
        for name, v in sorted(zip(RADIAL_LAWS, rmeans), key=lambda t: -t[1]):
            print(f"        radial  {name:14}{v:10.1f}")
        for name, v in sorted(zip(ANGULAR_LAWS, cmeans), key=lambda t: -t[1]):
            print(f"        angular {name:14}{v:10.1f}")

    # ---- the verdict ------------------------------------------------------
    a = anova["unique_cells"]
    print(f"\n  [3] verdict on UV coverage")
    dominant = max(("radial law", a["eta2_radial"]),
                   ("angular law", a["eta2_angular"]),
                   ("their interaction", a["eta2_interaction"]),
                   key=lambda t: t[1])
    print(f"      The {dominant[0]} explains {100*dominant[1]:.1f} % of the variance.")
    if a["eta2_interaction"] > max(a["eta2_radial"], a["eta2_angular"]):
        print("      Interaction dominates: the two laws cannot be chosen")
        print("      independently, and comparing composite families was right.")
    else:
        print("      The main effects dominate the interaction, so the two laws")
        print("      can be chosen largely independently -- and a ranking of")
        print("      composite families conflates two separable decisions.")

    save_csv(rows, RESULTS / "exp16_factorial.csv")
    save_json({"n_antennas": N_ANT, "r_min_m": R_MIN, "r_max_m": R_MAX,
               "replicates": N_REPLICATES, "jitter_m": JITTER_M,
               "radial_laws": list(RADIAL_LAWS),
               "angular_laws": list(ANGULAR_LAWS),
               "anova": anova, "rows": rows},
              RESULTS / "exp16_summary.json")
    _figures(anova, ps)
    print(f"\nWrote {RESULTS/'exp16_factorial.csv'} and exp16_summary.json")


def _figures(anova, ps) -> None:
    fig, ax = new_figure(1, 2, figsize=(12.0, 4.6))
    for k, (metric, cmap, label) in enumerate(
            (("unique_cells", "viridis", "occupied UV cells"),
             ("peak_sidelobe", "magma_r", "peak sidelobe"))):
        cm = np.array(anova[metric]["cell_mean"])
        im = ax[k].imshow(cm, cmap=cmap, aspect="auto")
        ax[k].set_xticks(range(len(ANGULAR_LAWS)))
        ax[k].set_xticklabels(ANGULAR_LAWS, rotation=30, ha="right", fontsize=7)
        ax[k].set_yticks(range(len(RADIAL_LAWS)))
        ax[k].set_yticklabels(RADIAL_LAWS, fontsize=7)
        ax[k].set_xlabel("angular law")
        ax[k].set_ylabel("radial law")
        ax[k].set_title(f"{label}\nmean over replicates", fontsize=9)
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax[k].text(j, i, f"{cm[i, j]:.3g}", ha="center", va="center",
                           fontsize=6,
                           color="white" if cmap == "viridis" else "black")
        fig.colorbar(im, ax=ax[k], fraction=0.046)
    stamp(fig, ps)
    save_figure(fig, "fig33_factorial_heatmap")

    fig2, ax2 = new_figure(figsize=(7.2, 4.3))
    metrics = ["unique_cells", "peak_sidelobe"]
    parts = ["eta2_radial", "eta2_angular", "eta2_interaction", "eta2_residual"]
    labels = ["radial law", "angular law", "interaction", "jitter noise"]
    colours = ["#2b6cb0", "#c53030", "#dd6b20", "#a0aec0"]
    bottom = np.zeros(len(metrics))
    for p, lab, col in zip(parts, labels, colours):
        vals = np.array([100 * anova[m][p] for m in metrics])
        ax2.bar(metrics, vals, bottom=bottom, color=col, label=lab)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 4:
                ax2.text(i, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                         fontsize=8, color="white")
        bottom += vals
    ax2.set_ylabel("% of total variance")
    ax2.set_title("What actually drives layout performance", fontsize=10)
    ax2.legend(fontsize=7, loc="lower right")
    stamp(fig2, ps)
    save_figure(fig2, "fig34_factorial_variance")


if __name__ == "__main__":
    main()
