"""Produce the Ladakh array design: pad coordinates, figures and reports.

Runs at both sites the project cares about -- Hanle, which has the
infrastructure, and the best-PWV candidate from the 2026 reanalysis study,
which does not -- and in two configurations at each, so the trade between
convenience and dryness, and between resolution and short spacings, is visible
rather than argued.

Outputs, all under ``outputs/``
-------------------------------
<site>_<config>_pads.csv     every pad: latitude, longitude, elevation, E/N
<site>_<config>_array.cfg    CASA observatory configuration, for simobserve
<site>_<config>_report.json  baselines, resolution, UV coverage, PSF
site_comparison.csv          every design side by side
fig_<site>_<config>_terrain.png   terrain, buildable ground, chosen pads
fig_<site>_<config>_uv.png        UV coverage and dirty beam
fig_site_comparison.png           the two sites against each other
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import design
import site_data
import terrain
from design import DesignSpec
from thz_opt.interferometry.earth_rotation import layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import uv_occupancy

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
SITES = ["hanle", "site_a"]

# Two configurations, not one.
#
# A single 16-antenna array cannot do both jobs. Spread over 3 km it reaches
# 0.09 arcsec, but its shortest baseline is around 250 m, so it is blind to
# anything wider than about 0.6 arcsec: the first design run produced exactly
# one baseline under 300 m out of 120. Concentrated enough to see extended
# emission, it loses the resolution entirely. This is not a defect of the
# optimiser, it is what sixteen antennas buy, and every real observatory
# answers it the same way -- by moving antennas between configurations on a
# shared pad field. ALMA does it, the VLA does it, and Section 4.5 of the
# manuscript measures what reconfiguration is worth.
#
# The compact configuration uses a finer candidate lattice because its whole
# purpose is short baselines, and a 120 m grid cannot express a 30 m spacing.
CONFIGS = [
    ("compact", dict(max_baseline_m=400.0, candidate_spacing_m=30.0)),
    ("extended", dict(max_baseline_m=3000.0, candidate_spacing_m=120.0)),
]


def write_pads_csv(path: Path, result) -> Path:
    """The deliverable: one row per pad, in coordinates a surveyor can use.

    Shares :func:`design.pad_table` with the CFG writer so the two files cannot
    disagree about which pad is which.
    """
    rows = design.pad_table(result["xy"], result["elev"], result["site"])
    keys = ["pad_id", "latitude_deg", "longitude_deg", "elevation_m",
            "east_m", "north_m", "radius_from_centre_m"]
    fmt = {"latitude_deg": "{:.7f}", "longitude_deg": "{:.7f}",
           "elevation_m": "{:.1f}", "east_m": "{:.2f}", "north_m": "{:.2f}",
           "radius_from_centre_m": "{:.1f}"}
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([fmt.get(k, "{}").format(r[k]) for k in keys])
    return path


def figure_terrain(result, name: str) -> None:
    dem, spec = result["dem"], result["spec"]
    cand, xy = result["candidates"], result["xy"]
    mask = terrain.buildable_mask(dem, spec.max_slope_deg, spec.array_radius_m)

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.4))
    ext = [dem.east.min() / 1e3, dem.east.max() / 1e3,
           dem.north.min() / 1e3, dem.north.max() / 1e3]
    im = ax[0].imshow(dem.elevation, origin="lower", extent=ext, cmap="terrain")
    ax[0].set_title(f"{result['site']['name']} -- terrain, "
                    f"{dem.pixel_m:.0f} m/pixel (SRTM/ASTER)", fontsize=9)
    fig.colorbar(im, ax=ax[0], label="elevation (m)", fraction=0.046)

    ax[1].imshow(mask, origin="lower", extent=ext, cmap="Greys_r", alpha=0.85)
    ax[1].scatter(cand[:, 0] / 1e3, cand[:, 1] / 1e3, s=4, color="#3182ce",
                  label=f"candidate pads ({len(cand)})")
    ax[1].scatter(xy[:, 0] / 1e3, xy[:, 1] / 1e3, s=70, marker="*",
                  color="#c53030", edgecolors="black", linewidths=0.5,
                  label=f"chosen pads ({len(xy)})")
    ax[1].set_title(f"buildable ground (slope < {spec.max_slope_deg:.0f} deg) "
                    f"and the selected layout", fontsize=9)
    ax[1].legend(fontsize=7, loc="upper right")
    for a in ax:
        a.set_xlabel("east (km)")
        a.set_aspect("equal")
    ax[0].set_ylabel("north (km)")
    fig.tight_layout()
    fig.savefig(OUT / f"fig_{name}_terrain.png", dpi=130)
    plt.close(fig)


def figure_uv(result, name: str) -> None:
    spec, grid, obs, xy = result["spec"], result["grid"], result["obs"], result["xy"]
    uv = layout_to_uv_tracks(xy, obs)
    occ, _ = uv_occupancy(uv.reshape(-1, 2), grid)
    beam, _ = dirty_beam(uv.reshape(-1, 2), grid)
    r = result["report"]

    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.4))
    # layout_to_uv_tracks already includes the Hermitian conjugate points
    ax[0].plot(uv[..., 0].ravel() / 1e3, uv[..., 1].ravel() / 1e3, ",",
               color="#2b6cb0")
    ax[0].set_xlabel("u (kilo-lambda)")
    ax[0].set_ylabel("v (kilo-lambda)")
    ax[0].set_aspect("equal")
    ax[0].set_title(f"UV coverage, {spec.frequency_hz/1e9:.0f} GHz, "
                    f"{2*spec.hour_angle_h:.0f} h track", fontsize=9)

    ax[1].imshow(occ.T, origin="lower", cmap="magma")
    ax[1].set_title(f"gridded sampling -- {r['uv_cells_occupied']} cells "
                    f"({100*r['uv_fill_fraction']:.1f} % filled)", fontsize=9)
    ax[1].set_xticks([])
    ax[1].set_yticks([])

    n = beam.shape[0]
    c = slice(n // 2 - n // 8, n // 2 + n // 8)
    ax[2].imshow(beam[c, c].T, origin="lower", cmap="viridis")
    ax[2].set_title(f"dirty beam (central eighth) -- "
                    f"{r['resolution_arcsec']:.3f} arcsec, "
                    f"peak sidelobe {r['peak_sidelobe']:.3f}", fontsize=9)
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / f"fig_{name}_uv.png", dpi=130)
    plt.close(fig)


def figure_comparison(results) -> None:
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
    names = [r["report"]["site"].split(",")[0][:22] for r in results]
    colours = ["#c53030", "#2f855a"]

    pwv = [100 * r["report"]["pwv_below_1mm_fraction"] for r in results]
    alma = 100 * site_data.site("alma")["pwv_below_1mm_fraction"]
    ax[0].bar(names, pwv, color=colours)
    ax[0].axhline(alma, ls="--", color="#2b6cb0", lw=1.3)
    ax[0].text(0.98, alma + 1, "ALMA", ha="right", fontsize=8, color="#2b6cb0",
               transform=ax[0].get_yaxis_transform())
    ax[0].set_ylabel("% of time PWV <= 1 mm")
    ax[0].set_title("Measured dryness (ERA5, 2010-2025)", fontsize=9)
    ax[0].tick_params(axis="x", labelsize=7)

    for r, c, nm in zip(results, colours, names):
        xy = r["xy"]
        ax[1].scatter(xy[:, 0] / 1e3, xy[:, 1] / 1e3, s=45, color=c, label=nm,
                      alpha=0.85, edgecolors="black", linewidths=0.4)
    ax[1].set_aspect("equal")
    ax[1].set_xlabel("east (km)")
    ax[1].set_ylabel("north (km)")
    ax[1].set_title("Extended configurations, terrain-constrained", fontsize=9)
    ax[1].legend(fontsize=7)

    relief = [r["report"]["elevation_spread_m"] for r in results]
    ax[2].bar(names, relief, color=colours)
    ax[2].set_ylabel("ground relief across the array (m)")
    ax[2].set_title("Terrain penalty (lower is better)", fontsize=9)
    ax[2].tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_site_comparison.png", dpi=130)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(exist_ok=True)

    print("=" * 76)
    print("Ladakh submillimetre array design")
    print("=" * 76)
    base = DesignSpec().as_dict()
    for k in ("n_antennas", "dish_diameter_m", "frequency_ghz", "wavelength_mm",
              "min_pad_separation_m", "declination_deg", "track_hours",
              "max_slope_deg", "objective"):
        print(f"  {k:28} {base[k]}")

    print("")
    print("Measured dryness (fraction of time PWV <= 1 mm):")
    for row in site_data.pwv_comparison():
        print(f"  {row['name'][:46]:48} {100*row['fraction_below_1mm']:5.0f} %"
              f"   {row['infrastructure']}")

    rows, results = [], {}
    for key in SITES:
        for cname, overrides in CONFIGS:
            spec = DesignSpec(**overrides)
            tag = f"{key}_{cname}"
            print("")
            print(f"--- {key} / {cname} configuration " + "-" * 34)
            r = design.design_site(key, spec)
            rep = r["report"]
            rep["configuration"] = cname
            rep["site_key"] = key
            print(f"  candidates   {rep['candidates']} buildable positions "
                  f"({spec.candidate_spacing_m:.0f} m lattice, "
                  f"{spec.array_radius_m:.0f} m radius)")
            print(f"  baselines    {rep['baseline_min_m']:.0f} to "
                  f"{rep['baseline_max_m']:.0f} m, median "
                  f"{rep['baseline_median_m']:.0f} m")
            print(f"  resolution   {rep['resolution_arcsec']:.3f} arcsec")
            print(f"  largest scale {rep['largest_angular_scale_arcsec']:.2f} arcsec")
            print(f"  UV coverage  {rep['uv_cells_occupied']} cells "
                  f"({100*rep['uv_fill_fraction']:.1f} %)")
            print(f"  ground       {rep['elevation_min_m']:.0f}-"
                  f"{rep['elevation_max_m']:.0f} m, relief "
                  f"{rep['elevation_spread_m']:.0f} m")

            write_pads_csv(OUT / f"{tag}_pads.csv", r)
            design.write_cfg(OUT / f"{tag}_array.cfg", r["xy"], r["elev"],
                             spec, r["site"])
            (OUT / f"{tag}_report.json").write_text(
                json.dumps({"spec": spec.as_dict(), "report": rep},
                           indent=2, default=str), encoding="utf-8")
            figure_terrain(r, tag)
            figure_uv(r, tag)
            rows.append(rep)
            results[tag] = r

    print("")
    print("=" * 76)
    print("What the two configurations achieve together")
    print("=" * 76)
    for key in SITES:
        c = results[f"{key}_compact"]["report"]
        e = results[f"{key}_extended"]["report"]
        span = c["largest_angular_scale_arcsec"] / e["resolution_arcsec"]
        print(f"  {c['site'][:46]}")
        print(f"    finest detail (extended)   {e['resolution_arcsec']:.3f} arcsec")
        print(f"    widest structure (compact) {c['largest_angular_scale_arcsec']:.2f} arcsec")
        print(f"    range of angular scales    {span:.0f}x")
        print(f"    worst ground relief        "
              f"{max(c['elevation_spread_m'], e['elevation_spread_m']):.0f} m")

    figure_comparison([results[f"{k}_extended"] for k in SITES])
    with (OUT / "site_comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        keys = ["site", "site_key", "configuration", "site_latitude_deg",
                "pwv_below_1mm_fraction", "n_antennas", "baselines",
                "baseline_min_m", "baseline_max_m", "resolution_arcsec",
                "largest_angular_scale_arcsec", "uv_cells_occupied",
                "peak_sidelobe", "elevation_spread_m", "candidates"]
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("")
    print(f"All outputs in {OUT}")


if __name__ == "__main__":
    main()
