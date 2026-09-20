"""Command-line array designer.

    python -m thz_opt.cli --site hanle --n 16 --dish 8 --max-baseline 3000
    python -m thz_opt.cli --site alma --n 20 --mode select
    python -m thz_opt.cli --config runs/my_design.yaml
    python -m thz_opt.cli --list-sites

Two modes, because there are two different problems
---------------------------------------------------
``select``
    A pad field already exists and the question is which N of its M pads to
    occupy. This is ALMA's problem, and it is the one the QUBO formulation in
    this repository addresses.

``place``
    Nothing is built and the question is where pads should go. Candidate
    positions are generated on real terrain, filtered by slope and confined to
    one landform, and N are chosen from them. This is the Ladakh problem.

The mode is inferred from whether the site has a published pad list, and can be
overridden.

On penalty weights
------------------
The classical search has none, and the interface deliberately does not offer
any. Its move set swaps one pad in and one out, which preserves the antenna
count exactly, so the search never leaves the feasible set and there is no
cardinality constraint to penalise. Penalty weights exist only on the QUBO
route, and there they are *derived* rather than chosen: ``--show-penalties``
prints the values ``selection_penalty_floor`` and ``rosenberg_penalty_floor``
compute, which are the smallest that provably cannot be cheated. Offering them
as free inputs would invite a number that silently makes the formulation wrong.

Every run writes ``run.yaml`` recording every parameter actually used, so a
result can be reproduced without reconstructing the command line.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

from .arrays.real_arrays import REAL_SITES, load_cfg
from .constraints.separation import SeparationConfig, shadow_pairs
from .interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from .interferometry.psf import dirty_beam
from .interferometry.uv import UVGrid, uv_occupancy
from .metrics.psf_metrics import psf_metrics
from .optimize import (
    UVState,
    build_pair_tables,
    greedy_addition,
    greedy_removal,
    local_search,
    max_unique_cells,
    simulated_annealing,
)
from .qubo.baseline_qubo import track_cell_multiplicities
from .sites import site_data

ARCSEC = np.pi / (180.0 * 3600.0)
C = 299_792_458.0
REPO = Path(__file__).resolve().parents[2]

#: Sites with a published pad list, so `select` mode is possible.
PAD_FILES = {
    "alma": REPO / "data" / "external" / "alma.all.cfg",
    "aca": REPO / "data" / "external" / "aca.all.cfg",
    "vla": REPO / "data" / "external" / "vla.a.cfg",
}


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def describe_sites() -> str:
    lines = ["Sites with a published pad list (mode 'select'):"]
    for key, path in PAD_FILES.items():
        s = REAL_SITES.get(key, {})
        ok = "present" if path.exists() else "MISSING"
        lines.append(f"  {key:8} {s.get('name', key)[:44]:46} {ok}")
    lines.append("")
    lines.append("Sites with terrain only, nothing built (mode 'place'):")
    for key, s in site_data.LADAKH_SITES.items():
        pwv = s.get("pwv_below_1mm_fraction")
        lines.append(f"  {key:8} {s['name'][:44]:46} "
                     f"PWV<1mm {100*pwv:.0f}% of the time"
                     if pwv else f"  {key:8} {s['name']}")
    return "\n".join(lines)


def layout_metrics(xy, elev, lam, obs, grid, dish_m) -> dict:
    """Everything a design report should carry, computed once."""
    n = len(xy)
    dz = (np.asarray(elev)[:, None] - np.asarray(elev)[None, :]
          if elev is not None else np.zeros((n, n)))
    d = np.sqrt((xy[:, None, 0] - xy[None, :, 0]) ** 2
                + (xy[:, None, 1] - xy[None, :, 1]) ** 2 + dz ** 2)
    b = d[np.triu_indices(n, 1)]

    uv = layout_to_uv_tracks(xy, obs, heights=elev)
    occ, outside = uv_occupancy(uv.reshape(-1, 2), grid)
    beam, binfo = dirty_beam(uv.reshape(-1, 2), grid)
    m = psf_metrics(beam, binfo["pixel_scale_arcsec"])

    from .imaging.instrument import primary_beam_fwhm_arcsec

    return {
        "n_antennas": n,
        "n_baselines": n * (n - 1) // 2,
        "baseline_min_m": float(b.min()),
        "baseline_max_m": float(b.max()),
        "baseline_median_m": float(np.median(b)),
        "resolution_arcsec": float(lam / b.max() / ARCSEC),
        "largest_angular_scale_arcsec": float(0.6 * lam / b.min() / ARCSEC),
        "primary_beam_arcsec": primary_beam_fwhm_arcsec(dish_m, lam),
        "uv_cells_occupied": int(np.count_nonzero(occ)),
        "uv_fill_fraction": float(np.count_nonzero(occ) / occ.size),
        "uv_samples_outside_grid": int(outside),
        "peak_sidelobe": float(m.get("peak_sidelobe_level", np.nan)),
        "psf_fwhm_arcsec": float(m.get("fwhm_arcsec", np.nan)),
    }


def write_outputs(out: Path, xy, elev, args, site_info, metrics, extra=None) -> None:
    out.mkdir(parents=True, exist_ok=True)

    if site_info and site_info.get("latitude_deg") is not None:
        from .sites import terrain as _t
        lat, lon = _t.enu_to_latlon(xy[:, 0], xy[:, 1],
                                    site_info["latitude_deg"],
                                    site_info["longitude_deg"])
    else:
        lat = lon = np.full(len(xy), np.nan)

    radius = np.hypot(xy[:, 0], xy[:, 1])
    order = np.lexsort((np.arctan2(xy[:, 0], xy[:, 1]), np.round(radius, 3)))
    el = np.asarray(elev) if elev is not None else np.full(len(xy), np.nan)

    with (out / "pads.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pad_id", "latitude_deg", "longitude_deg", "elevation_m",
                    "east_m", "north_m", "radius_from_centre_m"])
        for rank, k in enumerate(order, start=1):
            w.writerow([f"P{rank:03d}", f"{lat[k]:.7f}", f"{lon[k]:.7f}",
                        f"{el[k]:.1f}", f"{xy[k, 0]:.2f}", f"{xy[k, 1]:.2f}",
                        f"{radius[k]:.1f}"])

    # pairwise distances, which is what a cable or transporter plan needs
    n = len(xy)
    with (out / "baselines.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pad_i", "pad_j", "distance_m"])
        rank_of = {int(k): r for r, k in enumerate(order, start=1)}
        for i in range(n):
            for j in range(i + 1, n):
                dij = float(np.hypot(xy[i, 0] - xy[j, 0], xy[i, 1] - xy[j, 1]))
                w.writerow([f"P{rank_of[i]:03d}", f"P{rank_of[j]:03d}",
                            f"{dij:.2f}"])

    lam = C / (args.freq * 1e9)
    lines = [
        f"# observatory={site_info.get('name', args.site) if site_info else args.site}",
        "# coordsys=LOC (local tangent plane)",
        f"# frequency={args.freq:.1f} GHz  dish={args.dish:.1f} m  N={len(xy)}",
        "# PRELIMINARY GEOMETRIC DESIGN - not a construction plan.",
        "# x y z diam pad#",
    ]
    datum = float(np.nanmean(el)) if np.isfinite(el).any() else 0.0
    for rank, k in enumerate(order, start=1):
        z = (el[k] - datum) if np.isfinite(el[k]) else 0.0
        lines.append(f"{xy[k, 0]:.3f}\t{xy[k, 1]:.3f}\t{z:.3f}\t"
                     f"{args.dish:.1f}\tP{rank:03d}")
    (out / "array.cfg").write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = {"parameters": vars(args) | {"wavelength_mm": lam * 1e3},
              "site": site_info, "metrics": metrics, **(extra or {})}
    (out / "report.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8")

    import yaml
    (out / "run.yaml").write_text(
        yaml.safe_dump(vars(args), sort_keys=True), encoding="utf-8")


# --------------------------------------------------------------------------
# the two modes
# --------------------------------------------------------------------------

def run_select(args, obs, sep):
    """Choose N of an existing pad list."""
    path = Path(args.pads) if args.pads else PAD_FILES.get(args.site)
    if path is None or not path.exists():
        raise SystemExit(f"no pad file for site {args.site!r}; pass --pads")

    site_info = dict(REAL_SITES.get(args.site, {"name": args.site}))
    cfg = load_cfg(path, dish_diameter_m=args.dish if args.filter_dish else None)
    pads = np.asarray(cfg.xy)[:, :2]
    # the .cfg carries a local up coordinate; adding the site's published
    # elevation turns it into a real height, so the UV projection is
    # three-dimensional rather than coplanar
    base = site_info.get("elevation_m")
    heights = (np.asarray(cfg.z, dtype=float) + base
               if base is not None else np.asarray(cfg.z, dtype=float))
    order = np.argsort(np.hypot(pads[:, 0], pads[:, 1]))
    if args.m and args.m < len(pads):
        pads, heights = pads[order[:args.m]], heights[order[:args.m]]
    print(f"  {len(pads)} candidate pads from {path.name}, "
          f"relief {heights.max()-heights.min():.0f} m")

    lam = C / (args.freq * 1e9)
    forbidden = shadow_pairs(pads, sep)
    uv = layout_to_uv_tracks(pads, obs, heights=heights)
    grid = UVGrid(uv_max=1.02 * float(np.abs(uv).max()), n_cells=args.uv_cells)
    print(f"  precomputing UV tracks for {len(pads)*(len(pads)-1)//2} baselines...")
    t0 = time.time()
    mult = track_cell_multiplicities(pads, lam, grid, obs, heights=heights)
    cells, counts = build_pair_tables(mult)
    print(f"  done in {time.time()-t0:.1f} s")

    state = UVState(cells, counts, len(pads), grid.n_cells ** 2)
    sel = search(state, args, forbidden)
    return pads[sel], heights[sel], site_info, lam, grid


def run_place(args, obs, sep):
    """Decide where pads go, on real terrain."""
    sys.path.insert(0, str(REPO / "ladakh_array_design"))
    import design as _design

    spec = _design.DesignSpec(
        n_antennas=args.n, dish_m=args.dish,
        max_baseline_m=args.max_baseline, frequency_hz=args.freq * 1e9,
        declination_deg=args.dec, hour_angle_h=args.hours / 2,
        separation_factor=args.separation_factor,
        max_slope_deg=args.max_slope, candidate_spacing_m=args.spacing,
        uv_cells=args.uv_cells, uv_taper=args.uv_taper,
        max_elevation_spread_m=args.max_relief, seed=args.seed,
    )
    print(f"  fetching terrain for {args.site} (cached after the first run)...")
    r = _design.design_site(args.site, spec)
    ci = r["report"]["array_centre"]
    print(f"  {r['report']['candidates']} buildable candidates; array centre "
          f"{ci.get('moved_m', 0):.0f} m from the site marker")
    return r["xy"], r["elev"], r["site"], spec.wavelength_m, r["grid"]


def search(state, args, forbidden):
    """Run the requested heuristic and return the selected indices."""
    score = max_unique_cells()
    t0 = time.time()
    if args.method in ("greedy-add", "auto"):
        res = greedy_addition(state, args.n, score, forbidden)
    elif args.method == "greedy-remove":
        res = greedy_removal(state, args.n, score, forbidden)
    elif args.method == "annealing":
        res = simulated_annealing(state, args.n, score, forbidden,
                                  n_iterations=args.iterations, seed=args.seed)
    else:
        raise SystemExit(f"unknown method {args.method!r}")
    if args.method != "annealing" and not args.no_polish:
        state.set_selection(res.selection)
        res = local_search(state, score, forbidden)
    print(f"  search: {args.method}"
          f"{'' if args.no_polish else ' + 1-swap local search'}"
          f" in {time.time()-t0:.1f} s")
    return np.flatnonzero(res.selection)


# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="thz-design", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list-sites", action="store_true",
                   help="show the available sites and exit")
    p.add_argument("--config", help="YAML file of these same options")
    p.add_argument("--site", default="alma", help="site key, see --list-sites")
    p.add_argument("--mode", choices=("auto", "select", "place"), default="auto",
                   help="select from an existing pad list, or place new pads")
    p.add_argument("--pads", help="a CASA .cfg pad file, for --mode select")

    p.add_argument("--n", type=int, default=20, help="antennas to place or select")
    p.add_argument("--m", type=int, default=None,
                   help="candidate pads to consider (innermost M), select mode")
    p.add_argument("--dish", type=float, default=12.0, help="dish diameter, m")
    p.add_argument("--max-baseline", type=float, default=3000.0,
                   help="longest baseline, m (place mode)")
    p.add_argument("--freq", type=float, default=230.0, help="frequency, GHz")
    p.add_argument("--dec", type=float, default=-30.0,
                   help="source declination, degrees")
    p.add_argument("--hours", type=float, default=4.0, help="track length, h")
    p.add_argument("--times", type=int, default=41, help="samples along the track")
    p.add_argument("--separation-factor", type=float, default=1.5,
                   help="minimum pad separation as a multiple of the dish")
    p.add_argument("--uv-cells", type=int, default=64, help="UV grid size")

    p.add_argument("--max-slope", type=float, default=10.0,
                   help="buildable slope limit, degrees (place mode)")
    p.add_argument("--spacing", type=float, default=120.0,
                   help="candidate lattice pitch, m (place mode)")
    p.add_argument("--max-relief", type=float, default=60.0,
                   help="elevation band the array must fit in, m (place mode)")
    p.add_argument("--uv-taper", type=float, default=0.5,
                   help="Gaussian UV density width as a fraction of q_max; "
                        "0 disables it and maximises bare cell count")

    p.add_argument("--method", default="auto",
                   choices=("auto", "greedy-add", "greedy-remove", "annealing"))
    p.add_argument("--iterations", type=int, default=20000,
                   help="annealing iterations")
    p.add_argument("--no-polish", action="store_true",
                   help="skip the 1-swap local search after a greedy pass")
    p.add_argument("--filter-dish", action="store_true",
                   help="keep only pads matching --dish, select mode")
    p.add_argument("--seed", type=int, default=20260919)
    p.add_argument("--quick", action="store_true",
                   help="snapshot instead of Earth rotation, and a coarse grid; "
                        "seconds instead of minutes, for exploring")
    p.add_argument("--show-penalties", action="store_true",
                   help="print the derived QUBO penalty floors and exit")
    p.add_argument("-o", "--out", default="runs/design",
                   help="output directory")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_sites:
        print(describe_sites())
        return 0

    if args.config:
        import yaml
        cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
        for k, v in cfg.items():
            if hasattr(args, k):
                setattr(args, k, v)
            else:
                print(f"  warning: ignoring unknown option {k!r} in {args.config}")

    if args.quick:
        args.times = 1
        args.uv_cells = min(args.uv_cells, 32)
        print("  quick mode: snapshot UV, coarse grid -- for exploring, not for "
              "a final design")

    if args.mode == "auto":
        args.mode = "select" if args.site in PAD_FILES else "place"

    print(f"\nthz-design: {args.mode} {args.n} antennas at {args.site}, "
          f"{args.dish:.0f} m dishes, {args.freq:.0f} GHz")

    lat = (REAL_SITES.get(args.site) or site_data.LADAKH_SITES.get(args.site)
           or {}).get("latitude_deg")
    if lat is None:
        raise SystemExit(f"unknown site {args.site!r}; try --list-sites")

    obs = ObservationConfig(
        frequency_hz=args.freq * 1e9, declination_deg=args.dec,
        latitude_deg=lat, hour_angle_start_h=-args.hours / 2,
        hour_angle_end_h=args.hours / 2, n_times=max(args.times, 1),
    )
    sep = SeparationConfig(args.dish, args.separation_factor)

    if args.mode == "select":
        xy, elev, site_info, lam, grid = run_select(args, obs, sep)
    else:
        xy, elev, site_info, lam, grid = run_place(args, obs, sep)

    metrics = layout_metrics(xy, elev, lam, obs, grid, args.dish)
    print(f"\n  {metrics['n_antennas']} antennas, {metrics['n_baselines']} baselines")
    print(f"  baselines      {metrics['baseline_min_m']:.0f} to "
          f"{metrics['baseline_max_m']:.0f} m "
          f"(median {metrics['baseline_median_m']:.0f})")
    print(f"  resolution     {metrics['resolution_arcsec']:.4f} arcsec")
    print(f"  largest scale  {metrics['largest_angular_scale_arcsec']:.2f} arcsec")
    print(f"  primary beam   {metrics['primary_beam_arcsec']:.1f} arcsec")
    print(f"  UV coverage    {metrics['uv_cells_occupied']} cells "
          f"({100*metrics['uv_fill_fraction']:.1f} %)")
    print(f"  peak sidelobe  {metrics['peak_sidelobe']:.3f}")

    extra = None
    if args.show_penalties:
        extra = report_penalties(xy, lam, grid, obs, args)

    out = Path(args.out).resolve()
    write_outputs(out, xy, np.asarray(elev) if elev is not None else None,
                  args, site_info, metrics, extra)
    print(f"\n  wrote {out}")
    for f in ("pads.csv", "baselines.csv", "array.cfg", "report.json", "run.yaml"):
        print(f"    {f}")
    return 0


def report_penalties(xy, lam, grid, obs, args) -> dict:
    """The derived QUBO penalty floors, which are computed and not chosen."""
    from .qubo.baseline_qubo import (
        bonferroni_coverage_terms,
        rosenberg_penalty_floor,
        selection_penalty_floor,
    )
    from .qubo.coefficients import baseline_cell_sets

    cs = baseline_cell_sets(xy, lam, grid, observation=obs)
    terms = bonferroni_coverage_terms(cs, len(xy))
    lam_r = rosenberg_penalty_floor(terms)
    lam_s = selection_penalty_floor(terms, args.n)
    n = len(xy)
    print(f"\n  derived QUBO penalty floors (not user-set):")
    print(f"    lambda_rosenberg  {lam_r:.4g}")
    print(f"    lambda_select     {lam_s:.4g}")
    print(f"    QUBO variables    {n + n*(n-1)//2} "
          f"({n} pad + {n*(n-1)//2} baseline)")
    print(f"    These are the smallest weights that provably cannot be cheated;")
    print(f"    a smaller value makes the formulation wrong, which is why they")
    print(f"    are computed rather than offered as an input.")
    return {"qubo_penalties": {"lambda_rosenberg": float(lam_r),
                               "lambda_select": float(lam_s),
                               "n_variables": n + n * (n - 1) // 2}}


if __name__ == "__main__":
    raise SystemExit(main())
