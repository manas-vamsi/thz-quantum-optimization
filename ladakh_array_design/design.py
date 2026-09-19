"""Design a submillimetre array for a Ladakh site: where each pad goes.

The difference from the rest of this repository
-----------------------------------------------
Everywhere else the problem is *selection*: ALMA has 174 concrete pads and the
question is which twenty to occupy. Here nothing is built, so the problem is
*placement*: decide where pads should go on real ground. That is the easier
direction mathematically, and it is the one that produces a deliverable a
surveyor can act on.

The method is still discrete selection, because that is what the validated code
in ``thz_opt`` does and because it is the honest way to respect terrain: the
buildable ground is sampled into a candidate grid, and the optimiser chooses N
of those candidates. A continuous optimiser would return coordinates that look
more precise while quietly ignoring whether the ground at that point is flat.

Design parameters and why
-------------------------
=======================  ==============================================
Antennas, N = 16         120 baselines, enough for real imaging, and a
                         credible first national array.
Dish diameter, 8 m       Buildable domestically; sets the minimum pad
                         separation at 1.5 x D = 12 m, the shadowing rule
                         from the project note.
Maximum baseline, 3 km   Gives 0.09 arcsec at 230 GHz.
Frequency, 230 GHz       The only band with measured on-site transparency
                         at Hanle: a 220 GHz tipping radiometer has run
                         there since 1999. Designing for a band with no
                         local measurement would be guesswork.
Declination, +30 deg     Transits within 3 degrees of zenith at these
                         latitudes, so the Earth-rotation track is the
                         best case the site can offer.
Track, 4 hours           +/- 2 h of hour angle about transit.
=======================  ==============================================

Honest status of the output
---------------------------
This is a preliminary geometric design: pad coordinates, baselines, UV
coverage, point-spread function and resolution. It is **not** a construction
plan. There is no geotechnical survey, no land-rights check, no access-road or
power routing, and no environmental assessment. The terrain filter establishes
only that the ground at each pad is flat enough, from a 16 m digital elevation
model.

The atmospheric coherence figure is a **transfer** from ALMA's measured phase
structure function, scaled by nothing: no phase-stability measurement exists
for any Ladakh site, because that needs an operating interferometer. It is
reported as a model and must not be quoted as a site measurement.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from thz_opt.constraints.coherence import CoherenceConfig, coherence_summary
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.psf import dirty_beam
from thz_opt.interferometry.uv import UVGrid, uv_occupancy
from thz_opt.metrics.psf_metrics import psf_metrics
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_addition,
    local_search,
    max_unique_cells,
)
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities

import site_data
import terrain

__all__ = ["DesignSpec", "candidate_pads", "optimise_layout", "design_site",
           "write_cfg", "layout_report", "order_pads", "pad_table"]

C = 299_792_458.0
ARCSEC = np.pi / (180.0 * 3600.0)


class DesignSpec:
    """Everything that defines the array being designed."""

    def __init__(self, n_antennas=16, dish_m=8.0, max_baseline_m=3000.0,
                 frequency_hz=230e9, declination_deg=30.0,
                 hour_angle_h=2.0, n_times=41, separation_factor=1.5,
                 max_slope_deg=10.0, candidate_spacing_m=120.0,
                 uv_cells=64, uv_taper=0.5, max_elevation_spread_m=60.0,
                 seed=20260919):
        self.n_antennas = int(n_antennas)
        self.dish_m = float(dish_m)
        self.max_baseline_m = float(max_baseline_m)
        self.frequency_hz = float(frequency_hz)
        self.declination_deg = float(declination_deg)
        self.hour_angle_h = float(hour_angle_h)
        self.n_times = int(n_times)
        self.separation_factor = float(separation_factor)
        self.max_slope_deg = float(max_slope_deg)
        self.candidate_spacing_m = float(candidate_spacing_m)
        self.uv_cells = int(uv_cells)
        # None reproduces the pure cell-count objective, which yields a ring;
        # a float is the Gaussian target-density width as a fraction of q_max.
        self.uv_taper = None if uv_taper is None else float(uv_taper)
        # Candidates are confined to one landform of this vertical extent;
        # None disables the filter and allows an array split across a cliff.
        self.max_elevation_spread_m = (None if max_elevation_spread_m is None
                                       else float(max_elevation_spread_m))
        self.seed = int(seed)

    @property
    def wavelength_m(self) -> float:
        return C / self.frequency_hz

    @property
    def array_radius_m(self) -> float:
        """Longest baseline is a diameter, so the field radius is half of it."""
        return 0.5 * self.max_baseline_m

    @property
    def min_separation_m(self) -> float:
        return self.separation_factor * self.dish_m

    @property
    def resolution_arcsec(self) -> float:
        """lambda / b_max, the diffraction limit the design targets."""
        return self.wavelength_m / self.max_baseline_m / ARCSEC

    def as_dict(self) -> dict:
        return {
            "n_antennas": self.n_antennas, "dish_diameter_m": self.dish_m,
            "max_baseline_m": self.max_baseline_m,
            "frequency_ghz": self.frequency_hz / 1e9,
            "wavelength_mm": self.wavelength_m * 1e3,
            "target_resolution_arcsec": round(self.resolution_arcsec, 4),
            "min_pad_separation_m": self.min_separation_m,
            "declination_deg": self.declination_deg,
            "track_hours": 2 * self.hour_angle_h,
            "max_slope_deg": self.max_slope_deg,
            "candidate_spacing_m": self.candidate_spacing_m,
            "uv_taper_sigma_fraction": self.uv_taper,
            "max_elevation_spread_m": self.max_elevation_spread_m,
            "objective": ("gaussian target UV density"
                          if self.uv_taper else "maximum unique UV cells"),
        }


def candidate_pads(dem, spec: DesignSpec):
    """Buildable ground sampled into a candidate pad grid.

    Cells are taken on a coarse lattice so that candidates are already further
    apart than the shadowing limit, then filtered by slope and radius. The
    lattice is offset to the site origin so the array centre is always a
    candidate.
    """
    mask = terrain.buildable_mask(dem, spec.max_slope_deg, spec.array_radius_m)
    step = max(1, int(round(spec.candidate_spacing_m / dem.pixel_m)))

    # centre the lattice on the origin cell
    i0 = int(np.argmin(np.abs(dem.north)))
    j0 = int(np.argmin(np.abs(dem.east)))
    ii = np.arange(i0 % step, mask.shape[0], step)
    jj = np.arange(j0 % step, mask.shape[1], step)

    sel_i, sel_j = np.meshgrid(ii, jj, indexing="ij")
    keep = mask[sel_i, sel_j]
    east = dem.east[sel_j][keep]
    north = dem.north[sel_i][keep]
    elev = dem.elevation[sel_i, sel_j][keep]
    xy = np.column_stack([east, north])

    if spec.max_elevation_spread_m is not None and len(xy):
        xy, elev = _dominant_landform(xy, elev, spec.max_elevation_spread_m)
    return xy, elev


def _dominant_landform(xy, elev, band_m: float):
    """Keep only candidates on one landform, the one with the most of them.

    Flat ground is not one surface. At Hanle the observatory stands on an
    isolated summit about 240 m above a wide plain, and both are locally flat,
    so a pure slope filter accepts them equally. An earlier run then placed
    fourteen pads on the plain and one on the summit, and that single outlier
    accounted for the whole 235 m of reported relief -- an array that would
    need a 240 m climb between two of its antennas.

    Selecting the elevation band of width ``band_m`` containing the most
    candidates picks the dominant landform and discards the rest. It is a crude
    proxy for "one contiguous buildable area", but it removes the failure mode
    that matters: an array split across a cliff.
    """
    elev = np.asarray(elev, dtype=float)
    lo_edges = np.unique(elev)
    best_lo, best_n = lo_edges[0], -1
    for lo in lo_edges:
        n = int(np.sum((elev >= lo) & (elev <= lo + band_m)))
        if n > best_n:
            best_lo, best_n = lo, n
    keep = (elev >= best_lo) & (elev <= best_lo + band_m)
    return xy[keep], elev[keep]


def uv_taper_weights(grid, sigma_fraction: float) -> np.ndarray:
    """Gaussian target UV density, as a per-cell weight.

    Why a bare cell count is the wrong objective for a real array
    -------------------------------------------------------------
    Maximising the number of *distinct* occupied cells sounds neutral but is
    not: the outer annuli of a UV grid contain far more cells than the inner
    ones, because area grows as the radius. The optimiser therefore pushes
    every antenna to the boundary of the allowed field. The first design run
    here did exactly that -- all sixteen pads landed between 1064 m and 1488 m
    of a 1500 m field, none in the inner half -- producing a ring that resolves
    beautifully and cannot see anything wider than 0.72 arcsec.

    A ring is a real array design, but only for a source you already know is
    compact. A general-purpose instrument needs short spacings too, and the
    standard prescription is to match a Gaussian target density in UV radius
    (Boone 2001), which yields a well-behaved beam with low sidelobes rather
    than the sharp ring transform. Weighting each cell by

        w(q) = exp(-q^2 / (2 (f * q_max)^2))

    makes inner cells worth more, so the same optimiser now trades a little
    resolution for the short baselines an image actually needs. ``f`` is the
    only free number and is reported with the design.
    """
    n = grid.n_cells
    axis = (np.arange(n) + 0.5) * (2.0 * grid.uv_max / n) - grid.uv_max
    uu, vv = np.meshgrid(axis, axis, indexing="ij")
    q = np.hypot(uu, vv).ravel()
    sigma = max(float(sigma_fraction) * grid.uv_max, 1e-9)
    return np.exp(-0.5 * (q / sigma) ** 2)


def optimise_layout(pads_xy, spec: DesignSpec, latitude_deg: float):
    """Choose ``n_antennas`` candidate positions maximising UV coverage."""
    obs = ObservationConfig(
        frequency_hz=spec.frequency_hz,
        declination_deg=spec.declination_deg,
        latitude_deg=latitude_deg,
        hour_angle_start_h=-spec.hour_angle_h,
        hour_angle_end_h=spec.hour_angle_h,
        n_times=spec.n_times,
    )
    sep = SeparationConfig(spec.dish_m, spec.separation_factor)
    forbidden = shadow_pairs(pads_xy, sep)

    uv = layout_to_uv_tracks(pads_xy, obs)
    grid = UVGrid(uv_max=1.02 * float(np.abs(uv).max()), n_cells=spec.uv_cells)
    mult = track_cell_multiplicities(pads_xy, spec.wavelength_m, grid, obs)
    cells, counts = build_pair_tables(mult)

    if spec.uv_taper is None:
        state = UVState(cells, counts, len(pads_xy), grid.n_cells ** 2)
        score = max_unique_cells()
    else:
        w = uv_taper_weights(grid, spec.uv_taper)
        state = UVState(cells, counts, len(pads_xy), grid.n_cells ** 2,
                        cell_weights=w)

        def score(st) -> float:
            return -float(st.weighted_cells)
        score.name = "gaussian_density_match"
    t0 = time.time()
    # Addition, not removal. Greedy removal costs O((M-N)*M) evaluations and
    # each one touches every baseline of a nearly full array; at the candidate
    # counts a terrain grid produces (several hundred) that is minutes of work
    # for a worse starting point than local search can fix cheaply. Addition is
    # O(N*M) with a small working set, and the swap pass repairs its habit of
    # committing early.
    res = greedy_addition(state, spec.n_antennas, score, forbidden)
    state.set_selection(res.selection)
    res = local_search(state, score, forbidden)
    return np.asarray(res.selection, dtype=bool), grid, obs, time.time() - t0


def order_pads(xy, elev):
    """Sort pads once, by distance from the array centre, ties broken by bearing.

    Every output must label the same physical pad ``P001``. An earlier version
    sorted inside the CSV writer but not the CFG writer, so the two files gave
    the same identifier to different pads -- the sort of mistake that reaches a
    survey team. The ordering happens once, here, and both writers consume it.
    """
    xy = np.asarray(xy, dtype=float)
    elev = np.asarray(elev, dtype=float)
    radius = np.hypot(xy[:, 0], xy[:, 1])
    bearing = np.arctan2(xy[:, 0], xy[:, 1])
    order = np.lexsort((bearing, np.round(radius, 3)))
    return xy[order], elev[order], radius[order]


def pad_table(xy, elev, site_info) -> list:
    """The pad list, canonically ordered, in both coordinate systems."""
    xy, elev, radius = order_pads(xy, elev)
    lat, lon = terrain.enu_to_latlon(xy[:, 0], xy[:, 1],
                                     site_info["latitude_deg"],
                                     site_info["longitude_deg"])
    return [
        {"pad_id": f"P{k+1:03d}", "latitude_deg": float(lat[k]),
         "longitude_deg": float(lon[k]), "elevation_m": float(elev[k]),
         "east_m": float(xy[k, 0]), "north_m": float(xy[k, 1]),
         "radius_from_centre_m": float(radius[k])}
        for k in range(len(xy))
    ]


def layout_report(xy, elev, spec: DesignSpec, grid, obs, site_info) -> dict:
    """Everything an engineer or referee would ask about the chosen layout."""
    n = len(xy)
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    iu = np.triu_indices(n, 1)
    b = d[iu]

    uv = layout_to_uv_tracks(xy, obs)
    occ, n_outside = uv_occupancy(uv.reshape(-1, 2), grid)
    beam, beam_info = dirty_beam(uv.reshape(-1, 2), grid)
    pm = psf_metrics(beam, beam_info["pixel_scale_arcsec"])

    # coherence is a MODEL here: no Ladakh phase measurement exists
    coh = coherence_summary(xy, spec.wavelength_m,
                            CoherenceConfig.from_measured("wvr_corrected"))

    return {
        "n_antennas": n,
        "baselines": n * (n - 1) // 2,
        "baseline_min_m": float(b.min()),
        "baseline_max_m": float(b.max()),
        "baseline_median_m": float(np.median(b)),
        "resolution_arcsec": float(spec.wavelength_m / b.max() / ARCSEC),
        "largest_angular_scale_arcsec": float(0.6 * spec.wavelength_m / b.min() / ARCSEC),
        "psf_fwhm_arcsec": float(pm.get("fwhm_arcsec", np.nan)),
        "peak_sidelobe": float(pm.get("peak_sidelobe_level", np.nan)),
        "uv_cells_occupied": int(np.count_nonzero(occ)),
        "uv_fill_fraction": float(np.count_nonzero(occ) / occ.size),
        "uv_samples_outside_grid": int(n_outside),
        "elevation_min_m": float(elev.min()),
        "elevation_max_m": float(elev.max()),
        "elevation_spread_m": float(elev.max() - elev.min()),
        "pads_inside_half_radius": int(np.sum(np.hypot(xy[:, 0], xy[:, 1])
                                              < 0.5 * spec.array_radius_m)),
        "coherence_model_mean": coh["coherence_mean"],
        "coherence_model_note": ("transferred from ALMA Memo 624; no phase "
                                 "measurement exists for any Ladakh site"),
        "site": site_info["name"],
        "site_latitude_deg": site_info["latitude_deg"],
        "pwv_below_1mm_fraction": site_info.get("pwv_below_1mm_fraction"),
    }


def write_cfg(path: Path, xy, elev, spec: DesignSpec, site_info) -> Path:
    """A CASA observatory configuration file, same format as ``alma.all.cfg``.

    This is the deliverable that makes the design usable by anyone else: it can
    be dropped straight into CASA's ``simobserve`` to simulate observations
    with this array.
    """
    rows = pad_table(xy, elev, site_info)
    # Height datum is the array's own mean elevation, not the site's published
    # figure: at Hanle the observatory building stands about 230 m above the
    # flat ground the array occupies, so referencing to it would record every
    # pad at a spurious -230 m.
    datum = float(np.mean([r["elevation_m"] for r in rows]))
    lines = [
        f"# observatory={site_info['name']}",
        "# coordsys=LOC (local tangent plane)",
        f"# origin lat={site_info['latitude_deg']:.6f} "
        f"lon={site_info['longitude_deg']:.6f}",
        f"# height datum={datum:.1f} m (mean elevation of these pads)",
        f"# frequency={spec.frequency_hz/1e9:.0f} GHz  "
        f"dish={spec.dish_m:.1f} m  N={spec.n_antennas}",
        "# PRELIMINARY GEOMETRIC DESIGN - not a construction plan.",
        "# Terrain from AWS Terrain Tiles (SRTM/ASTER); slope filter only.",
        "# Pad ids match the accompanying _pads.csv row for row.",
        "# x y z diam pad#",
    ]
    for r in rows:
        lines.append(f"{r['east_m']:.3f}\t{r['north_m']:.3f}\t"
                     f"{r['elevation_m'] - datum:.3f}\t{spec.dish_m:.1f}\t"
                     f"{r['pad_id']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def design_site(site_key: str, spec: DesignSpec | None = None,
                half_size_m: float | None = None, verbose: bool = True):
    """Full pipeline for one site: terrain, candidates, optimise, report."""
    spec = spec or DesignSpec()
    info = site_data.site(site_key)
    half = half_size_m or (spec.array_radius_m + 1000.0)

    dem = terrain.fetch_dem(info["latitude_deg"], info["longitude_deg"], half)
    if info.get("elevation_m") is None:
        info["elevation_m"] = dem.at(0.0, 0.0)

    cand_xy, cand_elev = candidate_pads(dem, spec)
    if len(cand_xy) < spec.n_antennas:
        raise RuntimeError(f"only {len(cand_xy)} buildable candidates at "
                           f"{site_key}; loosen the slope or spacing limits")

    sel, grid, obs, secs = optimise_layout(cand_xy, spec, info["latitude_deg"])
    xy, elev = cand_xy[sel], cand_elev[sel]
    report = layout_report(xy, elev, spec, grid, obs, info)
    report["candidates"] = int(len(cand_xy))
    report["optimise_seconds"] = round(secs, 1)
    report["terrain"] = dem.summary()
    return {"site": info, "spec": spec, "dem": dem, "candidates": cand_xy,
            "candidate_elev": cand_elev, "xy": xy, "elev": elev,
            "grid": grid, "obs": obs, "report": report}
