"""Experiment 15 -- what the primary beam, smearing and real noise cost.

Three effects were missing from every earlier result in this repository, and
all three get worse exactly where this project pushes: long baselines at short
wavelengths. This experiment measures what each one does to the usable field of
view and to the observing time required, on the real ALMA pad geometry and with
the measured atmospheric opacity.

The questions
-------------
1. **How wide is the usable field?** Three separate limits apply -- the primary
   beam of one dish, bandwidth smearing, and time smearing -- and whichever is
   smallest is the answer. Which one binds depends on the array, so it has to
   be computed rather than assumed.
2. **What do smearing losses do to an image?** Measured by imaging the same sky
   with and without them.
3. **How long would this take?** Converted from arbitrary signal-to-noise into
   seconds using the measured zenith opacity, so the answer is quotable.

The headline tension
--------------------
Resolution wants long baselines. Field of view wants short ones, because both
smearing terms scale with the distance from the phase centre measured *in
synthesised beams*, and the beam shrinks as the array grows. Extending an array
therefore shrinks its usable field quadratically in a sense the coverage
metrics used elsewhere cannot express at all.

Outputs
-------
figures/fig31_field_of_view_limits.png
figures/fig32_smearing_cost.png
data/results/exp15_instrumental.csv
data/results/exp15_summary.json
"""

from __future__ import annotations

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

from thz_opt.arrays.real_arrays import REAL_SITES, load_cfg
from thz_opt.constraints.separation import SeparationConfig, shadow_pairs
from thz_opt.imaging import (
    apply_primary_beam,
    disc_with_gap_image,
    fidelity_metrics,
    gaussian_source,
    hogbom_clean,
    observe,
    point_sources,
    primary_beam_fwhm_arcsec,
    radiometer_noise,
    restore,
    sefd,
    smear_image,
    system_temperature,
)
from thz_opt.interferometry.earth_rotation import ObservationConfig, layout_to_uv_tracks
from thz_opt.interferometry.uv import UVGrid
from thz_opt.optimize import (
    UVState,
    build_pair_tables,
    greedy_addition,
    local_search,
    max_unique_cells,
)
from thz_opt.qubo.baseline_qubo import track_cell_multiplicities

ARCSEC = np.pi / (180.0 * 3600.0)
FREQ_HZ = 3.0e11
DISH_M = 12.0
N_SELECT = 20
GRID_CELLS = 128
SEED = 20260919
CFG_FILE = ROOT / "data" / "external" / "alma.all.cfg"

#: Observing setups, from a narrow spectral line to wide continuum.
SETUPS = [
    ("spectral line", dict(fractional_bandwidth=0.0005, integration_s=2.0)),
    ("narrow continuum", dict(fractional_bandwidth=0.01, integration_s=6.0)),
    ("wide continuum", dict(fractional_bandwidth=0.08, integration_s=30.0)),
]

#: Measured zenith opacities. The Hanle figure is from its on-site radiometer.
OPACITIES = [
    ("Hanle, best winter decile", 0.06),
    ("Hanle, typical", 0.15),
    ("Chajnantor, 1 mm PWV", 0.10),
]


def smearing_field_arcsec(b_max_m: float, wavelength_m: float,
                          fractional_bandwidth: float, integration_s: float,
                          declination_deg: float, loss: float = 0.1) -> dict:
    """Radius at which each smearing term costs a given fraction of the peak.

    Both smears displace a source by a distance proportional to its offset from
    the phase centre. Requiring that displacement to stay below ``loss`` of a
    synthesised beam gives a radius directly:

        bandwidth:  theta < loss * theta_beam / (delta_nu / nu)
        time:       theta < loss * theta_beam / (omega_E cos(dec) tau)

    The synthesised beam is ``lambda / b_max``, so both radii shrink as the
    array is extended even though the resolution improves.
    """
    from thz_opt.imaging.instrument import OMEGA_EARTH

    beam = wavelength_m / b_max_m / ARCSEC
    bw = (loss * beam / fractional_bandwidth
          if fractional_bandwidth > 0 else np.inf)
    omega = OMEGA_EARTH * abs(np.cos(np.deg2rad(declination_deg)))
    tm = loss * beam / (omega * integration_s) if omega * integration_s > 0 else np.inf
    return {"beam_arcsec": beam, "bandwidth_arcsec": bw, "time_arcsec": tm}


def main() -> None:
    ensure_dirs()
    cfg = load_config()
    site = REAL_SITES["alma"]
    if not CFG_FILE.exists():
        raise SystemExit(f"missing real pad file: {CFG_FILE}")

    lam = 299_792_458.0 / FREQ_HZ
    pads_cfg = load_cfg(CFG_FILE, dish_diameter_m=DISH_M)
    all_pads = pads_cfg.xy
    order = np.argsort(np.hypot(all_pads[:, 0], all_pads[:, 1]))
    pads = all_pads[order[:80]]

    obs = ObservationConfig(
        frequency_hz=FREQ_HZ, declination_deg=-30.0,
        latitude_deg=site["latitude_deg"],
        hour_angle_start_h=-2.0, hour_angle_end_h=2.0, n_times=41,
    )
    sep = SeparationConfig(DISH_M, cfg["constraints"]["minimum_separation_factor"])
    forbidden = shadow_pairs(pads, sep)

    ps = (f"real ALMA pads, N={N_SELECT}, {DISH_M:.0f} m dishes, "
          f"{FREQ_HZ/1e9:.0f} GHz, 4 h track")
    banner(f"Experiment 15 -- instrumental limits -- {ps}")

    # ---- an array to test with -------------------------------------------
    uv_all = layout_to_uv_tracks(pads, obs)
    grid = UVGrid(uv_max=1.02 * float(np.abs(uv_all).max()), n_cells=GRID_CELLS)
    mult = track_cell_multiplicities(pads, lam, grid, obs)
    cells_t, counts_t = build_pair_tables(mult)
    state = UVState(cells_t, counts_t, len(pads), grid.n_cells ** 2)
    score = max_unique_cells()
    res = greedy_addition(state, N_SELECT, score, forbidden)
    state.set_selection(res.selection)
    res = local_search(state, score, forbidden)
    sel = np.flatnonzero(res.selection)
    xy = pads[sel]
    uv = layout_to_uv_tracks(xy, obs)

    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    b_max = float(d[np.triu_indices(len(xy), 1)].max())
    pb = primary_beam_fwhm_arcsec(DISH_M, lam)
    pix_arcsec = 1.0 / (grid.n_cells * grid.cell_size) / ARCSEC

    print(f"\n  [1] the array")
    print(f"      {N_SELECT} pads, longest baseline {b_max/1e3:.1f} km, "
          f"{state.unique_cells} UV cells")
    print(f"      synthesised beam  {lam/b_max/ARCSEC:.4f} arcsec")
    print(f"      primary beam      {pb:.1f} arcsec (a {DISH_M:.0f} m dish)")
    print(f"      image grid        {GRID_CELLS}^2 of {pix_arcsec:.4f} arcsec "
          f"= {GRID_CELLS*pix_arcsec:.2f} arcsec")

    # ---- 2. which limit binds? -------------------------------------------
    print(f"\n  [2] usable field radius, at a 10 % peak loss")
    print(f"      {'setup':20}{'bandwidth':>12}{'time':>12}"
          f"{'primary/2':>12}{'binding':>18}")
    rows = []
    for name, kw in SETUPS:
        f = smearing_field_arcsec(b_max, lam, kw["fractional_bandwidth"],
                                  kw["integration_s"], obs.declination_deg)
        limits = {"bandwidth smearing": f["bandwidth_arcsec"],
                  "time smearing": f["time_arcsec"],
                  "primary beam": pb / 2}
        binding = min(limits, key=limits.get)
        rows.append({"setup": name, **kw, "beam_arcsec": f["beam_arcsec"],
                     "bandwidth_field_arcsec": f["bandwidth_arcsec"],
                     "time_field_arcsec": f["time_arcsec"],
                     "primary_half_arcsec": pb / 2,
                     "binding_limit": binding,
                     "usable_field_arcsec": limits[binding]})
        print(f"      {name:20}{f['bandwidth_arcsec']:>11.2f}\""
              f"{f['time_arcsec']:>11.2f}\"{pb/2:>11.2f}\"{binding:>18}")

    print(f"\n      Resolution and field of view pull against each other: both")
    print(f"      smearing radii are proportional to the synthesised beam, so")
    print(f"      extending the array shrinks the usable field in proportion.")

    # ---- 3. what smearing costs across the field -------------------------
    #
    # A compact source at the phase centre is immune to both effects by
    # construction, so imaging one measures nothing: an earlier version of this
    # section used the 0.5 arcsec disc and reported a 1 % change across every
    # setup. The effects are position-dependent, so the probe has to be spread
    # across the field. A ring of identical point sources at increasing radius
    # measures the response directly.
    print("")
    print("  [3] response across the field (identical sources, varying radius)")
    radii_arcsec = np.array([0.0, 1.0, 2.5, 5.0, 10.0, 15.0])
    offsets = [(r, 0.0) for r in radii_arcsec]
    probe = point_sources(grid, offsets, [1.0] * len(offsets))
    centre = grid.n_cells // 2
    px = [int(round(r / pix_arcsec)) for r in radii_arcsec]

    print("      " + f"{'setup':20}{'primary beam':>13}" +
          "".join(f"{r:>8.1f}\"" for r in radii_arcsec))
    img_rows = []
    for name, kw in [("none (idealised)", dict(fractional_bandwidth=0.0,
                                               integration_s=0.0))] + SETUPS:
        for use_pb in (False, True):
            src = (apply_primary_beam(probe, DISH_M, lam, grid)
                   if use_pb else probe)
            out = smear_image(src.image, declination_deg=obs.declination_deg,
                              **kw)
            peaks = [float(out[centre + k, centre]) for k in px]
            img_rows.append({"setup": name, "primary_beam": use_pb, **kw,
                             **{f"response_at_{r:g}as": v
                                for r, v in zip(radii_arcsec, peaks)}})
            print(f"      {name:20}{str(use_pb):>13}" +
                  "".join(f"{v:>9.3f}" for v in peaks))

    print("")
    print("      The phase centre is untouched in every case. At 15 arcsec the")
    print(f"      wide-continuum setup keeps a small fraction of the peak, and the")
    print(f"      primary beam of a {DISH_M:.0f} m dish removes most of what is left.")

    # ---- 4. how long would it take? --------------------------------------
    print(f"\n  [5] observing time, from measured zenith opacity")
    print(f"      {'site and condition':30}{'Tsys':>7}{'SEFD':>9}"
          f"{'1 h rms':>11}{'t for 0.1 mJy':>15}")
    noise_rows = []
    for label, tau in OPACITIES:
        t = system_temperature(tau, elevation_deg=60.0)
        s = sefd(t["t_sys_k"], DISH_M)
        n1 = radiometer_noise(s, N_SELECT, 8e9, 3600.0)
        target = 1e-4                       # 0.1 mJy
        hours = 3600.0 * (n1["sigma_jy"] / target) ** 2 / 3600.0
        noise_rows.append({"condition": label, "tau_zenith": tau,
                           "t_sys_k": t["t_sys_k"], "sefd_jy": s,
                           "sigma_1h_mjy": 1e3 * n1["sigma_jy"],
                           "hours_for_0p1mjy": hours})
        print(f"      {label:30}{t['t_sys_k']:>6.0f}K{s:>8.0f}Jy"
              f"{1e3*n1['sigma_jy']:>10.3f}mJy{hours:>13.2f} h")

    save_csv(rows + img_rows, RESULTS / "exp15_instrumental.csv")
    save_json({"array": {"n_select": N_SELECT, "b_max_m": b_max,
                         "uv_cells": int(state.unique_cells),
                         "dish_m": DISH_M, "frequency_ghz": FREQ_HZ / 1e9},
               "primary_beam_arcsec": pb,
               "field_limits": rows, "imaging": img_rows,
               "sensitivity": noise_rows},
              RESULTS / "exp15_summary.json")
    _figures(rows, img_rows, noise_rows, pb, b_max, lam, ps)
    print(f"\nWrote {RESULTS/'exp15_instrumental.csv'} and exp15_summary.json")


def _figures(rows, img_rows, noise_rows, pb, b_max, lam, ps) -> None:
    fig, ax = new_figure(1, 2, figsize=(11.5, 4.4))

    names = [r["setup"] for r in rows]
    x = np.arange(len(names))
    ax[0].bar(x - 0.25, [min(r["bandwidth_field_arcsec"], 1e4) for r in rows],
              0.25, label="bandwidth smearing", color="#2b6cb0")
    ax[0].bar(x, [min(r["time_field_arcsec"], 1e4) for r in rows], 0.25,
              label="time smearing", color="#c53030")
    ax[0].bar(x + 0.25, [r["primary_half_arcsec"] for r in rows], 0.25,
              label="primary beam (half power)", color="#2f855a")
    ax[0].set_yscale("log")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(names, rotation=15, ha="right", fontsize=8)
    ax[0].set_ylabel("usable field radius (arcsec)")
    ax[0].set_title("Which limit binds the field of view?\n"
                    "(smallest bar wins)", fontsize=9)
    ax[0].legend(fontsize=7)

    # field of view against array size, at fixed setup
    b = np.logspace(np.log10(200.0), np.log10(30000.0), 200)
    for (name, kw), col in zip(SETUPS, ("#2b6cb0", "#dd6b20", "#c53030")):
        f = [smearing_field_arcsec(bb, lam, kw["fractional_bandwidth"],
                                   kw["integration_s"], -30.0) for bb in b]
        usable = [min(v["bandwidth_arcsec"], v["time_arcsec"], pb / 2) for v in f]
        ax[1].loglog(b / 1e3, usable, lw=1.5, color=col, label=name)
    ax[1].axhline(pb / 2, ls="--", lw=1.0, color="grey")
    ax[1].text(0.02, pb / 2 * 1.1, "primary beam floor", fontsize=7, color="grey")
    ax[1].axvline(b_max / 1e3, ls=":", lw=1.0, color="black")
    ax[1].set_xlabel("longest baseline (km)")
    ax[1].set_ylabel("usable field radius (arcsec)")
    ax[1].set_title("Extending the array shrinks the field\n"
                    "(dotted: this array)", fontsize=9)
    ax[1].legend(fontsize=7)
    stamp(fig, ps)
    save_figure(fig, "fig31_field_of_view_limits")

    fig2, ax2 = new_figure(1, 2, figsize=(11.5, 4.3))
    # response against distance from the phase centre, which is where the
    # instrumental effects live: all of them are identically unity at the
    # centre and diverge only off-axis
    radii = sorted(float(k.split("_at_")[1][:-2])
                   for k in img_rows[0] if k.startswith("response_at_"))
    colours = {"none (idealised)": "#718096", "spectral line": "#2b6cb0",
               "narrow continuum": "#dd6b20", "wide continuum": "#c53030"}
    for r in img_rows:
        vals = [r[f"response_at_{x:g}as"] for x in radii]
        ax2[0].plot(radii, vals, "o-" if r["primary_beam"] else "o--",
                    lw=1.4, ms=4, color=colours.get(r["setup"], "black"),
                    alpha=1.0 if r["primary_beam"] else 0.55,
                    label=f"{r['setup']}"
                          f"{' + primary beam' if r['primary_beam'] else ''}")
    ax2[0].set_xlabel("distance from phase centre (arcsec)")
    ax2[0].set_ylabel("peak response")
    ax2[0].set_ylim(0, 1.05)
    ax2[0].set_title("Response across the field\n"
                     "(dashed: smearing only; solid: with primary beam)",
                     fontsize=9)
    ax2[0].legend(fontsize=6, loc="lower left", ncol=2)
    ax2[0].grid(alpha=0.3)

    labels = [r["condition"] for r in noise_rows]
    ax2[1].bar(labels, [r["sigma_1h_mjy"] for r in noise_rows],
               color=["#2f855a", "#dd6b20", "#2b6cb0"])
    ax2[1].set_ylabel("1 hour point-source rms (mJy)")
    ax2[1].set_title("Sensitivity from measured opacity", fontsize=9)
    ax2[1].tick_params(axis="x", labelsize=7, rotation=12)
    stamp(fig2, ps)
    save_figure(fig2, "fig32_smearing_cost")


if __name__ == "__main__":
    main()
