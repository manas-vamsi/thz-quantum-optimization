"""Tests for the Ladakh array design.

The coordinate conversions matter more here than anywhere else in the project.
Everywhere else a sign error in a transform shows up as a bad score; here it
shows up as a pad dug in the wrong place, so the round trips are checked
explicitly rather than assumed from the fact that the figures look plausible.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "ladakh_array_design"))

import design  # noqa: E402
from thz_opt.sites import site_data, terrain  # noqa: E402
from thz_opt.arrays.real_arrays import load_cfg  # noqa: E402

HANLE = site_data.site("hanle")


# --------------------------------------------------------------------------
# site data
# --------------------------------------------------------------------------

def test_no_site_claims_a_measured_phase_structure_function():
    """None exists for Ladakh, and inventing one is the failure mode here."""
    for key, entry in site_data.LADAKH_SITES.items():
        assert entry["phase_structure_function"] is None, key


def test_hanle_is_recorded_as_drier_than_it_is_not():
    """The project note calls Hanle ALMA-comparable; the data says otherwise.

    This is a guard against the claim drifting back in through a later edit.
    """
    hanle = site_data.site("hanle")["pwv_below_1mm_fraction"]
    alma = site_data.site("alma")["pwv_below_1mm_fraction"]
    assert hanle < alma / 5.0
    assert site_data.site("site_a")["pwv_below_1mm_fraction"] > 4 * hanle


def test_pwv_comparison_is_sorted_best_first():
    rows = site_data.pwv_comparison()
    fracs = [r["fraction_below_1mm"] for r in rows]
    assert fracs == sorted(fracs, reverse=True)


def test_unknown_site_is_rejected():
    with pytest.raises(ValueError):
        site_data.site("everest")


# --------------------------------------------------------------------------
# coordinate transforms -- the ones that put concrete in the wrong valley
# --------------------------------------------------------------------------

def test_enu_round_trip_returns_the_original_point():
    lat0, lon0 = HANLE["latitude_deg"], HANLE["longitude_deg"]
    east = np.array([0.0, 1500.0, -1500.0, 800.0])
    north = np.array([0.0, -1200.0, 900.0, 2000.0])
    lat, lon = terrain.enu_to_latlon(east, north, lat0, lon0)
    e2, n2 = terrain.local_enu(lat, lon, lat0, lon0)
    assert e2 == pytest.approx(east, abs=1e-6)
    assert n2 == pytest.approx(north, abs=1e-6)


def test_origin_maps_to_the_site_itself():
    lat0, lon0 = HANLE["latitude_deg"], HANLE["longitude_deg"]
    lat, lon = terrain.enu_to_latlon(0.0, 0.0, lat0, lon0)
    assert float(lat) == pytest.approx(lat0)
    assert float(lon) == pytest.approx(lon0)


def test_a_degree_of_latitude_is_about_111_km():
    """Catches a radians/degrees slip, which round-tripping alone would not."""
    _, north = terrain.local_enu(HANLE["latitude_deg"] + 1.0,
                                 HANLE["longitude_deg"],
                                 HANLE["latitude_deg"], HANLE["longitude_deg"])
    assert 110_000 < float(north) < 112_000


def test_longitude_scale_shrinks_with_latitude():
    """One degree of longitude is cos(lat) shorter than one of latitude."""
    east, _ = terrain.local_enu(HANLE["latitude_deg"], HANLE["longitude_deg"] + 1.0,
                                HANLE["latitude_deg"], HANLE["longitude_deg"])
    expected = 111_320.0 * math.cos(math.radians(HANLE["latitude_deg"]))
    assert float(east) == pytest.approx(expected, rel=0.01)


# --------------------------------------------------------------------------
# terrain
# --------------------------------------------------------------------------

def _fake_dem(slope_ew_deg=0.0, n=41, pitch=16.0):
    """A synthetic DEM with a known constant slope, for checking the gradient."""
    east = (np.arange(n) - n // 2) * pitch
    north = (np.arange(n) - n // 2) * pitch
    ee, _ = np.meshgrid(east, north)
    z = ee * math.tan(math.radians(slope_ew_deg))
    return terrain.DEM(32.0, 79.0, east, north, z, 13, 0)


def test_slope_recovers_a_known_plane():
    for want in (0.0, 5.0, 12.5, 30.0):
        s = terrain.slope_deg(_fake_dem(want))
        assert float(np.median(s)) == pytest.approx(want, abs=1e-6)


def test_buildable_mask_applies_both_slope_and_radius():
    dem = _fake_dem(3.0, n=81)
    flat_all = terrain.buildable_mask(dem, max_slope_deg=10.0)
    assert flat_all.all()                      # 3 degrees passes a 10 degree limit
    assert not terrain.buildable_mask(dem, max_slope_deg=1.0).any()

    limited = terrain.buildable_mask(dem, 10.0, max_radius_m=200.0)
    ee, nn = np.meshgrid(dem.east, dem.north)
    assert np.all(np.hypot(ee, nn)[limited] <= 200.0 + 1e-9)
    assert limited.sum() < flat_all.sum()


# --------------------------------------------------------------------------
# design specification
# --------------------------------------------------------------------------

def test_resolution_follows_lambda_over_baseline():
    spec = design.DesignSpec(max_baseline_m=3000.0, frequency_hz=230e9)
    lam = 299_792_458.0 / 230e9
    expected = lam / 3000.0 / (np.pi / (180 * 3600))
    assert spec.resolution_arcsec == pytest.approx(expected)
    # doubling the array halves the beam
    assert design.DesignSpec(max_baseline_m=6000.0).resolution_arcsec == \
        pytest.approx(spec.resolution_arcsec / 2)


def test_minimum_separation_is_the_shadowing_rule():
    spec = design.DesignSpec(dish_m=8.0, separation_factor=1.5)
    assert spec.min_separation_m == pytest.approx(12.0)


def test_array_radius_is_half_the_longest_baseline():
    assert design.DesignSpec(max_baseline_m=3000.0).array_radius_m == 1500.0


# --------------------------------------------------------------------------
# the deliverable
# --------------------------------------------------------------------------

def test_generated_cfg_reads_back_through_the_alma_loader(tmp_path):
    """The design must be readable by the same parser that reads real arrays.

    If `load_cfg` can open our output, so can CASA, and the pad list is a real
    interchange artefact rather than a private format.
    """
    spec = design.DesignSpec(n_antennas=6, dish_m=8.0)
    rng = np.random.default_rng(0)
    xy = rng.uniform(-1200, 1200, size=(6, 2))
    elev = np.full(6, 4500.0)
    path = design.write_cfg(tmp_path / "test_array.cfg", xy, elev, spec, HANLE)

    cfg = load_cfg(path, dish_diameter_m=8.0)
    assert len(cfg) == 6
    # write_cfg emits pads in the canonical order, not the order handed in, so
    # the comparison is against that ordering rather than the raw input
    expected = np.array([[r["east_m"], r["north_m"]]
                         for r in design.pad_table(xy, elev, HANLE)])
    assert cfg.xy == pytest.approx(expected, abs=1e-3)


def test_cfg_carries_the_not_a_construction_plan_warning(tmp_path):
    """An engineer opening the file must see the status before the numbers."""
    spec = design.DesignSpec(n_antennas=4)
    xy = np.zeros((4, 2))
    xy[:, 0] = [0.0, 100.0, 200.0, 300.0]
    path = design.write_cfg(tmp_path / "a.cfg", xy, np.full(4, 4500.0), spec, HANLE)
    header = path.read_text(encoding="utf-8")
    assert "PRELIMINARY" in header
    assert "not a construction plan" in header.lower()


def test_candidate_lattice_respects_slope_and_spacing():
    dem = _fake_dem(2.0, n=201, pitch=16.0)
    spec = design.DesignSpec(candidate_spacing_m=120.0, max_baseline_m=2000.0,
                             max_slope_deg=10.0)
    xy, elev = design.candidate_pads(dem, spec)
    assert len(xy) == len(elev) > 0
    assert np.all(np.hypot(xy[:, 0], xy[:, 1]) <= spec.array_radius_m + 1e-6)
    # every candidate is far enough apart to be a legal pad pair
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    off = d[~np.eye(len(xy), dtype=bool)]
    assert off.min() >= spec.min_separation_m


def test_candidates_vanish_on_impossible_ground():
    """A cliff yields no pads rather than a layout drawn on a vertical face."""
    dem = _fake_dem(45.0, n=101)
    spec = design.DesignSpec(max_slope_deg=10.0, max_baseline_m=1000.0)
    xy, _ = design.candidate_pads(dem, spec)
    assert len(xy) == 0


# --------------------------------------------------------------------------
# cross-file consistency -- the bug that would reach a survey team
# --------------------------------------------------------------------------

def test_csv_and_cfg_name_the_same_physical_pads(tmp_path):
    """`P001` must be one pad, not two.

    An earlier version sorted inside the CSV writer but not the CFG writer, so
    the two deliverables gave the same identifier to different coordinates.
    Everything downstream -- surveying, cable schedules, antenna assignment --
    reads one or the other.
    """
    import run_design

    rng = np.random.default_rng(7)
    xy = rng.uniform(-1400, 1400, size=(12, 2))
    elev = 4500.0 + rng.uniform(-40, 40, size=12)
    spec = design.DesignSpec(n_antennas=12)
    result = {"xy": xy, "elev": elev, "site": HANLE}

    csv_path = run_design.write_pads_csv(tmp_path / "p.csv", result)
    cfg_path = design.write_cfg(tmp_path / "p.cfg", xy, elev, spec, HANLE)

    import csv as _csv
    rows = list(_csv.DictReader(csv_path.open(encoding="utf-8")))
    cfg = [ln.split("\t") for ln in cfg_path.read_text(encoding="utf-8").splitlines()
           if ln and not ln.startswith("#")]

    assert len(rows) == len(cfg) == 12
    for r, c in zip(rows, cfg):
        assert r["pad_id"] == c[4].strip()
        assert float(r["east_m"]) == pytest.approx(float(c[0]), abs=0.02)
        assert float(r["north_m"]) == pytest.approx(float(c[1]), abs=0.02)


def test_pad_ordering_is_deterministic():
    """Two calls on the same layout must not permute the identifiers."""
    rng = np.random.default_rng(11)
    xy = rng.uniform(-1000, 1000, size=(20, 2))
    elev = np.full(20, 4500.0)
    a = design.pad_table(xy, elev, HANLE)
    b = design.pad_table(xy.copy(), elev.copy(), HANLE)
    assert [r["pad_id"] for r in a] == [r["pad_id"] for r in b]
    assert [r["east_m"] for r in a] == [r["east_m"] for r in b]


def test_cfg_height_datum_is_the_array_not_the_published_site_elevation():
    """Hanle's published 4500 m is the ridge-top observatory, not the pad field.

    Referencing to it put every pad at a spurious -225 m in an earlier run.
    """
    import re
    import tempfile

    elev = np.array([4260.0, 4270.0, 4280.0, 4290.0])
    xy = np.column_stack([np.arange(4) * 500.0, np.zeros(4)])
    spec = design.DesignSpec(n_antennas=4)
    with tempfile.TemporaryDirectory() as d:
        p = design.write_cfg(Path(d) / "a.cfg", xy, elev, spec, HANLE)
        text = p.read_text(encoding="utf-8")
        datum = float(re.search(r"height datum=([\d.]+)", text).group(1))
        assert datum == pytest.approx(elev.mean())
        z = [float(ln.split("\t")[2]) for ln in text.splitlines()
             if ln and not ln.startswith("#")]
        assert max(abs(v) for v in z) < 20.0      # not ~225
