"""Tests for real, published array configurations.

These check the loader against values published independently of the files, so a
silently corrupted or misparsed coordinate file cannot pass.
"""

from pathlib import Path

import numpy as np
import pytest

from thz_opt.arrays.real_arrays import (
    REAL_SITES,
    available_arrays,
    itrf_to_local,
    load_cfg,
)

DATA = Path(__file__).resolve().parents[1] / "data" / "external"
ALMA = DATA / "alma.all.cfg"
VLA = DATA / "vla.a.cfg"

pytestmark = pytest.mark.skipif(
    not ALMA.exists(),
    reason="real array files absent; see data/external/SOURCES.md",
)


def test_alma_pad_count_and_dish():
    c = load_cfg(ALMA, dish_diameter_m=12.0)
    assert len(c) == 174
    assert np.allclose(c.diameters, 12.0)
    assert c.coordsys.startswith("LOC")


def test_alma_survey_marker_is_not_a_pad():
    """alma.all.cfg ends with a MASTER0 reference point, which is not a station."""
    kept = load_cfg(ALMA, drop_markers=True)
    dropped = load_cfg(ALMA, drop_markers=False)
    assert len(dropped) == len(kept) + 1
    assert any(n.upper().startswith("MASTER") for n in dropped.names)
    assert not any(n.upper().startswith("MASTER") for n in kept.names)


def test_alma_maximum_baseline_matches_published():
    """ALMA's longest baseline is about 16 km. If the parse were wrong -- wrong
    columns, wrong units -- this would not land near the published value."""
    c = load_cfg(ALMA, dish_diameter_m=12.0)
    d = np.hypot(c.xy[:, None, 0] - c.xy[None, :, 0],
                 c.xy[:, None, 1] - c.xy[None, :, 1])
    assert 15_000 < d.max() < 17_000


def test_alma_seven_metre_pads_are_separable():
    """The file mixes 12 m and 7 m pads; one shadowing limit must not cover both."""
    twelve = load_cfg(ALMA, dish_diameter_m=12.0)
    seven = load_cfg(ALMA, dish_diameter_m=7.0)
    assert len(seven) == 18
    assert set(twelve.names).isdisjoint(seven.names)


@pytest.mark.skipif(not VLA.exists(), reason="vla.a.cfg absent")
def test_vla_itrf_conversion_recovers_published_extent():
    """VLA A-configuration spans about 36 km. These coordinates are geocentric,
    so recovering that number also validates the ITRF to local rotation."""
    c = load_cfg(VLA)
    assert c.coordsys.startswith("XYZ")
    assert len(c) == 27
    d = np.hypot(c.xy[:, None, 0] - c.xy[None, :, 0],
                 c.xy[:, None, 1] - c.xy[None, :, 1])
    assert 34_000 < d.max() < 39_000


def test_itrf_conversion_is_centred_and_isometric():
    """The rotation must preserve distances and put the centroid at the origin."""
    rng = np.random.default_rng(0)
    centre = np.array([1.6e6, -5.0e6, 3.5e6])
    p = centre + rng.normal(0, 3000, size=(30, 3))
    local = itrf_to_local(p)
    assert np.allclose(local.mean(axis=0), 0.0, atol=1e-6)

    def pdist(x):
        return np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)

    assert np.allclose(pdist(p), pdist(local), atol=1e-6)


def test_site_metadata_is_consistent_with_the_files():
    c = load_cfg(ALMA, dish_diameter_m=12.0)
    site = REAL_SITES["alma"]
    assert site["n_pads"] == len(c)
    assert site["dish_diameter_m"] == pytest.approx(float(np.median(c.diameters)))
    assert -24 < site["latitude_deg"] < -22        # Chajnantor
    assert 4000 < site["elevation_m"] < 5500


def test_summary_fields_are_finite():
    s = load_cfg(ALMA, dish_diameter_m=12.0).summary()
    assert all(np.isfinite(v) for v in s.values() if isinstance(v, float))
    assert s["separation_min_m"] > 0
    assert s["baseline_max_m"] > s["radius_max_m"]


def test_available_arrays_discovers_what_is_present():
    found = available_arrays(DATA)
    assert "alma" in found and found["alma"].exists()


def test_filtering_everything_out_is_an_error():
    with pytest.raises(ValueError):
        load_cfg(ALMA, dish_diameter_m=999.0)
