"""Tests for the command-line designer.

The CLI is the surface other people touch, so the tests here are about the
contract it offers rather than the physics underneath: that the outputs are
complete and mutually consistent, that a run records enough to be reproduced,
and that penalty weights cannot be supplied by hand.
"""

import csv
import json
import sys

import numpy as np
import pytest

from thz_opt import cli


# --------------------------------------------------------------------------
# the interface contract
# --------------------------------------------------------------------------

def test_list_sites_names_both_kinds_of_site():
    text = cli.describe_sites()
    assert "alma" in text and "hanle" in text
    assert "select" in text and "place" in text


def test_penalty_weights_are_not_accepted_as_input():
    """They are derived. Offering them would invite a wrong formulation.

    This is asserted rather than documented because a later convenience edit
    adding a --lambda flag would silently reintroduce the failure.
    """
    opts = {a for action in cli.build_parser()._actions
            for a in action.option_strings}
    for forbidden in ("--penalty", "--lambda", "--lambda-select",
                      "--lambda-rosenberg", "--penalty-weight"):
        assert forbidden not in opts, f"{forbidden} must not be a user input"
    assert "--show-penalties" in opts      # readable, not writable


def test_mode_is_inferred_from_whether_a_pad_list_exists():
    p = cli.build_parser()
    assert p.parse_args(["--site", "alma"]).mode == "auto"
    assert "alma" in cli.PAD_FILES
    assert "hanle" not in cli.PAD_FILES


def test_defaults_are_a_runnable_design():
    a = cli.build_parser().parse_args([])
    assert a.n > 1 and a.dish > 0 and a.freq > 0
    assert a.max_baseline > a.dish
    assert 0 < a.separation_factor <= 5


# --------------------------------------------------------------------------
# a real end-to-end run, small enough for the suite
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("cli_run") / "design"
    rc = cli.main(["--site", "alma", "--n", "8", "--m", "24", "--dish", "12",
                   "--quick", "--show-penalties", "-o", str(out)])
    assert rc == 0
    return out


def test_run_writes_every_promised_file(run):
    for name in ("pads.csv", "baselines.csv", "array.cfg", "report.json",
                 "run.yaml"):
        assert (run / name).exists(), f"{name} was not written"


def test_pad_ids_agree_between_the_csv_and_the_cfg(run):
    """The bug that would reach a survey team: same id, different pad."""
    rows = list(csv.DictReader((run / "pads.csv").open(encoding="utf-8")))
    cfg = [ln.split("\t") for ln in
           (run / "array.cfg").read_text(encoding="utf-8").splitlines()
           if ln and not ln.startswith("#")]
    assert len(rows) == len(cfg) == 8
    for r, c in zip(rows, cfg):
        assert r["pad_id"] == c[4].strip()
        assert float(r["east_m"]) == pytest.approx(float(c[0]), abs=0.02)
        assert float(r["north_m"]) == pytest.approx(float(c[1]), abs=0.02)


def test_baseline_table_is_complete_and_matches_the_pads(run):
    rows = list(csv.DictReader((run / "pads.csv").open(encoding="utf-8")))
    pos = {r["pad_id"]: (float(r["east_m"]), float(r["north_m"])) for r in rows}
    pairs = list(csv.DictReader((run / "baselines.csv").open(encoding="utf-8")))
    n = len(rows)
    assert len(pairs) == n * (n - 1) // 2
    for p in pairs:
        a, b = pos[p["pad_i"]], pos[p["pad_j"]]
        assert float(p["distance_m"]) == pytest.approx(
            float(np.hypot(a[0] - b[0], a[1] - b[1])), abs=0.02)


def test_reported_metrics_match_the_pad_table(run):
    """The summary must describe the layout that was actually written."""
    report = json.loads((run / "report.json").read_text(encoding="utf-8"))
    m = report["metrics"]
    rows = list(csv.DictReader((run / "pads.csv").open(encoding="utf-8")))
    xy = np.array([[float(r["east_m"]), float(r["north_m"])] for r in rows])
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    b = d[np.triu_indices(len(xy), 1)]
    # the report uses 3-D separations, so it can only be longer
    assert m["baseline_max_m"] >= b.max() - 1e-6
    assert m["n_antennas"] == len(rows)
    assert m["n_baselines"] == len(rows) * (len(rows) - 1) // 2


def test_elevations_are_real_not_placeholder(run):
    """The .cfg z column plus the site elevation, not NaN."""
    rows = list(csv.DictReader((run / "pads.csv").open(encoding="utf-8")))
    el = np.array([float(r["elevation_m"]) for r in rows])
    assert np.all(np.isfinite(el))
    assert 4000 < el.mean() < 6000          # Chajnantor, roughly 5000 m


def test_run_yaml_records_enough_to_reproduce(run):
    import yaml

    rec = yaml.safe_load((run / "run.yaml").read_text(encoding="utf-8"))
    for key in ("site", "n", "dish", "freq", "dec", "hours", "seed",
                "method", "uv_cells"):
        assert key in rec, f"run.yaml does not record {key!r}"
    assert rec["n"] == 8 and rec["site"] == "alma"


def test_derived_penalties_are_reported_and_positive(run):
    report = json.loads((run / "report.json").read_text(encoding="utf-8"))
    pen = report["qubo_penalties"]
    assert pen["lambda_rosenberg"] > 0
    assert pen["lambda_select"] > pen["lambda_rosenberg"]
    assert pen["n_variables"] == 8 + 8 * 7 // 2


def test_a_config_file_reproduces_a_command_line(tmp_path):
    import yaml

    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    assert cli.main(["--site", "alma", "--n", "8", "--m", "24", "--quick",
                     "-o", str(out_a)]) == 0

    rec = yaml.safe_load((out_a / "run.yaml").read_text(encoding="utf-8"))
    rec["out"] = str(out_b)
    cfg = tmp_path / "rerun.yaml"
    cfg.write_text(yaml.safe_dump(rec), encoding="utf-8")
    assert cli.main(["--config", str(cfg)]) == 0

    assert ((out_a / "pads.csv").read_text(encoding="utf-8")
            == (out_b / "pads.csv").read_text(encoding="utf-8"))


def test_an_unknown_site_fails_loudly(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["--site", "mount-doom", "-o", str(tmp_path / "x")])
