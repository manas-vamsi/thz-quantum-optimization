"""Shared plumbing for the experiment scripts.

Responsibilities kept here so the experiments stay readable:

* load ``configs/default.yaml`` and ``configs/experiments.yaml``
* build any of the layout models from one config block, so that every
  comparison is guaranteed to use identical N / radial extent / seed
* write results in machine-readable form (JSON for metrics, CSV for tables,
  NPZ for arrays) under ``data/``
* stamp every figure with the parameters that produced it
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from thz_opt.arrays.fibonacci import fibonacci_layout  # noqa: E402
from thz_opt.arrays.golden_spiral import golden_spiral_layout  # noqa: E402
from thz_opt.arrays.random_array import random_layout, uniform_ring_layout  # noqa: E402
from thz_opt.interferometry.earth_rotation import ObservationConfig  # noqa: E402
from thz_opt.interferometry.uv import wavelength_from_frequency  # noqa: E402

FIGURES = ROOT / "figures"
GENERATED = ROOT / "data" / "generated"
RESULTS = ROOT / "data" / "results"

LAYOUT_LABELS = {
    "golden": "Golden spiral",
    "fibonacci_radial": "Fibonacci (radial)",
    "fibonacci_golden_angle": "Fibonacci (golden angle)",
    "fibonacci_rational": "Fibonacci (rational angle)",
    "random": "Random (area-uniform)",
    "rings": "Uniform rings",
}


def load_config() -> dict:
    """Merge ``configs/default.yaml`` with ``configs/experiments.yaml``."""
    with open(ROOT / "configs" / "default.yaml", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    with open(ROOT / "configs" / "experiments.yaml", encoding="utf-8") as fh:
        cfg.update(yaml.safe_load(fh))
    return cfg


def merge(base: dict, override: dict) -> dict:
    """Recursive dict merge used to apply an experiment's overrides."""
    out = {k: (v.copy() if isinstance(v, dict) else v) for k, v in base.items()}
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = v
    return out


def build_layout(name: str, cfg: dict) -> np.ndarray:
    """Construct one layout model from the ``array`` block of a config."""
    a = cfg["array"]
    common = dict(n=a["n_antennas"], r_min=a["r_min_m"], r_max=a["r_max_m"])
    if name == "golden":
        return golden_spiral_layout(
            n_turns=a["n_turns"], theta_0=a["theta_0_rad"],
            orientation=a["orientation_rad"], **common
        )
    if name == "fibonacci_radial":
        return fibonacci_layout(variant="radial", skip=a["fibonacci_skip"], **common)
    if name == "fibonacci_golden_angle":
        return fibonacci_layout(variant="golden_angle", **common)
    if name == "fibonacci_rational":
        return fibonacci_layout(
            variant="rational_angle", skip=a["fibonacci_skip"],
            rational_index=a["fibonacci_rational_index"], **common
        )
    if name == "random":
        return random_layout(seed=a["random_seed"], **common)
    if name == "rings":
        return uniform_ring_layout(**common)
    raise ValueError(f"unknown layout model {name!r}")


def observation_config(cfg: dict) -> ObservationConfig:
    o = cfg["observation"]
    return ObservationConfig(
        frequency_hz=float(o["frequency_hz"]),
        declination_deg=o["declination_deg"],
        latitude_deg=o["latitude_deg"],
        hour_angle_start_h=o["hour_angle_start_h"],
        hour_angle_end_h=o["hour_angle_end_h"],
        n_times=o["n_times"],
    )


def wavelength(cfg: dict) -> float:
    return wavelength_from_frequency(float(cfg["observation"]["frequency_hz"]))


# --------------------------------------------------------------------------
# output helpers
# --------------------------------------------------------------------------

def ensure_dirs() -> None:
    for d in (FIGURES, GENERATED, RESULTS):
        d.mkdir(parents=True, exist_ok=True)


def save_json(obj, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=_jsonable)
    return path


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


def save_csv(rows: list, path: Path) -> Path:
    """Write a list of dicts as CSV (pandas is used only for this convenience)."""
    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def save_npz(path: Path, **arrays) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)
    return path


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def new_figure(*args, **kwargs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 160,
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
    })
    return plt.subplots(*args, **kwargs)


def stamp(fig, params: str) -> None:
    """Print the generating parameters on the figure itself.

    Every figure in this repository carries its parameters; a figure whose
    settings are invisible cannot be reproduced or criticised.
    """
    fig.text(0.005, 0.005, params, fontsize=6, va="bottom", ha="left", color="0.35", wrap=True)


def save_figure(fig, name: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)
    print(f"  figure -> {path.relative_to(ROOT)}")
    return path


def param_string(cfg: dict, extra: str = "") -> str:
    a, o, g = cfg["array"], cfg["observation"], cfg["uv_grid"]
    s = (
        f"N={a['n_antennas']}, r=[{a['r_min_m']:.0f}, {a['r_max_m']:.0f}] m, "
        f"nu={float(o['frequency_hz'])/1e9:.0f} GHz "
        f"(lambda={wavelength(cfg)*1e3:.3f} mm), UV grid {g['n_cells']}^2 cells"
    )
    return s + (f", {extra}" if extra else "")


def banner(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)
