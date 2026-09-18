"""Measured atmospheric parameters for the ALMA site.

These are published measurements, not assumed constants. Every value below is
quoted from a paper or memo and carries its source in the table that holds it,
so a reader can check the number rather than trust it.

Phase stability
---------------
From ALMA Memo 624 (Maud et al. 2023, arXiv:2304.08318), an analysis of over
17 000 ALMA observations from 2015-2019, reporting phase RMS as a
frequency-independent path-length variation on summary baselines of 500, 1000,
5000 and 10 000 m:

* median path RMS **200 um** on a 1000 m baseline over a 120 s timescale,
  without water-vapour-radiometer correction;
* **115 um** for the subset of observations below the median PWV of 1.24 mm;
* with WVR correction applied, **70 um** for wind speeds below 10 m/s and
  **140 um** above it.

The spatial structure function exponents are those of Matsushita et al. (2017),
as adopted in that memo: path RMS scales as ``b**0.65`` below a 1 km baseline
and ``b**0.22`` above it without WVR correction, and as ``b**0.60`` and
``b**0.29`` with it. Note these are *measured* exponents and differ from the
textbook Kolmogorov values of 5/6 and 1/3 -- the real atmosphere at this site is
shallower than the idealised model in both regimes.

The memo's own consistency check: scaling the WVR-corrected 500 m value by 1.51
gives 1000 m, then 1.59 gives 5000 m, then 1.22 gives 10 000 m. Those factors
imply exponents of 0.59, 0.29 and 0.29, matching the quoted slopes.

Precipitable water vapour
-------------------------
Two independent sources agree. Cortés et al. (2020, A&A 640, A126) analyse
twenty years of measurements and report a year-round median near 1 mm on the
Chajnantor plateau, with the extreme months January and August at 2.56 mm and
0.72 mm. ALMA Memo 624 Table 3 gives monthly percentiles from the observation
metadata, reproduced below.

What is still not measured here
-------------------------------
The phase values are medians over all ALMA observations, so they represent the
site rather than any particular night, and the memo notes that its sample misses
the very worst conditions in which no observation was attempted. Elevation is
not normalised: the quoted values are as measured across a 21-88 degree
elevation range with a median of 58 degrees.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "ALMA_PHASE",
    "ALMA_PWV_MONTHLY",
    "ALMA_PWV_SUMMARY",
    "phase_conditions",
    "pwv_percentile",
    "monthly_pwv",
]

#: Measured phase-stability conditions on the ALMA site.
#: ``sigma_1km_um`` is the median path-length RMS on a 1000 m baseline over a
#: 120 s timescale; ``alpha_short``/``alpha_long`` are the structure-function
#: exponents below and above the 1 km break.
ALMA_PHASE = {
    "uncorrected": {
        "sigma_1km_um": 200.0,
        "alpha_short": 0.65,
        "alpha_long": 0.22,
        "description": "median of all observations, no WVR correction",
        "source": "ALMA Memo 624 (Maud et al. 2023), arXiv:2304.08318",
    },
    "uncorrected_low_pwv": {
        "sigma_1km_um": 115.0,
        "alpha_short": 0.65,
        "alpha_long": 0.22,
        "description": "observations below the median PWV of 1.24 mm, no WVR",
        "source": "ALMA Memo 624 (Maud et al. 2023), arXiv:2304.08318",
    },
    "wvr_corrected": {
        "sigma_1km_um": 70.0,
        "alpha_short": 0.60,
        "alpha_long": 0.29,
        "description": "WVR corrected, wind speed below 10 m/s",
        "source": "ALMA Memo 624 (Maud et al. 2023), arXiv:2304.08318",
    },
    "wvr_corrected_windy": {
        "sigma_1km_um": 140.0,
        "alpha_short": 0.60,
        "alpha_long": 0.29,
        "description": "WVR corrected, wind speed above 10 m/s",
        "source": "ALMA Memo 624 (Maud et al. 2023), arXiv:2304.08318",
    },
}

BREAK_SCALE_M = 1000.0   #: structure-function break, measured at 1 km

#: PWV in mm by month and percentile, from ALMA Memo 624 Table 3.
#: February is absent from the source data and is omitted here rather than
#: interpolated.
ALMA_PWV_MONTHLY = {
    "percentiles": [10, 25, 50, 75],
    "months": ["Jan", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"],
    "pwv_mm": {
        75: [4.65, 3.27, 2.52, 1.45, 1.30, 1.53, 1.06, 1.56, 1.54, 2.14, 2.97],
        50: [3.05, 2.03, 1.51, 0.98, 0.88, 0.82, 0.70, 0.92, 0.94, 1.29, 1.64],
        25: [1.86, 1.23, 0.90, 0.64, 0.64, 0.58, 0.48, 0.58, 0.63, 0.69, 0.94],
        10: [1.12, 0.75, 0.60, 0.43, 0.51, 0.43, 0.37, 0.45, 0.47, 0.46, 0.62],
    },
    "source": "ALMA Memo 624 (Maud et al. 2023) Table 3",
}

#: Long-baseline summary from twenty years of measurements.
ALMA_PWV_SUMMARY = {
    "median_year_round_mm": 1.0,
    "median_january_mm": 2.56,
    "median_august_mm": 0.72,
    "site": "Llano de Chajnantor plateau, 5059 m",
    "span_years": 20,
    "source": "Cortes et al. (2020), A&A 640, A126, doi:10.1051/0004-6361/201937784",
}


def phase_conditions(condition: str = "wvr_corrected") -> dict:
    """Measured phase parameters for one observing condition.

    ``condition`` is a key of :data:`ALMA_PHASE`. The returned dictionary adds
    ``sigma_1km_m`` in metres and the 1 km break scale, so it can be handed
    straight to :class:`thz_opt.constraints.coherence.CoherenceConfig`.
    """
    if condition not in ALMA_PHASE:
        raise ValueError(f"unknown condition {condition!r}; "
                         f"choose from {tuple(ALMA_PHASE)}")
    entry = dict(ALMA_PHASE[condition])
    entry["sigma_1km_m"] = entry["sigma_1km_um"] * 1e-6
    entry["break_scale_m"] = BREAK_SCALE_M
    entry["condition"] = condition
    entry["measured"] = True
    return entry


def monthly_pwv(month: str, percentile: int = 50) -> float:
    """Measured PWV in mm for a month and percentile.

    ``month`` is a three-letter name as listed in :data:`ALMA_PWV_MONTHLY`;
    February is not in the source data and raises.
    """
    months = ALMA_PWV_MONTHLY["months"]
    if month not in months:
        raise ValueError(f"no measured data for {month!r}; "
                         f"available: {months} (February absent from the source)")
    if percentile not in ALMA_PWV_MONTHLY["pwv_mm"]:
        raise ValueError(f"percentile must be one of "
                         f"{sorted(ALMA_PWV_MONTHLY['pwv_mm'])}")
    return float(ALMA_PWV_MONTHLY["pwv_mm"][percentile][months.index(month)])


def pwv_percentile(percentile: int = 50) -> np.ndarray:
    """Measured PWV in mm across all months at one percentile."""
    if percentile not in ALMA_PWV_MONTHLY["pwv_mm"]:
        raise ValueError(f"percentile must be one of "
                         f"{sorted(ALMA_PWV_MONTHLY['pwv_mm'])}")
    return np.asarray(ALMA_PWV_MONTHLY["pwv_mm"][percentile], dtype=float)
