"""Measured site parameters for the Ladakh submillimetre candidate sites.

Every number here is published and carries its source. Nothing is estimated,
and where a quantity has not been measured the entry says so rather than
carrying a plausible-looking placeholder.

Why these four sites
--------------------
The project note names Hanle, describing it as "a comparable high-altitude,
low-PWV, radio-quiet site" to ALMA's Chajnantor plateau. The measured record
does not support the comparison. Raghunath et al. (2026) analyse 184 months of
ERA5 reanalysis (January 2010 to April 2025) over the Ladakh plateau and find
the fraction of time with precipitable water vapour at or below 1 mm -- the
threshold the project note itself sets -- to be about 5 % at Hanle and about
8 % at Merak, against about 23 % and 19 % at two unbuilt sites further north.
For reference ALMA's twenty-year year-round median PWV is near 1 mm, so
Chajnantor sits at roughly half the time below that threshold.

Hanle is therefore an order of magnitude worse than ALMA on the quantity that
matters most at these frequencies, while remaining the only Ladakh site with
road access, power and an existing observatory. That trade -- infrastructure
against dryness -- is a real design decision, so both are carried here rather
than one being quietly preferred.

What is measured and what is not
--------------------------------
Measured: geodetic position, elevation, precipitable water vapour, and at Hanle
the 220 GHz zenith opacity from an on-site tipping radiometer running since
1999.

Not measured anywhere in Ladakh: the atmospheric **phase** structure function.
That quantity requires an operating interferometer recording phase over years,
as ALMA Memo 624 did over 17 000 observations. No interferometer exists at any
of these sites, so no equivalent number exists, and none is invented here. Any
coherence figure computed for a Ladakh array is therefore a transfer from
Chajnantor and is labelled as a model, not a measurement.
"""

from __future__ import annotations

__all__ = ["LADAKH_SITES", "REFERENCE_SITES", "site", "pwv_comparison"]

#: Candidate sites on the Ladakh plateau.
LADAKH_SITES = {
    "hanle": {
        "name": "Indian Astronomical Observatory, Hanle",
        "description": "Mt Saraswati, Digpa-ratsa Ri. Existing 2 m Himalayan "
                       "Chandra Telescope and MACE gamma-ray telescope; road, "
                       "power and staff accommodation already present.",
        "latitude_deg": 32.7794,
        "longitude_deg": 78.9642,
        "elevation_m": 4500.0,
        "pwv_below_1mm_fraction": 0.05,
        "pwv_source": "Raghunath et al. 2026, arXiv:2604.13487 (ERA5, 184 months)",
        "tau_220ghz_note": "zenith opacity below 0.06 for about 40 % of winter",
        "tau_220ghz_source": "Ananthasayanam et al. 2001, arXiv:astro-ph/0111492",
        "on_site_radiometer": "220 GHz tipping radiometer, operating since 1999",
        "infrastructure": "existing observatory",
        "phase_structure_function": None,     # never measured, not invented
    },
    "site_a": {
        "name": "Ladakh candidate site A (unbuilt)",
        "description": "Highest-ranked site in the 2026 reanalysis study. No "
                       "road, power or buildings; coordinates are approximate, "
                       "read from that study's site table.",
        "latitude_deg": 34.25,
        "longitude_deg": 78.75,
        "elevation_m": None,                  # taken from the DEM at runtime
        "pwv_below_1mm_fraction": 0.23,
        "pwv_source": "Raghunath et al. 2026, arXiv:2604.13487 (ERA5, 184 months)",
        "infrastructure": "none",
        "phase_structure_function": None,
    },
    "site_b": {
        "name": "Ladakh candidate site B (unbuilt)",
        "description": "Second-ranked site in the same study.",
        "latitude_deg": 32.50,
        "longitude_deg": 79.00,
        "elevation_m": None,
        "pwv_below_1mm_fraction": 0.19,
        "pwv_source": "Raghunath et al. 2026, arXiv:2604.13487 (ERA5, 184 months)",
        "infrastructure": "none",
        "phase_structure_function": None,
    },
    "merak": {
        "name": "Merak, Pangong Tso",
        "description": "Proposed site near Pangong Tso with partial access. "
                       "Coordinates approximate.",
        "latitude_deg": 33.79,
        "longitude_deg": 78.62,
        "elevation_m": None,
        "pwv_below_1mm_fraction": 0.08,
        "pwv_source": "Raghunath et al. 2026, arXiv:2604.13487 (ERA5, 184 months)",
        "infrastructure": "partial road access",
        "phase_structure_function": None,
    },
}

#: Operating arrays, for scale. Used only for comparison, never as design input.
REFERENCE_SITES = {
    "alma": {
        "name": "ALMA, Llano de Chajnantor",
        "latitude_deg": -23.0229,
        "longitude_deg": -67.7551,
        "elevation_m": 5058.7,
        "pwv_median_mm": 1.0,
        "pwv_below_1mm_fraction": 0.50,    # implied by the 20-year median
        "pwv_source": "Cortes et al. 2020, doi:10.1051/0004-6361/202037784",
        "phase_structure_function": "ALMA Memo 624 (Maud et al. 2023)",
    },
}


def site(key: str) -> dict:
    """One site's published parameters."""
    table = {**LADAKH_SITES, **REFERENCE_SITES}
    if key not in table:
        raise ValueError(f"unknown site {key!r}; choose from {tuple(table)}")
    return dict(table[key])


def pwv_comparison() -> list:
    """Fraction of time below the 1 mm PWV threshold, best site first."""
    rows = [
        {"key": k, "name": v["name"],
         "fraction_below_1mm": v["pwv_below_1mm_fraction"],
         "infrastructure": v.get("infrastructure", "operating array")}
        for k, v in {**LADAKH_SITES, **REFERENCE_SITES}.items()
    ]
    return sorted(rows, key=lambda r: -r["fraction_below_1mm"])
