# Ladakh submillimetre array design

A preliminary array design for a submillimetre interferometer in Ladakh: where
each antenna pad would go, in coordinates a survey team could take to the site.

This folder is self-contained. It uses the optimiser from `src/thz_opt/` but
answers a different question from the rest of the repository. There, ALMA has
174 concrete pads and the problem is *which twenty to occupy*. Here nothing is
built, so the problem is *where to put them* — placement on real ground rather
than selection from a list.

> **Status: preliminary geometric design. Not a construction plan.**
> No geotechnical survey, no land-rights check, no access roads, no power
> routing, no environmental assessment. A pad position here means "the ground
> at this point is flat enough according to a 16 m elevation model", nothing
> more. Do not present this as site-ready.

---

## Run it

```bash
cd ladakh_array_design
python run_design.py          # ~8 minutes, downloads terrain on first run
```

Outputs land in `outputs/`. Elevation tiles are cached in `data/dem_cache/`
(git-ignored, re-downloaded on demand).

---

## 1. Which site, and the correction that came out of asking

The project onboarding note names Hanle and describes it as *"a comparable
high-altitude, low-PWV, radio-quiet site"* to ALMA's Chajnantor plateau.

**The measured record does not support that comparison.** Raghunath et al.
(2026) analysed 184 months of ERA5 reanalysis over the Ladakh plateau and
report the fraction of time with precipitable water vapour at or below 1 mm —
the threshold the project note itself sets:

| Site | Time with PWV ≤ 1 mm | Infrastructure |
|---|---:|---|
| ALMA, Chajnantor | ~50 % | operating array |
| **Ladakh site A** (34.25 °N, 78.75 °E) | **23 %** | none |
| Ladakh site B (32.50 °N, 79.00 °E) | 19 % | none |
| Merak, Pangong Tso | 8 % | partial road access |
| **Hanle (IAO)** | **5 %** | existing observatory |

Hanle is roughly **ten times worse than ALMA** on the quantity that matters
most at these frequencies, and four times worse than the best unbuilt Ladakh
site. It remains the only candidate with road, power and staff already there.

That trade — infrastructure against dryness — is a genuine decision, so the
design is produced at **both** Hanle and site A rather than one being quietly
preferred.

**A second trade the water-vapour numbers hide.** Once the array centre is
free to move onto the best available ground (§3), Hanle holds all three
configurations within 15 m of relief; site A needs 57 m and its extended
configuration finds only 179 buildable positions against Hanle's 544. Hanle's
plain is large and genuinely flat; site A's is broken ground that happens to
contain flat patches. The terrain penalty at site A is real and partly offsets
its being four times drier.

---

## 2. Design parameters, and why each one

| Parameter | Value | Reason |
|---|---|---|
| Antennas | 16 | 120 baselines; enough for real imaging, credible as a first national array |
| Dish diameter | 8 m | domestically buildable; sets minimum pad separation at 1.5 × D = **12 m**, the shadowing rule from the project note |
| Frequency | **230 GHz** (1.303 mm) | the only band with *measured on-site transparency* at Hanle — a 220 GHz tipping radiometer has run there since 1999. Designing for a band with no local measurement would be guesswork |
| Declination | +30° | transits within 3° of zenith at these latitudes: the best case the site offers |
| Track | 4 h | ±2 h of hour angle about transit |
| Slope limit | 10° | admits ground needing cut-and-fill, excludes ground that is impossible. An engineering assumption, not a measurement |

---

## 3. Why there are three configurations

The first design run maximised occupied UV cells and produced a **ring**: all
sixteen pads on the outer rim of the allowed field, none in the inner half,
exactly **one baseline under 300 m out of 120**.

That is not an optimiser bug. The outer annuli of a UV grid contain far more
cells than the inner ones, because area grows with radius, so "maximise
distinct cells" is quietly an instruction to flee to the boundary. The array
that results resolves beautifully and is blind to anything wider than about
half an arcsecond.

Four changes follow, each forced by the previous one.

**A Gaussian target UV density** (Boone 2001) replaces the bare cell count, so
inner cells are worth more and the optimiser buys short baselines. Necessary,
not sufficient.

**A landform filter.** Flat ground is not one surface. Hanle's observatory
stands on an isolated summit about 220 m above a wide plain, and both are
locally flat, so a slope test accepts them equally. An early run placed
fourteen pads on the plain and one on the summit, and that single outlier
accounted for the entire 235 m of reported relief — an array with a 240 m climb
between two of its antennas. Candidates are now confined to the 60 m elevation
band holding the most buildable ground.

**An array centre that is not the site marker.** A published coordinate names a
building. Centring the design on Hanle's observatory put the compact array on
the summit and the extended array on the plain — two installations 220 m apart
vertically, not a shared pad field. The centre is now chosen by convolving the
landform mask with a disc of one array radius and taking the maximum: *where
does a full disc of usable ground actually fit?* It moves 2049 m at Hanle and
1954 m at site A. The obvious alternative, the centroid of the landform, fails
precisely here — Hanle's plain is an annulus and the centroid of an annulus
sits in its hole, back on the summit.

**Three configurations, not one, and not two.** Sixteen antennas over 3 km
cannot also provide 32 m baselines. Two configurations only give continuous
coverage if they overlap — the smaller array must *resolve* finer than the
larger one can *see*:

    lambda / b_max(small)  <  0.6 lambda / b_min(large)

A 400 m compact array failed that test at Hanle, leaving structures between
0.48″ and 0.69″ sampled well by neither. Widening it to 800 m closed the gap
but spread its own pads, pushing its shortest baseline from 32 m to 90 m and
collapsing its largest recoverable scale from 5.0″ to 2.5″. One array cannot be
both the short-spacing array and the bridge.

Three rungs resolve it, which is why real arrays have several — ALMA has ten.
Compact keeps the 32 m baselines and the 5″ sensitivity, extended keeps the
resolution, and intermediate exists only to join them. **The overlap at each
joint is computed and printed on every run**, so a terrain change that reopens
a gap is reported rather than hidden.

## 4. The design

All quantities below are computed in three dimensions from the elevation model,
not under the coplanar approximation.

### Hanle (32.7794 °N, 78.9642 °E) — array centre 2049 m from the site marker

| | Compact | Intermediate | Extended |
|---|---:|---:|---:|
| Candidate positions | 43 | 585 | 544 |
| Baselines | 32 – 232 m | 64 – 878 m | 338 – 2788 m |
| **Resolves to** | 1.158″ | 0.306″ | **0.096″** |
| **Sees up to** | **5.01″** | 2.50″ | 0.48″ |
| Occupied UV cells | 1772 | 2038 | 2108 |
| Peak sidelobe | 0.102 | 0.105 | 0.100 |
| Ground relief | 11 m | 9 m | 15 m |

### Ladakh site A (34.25 °N, 78.75 °E) — array centre 1954 m from the site marker

| | Compact | Intermediate | Extended |
|---|---:|---:|---:|
| Candidate positions | 45 | 244 | 179 |
| Baselines | 32 – 228 m | 47 – 857 m | 179 – 2686 m |
| **Resolves to** | 1.178″ | 0.314″ | **0.100″** |
| **Sees up to** | **5.10″** | 3.40″ | 0.90″ |
| Occupied UV cells | 1920 | 2088 | 1600 |
| Peak sidelobe | 0.088 | 0.098 | 0.142 |
| Ground relief | 12 m | 57 m | 53 m |

### Combined

| | Hanle | Site A |
|---|---:|---:|
| Finest detail | 0.096″ | 0.100″ |
| Widest structure | 5.01″ | 5.10″ |
| **Range of angular scales** | **52×** | **51×** |
| compact / intermediate joint | overlap 1.35″ | overlap 2.22″ |
| intermediate / extended joint | overlap 0.17″ | overlap 0.59″ |
| Coverage | **continuous** | **continuous** |
| Worst ground relief | **15 m** | 57 m |

Hanle holds every configuration within 15 m of relief; site A needs 57 m. With
the array centre free to move, Hanle's plain is the better ground by a clear
margin, which partly offsets site A being four times drier.

## 5. Outputs

For each of `hanle_compact`, `hanle_intermediate`, `hanle_extended`,
`site_a_compact`, `site_a_intermediate`, `site_a_extended`:

| File | Contents |
|---|---|
| `<tag>_pads.csv` | **the deliverable** — pad id, latitude, longitude, elevation, local east/north, radius from centre |
| `<tag>_array.cfg` | CASA observatory configuration, same format as `alma.all.cfg`; drop into `simobserve` |
| `<tag>_report.json` | baselines, resolution, UV coverage, PSF, and the full parameter set |
| `fig_<tag>_terrain.png` | terrain, buildable ground, candidates, chosen pads |
| `fig_<tag>_uv.png` | UV coverage, gridded sampling, dirty beam |

Plus `site_comparison.csv` and `fig_site_comparison.png`.

Pad identifiers are consistent between the CSV and the `.cfg` — both are
generated from one ordering function, and a test enforces it. An earlier
version sorted in one writer but not the other, so `P001` named different
physical pads in the two files.

---

## 6. Data sources

| Quantity | Source | Status |
|---|---|---|
| Terrain elevation | [AWS Terrain Tiles](https://registry.opendata.aws/terrain-tiles/) (Tilezen terrarium, SRTM/ASTER derived), 16 m/pixel | measured |
| PWV, all Ladakh sites | Raghunath et al. 2026, [arXiv:2604.13487](https://arxiv.org/abs/2604.13487) — ERA5, 184 months, Jan 2010 – Apr 2025 | measured |
| 220 GHz opacity at Hanle | Ananthasayanam et al. 2001, [arXiv:astro-ph/0111492](https://arxiv.org/abs/astro-ph/0111492) — on-site tipping radiometer | measured |
| ALMA comparison PWV | Cortés et al. 2020, [doi:10.1051/0004-6361/202037784](https://doi.org/10.1051/0004-6361/202037784) | measured |
| Site coordinates | IAO published position | measured |
| **Phase structure function** | — | **does not exist for any Ladakh site** |

### Terrain pipeline validation

Checked against two independently published elevations before being trusted:

| Site | This pipeline | Published |
|---|---:|---:|
| IAO Hanle | 4491 m | 4500 m |
| ALMA array operations site | 5036 m | 5058.7 m |

The ALMA check matters because it is on the other side of the planet and
exercises the same code path.

---

## 7. What is assumed, stated plainly

**No phase-stability measurement exists for Ladakh.** ALMA's structure function
comes from over 17 000 observations by an operating interferometer; nothing
comparable can exist at a site with no interferometer. The coherence figures in
the reports (0.93 extended, 0.99 compact) are a **transfer from Chajnantor**
and are labelled as a model. They must not be quoted as a site measurement, and
they are the single largest unknown in this design.

**Elevation is now included.** Station heights enter the baseline through the
``Up`` term, so ``(u, v)`` and ``w`` both carry the real relief; the coplanar
approximation used by the controlled studies elsewhere in this repository is
not used here. It matters: 200 m of relief across 3 km shifts the UV track by
about 2.3 % of its extent, while 20 m shifts it by under 0.5 %. Fixing this
also *exposed* the coverage gap described in §3, which the wrong physics had
been masking — the constraints interact.

**The terrain model is 16 m/pixel.** It establishes that a region is broadly
flat. It cannot see a boulder field, a gully, or unstable ground.

**Not modelled at all:** land ownership, protected-area boundaries, existing
structures, access roads, power routing, water courses, permafrost, avalanche
and rockfall exposure, geotechnical bearing capacity, snow loading, wind.

---

## 8. Files

| File | Role |
|---|---|
| `site_data.py` | published coordinates and measured PWV for every candidate site, each with its source |
| `terrain.py` | elevation download and caching, slope, buildability, coordinate transforms |
| `design.py` | candidate generation, the optimiser call, reporting, CASA output |
| `run_design.py` | runs both sites in both configurations, writes everything |

Tests are in `tests/test_ladakh_design.py` at the repository root.
