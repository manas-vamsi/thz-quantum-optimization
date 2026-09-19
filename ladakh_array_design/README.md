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

**A second trade the water-vapour numbers hide.** Once pads are confined to a
single landform (§3), Hanle keeps 313 of 327 candidate positions and its
extended array spans only 31 m of relief. Site A keeps 118 of 262 and still
spans 57 m, losing 0.2 km of maximum baseline and 9 % of its UV coverage in the
process. Hanle's plain is large and genuinely flat; site A's is broken ground
that happens to contain flat patches. The terrain penalty at site A is real,
and it partly offsets being four times drier.

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

## 3. Why there are two configurations

The first design run maximised occupied UV cells, and produced a **ring**: all
sixteen pads between 1064 m and 1488 m of a 1500 m field, none in the inner
half, exactly **one baseline under 300 m out of 120**.

That is not an optimiser bug. Outer annuli of a UV grid contain far more cells
than inner ones, because area grows with radius, so "maximise distinct cells"
is quietly an instruction to flee to the boundary. The array that results
resolves beautifully and is blind to anything wider than 0.6 arcsec.

Two changes follow.

**A Gaussian target UV density** (Boone 2001) replaces the bare cell count, so
inner cells are worth more and the optimiser buys short baselines. This helps
but does not solve it.

**A landform filter.** Flat ground is not one surface. Hanle's observatory
stands on an isolated summit about 220 m above a wide plain, and both are
locally flat, so a slope test accepts them equally. An early run placed
fourteen pads on the plain and one on the summit — and that single outlier
accounted for the entire 235 m of reported relief, describing an array with a
240 m climb between two of its antennas. Candidates are now confined to the
60 m elevation band holding the most buildable ground. At Hanle this costs
almost nothing (313 of 327 candidates survive, relief falls from 235 m to
31 m); at site A it removes more than half (262 to 118), which is itself the
measurement that site A's ground is broken.

**Two configurations, not one.** Sixteen antennas spread over 3 km cannot also
provide 30 m baselines — seeing a 5 arcsec source needs one, and no
rearrangement of sixteen pads over 3 km delivers both. Every real observatory
answers this the same way, by moving antennas between configurations on a
shared pad field; ALMA and the VLA both do it, and §4.5 of the manuscript
measures what reconfiguration is worth.

---

## 4. The design

### Hanle (32.7794 °N, 78.9642 °E)

| | Compact | Extended |
|---|---:|---:|
| Candidate positions | 36 | 313 |
| Baselines | 32 – 328 m | 252 – 2950 m |
| **Resolution** | 0.819″ | **0.091″** |
| **Largest angular scale** | **5.01″** | 0.64″ |
| Occupied UV cells | 1342 (32.8 %) | 2164 (52.8 %) |
| Peak sidelobe | 0.133 | 0.096 |
| Ground relief | 30 m | 31 m |

### Ladakh site A (34.25 °N, 78.75 °E)

| | Compact | Extended |
|---|---:|---:|
| Candidate positions | 123 | 118 |
| Baselines | 32 – 370 m | 127 – 2738 m |
| **Resolution** | 0.727″ | **0.098″** |
| **Largest angular scale** | **5.10″** | 1.27″ |
| Occupied UV cells | 1878 (45.8 %) | 1872 (45.7 %) |
| Peak sidelobe | 0.101 | 0.124 |
| Ground relief | 22 m | 57 m |

Hanle spans **55×** in angular scale, site A **52×**, from roughly 0.1″ detail
to 5″ structure.

Two caveats visible only in these tables.

**Hanle's compact configuration has just 36 candidate positions** within the
200 m radius, against 123 at site A: choosing sixteen antennas from thirty-six
options is a tight design space, because the observatory stands on a summit
whose plateau is small.

**The two Hanle configurations sit on different landforms.** The compact array
occupies the summit plateau at about 4490 m; the extended array occupies the
plain at about 4270 m, 220 m below. They are therefore *not* a shared,
reconfigurable pad field — moving antennas between them means a 220 m climb.
Reconfiguration as analysed in §4.5 of the manuscript assumes one field, and
that assumption does not hold at Hanle. Site A's two configurations are 60 m
apart vertically, which is far more workable.

---

## 5. Outputs

For each of `hanle_compact`, `hanle_extended`, `site_a_compact`,
`site_a_extended`:

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

**Elevation is ignored in the UV calculation.** The array is treated as
coplanar, with no *w* term. With the landform filter in place the relief is
22–57 m across every configuration, which is small against baselines of
hundreds of metres to 3 km, so the approximation is defensible here. It was not
defensible before that filter existed, and it would stop being defensible again
if the elevation band were widened.

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
