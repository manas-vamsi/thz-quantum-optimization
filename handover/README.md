# Handover: what was asked, what was built, and what to do next

**For:** Vishwas, Space Applications Centre, ISRO
**From:** Manasa Vamsi
**Repository:** <https://github.com/manas-vamsi/thz-quantum-optimization>

This folder exists so you do not have to read 56 modules to find out what
happened. It is the whole story in about ten minutes, with pointers to the
evidence for anything you want to check.

| Document | What it covers |
|---|---|
| **README.md** (this file) | The request, the answer, and the three results that change the plan |
| [`RESULTS.md`](RESULTS.md) | Every number, with the script that produced it |
| [`HOW_TO_RUN.md`](HOW_TO_RUN.md) | Reproduce anything, or design an array yourself, in five commands |
| [`WHAT_REMAINS.md`](WHAT_REMAINS.md) | Open questions, what is deliberately unfinished, and how to continue |

---

## 1. What you asked for

> *the mathematical formulation of the objective function for QUBO application*

**Delivered.** The derivation is in
[`docs/objective_function.md`](../docs/objective_function.md), the implementation
in `src/thz_opt/qubo/`, and the write-up in Section 2 of
[`paper/manuscript.pdf`](../paper/manuscript.pdf).

Two things about the onboarding note fell out of doing it properly.

### Your §4 claim is correct, and now proven

The note states that minimising UV auto-correlation minimises PSF sidelobes.
That is true, and it follows from Parseval's theorem: the energy in the dirty
beam equals the sum of squared gridded sampling counts. Verified numerically —
the ratio comes out **1.000000000000**.

This matters because it licenses the whole objective. Without it, "minimise
sidelobes" and "minimise UV overlap" are two different goals that happen to
sound similar.

### Your `Q_ij` recipe cannot be written as stated

The note's recipe puts baseline-overlap terms into `Q_ij`. It cannot hold them.
Whether two *baselines* overlap in the UV plane depends on four antennas —
i, j, k and l — and `Q_ij` has two index slots. Written directly, the objective
is **quartic** in the pad variables, not quadratic, so it is not a QUBO.

**The fix, and the core contribution of the work:** introduce one activation
variable `y_k` per antenna pair and enforce `y_k = x_i · x_j` through Rosenberg
quadratization. Overlap terms then become plain `y_k y_l` couplings, and the
whole thing is quadratic again.

The cost is honest and stated: `M + M(M−1)/2` variables instead of `M`. For the
full 174-pad ALMA field that is about 15,000 variables.

The penalty weights this needs are **derived, not chosen** —
`rosenberg_penalty_floor()` and `selection_penalty_floor()` compute the smallest
values that provably cannot be cheated, and there is a test confirming that
below them the constraint actually breaks.

---

## 2. Three results that should change the plan

These are the reason this handover exists. Two of them are uncomfortable, and
all three are backed by evidence you can re-run.

### 2.1 The QUBO is correct — and it loses to a simple classical search

The formulation is sound. On real ALMA geometry I enumerated **all 2,324
feasible selections** of a 16-pad instance: the QUBO's own optimum *is* the true
optimum. The **formulation gap is zero.**

But solving it is another matter:

| Instance | Proven optimum | QUBO annealing | Direct classical search |
|---|---:|---:|---:|
| 16 pads, pick 6 | 87 cells | 79 (90.8%) | **85 (97.7%)** |
| 24 pads, pick 8 | 162 cells | 149 (92.0%) | **160 (98.8%)** |
| 40 pads, pick 10 | 198 cells | 170 (85.9%) | **194 (98.0%)** |

The classical heuristic is **more accurate and 35 to 72 times faster**, and the
gap widens with problem size. Every annealing sample was feasible, so this is
not a tuning failure — the formulation works, the route does not pay.

**What this means for the project:** a quantum result here has to beat a
0.35-second classical method that is already within 2% of proven optimal. That
is a demanding bar, and it should be stated before anyone commits hardware time.

*Evidence:* `experiments/12_qubo_solver.py`, manuscript §4.2.2

### 2.2 The golden ratio is not where the performance comes from

The project started from the golden-spiral paper. Named layout families change
two things at once — how radius grows with index, and how bearing advances — so
a ranking of families cannot say which one matters.

I crossed five radial laws with five angular laws, held everything else fixed,
and decomposed the variance:

| Source | Share of variance in UV coverage |
|---|---:|
| **Radial law** | **84.6%** |
| Angular law | 11.2% |
| Interaction | 3.5% |

Concretely: uniform-in-area beats a Gaussian radial profile by a factor of
**2.7**. The golden angle beats a *deliberately degenerate* rational angle —
one that only ever uses five bearings — by **2.5%**.

**The radial distribution is what matters. The angular law is a second-order
correction.** Sidelobes behave differently (35.7% / 31.9%, with a large
interaction), so this does not transfer between objectives.

*Evidence:* `experiments/16_factorial_layout.py`, manuscript §4.6.1

### 2.3 Hanle is not comparable to Chajnantor on water vapour

The onboarding note describes Hanle as *"a comparable high-altitude, low-PWV,
radio-quiet site"* to ALMA's plateau. Against the note's own threshold of
precipitable water vapour below 1 mm:

| Site | Time below 1 mm PWV | Infrastructure |
|---|---:|---|
| ALMA, Chajnantor | ~50% | operating array |
| Ladakh site A (34.25°N, 78.75°E) | **23%** | none |
| Ladakh site B (32.50°N, 79.00°E) | 19% | none |
| Merak | 8% | partial road |
| **Hanle** | **5%** | existing observatory |

Source: Raghunath et al. 2026 (arXiv:2604.13487), 184 months of ERA5 reanalysis.

Hanle is roughly an order of magnitude drier-limited than ALMA, and four times
worse than the best unbuilt Ladakh site.

**This does not rule Hanle out.** It is the only candidate with road, power and
a working observatory, and our terrain analysis found something in its favour
that the meteorology hides: Hanle holds all three array configurations within
**15 m** of ground relief, against 57 m at site A. Its plain is genuinely flat.

**The recommendation is about framing, not about the site.** Argue Hanle as
*infrastructure against dryness*, with the numbers stated. If we claim parity
with ALMA, a reviewer will check and we lose credibility on everything else in
the paper.

*Evidence:* `ladakh_array_design/README.md`

---

## 3. What else got built

Beyond the formulation, because the formulation alone could not be evaluated:

**A working optimiser.** The code could score a layout but never search for one.
Now: greedy removal, greedy addition, 1-swap local search, simulated annealing,
and mixed-integer programming that returns **proven** optima with a dual bound,
not merely the best found.

**Real data throughout.** 174 published ALMA pads (verified: 16,195 m longest
baseline against a published ~16 km). Measured atmospheric phase from 17,000
ALMA observations. Twenty years of water-vapour statistics. Real terrain from
SRTM for the Ladakh sites.

**A science-derived objective.** Previously the score was "occupied UV cells" —
a number invented for convenience, which makes optimising it circular. It is now
Fisher information on a real source model, so two arrays compare by the error
bar they deliver on a physical quantity.

**End-to-end imaging.** Sky model → visibilities → noise → CLEAN → compare with
truth. This was built to test whether our cheap metric means anything. It does
for structure, and **not at all for flux** — see [`RESULTS.md`](RESULTS.md).

**Instrumental physics.** Primary beam, bandwidth smearing, time smearing, and
sensitivity in mJy from measured opacity.

**A terrain-constrained Ladakh array design.** Actual pad coordinates with
latitude, longitude and elevation, in CASA format. Three configurations giving
continuous coverage from 0.096 to 5.01 arcseconds.

**A command-line tool.** `python design_array.py --site alma --n 20 --dish 12`
returns pad locations and distances in seconds. See
[`HOW_TO_RUN.md`](HOW_TO_RUN.md).

---

## 4. Scale and state

| | |
|---|---:|
| Automated tests, all passing | **272** |
| Experiments | 17 |
| Python modules | 56 |
| Manuscript | 13 pages, 17 figures |
| References, none broken | 25 |

Every figure and table in the manuscript is generated by a script in the
repository. Nothing is hand-entered.

---

## 5. What I need from you

1. **Do you want to be a coauthor?** Your entry is prepared in
   `paper/metadata.yaml` but commented out. Nobody goes on a paper without
   having read it and agreed.
2. **Seven metadata fields** block submission: affiliation, correspondence
   email, funding, acknowledgements, competing interests, target venue,
   repository DOI. `python paper/build_pdf.py --check-metadata` lists them.
3. **Which story do we tell?** "We formulated it properly, solved it, and the
   quantum route does not currently pay" is honest, defensible and publishable.
   "Quantum optimization of terahertz arrays" is not supportable by what we
   found. That is a strategy decision, and it is yours and mine, not the code's.
4. **What should this array actually observe?** This is the one input I could
   not derive. Everything downstream — which baselines matter, which
   configuration to build first — follows from it.

See [`WHAT_REMAINS.md`](WHAT_REMAINS.md) for the technical continuation.
