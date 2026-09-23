# Every result, and the script that produced it

Nothing here is hand-entered. Each row names the experiment that generates it,
so any number can be re-derived rather than trusted.

---

## 1. The data is real, and checked against published values

| Quantity | Computed here | Published | Source |
|---|---:|---:|---|
| ALMA 12 m pads | 174 | 174 | CASA `alma.all.cfg` |
| ALMA longest baseline | **16,195 m** | ~16,000 m | ALMA extended configuration |
| VLA A-config longest baseline | **36,623 m** | ~36,000 m | VLA A configuration |
| Hanle elevation from terrain model | **4,491 m** | 4,500 m | IAO published |
| ALMA site elevation from terrain model | **5,036 m** | 5,058.7 m | ALMA published |

The VLA and ALMA-elevation rows are the load-bearing ones. VLA coordinates are
geocentric and had to be rotated into a local frame; recovering 36 km
independently proves the transform rather than the file. The ALMA elevation is
on the far side of the planet from Hanle and exercises the same terrain code.

**Atmosphere.** ALMA Memo 624 publishes its own scaling factors between summary
baselines. Our implementation reproduces them:

| Baselines | Our model | Memo 624 |
|---|---:|---:|
| 500 → 1000 m | 1.516 | 1.51 |
| 1000 → 5000 m | 1.595 | 1.59 |
| 5000 → 10000 m | 1.223 | 1.22 |

*Scripts:* `experiments/11_real_alma.py`, `tests/test_real_arrays.py`,
`tests/test_coherence.py`

---

## 2. Searching beats hand-drawn curves — on real pads

Selecting 20 of the real ALMA pads, every option constrained to be actually
buildable:

| Selection | Occupied UV cells | Feasible? |
|---|---:|:---:|
| Reuleaux curve, idealised | 916 | **no** |
| Reuleaux, snapped to real pads | 681 | yes |
| Golden spiral, idealised | 693 | **no** |
| Golden spiral, snapped to real pads | 463 | yes |
| Best of 16 random feasible selections | 301 | yes |
| **Searched (greedy + local search)** | **750** | yes |

Search improves on the best feasible analytic layout by **10.1%**, and on random
feasible selection by **149%**.

The idealised Reuleaux scores higher at 916 — but it places antennas where no
pad exists. The gap between 916 and 750 is the price the real pad field charges,
and no better curve recovers it.

**Qualification (see §5):** this is a gain in *structural* accuracy, not in
flux recovery.

*Script:* `experiments/11_real_alma.py`

---

## 3. The QUBO: correct formulation, losing route

**Formulation gap — zero.** On a 16-pad real instance, all **2,324** feasible
selections were enumerated. The QUBO's own optimum coincides with the true
coverage optimum. Correlation between QUBO score and true coverage: **0.929**.

**Solver gap — large, and widening:**

| Instance | QUBO variables | Couplings | Proven optimum | Annealing (500 reads) | Direct classical | Speed ratio |
|---|---:|---:|---:|---:|---:|---:|
| M=16, N=6 | 136 | 1,506 | 87 | 79 (90.8%) | **85 (97.7%)** | 35× |
| M=24, N=8 | 300 | 5,182 | 162 | 149 (92.0%) | **160 (98.8%)** | 63× |
| M=40, N=10 | 820 | 33,689 | 198 | 170 (85.9%) | **194 (98.0%)** | 72× |

**Every returned sample was feasible** at every size — correct cardinality, no
violated Rosenberg relation. The derived penalty floors hold against a real
solver, and a test confirms that lowering them below the floor does break
feasibility, so the check is informative rather than vacuous.

**Hardware implication.** A 40-pad instance already needs 820 variables and
33,689 couplings on a dense graph. The full 174-pad ALMA field would need about
15,000 variables. Current annealing hardware has sparse fixed connectivity, so
minor-embedding would multiply the physical qubit count substantially.

*Script:* `experiments/12_qubo_solver.py`

---

## 4. Radial law versus angular law

Five radial laws crossed with five angular laws, everything else held fixed,
five jittered replicates per cell so the decomposition has a real error term.

| Source | UV coverage | F | Peak sidelobe | F |
|---|---:|---:|---:|---:|
| **Radial law** | **84.6%** | 3049 | 35.7% | 82.9 |
| Angular law | 11.2% | 404 | 31.9% | 73.9 |
| Interaction | 3.5% | 31.5 | 21.6% | 12.5 |
| Replication noise | 0.7% | — | 10.8% | — |

Marginal means, averaged over the other factor:

| Radial law | Mean cells | | Angular law | Mean cells |
|---|---:|---|---|---:|
| **Uniform in area** | **2054** | | **Golden angle** | **1497** |
| Uniform in radius | 1635 | | Rational angle (5 bearings only) | 1461 |
| Power law | 1191 | | Multi-arm | 1393 |
| Exponential | 1104 | | Uniform random | 1358 |
| Gaussian | 761 | | Equal spacing | 1037 |

Choosing uniform-in-area over Gaussian is a factor of **2.7**. The golden angle
beats a deliberately degenerate rational angle by **2.5%**.

The ANOVA is verified against tables with a pure row effect, a pure column
effect, and a pure interaction with zero marginals — a decomposition that
misattributed variance would give a confident and entirely wrong answer about
which design choice matters.

*Script:* `experiments/16_factorial_layout.py`

---

## 5. Does UV coverage actually predict image quality?

This tests the assumption the entire framework rests on. Sky model →
visibilities → noise → Högbom CLEAN → compare with truth.

**First, the transform pair is verified.** Complete sampling returns the input
sky to **2×10⁻¹⁶** with a single-pixel delta beam. An inconsistent pair would
bias every fidelity number by the same silent factor rather than failing.

| Sky model | Cells predict structure | Cells predict flux |
|---|---:|---:|
| Point source | **+1.000** | +1.000 |
| Close double | **+1.000** | +0.600 |
| Smooth Gaussian | **+1.000** | **0.000** |
| Disc with gap | **+1.000** | +0.300 |

**Occupied UV cells predict structural accuracy perfectly.** Rank correlation of
exactly +1.000 against relative RMS image error, for every sky model. The cheap
metric used throughout the project is validated for that purpose.

**It says nothing about flux on extended sources:**

| Array | Cells | Gaussian flux recovered |
|---|---:|---:|
| Searched, maximum cells | 5204 | **42%** |
| Golden spiral, snapped to pads | 3350 | **91%** |

The cell-maximising array has 55% more cells and recovers **less than half** the
source. Maximising cells drives antennas to long baselines, which resolves
extended emission out — and no cell count, nor any objective built on one, can
see that happening.

**This directly qualifies §2.** The 10.1% improvement is real and is
structural; it is not a gain in flux.

*Script:* `experiments/14_imaging_validation.py`

---

## 6. Instrumental limits

Which effect binds the usable field of view is not fixed and cannot be assumed:

| Observing setup | Bandwidth limit | Time limit | Primary beam | Binds |
|---|---:|---:|---:|---|
| Spectral line (0.05%, 2 s) | 108″ | 427″ | 9.7″ | primary beam |
| Narrow continuum (1%, 6 s) | 5.4″ | 142″ | 9.7″ | **bandwidth** |
| Wide continuum (8%, 30 s) | **0.67″** | 28″ | 9.7″ | **bandwidth** |

Response measured across the field with identical sources, 12 m dishes,
300 GHz:

| Setup | 0″ | 2.5″ | 5″ | 10″ | 15″ |
|---|---:|---:|---:|---:|---:|
| Idealised | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Narrow continuum | 0.997 | 0.972 | 0.944 | 0.889 | 0.838 |
| Wide continuum | 0.980 | 0.781 | 0.559 | 0.293 | 0.202 |
| Wide continuum + primary beam | 0.980 | 0.749 | 0.463 | **0.138** | **0.039** |

**The design consequence:** both smearing radii are proportional to the
synthesised beam, and the beam shrinks as the array grows. **Extending an array
for resolution shrinks its usable field in exact proportion.** None of the
coverage objectives can express that cost.

**Sensitivity, from measured opacity.** Twenty 12 m antennas, 8 GHz bandwidth,
60° elevation:

| Condition | zenith τ | T_sys | SEFD | 1 h rms |
|---|---:|---:|---:|---:|
| Hanle, best winter decile | 0.06 | 86 K | 3001 Jy | 0.029 mJy |
| Chajnantor, 1 mm PWV | 0.10 | 103 K | 3575 Jy | 0.034 mJy |
| Hanle, typical | 0.15 | 124 K | 4332 Jy | 0.041 mJy |

*Script:* `experiments/15_instrumental_limits.py`

---

## 7. Atmospheric coherence at 300 GHz, measured

Applied to the searched 20-pad ALMA layout using the four measured regimes:

| Condition | Path RMS at 1 km | Mean coherence | Baselines below 0.5 |
|---|---:|---:|---:|
| No WVR correction | 200 μm | 0.174 | 98% |
| No WVR, below median PWV | 115 μm | 0.544 | 39% |
| **WVR corrected, wind < 10 m/s** | **70 μm** | **0.738** | **0%** |
| WVR corrected, wind > 10 m/s | 140 μm | 0.316 | 90% |

Water-vapour-radiometer correction is the difference between an array where
almost nothing is coherent and one where everything is, and **wind speed alone**
moves the result across most of that range.

**A correction worth recording.** An earlier version of this analysis assumed a
1 mm path RMS at 1 km and concluded that long baselines were unusable at
300 GHz. The measured value is 70 μm — fourteen times smaller. That conclusion
was an artefact of the assumption, and is why the atmospheric constants had to
be measured rather than estimated.

*Script:* `experiments/11_real_alma.py`

---

## 8. The Ladakh array design

Three configurations on a shared pad field at Hanle, array centre placed on the
plain 2,049 m from the observatory marker:

| | Compact | Intermediate | Extended |
|---|---:|---:|---:|
| Baselines | 32–232 m | 64–878 m | 338–2788 m |
| Resolves to | 1.158″ | 0.306″ | **0.096″** |
| Sees up to | **5.01″** | 2.50″ | 0.48″ |
| Ground relief | 11 m | 9 m | 15 m |

Continuous coverage from **0.096″ to 5.01″** — a 52× range — with measured
overlap at both joints, so nothing falls between configurations.

**Three configurations, not one, and not two.** A single 16-antenna array over
3 km cannot also provide 32 m baselines. Two configurations only give continuous
coverage if the smaller *resolves* finer than the larger can *see*; at 400 m the
compact array failed that test, and widening it to 800 m closed the gap but
pushed its own shortest baseline from 32 m to 90 m, collapsing its reach from
5.0″ to 2.5″. One array cannot be both the short-spacing array and the bridge.

*Script:* `ladakh_array_design/run_design.py`
