# What to optimise, and how: the array-configuration design space

This document answers a question that experiments 01–05 could not: *if not the
golden spiral, and not Fibonacci, then what?*

It is organised as four questions, in the order they have to be answered:

1. What does the source paper actually establish?
2. What is the objective? (the question almost everybody skips)
3. What geometries and optimisers does the field already use?
4. What of that can a QUBO hold exactly?

Everything marked **[measured]** is a number produced by a script in
`experiments/` in this repository and reproducible from the committed config.

---

## 1. What the source paper establishes, and what it does not

Quiroga Rodríguez, *JINST* **20** (2025) P05026, proposes a golden-ratio
logarithmic spiral with antennas on rails that rotate the spiral during an
observation. Reading the paper together with its published appendix code:

**What is there.** A single-arm logarithmic spiral,
`r = a * PHI**(theta/2pi)` with `a = 10 m`, 100 antennas, `theta` linearly
spaced to reach `r = 5 km`; all `N(N-1)/2` baselines with conjugates; a
plot of the resulting UV points; and a sky-image simulation.

**What is not there.**

- *No UV-coverage computation.* The appendix code plots UV points and never
  bins, grids or counts them. The 95 % and 60 % figures in the paper's Table 2
  are asserted, not computed by the published code, and no cell size, grid
  extent or definition of "coverage" is stated. Without those three numbers a
  filled fraction is not a measurable quantity (see
  [`uv_metrics.md`](uv_metrics.md)).
- *The resolution claim does not follow.* The paper states that rotating the
  spiral around the array centre "would increase the effective D_max, reducing
  theta", and Table 2 reports resolution improving from 0.001 to 0.0005 arcsec.
  A rigid rotation is an isometry: it maps the baseline set onto a rotated copy
  of itself and leaves every baseline *length* unchanged, so `b_max`, and hence
  `lambda / b_max`, is invariant. Rotation fills in position angles; it cannot
  improve resolution. Only radial motion of the antennas changes `b_max`.
- *The sensitivity claim is not configuration-specific.* A +41 % SNR from
  `k = 2` positions is the ordinary `sqrt(k)` radiometric gain from collecting
  twice as much data. Any array integrating twice as long gets it.
- *The imaging simulation samples the image, not its transform.* In the
  appendix code for figures 3–4 the visibility grid is filled with
  `source_model[v, u]` — values of the sky *image* array read at UV pixel
  indices — rather than with the Fourier transform of the sky sampled at
  `(u, v)`. The "dirty image" is therefore not the interferometric response of
  the array.

**What survives, and deserves to.** Two things, both worth keeping:

1. *A logarithmic radial law is a sound instinct.* It gives baselines spread
   over decades of length, which is what multi-scale imaging needs, and real
   arrays use it: ALMA's configuration design is built on Conway's "zoom
   spiral" (ALMA Memos 283, 292), SKA-Mid and the VLA both use spiral or
   quasi-spiral arms.
2. *Reconfiguration genuinely adds UV coverage.* Moving antennas between
   observations is exactly what ALMA, the SMA and IRAM do, and the
   pad-selection problem that creates is the problem this repository is about.

The specific choice of `PHI` as the growth rate, however, is decoration —
see §3.3 **[measured]**.

---

## 2. The objective is the whole problem

The field's central result is not a shape, it is a relation: **the synthesised
beam is the Fourier transform of the UV sampling density**, so the *shape* of
that density determines the beam.

- Uniformly filled UV disc → `J1(r)/r` beam: the sharpest main lobe available
  for a given `b_max`, and slowly decaying sidelobe rings.
- Gaussian UV density → Gaussian beam: **no sidelobe rings at all**, at the
  price of a broader main lobe (Boone 2002, A&A 386, 1160).
- Sharp outer cutoff → high sidelobes; apodised (smoothed) outer edge → low
  sidelobes, lower resolution (Keto 2012, JAI 1, 1250008).

So "maximise UV coverage" is not a goal, it is one corner of a trade-off
triangle — **resolution ↔ sidelobe level ↔ completeness** — and you cannot have
all three. Any project that optimises a single coverage percentage has silently
picked a corner without saying so.

The objectives that are actually used:

| objective | reference | what it favours | implemented here |
|---|---|---|---|
| unique UV cells filled | folklore | completeness | `metrics/uv_coverage.py` |
| UV-point repulsion energy `F = R² Σ 1/d²` | Cornwell 1988, IEEE AP-36 1165; Karastergiou+ 2006, ApJS 164 552 | uniform spread | `metrics/density_matching.cornwell_energy` |
| match a target radial UV density | Boone 2001/2002 | beam *shape* | `metrics/density_matching.density_match_chi2` |
| PSF peak sidelobe / RMS | Kogan 2000, IEEE AP-48 1075 | image dynamic range | `metrics/psf_metrics.py` |
| multi-objective (imaging + cable length) | Cohanim, Hewitt & de Weck 2004, ApJS 154 705 | buildability | not implemented |
| minimum variance to a target distribution | Panduranga Rao et al. 2009, arXiv:0901.4901 | pad selection from a larger set | closest to our QUBO |

**These objectives disagree, and the disagreement is measurable.**
At N = 20, same extent, same grid, same observing geometry **[measured,
experiment 06]**:

| objective | best layout of the twelve tested |
|---|---|
| unique UV cells (max) | Reuleaux, single ring |
| Cornwell F (min) | Reuleaux, single ring |
| Gaussian density χ² (min) | golden spiral |
| PSF peak sidelobe (min) | golden spiral |
| PSF FWHM (min) | Reuleaux, single ring |

The golden spiral is not bad — it wins two of the five. It simply is not the
answer to a question nobody has posed yet.

---

## 3. The geometries the field actually uses

### 3.1 Ranked on this repository's own runs

N = 20, r ∈ [20, 1000] m, 300 GHz, shared 256² grid, 4 h of Earth rotation
**[measured, experiment 06]**:

| layout | unique cells | Cornwell F ↓ | Gauss χ² ↓ | PSF sidelobe ↓ | min sep (m) |
|---|---|---|---|---|---|
| **Reuleaux (1 ring)** | **10226** | **30.5** | 1.05 | 0.132 | 230 |
| hierarchical 3³ (Keto) | 9596 | 208.9 | 0.93 | 0.201 | 144 |
| random (area-uniform) | 9266 | 40.0 | 0.067 | 0.178 | 65 |
| Vogel golden angle | 9018 | 34.6 | 0.251 | 0.124 | 172 |
| Gaussian (Boone target) | 7844 | 46.3 | **0.096** | 0.151 | 118 |
| log spiral, 3 arms | 7736 | 70.1 | 0.214 | 0.098 | 35 |
| Reuleaux (3 rings) | 7708 | 50.7 | 0.224 | 0.105 | 179 |
| log spiral, 5 arms | 7644 | 111.0 | 1.539 | **0.090** | 24 |
| **golden spiral (1 arm)** | 6872 | 74.9 | **0.023** | **0.090** | 37 |
| uniform rings | 5116 | 308.1 | 5.75 | 0.159 | 17 |
| Fibonacci radial | 1813 | 1322.6 | 14.4 | 0.386 | 9 |

Headline: **the Reuleaux triangle covers 49 % more unique UV cells than the
golden spiral and has less than half its Cornwell energy**, which is exactly
what Keto (1997, ApJ 475, 843) predicted when he showed that curves of constant
width give the most uniform snapshot sampling — the reason the Submillimeter
Array's pads are four nested Reuleaux triangles.

### 3.2 The families, and when each is right

- **Curves of constant width (Reuleaux)** — most uniform UV sampling, circular
  beam, best when antennas are few and every baseline is precious. Weakness: it
  is a ring, so there are no short baselines (min separation 230 m here) and no
  compact core. Bad for extended emission.
- **Multi-arm logarithmic spirals** — the ALMA/VLA/SKA answer. Log radial law
  gives baselines across decades; several arms fix the angular coverage that a
  single arm cannot provide. Naturally supports "zoom" reconfiguration.
- **Gaussian (Boone target)** — draw antennas from a 2D Gaussian and the UV
  density is *exactly* Gaussian (difference of Gaussians is Gaussian, with
  `sigma_uv = sqrt(2) sigma_ant`). Lowest-sidelobe beam by construction, no
  optimiser needed. Best χ² of any structured layout here.
- **Hierarchical / self-similar (Keto 2012)** — repeat one good sub-array at
  several scales. Smooth multi-scale coverage without numerical optimisation,
  and a tunable resolution-vs-sidelobe knob (the inter-level scale factor).
  Weakness, visible above: repeating a pattern repeats baselines, giving 67 %
  redundancy and a Cornwell F of 209.
- **Minimum-redundancy / perfect non-redundant arrays** — Golay (1971); the
  largest *perfect* non-redundant 2D array has just six elements (the "Manx"
  array, McKay et al. 2022, Radio Science 57). Exact non-redundancy does not
  scale, so it matters for small N only.
- **Random / quasi-random** — a strong baseline in every sense. It beat the
  golden spiral on cell count at every N tested **[measured, experiment 01]**
  and comes third here. Incoherent sampling is also what compressed-sensing
  reconstruction theory prefers.

### 3.3 Is the golden ratio special? (direct test)

The usual argument for `PHI` is that it is the *most badly approximable*
irrational, so a golden-angle step avoids commensurate (repeating) position
angles for as long as possible. That argument is about irrationality, not about
`PHI`. Scanning the phyllotaxis angle `alpha` at fixed N, extent and grid
**[measured, experiment 06]**:

| alpha | unique cells | Cornwell F |
|---|---|---|
| √3 − 1 = 0.7321 | 9382 | 37.8 |
| 1/e = 0.3679 | 9338 | 33.6 |
| 5/13 = 0.3846 (Fibonacci ratio) | 9328 | 33.4 |
| √2 − 1 = 0.4142 | 9240 | 34.1 |
| 1/π = 0.3183 | 9180 | 32.4 |
| **1/Φ² = 0.3820 (golden angle)** | **9018** | **34.6** |
| 21/55 = 0.3818 (Fibonacci ratio) | 8996 | 34.2 |
| 1/3 (rational) | 8672 | 46.8 |
| 1/2 (rational) | 5252 | 238.0 |

The golden angle ranks **sixth of ten**. Every tested irrational is within a
few per cent of it; only low-order rationals are clearly bad, and they are bad
because they collapse the position angles onto `q` directions. What matters is
"not a low-order rational", not "golden". Note also that 21/55 — a Fibonacci
ratio, and therefore rational — scores essentially the same as `PHI` itself,
because at N = 20 the array cannot tell 0.3818 from 0.3820.

**Conclusion for the project.** Neither the golden ratio nor Fibonacci has a
demonstrated interferometric virtue. They are reasonable *initialisations*.
Treating them as the design is what should be dropped.

---

## 4. Optimisation methods

Closed-form curves are initial conditions. The literature optimises:

1. **Cornwell (1988)** — treat UV points as repelling charges, minimise
   `F = R² Σ 1/|V_i − V_j|²`. Still the standard scalar objective.
2. **Boone (2001, 2002)** — move pads to match a target UV density; the ALMA
   lineage.
3. **Karastergiou, Neri & Gurwell (2006)** — the *pad-selection* version:
   `N` antennas, `N_p > N` pads, search which pads to occupy by minimising
   Cornwell's `F`. This is our problem, solved by exhaustive/greedy search.
4. **Kogan (1997–2000)** — direct sidelobe minimisation with terrain
   constraints (ALMA memos).
5. **Cohanim, Hewitt & de Weck (2004)** — multi-objective GA + simulated
   annealing, imaging performance against cable length, for 27–160 stations.
6. **Panduranga Rao et al. (2009)** — greedy "minimum variance" removal from
   `M` pads down to `N`, `O((M−N)M³)`, targeting a Gaussian UV distribution.
7. **de Villiers (2007)** — tomographic projection: compare 1-D projections of
   achieved and ideal UV distributions, map the discrepancy back to positions.

None of these is a quantum method, and no prior application of QUBO or quantum
annealing to radio-interferometer configuration surfaced in the searches done
for this document. That is the opening this project has — but it also means the
classical baselines above are the bar to beat, and they are strong.

**What is missing from this repository:** a classical optimiser. Greedy removal
(method 6) is ~50 lines and would immediately show how much of the gap between
"golden spiral" and "best found" is real. It should exist before any solver
comparison is claimed.

---

## 5. What a QUBO can hold exactly

This is the part where the project can make a genuine contribution, and the
answer changed once the problem was looked at in the right variables.

### 5.1 The structural fact

In pad-selection variables `x`, the only exactly-quadratic objectives are those
that decompose over *single* baselines, `Σ w_ij x_i x_j`. Every objective in §2
couples *different baselines to each other* and is therefore quartic in `x`.

But introduce the baseline-activation variables

```
y_k = x_i x_j          (one per candidate pair, M(M-1)/2 of them)
```

and all three of the main objectives become **exactly quadratic in `y`**:

| objective | form in `y` |
|---|---|
| Cornwell repulsion | `Σ_k E_kk y_k + Σ_{k<l} E_kl y_k y_l` |
| target-density χ² | expand `Σ_b w_b (Σ_k a_bk y_k − t_b)²`, using `y² = y` |
| unique-cell coverage, to 2nd order | `Σ_k |C_k| y_k − Σ_{k<l} |C_k ∩ C_l| y_k y_l` |

and `y_k = x_i x_j` is itself enforced exactly by the Rosenberg penalty
`x_i x_j − 2 x_i y_k − 2 x_j y_k + 3 y_k`, which is 0 iff the identity holds
and ≥ 1 otherwise. So the model is a genuine QUBO over `M + M(M−1)/2`
variables **with no approximation of the objective**.

### 5.2 The coverage term is a Bonferroni bound

The third row is worth stating properly. By the second Bonferroni inequality,
for any sets `C_k`,

```
|∪_k C_k|  ≥  Σ_k |C_k| − Σ_{k<l} |C_k ∩ C_l|
```

so the quadratic form is a **provable lower bound** on unique UV-cell
coverage, and it is *exact* whenever no cell is covered three or more times.
This replaces the ad-hoc "discount contested cells" weighting of
`qubo/coefficients.py` with something that has a known direction of error and a
known exactness condition.

### 5.3 Measured cost and accuracy

Toy problem M = 8, N = 4, Earth-rotation sampling, all 70 feasible selections
enumerated **[measured, experiment 07]**:

| cells/axis | A: pairwise ρ | B: Bonferroni ρ | A finds optimum | B finds optimum | bound gap | B fill |
|---|---|---|---|---|---|---|
| 8 | 0.896 | 0.910 | no | **yes** | 24.9 % | 29 % |
| 16 | 0.964 | 0.994 | no | **yes** | 1.9 % | 10 % |
| 32 | 0.981 | 0.996 | no | **yes** | 0.3 % | 7 % |
| 64 | 0.991 | **1.000** | yes | **yes** | 0.05 % | 4 % |
| 128 | 0.993 | **1.000** | yes | **yes** | 0.01 % | 2 % |

B finds the true optimum at **every** resolution tested; A misses it at four of
eight. Variable cost at 16 cells/axis:

| formulation | variables |
|---|---|
| A pairwise surrogate (pad space) | 8 |
| **B Bonferroni bound (baseline space)** | **36** (37 couplings, 10 % fill) |
| C exact coverage (cell + slack variables) | 280 |

and B scales as `M + M(M−1)/2`: 210 variables at M = 20, 1275 at M = 50, 5050
at M = 100, 20 100 at M = 200. The coupling matrix stays sparse because two
baselines interact only if their UV tracks actually overlap.

**Recommendation.** Make B the primary formulation. Keep A as a cheap
baseline and C as ground truth for validation. Note that B's Cornwell variant
is *dense* (every baseline pair repels), so for that objective the distance
cutoff in `cornwell_terms(coupling_cutoff=...)` is not optional at scale.

### 5.4 Why this problem is worth a solver at all

With the cardinality constraint `Σ x_i = N` and a pairwise reward, the pad
selection problem is exactly **densest-k-subgraph** on a geometric graph: pick
`k` vertices maximising the weight of the induced edges. DkS is NP-hard, has no
known constant-factor approximation, and cardinality-constrained QUBO instances
of it are an active target for quantum and quantum-inspired solvers. So the
motivation for a QUBO here is structural, not fashionable — which is worth
stating explicitly in any write-up, because it is the part a referee will test.

---

## 6. Recommended order of work

1. **Add a classical optimiser** — greedy removal against the Bonferroni
   objective, then simulated annealing. Without it there is no baseline and no
   way to know whether any configuration is good.
2. **Fix a science case and derive a target UV density from it.** Until then,
   every "best layout" here is best only under an arbitrary metric.
3. **Factorial layout experiment** — radial law × angular law, separately, to
   attribute the differences seen in §3.1.
4. **Use the baseline-space (Bonferroni) QUBO** as the model, at M ≈ 50–100
   pads, benchmarked against greedy and SA before any quantum solver is
   involved.
5. **Then** quantum: `dimod`/`dwave-neal` exact on the 36-variable toy,
   hybrid at M = 100 (5050 variables is within reach of hybrid solvers, not of
   a bare Advantage QPU at this density).
6. **Replace the phase and PWV placeholders** with something physically
   defensible before any THz-specific claim is made.

## References

- Boone, F. 2001, A&A 377, 368 — *Interferometric array design: optimizing the locations of the antenna pads*
- Boone, F. 2002, A&A 386, 1160 — *Interferometric array design: distributions of Fourier samples for imaging*
- Cohanim, B. E., Hewitt, J. N. & de Weck, O. 2004, ApJS 154, 705 — *The design of radio telescope array configurations using multiobjective optimization*
- Conway, J. 1998/2000, ALMA Memos 283 & 292 — spiral / "zoom" configurations
- Cornwell, T. J. 1988, IEEE Trans. Ant. Prop. 36, 1165 — *A novel principle for optimization of the instantaneous Fourier plane coverage of correlation arrays*
- de Villiers, M. 2007, A&A 469, 793 — *Interferometric array layout design by tomographic projections*
- Golay, M. J. E. 1971, JOSA 61, 272 — *Point arrays having compact, nonredundant autocorrelations*
- Karastergiou, A., Neri, R. & Gurwell, M. A. 2006, ApJS 164, 552 — *Adapting and expanding interferometric arrays*
- Keto, E. 1997, ApJ 475, 843 — *The shapes of cross-correlation interferometers*
- Keto, E. 2012, JAI 1, 1250008 (arXiv:1209.0692) — *Hierarchical configurations for cross-correlation interferometers with many elements*
- Kogan, L. 2000, IEEE Trans. Ant. Prop. 48, 1075 — *Optimizing a large array configuration to minimize the sidelobes*
- McKay, J. et al. 2022, Radio Science 57, e2022RS007500 — *Manx arrays: perfect non-redundant interferometric geometries*
- Panduranga Rao, M. V. et al. 2009, arXiv:0901.4901 — *A minimum variance method for problems in radio antenna placement*
- Quiroga Rodríguez, E. 2025, JINST 20, P05026 — *Golden spiral radio interferometry*
- Su, Y., Nan, R. D. & Peng, B. 2004, ChJAA 4, 198 — uniform-weight pad selection
- Treloar, N. C. 1989, JRASC 83, 92 — early comparison of array configurations
