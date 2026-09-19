# Project brief — full context for review

**Purpose of this document.** A self-contained handoff: what the project is,
what has been built, what is known, what is unresolved, and what the plausible
next methods are. Written so that a reviewer with no prior context can read it
once and give useful technical advice. Weak points are stated rather than
hidden — the point is to get them attacked.

Repository: `manas-vamsi/thz-quantum-optimization` (private).
175 automated tests, all passing. The QUBO is now solved with an annealing
sampler (`dwave-samplers`, classical) and benchmarked against MILP-certified
optima; no quantum hardware has been used.

---

## 1. The project

**Goal.** Design the antenna layout of a terahertz / sub-millimetre
interferometric array by formulating the layout choice as a QUBO (Quadratic
Unconstrained Binary Optimization) problem, so that it can eventually be solved
on quantum hardware (annealing or QAOA).

**Setting.** Informal research collaboration with a scientist at ISRO Space
Applications Centre. The output is intended as a paper; candidate venues named
in the project's onboarding guide are A&A / MNRAS (astronomy impact), IEEE
Transactions on Antennas and Propagation (engineering), and QIP / npj Quantum
Information (algorithmic).

**The immediate deliverable that was requested:** *"come up with an idea for a
mathematical formulation of the objective function for QUBO application."*
Status: done, derived and numerically validated — see §5.

---

## 2. The physical problem

An interferometer does not have one big dish; it has `N` small antennas whose
**pairs** each measure one Fourier component of the sky. For antennas at
positions `r_i`, the pair `(i,j)` measures spatial frequency

```
(u, v) = (r_j − r_i) / λ          [units: wavelengths]
```

The set of all sampled `(u,v)` points is the **UV coverage**. The image the
array produces is the Fourier transform of what it sampled, so:

- **UV coverage → image quality.** Gaps in UV coverage become artefacts.
- The synthesised beam (**PSF**, "dirty beam") is the Fourier transform of the
  UV sampling density.
- Earth rotation sweeps each baseline through an arc over hours, filling more
  of the UV plane ("aperture synthesis").

**The design problem.** There are `M` possible concrete pads. Only `N` antennas
are affordable. Which `N` pads?

`M choose N` possibilities. For `M=100, N=50` that is ~10²⁹. It is a
combinatorial selection problem whose objective depends on *pairs*, which is
why it maps naturally onto a QUBO.

**Terahertz makes it harder** (this is the project's specific angle):

1. **Water vapour** — the atmosphere is opaque at THz unless PWV < ~1 mm.
   `τ₂₂₅ ≈ 0.049·PWV[mm] + 0.018`.
2. **Phase coherence** — turbulent water vapour scrambles phase; coherence
   times drop below a second on long baselines. This is the dominant THz
   constraint (see §6.2).
3. **Shadowing** — dishes block each other at low elevation, so
   `d_min ≳ 1.5 × dish diameter`.
4. **Discrete pads** — antenna transporters can only reach prepared pads.

---

## 3. The QUBO formulation

Binary variable per pad: `x_i = 1` if pad `i` is occupied. A baseline exists
iff `x_i x_j = 1`.

```
H_total = H_select + H_shadow + H_pwv + H_uv + H_phase
```

Convention used (upper triangular, one entry per unordered pair — this is also
`dimod`'s convention):

```
H(v) = Σ_a Q_aa v_a + Σ_{a<b} Q_ab v_a v_b + const
```

| term | form | degree in `x` |
|---|---|---|
| `H_select` | `λ_sel (Σx_i − N)²` → `λ_sel[(1−2N)Σx_i + 2Σ_{i<j}x_i x_j + N²]` | 2 ✓ |
| `H_shadow` | `λ_sh Σ_{d_ij < d_min} x_i x_j` | 2 ✓ |
| `H_phase` | `λ_ph Σ (1 − γ_ij) x_i x_j`, γ = coherence factor | 2 ✓ |
| `H_pwv` | `λ_pwv Σ c_i x_i` | 1 ✓ |
| `H_uv` | UV coverage / sidelobe quality | **4** ✗ |

**`H_uv` is the whole difficulty** and is discussed in §5.

---

## 4. What has been built

Full classical pipeline, all validated:

```
layout generator → baselines → (u,v) → Earth rotation → UV cell grid
                 → coverage metrics → dirty beam (PSF) → QUBO → validation
```

**Modules**
- `arrays/` — golden spiral, 3 Fibonacci variants, multi-arm log spiral,
  Reuleaux triangle, hierarchical (Keto), Gaussian, phyllotaxis, random, rings
- `interferometry/` — baselines, UV mapping, Earth-rotation synthesis
  (Thompson–Moran–Swenson convention), dirty beam via FFT
- `metrics/` — unique cells, filled fraction, redundancy, density/angular
  uniformity, baseline-length statistics, power-law fit, PSF sidelobes/FWHM/
  ellipticity, Cornwell repulsion energy, Gaussian-density matching
- `constraints/` — minimum separation, PWV (synthetic), phase coherence
  (Kolmogorov + decorrelation), graph stepping-stone placeholder
- `qubo/` — pad-space QUBO, baseline-space QUBO, exhaustive enumeration,
  validation

**Validation already passing**

| check | result |
|---|---|
| Matrix QUBO vs direct term-by-term objective, all 2⁸ configurations | difference **0.0** |
| Upper-triangular vs symmetric convention | difference **0.0** |
| Earth rotation reduces to snapshot at zenith | 2×10⁻¹⁰ λ |
| Dirty beam is real (Hermitian sampling) | imaginary part 5×10⁻¹⁷ |
| Sidelobe-energy quadratic form vs direct gridding | difference **0.0** |
| Parseval identity | ratio **1.000000000000** |
| Bonferroni bound ≤ exact coverage | holds always |

---

## 5. The central mathematical result

### 5.1 The obstruction

The guide specifies `H_uv` as: *score each baseline by the unique UV cells it
samples, penalise it if its track overlaps another pair's track, load the
result into `Q_ij`.*

The first half works. The second does not: "pair `(i,j)` overlaps pair `(k,l)`"
involves **four** pad indices. `Q_ij` holds two. Written in `x`, the overlap
term is degree-4 — not a QUBO term.

Formally: with `C_k` the set of UV cells touched by baseline `k`, coverage is
`|∪_k C_k|` over active baselines. The union is an OR of products, and
`z_c = 1 − Π(1 − y_k)` expands to arbitrarily high degree.

### 5.2 The resolution

Introduce baseline-activation variables `y_k = x_i x_j` as genuine variables,
enforced exactly by the **Rosenberg quadratization**:

```
R_k = λ_R (x_i x_j − 2 x_i y_k − 2 x_j y_k + 3 y_k)
```

which is 0 iff `y_k = x_i x_j` and ≥ λ_R otherwise. Every term is quadratic.
In `y`, several objectives are **exactly quadratic**. Note the distinction in
the first row: what is exact there is the *bound*, not the coverage.

| objective | form in `y` | source |
|---|---|---|
| unique-cell coverage (2nd order) | `Σ|C_k| y_k − Σ|C_k∩C_l| y_k y_l` | Bonferroni inequality — a provable lower bound |
| PSF sidelobe energy | `Σ A_k y_k + 2Σ B_kl y_k y_l`, `A_k=Σ_c a_ck²`, `B_kl=Σ_c a_ck a_cl` | Parseval — proved, see below |
| UV-point repulsion | `Σ E_kk y_k + Σ E_kl y_k y_l` | Cornwell 1988 |
| target UV density χ² | expand `Σ_b w_b(Σ_k a_bk y_k − t_b)²` | Boone 2002 |

Cost: `M + M(M−1)/2` variables → 210 at M=20, 1275 at M=50, **5050 at M=100**.

### 5.3 The Parseval proof

With `n_c = Σ_k a_ck y_k` the gridded sampling function, `PSF = F⁻¹[n]`:

```
Σ_{l,m} |PSF|²  =  (1/N²) Σ_c n_c²           (Parseval)
PSF_peak        =  (1/N²) Σ_c n_c
```

So with the peak normalised to 1, total beam energy is `N² Σn_c²/(Σn_c)²`.
Minimising `Σ_c n_c²` at fixed sample count **minimises sidelobe energy** —
which is exactly the guide's claim, now proved and numerically verified.

### 5.4 Measured comparison of formulations

Toy problem M=8, N=4, Earth-rotation sampling, all 70 feasible selections
enumerated. Spearman ρ against exact unique-cell coverage:

| UV cells/axis | pairwise (M vars) | Bonferroni (M+pairs) | pairwise finds optimum | Bonferroni finds optimum |
|---|---|---|---|---|
| 8 | 0.896 | 0.910 | no | **yes** |
| 16 | 0.964 | 0.994 | no | **yes** |
| 32 | 0.981 | 0.996 | no | **yes** |
| 64 | 0.991 | **1.000** | yes | **yes** |
| 128 | 0.993 | **1.000** | yes | **yes** |

Variable cost at 16 cells/axis: 8 (pairwise, approximate) vs **36**
(baseline-space, exact) vs 280 (cell-variable construction, exact).

Coupling matrix fill for the Bonferroni objective: 2–10 % — two baselines
interact only if their tracks actually overlap.

### 5.5 Complexity

With `Σx_i = N` and a pairwise reward, this is exactly **densest-k-subgraph**:
pick `k` vertices maximising induced edge weight. NP-hard, no known
constant-factor approximation. That is the structural argument for using a
solver at all, and it belongs in the paper.

---

## 6. Other findings

### 6.1 Layout geometry — the golden ratio is not special

Twelve layout families, N=20, identical extent / grid / observing geometry:

| layout | unique UV cells | Cornwell F ↓ | Gaussian χ² ↓ | PSF sidelobe ↓ |
|---|---|---|---|---|
| **Reuleaux (1 ring)** | **10226** | **30.5** | 1.05 | 0.132 |
| hierarchical 3³ | 9596 | 208.9 | 0.93 | 0.201 |
| random (area-uniform) | 9266 | 40.0 | 0.067 | 0.178 |
| Vogel golden angle | 9018 | 34.6 | 0.251 | 0.124 |
| Gaussian | 7844 | 46.3 | **0.096** | 0.151 |
| log spiral, 3 arms | 7736 | 70.1 | 0.214 | 0.098 |
| **golden spiral (1 arm)** | 6872 | 74.9 | **0.023** | **0.090** |
| Fibonacci radial | 1813 | 1322.6 | 14.4 | 0.386 |

1. **Reuleaux triangle covers 49 % more unique cells than the golden spiral**
   with less than half the Cornwell energy — as predicted by Keto (1997); the
   SMA is built this way.
2. **The objectives disagree.** Reuleaux wins cells/Cornwell/FWHM; the golden
   spiral wins density-match and ties on peak sidelobe. "Best layout" is
   undefined until the science case fixes the objective.
3. **Angular-constant sweep:** golden angle ranks **6th of 10**. √3−1, 1/e,
   5/13, √2−1, 1/π all score higher on cell count. Only low-order rationals
   (1/2, 1/3) are clearly bad, because they collapse position angles onto `q`
   directions. The rule is "not a low-order rational", not "golden".
4. **Fibonacci radial is disqualified, not merely worse**: normalised Fibonacci
   radii put ~93 of 100 antennas in the inner 5 % of the array, producing 532
   minimum-separation violations at N=100.

### 6.2 Phase coherence dominates at THz

Kolmogorov broken power law for RMS excess path:
`σ_path(b) = κ·σ₁·(b/1km)^α`, α = 5/6 (b ≲ 1 km), 1/3 (longer), saturating past
the outer scale ~6 km. Then `σ_φ = 2πσ_path/λ` and the decorrelation factor is
`γ = exp(−σ_φ²/2)`.

Surviving coherence at λ = 1 mm:

| κ (residual after phase referencing) | b=100 m | b=500 m | b=1 km | b=5 km |
|---|---|---|---|---|
| 1.0 (uncorrected) | 0.65 | 0.002 | ~0 | ~0 |
| 0.2 | 0.98 | 0.78 | 0.45 | 0.10 |
| 0.1 | 0.996 | 0.94 | 0.82 | 0.56 |

Uncorrected, **96 % of baselines lose more than half their coherence**. So `κ`
— how well phase referencing / water-vapour radiometry works at the site — is
the single most influential unknown in the whole problem. If κ ≈ 1 the
optimiser simply builds a compact array and there is no interesting
long-baseline trade-off left.

Because `γ` depends only on the pair, `H_phase` is exactly quadratic in `x`.

---

## 7. What is missing (honest list)

1. **No optimiser at all.** The code can *evaluate* a configuration but never
   *searches*. No greedy, no annealing, no solver. This is the single biggest
   gap: every "best layout" above is the best of a hand-made shortlist, not an
   optimised result.
2. **No science case**, therefore no principled choice among the objectives in
   §5.2. They disagree.
3. **No real data** — pad coordinates, site PWV statistics, and the phase
   constants (σ₁, κ, turbulence breakpoints) are all synthetic or typical
   values.
4. **Physics simplifications**: coplanar array (`w = 0`), no primary beam, no
   bandwidth/time smearing, no noise model, no deconvolution (CLEAN), single
   declination, single frequency.
5. **No classical baseline**, so no quantum claim can currently be justified.
6. **Not factorially controlled**: the layout families differ in radial law
   *and* angular law simultaneously, so differences cannot be attributed.

---

## 8. Methods not yet tried — the menu

### 8.1 Layout generators worth adding

| method | rationale |
|---|---|
| Minimum-redundancy / Golay arrays | maximise distinct baselines for small N; classical (Golay 1971) |
| Perfect difference sets / Costas arrays | non-redundant by construction; the largest *perfect* 2D non-redundant array has only 6 elements ("Manx", McKay 2022), so relevant at small N only |
| Low-discrepancy sequences (Halton, Sobol) | quasi-random but more uniform than random; cheap, and random already performs well here |
| Poisson-disk / blue-noise sampling | enforces `d_min` **by construction** — removes shadowing violations for free |
| Nested/"zoom" configurations | ALMA's actual approach; supports reconfiguration between compact and extended |
| Terrain-constrained variants | needed once real pads exist |

### 8.2 Classical optimisers (none implemented — highest priority)

| method | cost | note |
|---|---|---|
| Greedy removal from M pads down to N | `O((M−N)M³)` | Panduranga Rao et al. 2009 did exactly this problem |
| Greedy addition + local swap | cheap | strong practical baseline |
| Simulated annealing | cheap | the standard comparator for annealing claims |
| Tabu search | cheap | often beats SA on QUBO |
| Genetic / NSGA-II multi-objective | moderate | Cohanim et al. 2004 used this for imaging-vs-cable-length |
| **Exact MILP / branch-and-bound (Gurobi, CPLEX)** | expensive | **could give provably optimal answers at moderate M — the strongest possible baseline, and rarely done in this literature** |
| Semidefinite relaxation + rounding | moderate | standard for densest-k-subgraph; gives an optimality *bound* |
| **Differentiable relaxation** (`x ∈ [0,1]`, soft UV binning with a Gaussian kernel, autodiff in JAX/PyTorch, then round) | cheap on GPU | modern, scales to large M, and optimises antenna *positions* continuously rather than selecting from pads |
| Bayesian optimisation | expensive per eval | only if the objective becomes costly (e.g. full imaging simulation) |
| Reinforcement learning | expensive | probably premature here |

### 8.3 Quantum / quantum-inspired

| method | practicality now |
|---|---|
| `dwave-neal` simulated annealing on the QUBO | immediate — the honest first "solver" |
| D-Wave hybrid (`LeapHybridSampler`) | handles 10⁵–10⁶ variables; the realistic route for M≈100 (5050 vars) |
| D-Wave QPU direct | limited — dense coupling means minor-embedding blows up; ~180 fully connected logical variables max |
| QAOA (Qiskit) | toy sizes only (tens of qubits); the guide names it, so a small demonstration has narrative value |
| Quantum-inspired (Fujitsu Digital Annealer, Toshiba SBM, SimCIM) | strong classical competitors — worth including precisely because they are the fair comparison |
| Grover adaptive search | theoretical interest only |

### 8.4 Objective functions worth adding

| objective | why |
|---|---|
| **D-optimality / A-optimality** (maximise `log det` of the Fisher information of the measurement operator) | array design *is* optimal experimental design; this is the statistically principled objective and it is essentially absent from the radio-array literature — **possibly the most novel direction available** |
| Mutual coherence of the sensing matrix | compressed-sensing reconstruction guarantees; explains why random layouts do well |
| Imaging fidelity on test sky models (with CLEAN) | the only objective a referee will fully trust |
| Direct PSF peak-sidelobe minimisation | Kogan 2000 |
| Multi-objective Pareto front (coverage vs sidelobes vs cost vs coherence) | avoids having to pick one objective at all — arguably the correct framing |

---

## 9. Where a second opinion would help most

Ranked by how much the answer would change the work:

1. **Is the baseline-variable (`y_k = x_i x_j`) formulation the right call?**
   It is exact but costs `M(M−1)/2` extra variables. Alternatives: accept the
   pairwise approximation; use a smarter quadratization (Ishikawa, Boros–Gruber
   aggregative splitting); or reformulate so the quartic term never appears.
   Is there a formulation that is both exact and `O(M)` in variables?

2. **Is there a better objective than cell-counting?** D-optimality and
   mutual-coherence framings are argued for in §8.4 but not implemented. Which
   would a referee find most defensible?

3. **Should the problem be selection at all?** Continuous position optimisation
   with a differentiable relaxation may dominate discrete pad selection —
   except that the physical problem genuinely *is* discrete (pads exist).
   Is a hybrid (continuous relax → snap to nearest pad → local repair) better?

4. **Penalty weights.** Current rule: `λ_R = s·W`, `λ_sel = s·N·W` where `W` is
   the largest single-variable objective swing. Sufficient but conservative;
   too-large penalties flatten the spectrum. Is there a better principled
   choice for this problem structure?

5. **Does the quantum part earn its place?** Densest-k-subgraph is NP-hard, but
   at M ≈ 100 classical heuristics are excellent. What would make a quantum
   result genuinely interesting rather than decorative?

6. **What is being missed physically?** Candidates: `w`-term / non-coplanarity,
   mosaicking, frequency synthesis (multi-frequency UV filling), polarisation,
   the fact that reconfiguration between observations changes the whole
   problem.

---

## 10. Recommended order of work

1. **Implement a classical optimiser** (greedy removal + simulated annealing,
   plus exact MILP at small M). Nothing else can be evaluated without it.
2. **Fix a science case** and derive the objective from it.
3. **Factorial layout experiment** — radial law × angular law separately.
4. **Get the real numbers**: M, N, pad coordinates, site latitude, frequency,
   κ, σ₁.
5. **Scale the baseline-space QUBO** to M ≈ 50–100 and benchmark against the
   classical baselines.
6. **Then** quantum: `dwave-neal` → hybrid → (optionally) a small QAOA
   demonstration.
7. Replace the PWV placeholder or drop the term.

---

## 11. Key references

- Cornwell 1988, IEEE AP-36 1165 — UV-point repulsion objective
- Keto 1997, ApJ 475 843 — curves of constant width / Reuleaux
- Keto 2012, JAI 1 1250008 — hierarchical arrays
- Boone 2001 A&A 377 368; 2002 A&A 386 1160 — pad optimisation, target UV density
- Kogan 2000, IEEE AP-48 1075 — sidelobe minimisation
- Karastergiou, Neri & Gurwell 2006, ApJS 164 552 — pad *selection* by Cornwell energy
- Cohanim, Hewitt & de Weck 2004, ApJS 154 705 — multi-objective GA/SA
- Panduranga Rao et al. 2009, arXiv:0901.4901 — greedy minimum-variance pad selection
- Golay 1971, JOSA 61 272 — minimum-redundancy arrays
- McKay et al. 2022, Radio Science 57 — perfect non-redundant ("Manx") arrays
- Carilli & Holdaway 1999, Radio Science 34 817 — tropospheric phase, Kolmogorov
- Thompson, Moran & Swenson — the interferometry reference (uvw convention, decorrelation)
- Rosenberg 1975; Boros & Hammer 2002; Boros & Gruber 2014 (arXiv:1404.6538) — quadratization
- Verma & Lewis 2022, arXiv:2206.11040 — QUBO penalty weights
- Quiroga Rodríguez 2025, JINST 20 P05026 — the golden-spiral paper this project started from
