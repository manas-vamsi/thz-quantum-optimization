# Quantum Optimization of Terahertz Interferometric Architectures

Research prototype. Classical foundation first: array layouts → baselines →
(u,v) → UV-cell occupancy → metrics → PSF → a small, **validated** QUBO.
No quantum solver is used yet, on purpose.

**Status:** preliminary. Everything below is a first milestone, not a result
ready for publication.

---

## Reading guide

Claims in this repository are tagged. Please keep the tags when quoting
anything from it.

| tag | meaning |
|---|---|
| **[SOURCE]** | stated in the cited paper or the project onboarding note |
| **[CHOICE]** | an implementation decision made here, with a stated reason |
| **[HYPOTHESIS]** | something to be tested, not yet tested |
| **[RESULT]** | produced by a script in `experiments/`, reproducible from a fixed seed and config |

---

## 1. Motivation

Terahertz interferometry constrains array design harder than centimetre-wave
interferometry does: atmospheric precipitable water vapour, short phase
coherence times, minimum antenna separation set by dish size, a discrete set of
buildable pads, and the UV sampling the science requires. **[SOURCE, project
note]**

Choosing which $N$ of $M$ candidate pads to occupy is a combinatorial
selection problem, which makes it a candidate for QUBO formulation and,
eventually, quantum annealing or QAOA. The project note is explicit that the
physics — UV calculations, forward simulation — must be done **offline** and
compiled into QUBO coefficients, never inside the optimisation loop.
**[SOURCE, project note]** This repository follows that split exactly:
`thz_opt.qubo.coefficients` does all the physics once; the optimiser sees only
numbers.

## 2. Problem statement

Binary variables $x_i \in \{0,1\}$, one per candidate pad. A baseline exists
iff both ends are selected. Minimise

$$H_{\rm total} = H_{\rm select} + H_{\rm shadow} + H_{\rm pwv}
  + H_{\rm uv} + H_{\rm phase}$$

Full definitions and expansions: [`docs/mathematical_formulation.md`](docs/mathematical_formulation.md).

## 3. The existing golden-spiral approach

Elio Quiroga Rodríguez, *"Golden Spiral Radio Interferometry: Dynamic Antenna
Arrays for Optimal UV Coverage and Sub-Arcsecond Resolution (A Proposed
Implementation)"*, **JINST 20 (2025) P05026**,
[doi:10.1088/1748-0221/20/05/P05026](https://doi.org/10.1088/1748-0221/20/05/P05026).

The paper proposes placing antennas on a golden-ratio logarithmic spiral and
moving/rotating them dynamically to improve UV coverage, and reports
approximately 95 % UV coverage for one simulated configuration. **[SOURCE]**

That 95 % is treated here strictly as a paper-specific reported number. It is
**not** used as a target, benchmark or validation criterion anywhere in this
repository, because a filled-cell fraction is meaningless without the grid and
extent that produced it — see
[`docs/uv_metrics.md`](docs/uv_metrics.md). **[CHOICE]**

The dynamic-reconfiguration aspect of the paper is **not** implemented here;
only static layouts are compared.

## 3b. What the paper's own appendix code does and does not do

Read together with its published appendix (analysis in
[`docs/design_space.md`](docs/design_space.md) §1): the code generates a
single-arm log spiral, forms all baselines with conjugates, and plots them. It
contains **no UV-coverage computation at all** - no grid, no cell size, no
count - so the 95 % / 60 % figures in the paper's Table 2 are not outputs of
the published code. Two further claims do not follow: a rigid rotation of an
array is an isometry and cannot change `b_max`, so it cannot improve
resolution (only radial motion can); and the +41 % SNR is the ordinary
`sqrt(k)` gain from integrating longer, which any array gets. **[CHOICE: this
repository therefore treats the paper's geometry as a hypothesis to test, and
its quantitative claims as unreproduced.]**

## 4. The Fibonacci alternative

Substituting a Fibonacci number for $\Phi$ in $r = a\Phi^{\theta/2\pi}$ is not
a defined operation — $\Phi$ is a growth *rate*, $F_n$ is an index-dependent
*integer*. Three explicit Fibonacci constructions are implemented instead
(`radial`, `golden_angle`, `rational_angle`), each mathematically stated.
**[CHOICE]** Reasoning and the sharp distinction (rational Fibonacci angular
steps close after $F_k$ antennas; the irrational golden angle never does):
[`docs/fibonacci_vs_golden.md`](docs/fibonacci_vs_golden.md).

## 5–8. Formulation, baselines, UV, metrics

- Baselines: $\mathbf b_{ij} = \mathbf r_j - \mathbf r_i$ for $i<j$, exactly
  $N(N-1)/2$ of them; $(-u,-v)$ conjugate points added explicitly.
- Snapshot UV: $u = \Delta x/\lambda$, $v = \Delta y/\lambda$ (zenith,
  coplanar).
- Earth rotation: Thompson–Moran–Swenson convention, in its own module, and
  numerically verified to reduce to the snapshot case at the zenith.
- UV grid: square, cell-**centred** on the origin so the binning is symmetric
  under $(u,v)\to(-u,-v)$ — otherwise the dirty beam picks up a ~30 %
  spurious imaginary part. **[CHOICE, found during development]**
- Metrics: unique cells, filled fraction (full grid *and* reachable annulus),
  redundancy, density/angular uniformity, radial and angular profiles,
  baseline-length histogram/CDF/power-law fit, PSF peak sidelobe level,
  sidelobe RMS, FWHM, ellipticity.

Details: [`docs/mathematical_formulation.md`](docs/mathematical_formulation.md),
[`docs/uv_metrics.md`](docs/uv_metrics.md).

## 9. QUBO formulation

$$H_{\rm select} = \lambda_{\rm sel}\Big(\textstyle\sum_i x_i - N\Big)^2,
\qquad
H_{\rm shadow} = \lambda_{\rm sh}\!\!\!\sum_{i<j,\,d_{ij}<d_{\min}}\!\!\! x_ix_j,
\qquad
H_{\rm uv} = -\sum_{i<j} w_{ij}x_ix_j$$

Convention: upper-triangular $Q$, each pair once,
$E(x)=\sum_i Q_{ii}x_i+\sum_{i<j}Q_{ij}x_ix_j+\text{offset}$. The symmetric
$x^{\mathsf T}Ax$ form is available and the factor of two between them is
**tested**, not assumed.

Penalty weights come from a stated rule, not from taste:
$\lambda = s\cdot U_{\max}$ where $U_{\max}$ is the largest UV reward any
feasible selection can collect. **[CHOICE]**

The hard part — whether global unique-cell coverage can be a quadratic form at
all — has its own document:
[`docs/qubo_mapping.md`](docs/qubo_mapping.md). Short version: it is exactly
quadratic **iff** no UV cell is touched by two candidate pairs, and in the toy
problem that regime is also the regime where the objective is constant and
useless. So the pairwise form is an approximation, and the repository measures
how good an approximation it is rather than asserting one.

## 9b. Beyond golden vs Fibonacci

Both are closed-form curves chosen for their mathematics rather than for an
imaging objective. [`docs/design_space.md`](docs/design_space.md) surveys what
the field uses instead - curves of constant width (Keto 1997), multi-arm log
spirals (Conway/ALMA), Gaussian-density layouts (Boone 2002), hierarchical
arrays (Keto 2012), minimum-redundancy arrays (Golay 1971) - and what it
optimises: Cornwell's UV-point repulsion energy, target UV-density matching,
PSF sidelobes, multi-objective imaging-vs-cable-length. All of these are
implemented in `arrays/geometries.py` and `metrics/density_matching.py` and
measured in experiment 06.

The key structural result, in [`docs/design_space.md`](docs/design_space.md)
§5: rewriting the objective in **baseline-activation variables**
`y_k = x_i x_j` makes Cornwell's energy, target-density matching, and a
second-order (Bonferroni) bound on unique-cell coverage **all exactly
quadratic**, with `y_k = x_i x_j` enforced exactly by a Rosenberg penalty. The
model becomes a true QUBO over `M + M(M-1)/2` variables with no approximation
of the objective.

## 10. Preliminary results

All numbers from `experiments/`, config `configs/default.yaml`,
$N=20$, $r\in[20,1000]$ m, $\nu = 300$ GHz, shared $256^2$ UV grid,
4 h hour-angle window, 41 time samples. **[RESULT]**

### Golden spiral vs Fibonacci (Earth-rotation sampling, N = 20)

| Metric | Golden | Fib. radial | Fib. golden-angle | Random |
|---|---|---|---|---|
| unique UV cells | 6872 | 1813 | 9018 | 9266 |
| redundancy fraction | 0.559 | 0.884 | 0.421 | 0.405 |
| UV density uniformity | 0.968 | 0.841 | — | — |
| shortest baseline (m) | 37.0 | 9.1 | 172.1 | 65.1 |
| baseline log-spacing uniformity | 0.855 | 0.929 | — | — |
| baseline power-law $\alpha$ ($R^2$) | 0.82 (0.61) | 0.07 (0.02) | — | — |
| PSF peak sidelobe level | 0.090 | 0.386 | 0.124 | 0.178 |
| PSF ellipticity | 0.065 | 0.479 | — | — |

Full table for every metric and every model: `data/results/exp01_metrics.csv`.

What this **does** support:

1. The three Fibonacci constructions behave completely differently from one
   another. `golden_angle` covers *more* unique cells than the golden spiral;
   `radial` covers far fewer. "Fibonacci vs golden" is not a well-posed
   comparison until the construction is named.
2. The `radial` Fibonacci construction concentrates antennas near $r_{\min}$
   so strongly that from $N=50$ it violates the minimum-separation constraint
   (91 violating pairs at $N=50$, 532 at $N=100$). It is not a valid array at
   that size, so its poor UV metrics there are a disqualification rather than
   a fair loss.
3. A plain area-uniform **random** layout covers more unique UV cells than the
   golden spiral at every $N$ tested. Whatever the golden spiral is good for
   in this metric set, it is not raw cell count.
4. A good spread of baseline *lengths* does not imply good UV coverage:
   `fibonacci_radial` has the flattest per-log-bin length distribution
   ($\alpha = 0.07$) and the worst coverage and sidelobes of any valid layout
   tested. This directly tests, and does not support, the onboarding note's
   assertion as stated. **[tests the HYPOTHESIS]**

What this does **not** support: any statement that one layout family is better.
Single declination, single extent, single wavelength, single grid, no noise,
no deconvolution, no science case. The constructions differ in radial *and*
angular law simultaneously, so no difference can be attributed to either.

### Layout shootout, twelve families, N = 20 **[RESULT, experiment 06]**

| layout | unique UV cells | Cornwell F (lower better) | Gaussian chi2 (lower better) | PSF sidelobe (lower better) |
|---|---|---|---|---|
| **Reuleaux (1 ring)** | **10226** | **30.5** | 1.05 | 0.132 |
| hierarchical 3^3 | 9596 | 208.9 | 0.93 | 0.201 |
| random (area-uniform) | 9266 | 40.0 | 0.067 | 0.178 |
| Vogel golden angle | 9018 | 34.6 | 0.251 | 0.124 |
| Gaussian (Boone target) | 7844 | 46.3 | **0.096** | 0.151 |
| log spiral, 3 arms | 7736 | 70.1 | 0.214 | 0.098 |
| **golden spiral (1 arm)** | 6872 | 74.9 | **0.023** | **0.090** |
| Fibonacci radial | 1813 | 1322.6 | 14.4 | 0.386 |

1. The Reuleaux triangle covers **49 % more unique UV cells** than the golden
   spiral with **less than half** its Cornwell energy - as Keto (1997)
   predicted, and the reason the SMA's pads are nested Reuleaux triangles.
2. **The objectives disagree.** Reuleaux wins cell count, Cornwell energy and
   FWHM; the golden spiral wins Gaussian-density match and ties on peak
   sidelobe. "Best layout" is undefined until the science case fixes the
   objective.
3. **The golden angle is not special.** Scanning the phyllotaxis angle at fixed
   N, extent and grid, it ranks **6th of 10**: sqrt(3)-1, 1/e, 5/13, sqrt(2)-1
   and 1/pi all score higher on cell count. Only low-order rationals (1/2, 1/3)
   are clearly bad, and they are bad because they collapse position angles onto
   q directions. What matters is "not a low-order rational", not "golden".

### Toy QUBO (M = 8 pads, N = 4, $\lambda$ = 1 mm) **[RESULT]**

- Direct objective vs matrix QUBO over all $2^8 = 256$ configurations:
  max difference **0.0** — exact agreement.
- Upper-triangular vs symmetric convention: max difference **0.0**.
- With $\lambda = 2U_{\max} = 24$: best feasible energy $-12.0$, best
  infeasible $+4.0$ — the penalty rule provably keeps the optimum feasible
  here.
- Pairwise surrogate vs exact coverage over the 70 feasible selections, under
  Earth-rotation sampling: Spearman $\rho$ rises 0.88 → 0.999 as the grid is
  refined; the surrogate's own optimum is still up to **16 % short** of true
  coverage at coarse grids, and 5–8 % short even at $\rho\approx0.96$.
- Exact auxiliary-variable construction for the same toy problem:
  **138 binary variables** (8 pads + 28 pairs + 46 cells + 56 slack) to select
  4 antennas from 8. Exactness is affordable only at toy scale.

### QUBO formulation comparison **[RESULT, experiment 07]**

Same toy problem under Earth-rotation sampling, all 70 feasible selections
enumerated:

| UV cells/axis | pairwise surrogate rho | Bonferroni rho | pairwise finds optimum | Bonferroni finds optimum |
|---|---|---|---|---|
| 8 | 0.896 | 0.910 | no | **yes** |
| 16 | 0.964 | 0.994 | no | **yes** |
| 32 | 0.981 | 0.996 | no | **yes** |
| 64 | 0.991 | **1.000** | yes | **yes** |
| 128 | 0.993 | **1.000** | yes | **yes** |

The baseline-space Bonferroni formulation finds the true optimum at **every**
resolution; the pairwise surrogate misses it at four of eight. Cost at 16
cells/axis: 8 variables (pairwise) vs **36** (Bonferroni, 10 % coupling fill)
vs 280 (exact cell variables). Bonferroni scales as `M + M(M-1)/2` - 1275
variables at M = 50, 5050 at M = 100.

With `sum_i x_i = N` and a pairwise reward this is exactly
**densest-k-subgraph**: NP-hard, no known constant-factor approximation. That
is the structural reason a QUBO solver is worth trying here at all.

## 11. Current limitations

- No quantum solver. Deliberate: the classical objective is not yet settled.
- PWV costs are **synthetic** (`constant`, `altitude`). No measured atmospheric
  data is shipped and none is invented. The term is an interface, not physics.
- The phase term is a **graph-connectivity placeholder**, not a coherence
  model. Its parameters (`coherence_length_m`, `core_radius_m`) are free
  numbers, not measurements.
- Arrays are strictly coplanar; no station heights, no $w$-term, no
  non-coplanar correction.
- PSF: no gridding kernel, no taper, no deconvolution, no noise. Absolute
  sidelobe levels are grid-dependent and comparable only within one figure
  (the same layout scores 0.090 and 0.079 in two scripts whose grids differ).
- No primary beam, no bandwidth/time smearing, no elevation cut applied
  (elevation is reported so the user can check).
- One observing geometry throughout. Declination dependence is untested.
- `float64` holds Fibonacci numbers exactly only to $F_{78}$; beyond that the
  sequence is approximate (irrelevant after normalisation, stated for honesty).

## 12. Experiments

| script | what it does | outputs |
|---|---|---|
| `01_golden_vs_fibonacci.py` | full pipeline, 6 layouts × N = 10/20/50/100, comparison table | figs 1, 2, 5, 9; `exp01_metrics.{csv,json}` |
| `02_snapshot_uv.py` | snapshot UV only, radial/angular densities | figs 3, 4, 6, 7; `exp02_snapshot_metrics.csv` |
| `03_earth_rotation.py` | UV tracks, coverage growth vs integration time | figs 11, 12; `exp03_*.csv` |
| `04_psf_comparison.py` | dirty beams and cuts | fig 8; `exp04_psf_metrics.csv` |
| `05_small_qubo.py` | toy QUBO: validate, enumerate, measure the surrogate | fig 10; `exp05_*.{json,csv}` |
| `06_layout_shootout.py` | twelve layout families on three objectives + phyllotaxis-angle sweep | figs 13, 14, 15; `exp06_*.csv` |
| `07_formulation_comparison.py` | pairwise vs Bonferroni vs exact QUBO formulations | fig 16; `exp07_*.{csv,json}` |

Every figure carries its generating parameters in the corner. Every run writes
coordinates, baselines, UV samples, occupancy grids and metrics to
`data/generated/*.npz` and `data/results/*.{csv,json}`.

## 13. Reproducibility

```bash
git clone <this repo>
cd thz-quantum-optimization
python -m pip install -e .          # or: pip install -r requirements.txt

PYTHONPATH=src python -m pytest -q  # 70 tests

cd experiments
python 01_golden_vs_fibonacci.py
python 02_snapshot_uv.py
python 03_earth_rotation.py
python 04_psf_comparison.py
python 05_small_qubo.py
```

Requires Python ≥ 3.10 and numpy / scipy / matplotlib / pandas / pyyaml.
Optional quantum stack (`dimod`, `dwave-neal`, `qiskit`) is **not** required
and is not used by any of the above. Random layouts are seeded from
`configs/default.yaml`; every other generator is deterministic.

## 14. What should be done next

Ordered by what unblocks the most, per
[`docs/design_space.md`](docs/design_space.md) §6:

1. **A classical optimiser** - greedy removal against the Bonferroni objective
   (Panduranga Rao et al. 2009 style), then simulated annealing. Without it
   there is no baseline, and no way to tell whether any configuration here is
   good.
2. **A science case, and a target UV density derived from it.** Every "best
   layout" above is best only under an arbitrary metric.
3. **Factorial layout experiment** - radial law x angular law separately, so
   the differences in experiment 06 can be attributed.
4. **Declination and hour-angle sweep** - every result here is at one geometry.
5. **Scale the baseline-space QUBO** to M = 50-100 pads and benchmark against
   greedy and SA *before* any quantum solver is involved.
6. **Then** quantum: `dimod`/`dwave-neal` exact on the 36-variable toy, hybrid
   at M = 100 (5050 variables is within reach of hybrid solvers, not of a bare
   Advantage QPU at this coupling density).
7. **Replace the phase and PWV placeholders** with something physically
   defensible before any THz-specific claim is made.

## Repository layout

```
src/thz_opt/         arrays/ interferometry/ metrics/ constraints/ qubo/
experiments/         01..05 + common.py
tests/               70 tests across arrays, baselines, uv, metrics, qubo
configs/             default.yaml, experiments.yaml
docs/                mathematical_formulation, fibonacci_vs_golden, uv_metrics,
                     qubo_mapping, design_space
data/                generated/ (npz)  results/ (csv, json)
figures/             fig01 .. fig16
```

## Licence

MIT — see [LICENSE](LICENSE).
