# The QUBO objective function for THz interferometric pad selection

**Deliverable for:** the objective-function formulation requested in the project
onboarding guide, §3–4.

**Scope.** This document derives every term of `H_total` from the underlying
physics, expands each into explicit QUBO matrix entries, states exactly which
terms are quadratic and which are not, and reports the numerical validation of
all of it. Nothing here is asserted without either a derivation or a measured
number; open questions are listed in §10 rather than papered over.

**Status of each term**

| term | physical basis | exactly quadratic? | status |
|---|---|---|---|
| `H_select` | cardinality constraint | yes, in `x` | complete |
| `H_shadow` | dish obscuration, `d_min ≥ 1.5 D` | yes, in `x` | complete |
| `H_phase` | Kolmogorov phase structure function + decorrelation | **yes, in `x`** | derived here, replaces the earlier graph placeholder |
| `H_pwv` | site opacity, `τ₂₂₅ ≈ 0.049·PWV + 0.018` | yes, linear in `x` | needs per-pad PWV data |
| `H_uv` | UV-cell occupancy from precomputed tracks | **no in `x`; yes in baseline variables `y`** (exactly for sidelobe energy / Cornwell / density; as a bound for union coverage) | complete, with the resolution in §6 |

---

## 1. Problem statement and variables

`M` candidate pads with known positions `r_i`, of which exactly `N` are to carry
antennas.

$$x_i \in \{0,1\}, \quad i = 1,\dots,M, \qquad x_i = 1 \iff \text{pad } i \text{ is occupied}$$

A baseline between pads `i` and `j` exists **iff both are occupied**, so it is
represented by the product `x_i x_j`. This is why interferometry is a natural
QUBO problem, exactly as the onboarding guide states: the pairwise structure is
physical, not imposed.

It is convenient to index the `P = M(M-1)/2` candidate pairs by a single index
`k ↔ (i,j)` with `i < j`, and to define

$$y_k \;=\; x_i x_j$$

`y_k` is **not** a new degree of freedom — it is determined by `x` — but making
it an explicit variable is what allows §6 to be exact.

## 2. QUBO convention

Following the guide's §2, the objective is

$$H(\mathbf v) \;=\; \sum_a Q_{aa} v_a \;+\; \sum_{a<b} Q_{ab} v_a v_b \;+\; \text{const}
\tag{2.1}$$

with `Q` **upper triangular**: each unordered pair appears exactly once and
carries its whole coefficient. This is also the convention of
`dimod.BinaryQuadraticModel.from_qubo`, so the matrix can be handed to a D-Wave
sampler with no rescaling.

The alternative symmetric convention `H = v^T A v` splits each pair coefficient
in half, `A_ab = A_ba = Q_ab/2`. Mixing the two is the classic factor-of-two
bug; the code provides conversions both ways and asserts their equality on
random binary vectors (`tests/test_qubo.py::test_conventions_agree`).

Throughout, binary variables satisfy `v_a² = v_a`, which is what allows every
square to collapse into linear + quadratic terms.

## 3. The offline precomputation pipeline

Per the guide's §4, the solver never runs a simulation. Everything below is
computed once, before optimisation:

1. **Baselines.** For every candidate pair, `b_ij = r_j − r_i` in local
   East–North–Up metres. There are exactly `M(M−1)/2` of them.

2. **UV tracks.** For an observation over hour angles
   `H ∈ [H_start, H_end]` towards declination `δ` from latitude `φ`, using the
   Thompson–Moran–Swenson convention:

   $$\begin{pmatrix} u \\ v \\ w \end{pmatrix} = \frac{1}{\lambda}
   \begin{pmatrix}
   \sin H & \cos H & 0 \\
   -\sin\delta\cos H & \sin\delta\sin H & \cos\delta \\
   \cos\delta\cos H & -\cos\delta\sin H & \sin\delta
   \end{pmatrix}
   \begin{pmatrix} X \\ Y \\ Z\end{pmatrix}$$

   with `X = −n sinφ + h cosφ`, `Y = e`, `Z = n cosφ + h sinφ`. Each baseline
   traces an elliptical arc. Both `(u,v)` and `(−u,−v)` are sampled.

   *Verified:* at `δ = φ`, `H = 0`, `h = 0` this reduces to the snapshot case
   `u = Δx/λ`, `v = Δy/λ` to within `2×10⁻¹⁰` wavelengths
   (`tests/test_uv.py::test_earth_rotation_reduces_to_snapshot_at_zenith`).

3. **Cell binning.** The UV plane is discretised into a square grid of
   `n_c` cells per axis over `[−u_max, u_max]`. Cells are **centred** on the
   origin, i.e. cell `k` spans a box centred at `(k − n_c/2)·Δu`. This matters:
   binning on cell *edges* breaks the `(u,v) → (−u,−v)` symmetry of the
   occupancy grid by half a cell and gives the dirty beam a spurious imaginary
   part of order 30 %. With centred cells the residual is `5×10⁻¹⁷`.

4. **Per-baseline occupancy.** For pair `k`, record

   $$a_{ck} \;=\; \text{number of UV samples that baseline } k \text{ places in cell } c$$

   The *set* of touched cells is `C_k = {c : a_{ck} > 0}`. Both the set and the
   multiplicities are needed: the coverage objective uses `C_k`, the
   sidelobe-energy objective uses `a_{ck}` (an Earth-rotation track can revisit
   a cell).

These `a_{ck}` are the entire physical input to the optimisation.

## 4. Terms that are quadratic in `x` directly

### 4.1 Selection

$$H_{\rm select} = \lambda_{\rm sel}\Big(\sum_i x_i - N\Big)^2$$

Expanding and using `x_i² = x_i`:

$$H_{\rm select} = \lambda_{\rm sel}\Big[(1-2N)\sum_i x_i \;+\; 2\sum_{i<j} x_i x_j \;+\; N^2\Big]
\tag{4.1}$$

so it contributes `λ_sel(1−2N)` to every diagonal entry, `2λ_sel` to every
off-diagonal entry among the `x` block, and a constant `λ_sel N²`.

Note a practical consequence flagged in the QUBO literature: this term makes
the `x` block **fully connected** regardless of how sparse the physics is —
`M(M−1)/2` couplers, which dominates the embedding cost on hardware with
bounded qubit degree.

### 4.2 Shadowing

Tightly packed dishes obscure each other at low elevation. The guide gives
`d_min ≈ 1.5 D` with `D` the dish diameter:

$$H_{\rm shadow} = \lambda_{\rm sh}\!\!\!\sum_{i<j\,:\,d_{ij} < d_{\min}}\!\!\! x_i x_j,
\qquad d_{\min} = f\,D,\ f = 1.5
\tag{4.2}$$

`f` is a parameter, not a constant — 1.5 is what the guide states and the true
value depends on dish design and the minimum working elevation.

### 4.3 Phase coherence — a physically grounded form

The earlier version of this project used a graph-connectivity placeholder
("is there a stepping-stone pad within a coherence distance?"). That can be
replaced by standard millimetre-interferometry physics, and the result turns
out to be **exactly quadratic in `x`**, which the graph version was not.

Tropospheric water-vapour fluctuations follow Kolmogorov turbulence, so the RMS
excess path length on a baseline of length `b` follows a broken power law:

$$\sigma_{\rm path}(b) \;=\; \sigma_1 \left(\frac{b}{1\,\text{km}}\right)^{\alpha},
\qquad
\alpha = \begin{cases}
5/6 & b \lesssim L_{\rm 3D} \quad (\text{3D turbulence, } L_{\rm 3D} \sim 0.5\text{–}2\ \text{km})\\
1/3 & L_{\rm 3D} \lesssim b \lesssim L_{\rm out}\quad (\text{2D turbulence})\\
0 & b \gtrsim L_{\rm out} \quad (\text{saturated, } L_{\rm out} \sim 5\text{–}10\ \text{km})
\end{cases}
\tag{4.3}$$

with `σ₁` the RMS path at a 1 km baseline (of order 1 mm at a good site, varying
by a factor of a few with weather). The corresponding RMS phase is

$$\sigma_\phi(b) = \frac{2\pi}{\lambda}\,\sigma_{\rm path}(b)$$

and for Gaussian phase errors the measured visibility amplitude is reduced by
the standard **decorrelation factor**

$$\gamma_{ij} \;=\; \exp\!\left(-\tfrac{1}{2}\,\sigma_\phi^2(d_{ij})\right)
\tag{4.4}$$

`γ_ij ∈ (0,1]` is the fraction of a baseline's coherence that survives. It
depends only on the pair, so

$$\boxed{\,H_{\rm phase} \;=\; \lambda_{\rm ph}\sum_{i<j}\bigl(1 - \gamma_{ij}\bigr)\,x_i x_j\,}
\tag{4.5}$$

is exactly quadratic in `x`, with no auxiliary variables.

This is a strict improvement on the stepping-stone heuristic in two ways: it is
derived from measured atmospheric behaviour rather than invented, and it is
continuous rather than a binary graph condition. A residual-correction factor
`κ ∈ (0,1]` multiplying `σ_path` represents phase referencing / water-vapour
radiometry, which at ALMA leaves >200 µm residual path on 10 km baselines at
high frequency — i.e. `κ` is not small, and that is a number to pin down with
Vishwas rather than guess.

**The same `γ_ij` should also weight the UV reward**: an incoherent baseline
contributes sampling but not signal. Multiplying the per-baseline UV
coefficients of §5–6 by `γ_ij` costs nothing structurally.

### 4.4 Precipitable water vapour

At THz the atmosphere is opaque wherever water vapour is not scarce. The
standard site relation between 225 GHz zenith opacity and PWV (Chajnantor
calibration) is

$$\tau_{225} \;\approx\; 0.049\,\mathrm{PWV[mm]} + 0.018$$

with transmission `exp(−τ/\sin(\text{elevation}))`. A per-pad cost `c_i`
derived from the pad's own PWV statistics enters linearly:

$$H_{\rm pwv} = \lambda_{\rm pwv}\sum_i c_i x_i
\tag{4.6}$$

**Honest caveat.** For a compact array on one plateau, PWV differences between
pads are negligible and this term is a constant offset under the fixed-`N`
constraint — it changes nothing. It becomes meaningful only if the candidate
pads span significant altitude or terrain differences (e.g. a valley vs a
ridge). *No measured PWV data is shipped in this repository, and none is
invented.* The implemented cost models are labelled synthetic.

## 5. The UV term: what it is supposed to do

The guide's §4 specifies the intent precisely:

> "We score the baseline based on which unique cells it samples, and apply
> redundancy penalties if it overlaps with another pair's track… minimising the
> QUBO energy landscape mathematically translates to minimising the
> auto-correlation of (u,v) density, which directly minimises the synthesised
> beam's PSF sidelobes."

**This claim is correct, and it can be proved.** Let the gridded sampling
function of a configuration be

$$n_c \;=\; \sum_k a_{ck}\, y_k$$

The dirty beam is `PSF = \mathcal F^{-1}[n]`. By Parseval's theorem, with the
`1/N²` convention of `numpy.fft.ifft2` on an `N × N` grid,

$$\sum_{l,m} \bigl|\mathrm{PSF}(l,m)\bigr|^2 \;=\; \frac{1}{N^2}\sum_c n_c^2 ,
\qquad
\mathrm{PSF}_{\rm peak} = \frac{1}{N^2}\sum_c n_c$$

so, normalising the beam peak to 1, the **total energy in the dirty beam** is

$$E_{\rm beam} \;=\; N^2\,\frac{\sum_c n_c^2}{\bigl(\sum_c n_c\bigr)^2}
\tag{5.1}$$

Since the main lobe is pinned at unit height, minimising `Σ_c n_c²` at a fixed
number of samples minimises the energy **outside** the main lobe — the
sidelobes. The guide's "auto-correlation of UV density" is exactly `Σ_c n_c²`.

*Verified numerically:* the Parseval ratio in (5.1) is `1.000000000000` for
four different array layouts, and the predicted `E_beam` matches the measured
beam energy to 12 significant figures.

Expanding the square, and using `y_k² = y_k`:

$$\boxed{\;\sum_c n_c^2 \;=\; \sum_k \underbrace{\Bigl(\sum_c a_{ck}^2\Bigr)}_{A_k} y_k
\;+\; 2\sum_{k<l} \underbrace{\Bigl(\sum_c a_{ck}a_{cl}\Bigr)}_{B_{kl}} y_k y_l \;}
\tag{5.2}$$

This is **exact** — verified to `0.000e+00` absolute error against direct
gridding across all 70 feasible selections of the toy problem.

`A_k` is the baseline's self-energy, `B_kl` is the guide's "redundancy penalty
for overlapping tracks". Both are pure offline numbers.

An alternative reading of the guide's step 3 — *maximise unique cells covered*
— gives a different but equally quadratic form. By the second **Bonferroni
inequality**, for any sets,

$$\Bigl|\bigcup_k C_k\Bigr| \;\ge\; \sum_k |C_k|\,y_k \;-\; \sum_{k<l} |C_k \cap C_l|\, y_k y_l
\tag{5.3}$$

a **provable lower bound** on unique-cell coverage, exact whenever no cell is
covered three or more times.

The two objectives are related but not identical: (5.2) wants the occupancy
spread evenly, (5.3) wants as many distinct cells as possible. They agree that
redundancy is bad and disagree about everything else, so one must be chosen
deliberately, or they must be combined with a stated weight:

$$H_{\rm uv} = -\,\alpha\Bigl[\sum_k |C_k| y_k - \sum_{k<l}|C_k\cap C_l| y_k y_l\Bigr]
\;+\; \beta\Bigl[\sum_k A_k y_k + 2\sum_{k<l} B_{kl}\, y_k y_l\Bigr]
\tag{5.4}$$

**Which of `α`, `β` matters is a science question, not a mathematical one.**
It is the first item in §10.

## 6. The one place the guide's recipe needs amending

The guide says these rewards and penalties are "loaded directly into the
off-diagonal coefficient `Q_ij`". The first part can be: `Σ_k A_k y_k` and
`Σ_k |C_k| y_k` are sums over single baselines, and

$$\sum_k w_k y_k \;=\; \sum_{i<j} w_{ij}\, x_i x_j$$

lands straight in `Q_ij`. **The overlap part cannot.** A term `y_k y_l`
involves pads `(i,j)` *and* `(k,l)` — **four** indices:

$$y_k y_l \;=\; x_i x_j x_k x_l$$

`Q_ij` has room for two. Written in the pad variables `x`, the redundancy term
is a **degree-4 pseudo-Boolean function**, not a QUBO term. This is not a flaw
in the physics — it is a statement about what a quadratic form can hold.

There are two consistent ways forward.

### Option A — stay in `x`, approximate

Fold the overlap into pairwise weights, e.g. by discounting each cell by how
many *candidate* pairs could supply it. Cost: `M` variables. Accuracy: measured
below — it picked the wrong optimum in 4 of 8 tested UV resolutions.

### Option B — promote `y` to real variables, stay exact

Introduce the `P = M(M−1)/2` variables `y_k` and enforce `y_k = x_i x_j` with
the standard **Rosenberg quadratization** (Rosenberg 1975; see Boros & Hammer
2002, Boros & Gruber 2014 for the general theory):

$$R_k \;=\; \lambda_R\bigl(x_i x_j - 2x_i y_k - 2x_j y_k + 3y_k\bigr)
\tag{6.1}$$

Check the four cases: `R_k = 0` whenever `y_k = x_i x_j`, and `R_k ≥ λ_R > 0`
otherwise. Every term in (6.1) is quadratic, so the whole model

$$H_{\rm total} = H_{\rm select} + H_{\rm shadow} + H_{\rm phase} + H_{\rm pwv}
  + H_{\rm uv}(y) + \sum_k R_k$$

is a **genuine QUBO over `M + M(M−1)/2` variables**.

**What "exactly" covers, precisely.** The quadratic form in `y` *equals* the
objective for PSF sidelobe energy (5.2), Cornwell repulsion, and target-density
matching. For unique-cell coverage it equals the **Bonferroni lower bound**
(5.3), not the coverage itself: the two coincide only when no UV cell is
touched by three or more active baselines. Exact union coverage requires cell
variables or a higher-order inclusion–exclusion construction. The size of that
gap is measured, not assumed — see the `bound gap` column in §6.

### Measured comparison

Toy problem `M = 8`, `N = 4`, Earth-rotation sampling, all 70 feasible
selections enumerated exhaustively; Spearman rank correlation against exact
unique-cell coverage:

| UV cells/axis | A: pairwise ρ | B: Bonferroni ρ | A finds optimum | B finds optimum | bound gap | B coupling fill |
|---|---|---|---|---|---|---|
| 8 | 0.896 | 0.910 | no | **yes** | 24.9 % | 29 % |
| 12 | 0.923 | 0.979 | yes | **yes** | 5.6 % | 18 % |
| 16 | 0.964 | 0.994 | no | **yes** | 1.9 % | 10 % |
| 24 | 0.967 | 0.997 | no | **yes** | 0.5 % | 8 % |
| 32 | 0.981 | 0.996 | no | **yes** | 0.3 % | 7 % |
| 48 | 0.993 | 1.000 | yes | **yes** | 0.03 % | 4 % |
| 64 | 0.991 | **1.000** | yes | **yes** | 0.05 % | 4 % |
| 128 | 0.993 | **1.000** | yes | **yes** | 0.01 % | 2 % |

Option B finds the true optimum at **every** resolution; Option A misses it at
four of eight. Variable cost at 16 cells/axis:

| formulation | variables |
|---|---|
| A — pairwise surrogate, pad space | 8 |
| **B — exact, baseline space** | **36** (37 couplings, 10 % fill) |
| C — exact via one variable per UV cell + slack | 280 |

(Option C is the textbook OR-linearisation of "cell is covered". It is exact
too, but needs a variable per cell plus slack bits, and is included only as
ground truth.)

## 7. Assembling `Q`

Variable ordering: `x_0 … x_{M−1}`, then `y_0 … y_{P−1}`.

| block | source | entries |
|---|---|---|
| `Q_ii`, `i < M` | (4.1) selection, (4.6) PWV | `λ_sel(1−2N) + λ_pwv c_i` |
| `Q_ij`, `i<j<M` | (4.1) selection, (4.2) shadow, (4.5) phase, (6.1) Rosenberg | `2λ_sel + λ_sh S_ij + λ_ph(1−γ_ij) + λ_R` |
| `Q_{i,M+k}` | (6.1) Rosenberg | `−2λ_R` for `i ∈ {i_k, j_k}` |
| `Q_{M+k,M+k}` | (6.1) + linear UV | `3λ_R + (\text{linear coefficient of } y_k)` |
| `Q_{M+k,M+l}` | overlap / repulsion | `(\text{quadratic coefficient of } y_k y_l)` |
| constant | (4.1) | `λ_sel N²` |

Scaling of the variable count:

| `M` pads | variables `M + M(M−1)/2` |
|---|---|
| 8 | 36 |
| 20 | 210 |
| 50 | 1 275 |
| 100 | **5 050** |
| 200 | 20 100 |

The objective's coupling matrix is sparse (2–10 % fill at useful resolutions,
since two baselines interact only if their tracks actually overlap). The
*constraint* is dense, however — §4.1. At `M = 100` the model is within reach
of hybrid quantum-classical solvers; it is **not** within reach of a bare
5000-qubit annealer at this connectivity, because minor-embedding a
densely-coupled logical graph costs many physical qubits per logical one.

## 8. Penalty weights

Penalty weights are derived, not chosen by taste. The standard rule is the
**upper-bound strategy**: make a constraint violation cost more than the most
the objective could possibly gain by committing it.

Let

$$W \;=\; \max_k \Bigl( |{\rm lin}_k| + \sum_l |{\rm quad}_{kl}| \Bigr)$$

be the largest change the objective can see from a single baseline variable
flipping. Then:

$$\lambda_R \;=\; s\,W, \qquad \lambda_{\rm sel} \;=\; s\,N\,W, \qquad s > 1
\tag{8.1}$$

The factor `N` in `λ_sel` is not cosmetic: adding one antenna beyond `N`
switches on up to `N` new baselines, so `N` small rewards can otherwise cancel
one over-count penalty **exactly**. That produced a genuine energy tie with an
infeasible configuration during development, which is why the rule is stated
this way.

**Caveat from the QUBO literature:** these are *sufficient* bounds and
deliberately conservative. Penalties that are too large flatten the low-energy
spectrum — the solver spends all its resolution on feasibility and returns
feasible-but-mediocre states. In practice one should start at (8.1) and reduce
`s` while the optimum stays feasible. The rule and its inputs are returned
alongside the numbers so this can be done knowingly.

## 9. Validation performed

All checks are automated (97 tests) and reproduce from a fixed seed.

| check | result |
|---|---|
| Matrix QUBO vs direct term-by-term objective, all `2^8 = 256` configurations | max difference **0.0** |
| Upper-triangular vs symmetric convention | max difference **0.0** |
| Sidelobe-energy quadratic form (5.2) vs direct gridding, all 70 selections | max difference **0.0** |
| Parseval identity (5.1), four layouts | ratio **1.000000000000** |
| Bonferroni bound (5.3) ≤ exact coverage | holds for every selection |
| Bonferroni exact when no triple coverage | holds |
| Rosenberg penalty: no auxiliary assignment beats the intended one | verified over random assignments |
| Penalty rule (8.1) keeps the global optimum feasible | verified by enumeration |
| Earth rotation reduces to snapshot at zenith | `2×10⁻¹⁰` wavelengths |
| Dirty beam is real (Hermitian sampling) | imaginary residual `5×10⁻¹⁷` |
| Baseline count `N(N−1)/2`, conjugate symmetry, no NaNs | verified |

## 10. Open questions — needed before the next step

These cannot be resolved from this side and are the natural agenda for the
meeting:

1. **Which UV objective?** Equation (5.4) has two terms. Maximum unique-cell
   coverage and minimum sidelobe energy are *different* goals and pick
   different arrays. The choice follows from the science case (an EHT-style
   shadow measurement wants different sampling from an extended-source survey).
2. **`M` and `N`?** The formulation is size-independent but the solver strategy
   is not: `M ≤ 20` is exhaustively verifiable, `M ≈ 100` needs hybrid methods.
3. **Real candidate pad coordinates**, or continue with synthetic geometry?
4. **Site parameters:** latitude, observing frequency, and whether pads differ
   enough in PWV for (4.6) to be non-trivial.
5. **Phase model constants:** `σ₁`, the turbulence breakpoints `L_3D`,
   `L_out`, and the residual factor `κ` after phase referencing. Equation (4.5)
   is only as good as these.
6. **Observation window:** hour-angle range and declination fix every `a_ck`.
7. **Access to the UV-track repository** for cross-validation of §3 step 2.

## 11. What is deliberately not claimed

- No quantum **hardware** has been run. The QUBO is solved with the classical
  simulated-annealing sampler in `dwave-samplers`, behind the same `dimod`
  interface a quantum annealer would use. The comparison it produced is
  unfavourable to the QUBO route: on real ALMA geometry the formulation gap is
  zero, but the sampler reaches 86-92 % of the proven optimum where a direct
  classical heuristic reaches 98 % roughly a hundred times faster. Any quantum
  claim must beat that, not a straw man.
- The PWV term is plumbing, not physics, until per-pad data exists.
- The phase term (4.5) now uses measured ALMA constants (ALMA Memo 624), not
  assumed ones, and reproduces that memo's own published scaling factors. What
  remains unpinned is time: the values are medians over thousands of
  observations, so they describe the site rather than a given night.
- Arrays are treated as coplanar (`w = 0`), with no primary beam, no bandwidth
  or time smearing, and no noise model.

## References

- Boros, E. & Hammer, P. L. 2002, *Pseudo-Boolean optimization*, Discrete Appl. Math. 123, 155
- Boros, E. & Gruber, A. 2014, *On quadratization of pseudo-Boolean functions*, arXiv:1404.6538
- Rosenberg, I. G. 1975 — reduction of nonlinear 0-1 optimisation to quadratic form
- Thompson, A. R., Moran, J. M. & Swenson, G. W., *Interferometry and Synthesis in Radio Astronomy* — `(u,v,w)` convention, decorrelation factor
- Cornwell, T. J. 1988, IEEE Trans. Ant. Prop. 36, 1165 — UV-point repulsion objective
- Boone, F. 2002, A&A 386, 1160 — target UV density and beam shape
- Keto, E. 1997, ApJ 475, 843 — uniform UV sampling and array shape
- Carilli, C. L. & Holdaway, M. A. 1999, Radio Science 34, 817 — tropospheric phase, Kolmogorov structure function
- Verma, A. & Lewis, M. 2022, arXiv:2206.11040 — penalty weights in QUBO formulations
