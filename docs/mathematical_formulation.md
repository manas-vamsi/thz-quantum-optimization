# Mathematical formulation

Every symbol used by the code, with the convention the code actually follows.
Where a convention could reasonably have gone the other way, the choice is
stated and a unit test pins it (test file named in brackets).

---

## 1. Fibonacci sequence

$$F_n = F_{n-1} + F_{n-2}, \qquad F_1 = F_2 = 1$$

so $1, 1, 2, 3, 5, 8, 13, 21, \dots$  Implemented in
`thz_opt.arrays.fibonacci.fibonacci_sequence` with configurable seeds
$F_1, F_2$. [`tests/test_arrays.py::test_fibonacci_sequence_recurrence`]

Binet's formula gives the asymptotics used later:

$$F_n = \frac{\Phi^n - (-\Phi)^{-n}}{\sqrt 5} \;\sim\; \frac{\Phi^n}{\sqrt 5}$$

## 2. Golden ratio

$$\Phi = \frac{1+\sqrt 5}{2} = 1.6180339887\ldots, \qquad \Phi^2 = \Phi + 1$$

$$\lim_{n\to\infty} \frac{F_{n+1}}{F_n} = \Phi$$

[`tests/test_arrays.py::test_phi_value`, `::test_fibonacci_ratio_converges_to_phi`]

The **golden angle** is the circle divided in the golden ratio:

$$\theta_g = 2\pi\left(1 - \Phi^{-1}\right) = 2.399963\ \mathrm{rad} = 137.5078^\circ$$

## 3. Golden-ratio logarithmic spiral (Model 1)

$$r(\theta) = a\,\Phi^{\theta/2\pi}$$

i.e. the radius multiplies by exactly $\Phi$ per full turn. $N$ antennas are
placed at

$$\theta_n = \theta_0 + \frac{2\pi T (n-1)}{N-1}, \qquad n = 1 \dots N$$

over $T$ turns, then

$$x_n = r(\theta_n)\cos(\theta_n + \psi), \qquad
  y_n = r(\theta_n)\sin(\theta_n + \psi)$$

with $\psi$ a rigid orientation.

**Implementation choice.** The radii are affinely rescaled onto
$[r_{\min}, r_{\max}]$:

$$r \mapsto r_{\min} + (r_{\max}-r_{\min})\,
  \frac{r - \min_k r_k}{\max_k r_k - \min_k r_k}$$

This is not in the source paper. It exists so that every layout in a
comparison occupies the same physical extent; otherwise the comparison would
confound shape with size. The growth-per-turn property of the *unnormalised*
spiral is tested separately.
[`tests/test_arrays.py::test_golden_spiral_growth_per_turn`]

## 4. Fibonacci radial construction (Model 2)

Three distinct constructions, because substituting $F_n$ for $\Phi$ in
$r = a\Phi^{\theta/2\pi}$ is not a defined operation (see
[fibonacci_vs_golden.md](fibonacci_vs_golden.md)).

**(a) `radial`** — Fibonacci numbers set the radii:

$$r_n = r_{\min} + (r_{\max}-r_{\min})\,\frac{F_n - F_{\min}}{F_{\max}-F_{\min}},
\qquad \theta_n = \theta_0 + n\,\Delta\theta$$

with $\Delta\theta$ configurable (default $\theta_g$). Since
$F_n \sim \Phi^n/\sqrt5$, the normalised radii are close to a geometric
progression — but with $\max_k F_k$ in the denominator, $r_n \to r_{\min}$ for
all but the last handful of $n$, which is a *strong* radial concentration and
shows up clearly in the experiments.

**(b) `golden_angle`** — Vogel phyllotaxis:

$$r_n = r_{\min} + (r_{\max}-r_{\min})\sqrt{\tfrac{n+\tfrac12}{N}} \ \text{(rescaled)},
\qquad \theta_n = \theta_0 + n\,\theta_g$$

Area-uniform, not logarithmic. The Fibonacci content is in the parastichy
counts, not in the coordinates.

**(c) `rational_angle`** — radii as in (a), but

$$\Delta\theta = 2\pi\,\frac{F_{k-1}}{F_k}$$

which is *rational*, so the pattern closes after $F_k$ steps:
$F_k \Delta\theta \equiv 0 \pmod{2\pi}$. The irrational golden angle never
closes. This is the cleanest mathematical statement of the
Fibonacci-vs-golden distinction.
[`tests/test_arrays.py::test_rational_angle_is_periodic_but_golden_angle_is_not`]

## 5. Antenna coordinates

$$\mathbf r_i = (x_i, y_i), \quad i = 1 \dots N$$

Local East–North–Up metres: $x$ East, $y$ North, $z = 0$ (strictly coplanar
array; station height is outside the current model).

## 6. Baselines

$$\mathbf b_{ij} = \mathbf r_j - \mathbf r_i, \qquad i < j$$

$$N_b = \binom{N}{2} = \frac{N(N-1)}{2}$$

Only $i<j$ is enumerated; $(j,i)$ carries $-\mathbf b_{ij}$ and is added as the
conjugate point where UV sampling is formed. Baseline length
$d_{ij} = \lVert \mathbf b_{ij}\rVert$.
[`tests/test_baselines.py`]

## 7. UV coordinates

**Snapshot (zenith, coplanar):**

$$u = \frac{\Delta x}{\lambda}, \qquad v = \frac{\Delta y}{\lambda}
\qquad [\text{wavelengths}]$$

with $\lambda = c/\nu$. Both $(u,v)$ and $(-u,-v)$ are sampled.

**Earth rotation** (Thompson, Moran & Swenson convention). Local ENU baseline
$(e, n, h)$ at latitude $\phi$ to equatorial components:

$$X = -n\sin\phi + h\cos\phi, \qquad Y = e, \qquad Z = n\cos\phi + h\sin\phi$$

and towards a source at hour angle $H$, declination $\delta$:

$$\begin{pmatrix} u \\ v \\ w \end{pmatrix} = \frac{1}{\lambda}
\begin{pmatrix}
\sin H & \cos H & 0 \\
-\sin\delta\cos H & \sin\delta\sin H & \cos\delta \\
\cos\delta\cos H & -\cos\delta\sin H & \sin\delta
\end{pmatrix}
\begin{pmatrix} X \\ Y \\ Z\end{pmatrix}$$

**Consistency.** For $h=0$, $\delta = \phi$, $H = 0$ this reduces exactly to
the snapshot formulae, which is asserted numerically rather than assumed.
[`tests/test_uv.py::test_earth_rotation_reduces_to_snapshot_at_zenith`]

## 8. UV-cell discretisation

Square grid, $n_c$ cells per axis over $[-u_{\max}, u_{\max}]$, cell size
$\Delta u = 2u_{\max}/n_c$. Cell $k$ along an axis is **centred** on

$$u_k = \left(k - \tfrac{n_c}{2}\right)\Delta u$$

so one cell is centred on the origin and the binning is exactly symmetric
under $(u,v)\mapsto(-u,-v)$. Binning on cell *edges* instead introduces a
half-cell offset which makes the gridded sampling function non-Hermitian and
gives the dirty beam a spurious imaginary part of order 30 % — this was
observed during development and is why the convention is what it is.
[`tests/test_uv.py::test_origin_sits_in_a_single_cell_and_binning_is_symmetric`]

Flat cell id: $c = i_u n_c + i_v$. Samples outside the extent are dropped and
counted, never clipped into edge cells.

## 9. UV-coverage metrics

Let $G$ be the set of considered cells, $n_c$ the number of samples in cell
$c$, $O = \{c : n_c > 0\}$, $S = \sum_c n_c$.

| metric | definition |
|---|---|
| unique cells | $\lvert O\rvert$ |
| filled fraction | $\lvert O\rvert / \lvert G\rvert$ — **grid-dependent**, quote $G$ |
| redundancy fraction | $(S - \lvert O\rvert)/S$ |
| mean multiplicity | $S/\lvert O\rvert$ |
| density uniformity | $-\sum_{c\in O} p_c\ln p_c / \ln\lvert O\rvert$, $p_c = n_c/S$ |
| angular uniformity | same entropy over position-angle bins, angles folded onto $[0,\pi)$ |
| radial density | counts per annulus divided by $\pi(r_{\rm out}^2 - r_{\rm in}^2)$ |
| log-spacing uniformity | entropy of baseline lengths over log-spaced bins |
| power-law exponent | least squares $\log_{10} N(d) = \alpha\log_{10} d + c$ |

No single one of these is "UV coverage". A filled fraction quoted without its
grid and extent is meaningless, since $\lvert G\rvert$ is a free parameter.

## 10. PSF relation

$$\mathrm{PSF}(l,m) = \mathcal F^{-1}\!\left[S(u,v)\right]$$

computed as `fftshift(ifft2(ifftshift(S)))`, normalised to unit peak.
$S$ is either the occupancy (natural weighting) or its indicator (uniform).
Flat sky, $w = 0$, no gridding kernel, no taper. Pixel scale
$1/(n_c\Delta u)$ radians; field of view $1/\Delta u$ radians.

Metrics: peak sidelobe level and sidelobe RMS outside the main lobe (main lobe
= first minimum of the azimuthally averaged profile), FWHM of that profile,
and the ellipticity of the half-power region. Absolute values depend on the
grid, so they compare only within one figure.

## 11. Binary selection variables

$$x_i \in \{0,1\}, \quad i = 1 \dots M \qquad
x_i = 1 \iff \text{candidate pad } i \text{ carries an antenna}$$

A baseline exists iff both ends are selected: $y_{ij} = x_i x_j$.

## 12. Selection penalty

$$H_{\rm select} = \lambda_{\rm sel}\Big(\sum_i x_i - N\Big)^2$$

Expanding with $x_i^2 = x_i$:

$$H_{\rm select} = \lambda_{\rm sel}\left[(1-2N)\sum_i x_i
  + 2\sum_{i<j} x_i x_j + N^2\right]$$

so it contributes $\lambda_{\rm sel}(1-2N)$ to every diagonal entry,
$2\lambda_{\rm sel}$ to every off-diagonal entry, and a constant
$\lambda_{\rm sel}N^2$.

## 13. Shadow penalty

$$H_{\rm shadow} = \lambda_{\rm sh}\!\!\sum_{i<j,\; d_{ij} < d_{\min}}\!\! x_i x_j,
\qquad d_{\min} = f\,D$$

with $D$ the dish diameter and $f$ the separation factor ($f = 1.5$ by
default, from the project note; not a universal constant).

## 14. UV objective

$$H_{\rm uv} = -\sum_{i<j} w_{ij}\,x_i x_j, \qquad w_{ij}\ge 0$$

The construction of $w_{ij}$, and the sense in which this can or cannot
represent unique UV-cell coverage, is the subject of
[qubo_mapping.md](qubo_mapping.md).

## 15. Complete QUBO

$$H_{\rm total} = H_{\rm select} + H_{\rm shadow} + H_{\rm pwv}
  + H_{\rm uv} + H_{\rm phase}$$

with the linear PWV term $H_{\rm pwv} = \lambda_{\rm pwv}\sum_i c_i x_i$
(synthetic cost model — see `thz_opt.constraints.pwv`) and the phase surrogate
$H_{\rm phase} = \lambda_{\rm ph}\sum_{i<j} P_{ij}x_i x_j$ (graph
connectivity, explicitly a placeholder).

## 16. Matrix form and the factor of two

The code stores an **upper-triangular** $Q$ and a scalar offset, with

$$E(x) = \sum_i Q_{ii}x_i + \sum_{i<j} Q_{ij}x_i x_j + \text{offset}
\tag{1}$$

Each unordered pair appears **once**, carrying the whole coefficient. This is
the convention of `dimod`/D-Wave QUBO dictionaries, so `to_qubo_dict` needs no
rescaling.

The symmetric convention

$$E(x) = x^{\mathsf T} A x + \text{offset} \tag{2}$$

splits each pair coefficient in half: $A_{ij} = A_{ji} = Q_{ij}/2$ for
$i \ne j$, $A_{ii} = Q_{ii}$. `to_symmetric` and `to_upper_triangular`
convert between them, and equality of (1) and (2) is asserted on random binary
vectors — this is the classic factor-of-two bug and it is tested, not trusted.
[`tests/test_qubo.py::test_conventions_agree`]

Note for (1): since $x_i^2 = x_i$, $x^{\mathsf T}\,\mathrm{triu}(Q)\,x$ already
evaluates (1) exactly for binary $x$; `qubo_energy` rejects non-binary input
rather than silently returning something else.

## 17. Penalty weights

Not chosen by hand. Let

$$U_{\max} = \sum \text{(the } \tbinom{N}{2} \text{ largest } w_{ij})$$

be the most any selection of $N$ pads can gain from $H_{\rm uv}$. A penalty
cannot be bought off by UV gain when one unit of violation costs more than
$U_{\max}$, so

$$\lambda_{\rm sel} = \lambda_{\rm sh} = s\,U_{\max}, \qquad s > 1$$

with $s = 2$ by default. This is *sufficient*, deliberately conservative, and
the inputs are returned with the numbers so it can be tightened knowingly.
That the rule actually keeps the global optimum feasible is verified by
exhaustive enumeration.
[`tests/test_qubo.py::test_penalty_rule_is_sufficient`]
