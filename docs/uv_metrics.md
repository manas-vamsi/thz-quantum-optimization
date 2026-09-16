# UV metrics: what each number means and what it does not

The single most common failure mode in this area is quoting one percentage as
"UV coverage". This document states what each implemented metric measures,
what it is sensitive to, and how it can mislead.

## The grid is a free parameter

Every cell-based metric depends on two arbitrary choices:

- **extent** $u_{\max}$ — how far out the grid reaches
- **resolution** $n_c$ — how many cells span it

Filled fraction is $\lvert O\rvert/\lvert G\rvert$. Doubling $n_c$ multiplies
$\lvert G\rvert$ by 4 while $\lvert O\rvert$ grows at most linearly in the
number of samples, so the *same array* can be reported at almost any
percentage by choosing the grid. A filled fraction without its grid is not a
result.

Consequences enforced in the code:

- `UVGrid` is a frozen dataclass whose parameters are written into every
  results file and stamped on every figure.
- Every comparison uses **one shared grid**, sized from the most extended model
  in that comparison.
- `coverage_summary` reports the filled fraction twice: over the full square
  grid, and over the annulus between the shortest and longest sampled UV
  radius. The square-grid number is always lower because the corners are
  unreachable by any circularly bounded array.

### On the 95 % figure

The source paper (Quiroga Rodríguez, *JINST* **20** (2025) P05026) reports
approximately 95 % UV coverage for one of its simulated configurations. That
is a result *of that paper's* configuration, grid and definition. It is not
adopted here as a target, a benchmark, or a validation criterion, and no
number produced by this repository should be compared against it without first
reproducing that paper's grid and definition — which this repository has not
done.

## Metric by metric

### Unique cells $\lvert O \rvert$

Count of distinct cells with at least one sample. Grid-dependent but not
extent-normalised, so it compares cleanly *within* a fixed grid. This is the
quantity the QUBO UV objective tries to approximate.

### Filled fraction

$\lvert O\rvert/\lvert G\rvert$. See above. Report the grid or do not report
the number.

### Redundancy fraction

$(S - \lvert O\rvert)/S$ where $S$ is the total number of samples. High
redundancy is **not** automatically bad: repeated baselines buy
signal-to-noise and enable redundancy calibration. It is bad only relative to
a science goal that values new spatial frequencies over sensitivity. The code
therefore measures it and does not sign it.

### Mean and max multiplicity

$S/\lvert O\rvert$ and $\max_c n_c$. Useful for spotting a layout whose
"coverage" is one cell hit a thousand times.

### Density uniformity

Normalised Shannon entropy of $\{n_c\}$ over occupied cells only:

$$H = -\frac{1}{\ln\lvert O\rvert}\sum_{c\in O} p_c\ln p_c,
\qquad p_c = n_c/S$$

1 means every occupied cell carries equally many samples. Deliberately
independent of the extent, so it complements filled fraction rather than
duplicating it. It says nothing about *which* cells are filled — a layout
covering only short baselines can score 1.

### Angular uniformity

Same entropy over position-angle bins, with angles folded onto $[0,\pi)$.
Folding is essential: $(u,v)$ and $(-u,-v)$ are the same physical baseline, so
an unfolded histogram is symmetric by construction and every array scores
near 1. Values near 1 mean no preferred baseline orientation — relevant to
beam ellipticity, which is reported alongside it.

### Radial profile

Counts per annulus and counts per unit UV *area*. The per-area version is the
one to read: a flat curve means uniform coverage per unit area, and the
characteristic falling curve of a centrally condensed array is visible
immediately.

### Baseline-length distribution

Linear histogram, log-spaced histogram, empirical CDF, plus:

- **log-spacing uniformity** — entropy over log-spaced length bins. 1 means
  baselines are spread evenly across decades, which is the operational reading
  of "logarithmic baseline distribution".
- **power-law exponent $\alpha$** and its $R^2$, from a least-squares fit of
  $\log_{10}N(d)$ against $\log_{10}d$ on log bins. A log bin at $d$ has width
  $\propto d$, so a sample uniform in $d$ gives $\alpha = 1$ and a sample
  uniform in $\log d$ gives $\alpha = 0$.
  [`tests/test_metrics.py::test_powerlaw_fit_recovers_a_known_exponent`]

The project note asserts that useful UV sampling needs a logarithmic or
power-law spread of baseline lengths. These metrics let that be *tested*
rather than asserted. What the runs here show is that the assertion is not
straightforwardly supported as stated: at $N = 20$ the `fibonacci_radial`
layout has the **highest** log-spacing uniformity (0.93 vs 0.86 for the golden
spiral) and the flattest per-log-bin distribution ($\alpha = 0.07$ vs 0.82),
yet it covers by far the **fewest** unique UV cells (1813 vs 6872) and has the
worst peak sidelobe level (0.39 vs 0.09). A good spread of baseline *lengths*
does not imply a good spread of baseline *vectors*; the angular distribution
and the actual pad separations matter at least as much.

### PSF metrics

See §10 of [mathematical_formulation.md](mathematical_formulation.md). Two
cautions:

- Peak sidelobe level and sidelobe RMS depend on the grid, the weighting and
  the absence of a taper, so they compare only within one figure produced by
  one script. The golden spiral scores PSL 0.090 in experiment 01 and 0.079 in
  experiment 04 — same layout, same weighting, different grid extent because
  the set of compared models differed. That spread is the size of the
  systematic.
- FWHM comes from an azimuthal average, which is a poor summary for an
  elliptical beam. The ellipticity is reported next to it precisely so that a
  reader can see when the FWHM should not be trusted
  (`fibonacci_radial` reaches 0.48 at $N=20$).

## What is not implemented

- No weighting schemes beyond natural and uniform (no Briggs/robust, no taper).
- No $w$-term, no non-coplanar correction, no primary beam.
- No noise model, so nothing here speaks to sensitivity.
- No deconvolution, so "imaging quality" is not measured — only the beam is.
