# Fibonacci is not the golden ratio

This document exists because the obvious move — "replace $\Phi$ with $F_n$ in
the spiral equation" — is not a mathematically defined operation, and doing it
anyway would produce a result that means nothing.

## The distinction

**Golden-ratio logarithmic spiral.** A continuous curve parameterised by a
single irrational constant:

$$r(\theta) = a\,\Phi^{\theta/2\pi}, \qquad \Phi = \tfrac{1+\sqrt5}{2}$$

$\Phi$ enters as a *growth rate*. It has units of nothing; it is the factor by
which the radius multiplies per turn.

**Fibonacci sequence.** A discrete integer sequence:

$$F_n = F_{n-1} + F_{n-2}, \qquad F_1 = F_2 = 1$$

$F_n$ is an *index-dependent integer*, not a rate. Writing
$r(\theta) = a F_n^{\theta/2\pi}$ requires choosing which $n$, and the answer
changes the curve completely: $F_3 = 2$ gives doubling per turn,
$F_{10} = 55$ gives a 55-fold expansion per turn. There is no principled
choice, so the substitution is meaningless.

The two are related only asymptotically:

$$\lim_{n\to\infty}\frac{F_{n+1}}{F_n} = \Phi,
\qquad F_n \sim \frac{\Phi^n}{\sqrt5}$$

Convergence is fast ($F_{11}/F_{10} = 1.61818$) but the *early* terms —
$1, 1, 2, 3$ — are exactly where the sequence is least golden, and they are
the ones that end up in the compact core of a small array.

## What this repository implements instead

Three constructions, each with a stated definition, selected by `variant=`:

### `radial` — Fibonacci numbers as radii

$$r_n = r_{\min} + (r_{\max}-r_{\min})\frac{F_n - F_{\min}}{F_{\max}-F_{\min}},
\qquad \theta_n = \theta_0 + n\,\Delta\theta$$

This is the honest reading of "a Fibonacci distribution of antennas": the
*sequence* sets the radial positions. The angular step is a separate, freely
configurable parameter (default: the golden angle), and that separation is the
point — the radial law and the angular law are independent design choices that
the golden-spiral formulation ties together.

**Consequence, measured in experiment 01.** Because $F_N$ dominates the
normalisation, $r_n \approx r_{\min}$ for all but the last few $n$. With
$r \in [20, 1000]$ m and `skip=2`:

| $N$ | antennas within the inner 5 % of the radial range | min. pad separation | pairs violating $d_{\min} = 9$ m |
|---|---|---|---|
| 20 | 13 | 9.13 m | 0 |
| 50 | 43 | 1.83 m | 91 |
| 100 | 93 | 1.02 m | 532 |

So beyond $N \approx 20$ the `radial` construction is not a physically valid
array at this extent at all — it packs almost every antenna into a crowded
core. This is a property of the construction, not a bug, and it is the main
reason the variant behaves the way it does in experiment 01. The
`rational_angle` variant is worse still: at $N = 100$ two pads land on top of
each other (separation $2.6\times10^{-11}$ m), because the rational angular
step repeats exactly while the radii have already saturated at $r_{\min}$.

### `golden_angle` — Vogel phyllotaxis (the "Fibonacci sunflower")

$$r_n \propto \sqrt{n}, \qquad \theta_n = \theta_0 + n\,\theta_g$$

This is what most sources actually mean by "Fibonacci pattern": the visible
spiral arm counts of a sunflower head are Fibonacci numbers, but the
*generating rule* uses the golden angle $\theta_g = 2\pi(1-\Phi^{-1})$, an
irrational constant, not the sequence. The radial law is area-uniform, not
logarithmic — a genuinely different family from the golden spiral, despite the
shared name.

Calling this "Fibonacci" and the log spiral "golden" and then comparing them
is comparing an angular rule plus a $\sqrt n$ radial law against an exponential
radial law. That is a legitimate experiment, but the labels hide what is
actually varying, which is why both are reported separately here.

### `rational_angle` — the sharpest statement of the difference

Radii as in `radial`, but the angular step is a **rational** Fibonacci
approximant:

$$\Delta\theta = 2\pi\frac{F_{k-1}}{F_k}$$

Then $F_k\,\Delta\theta \equiv 0 \pmod{2\pi}$: the pattern closes exactly
after $F_k$ antennas and the array acquires $F_k$-fold angular symmetry. The
golden angle, being irrational, never closes — no finite number of steps
returns to the starting angle.

**This is the operational difference between "Fibonacci" and "golden".**
Fibonacci ratios are rational approximants; they produce periodic, redundant
angular structure. $\Phi$ is irrational; it produces aperiodic structure. For
an interferometer, angular periodicity means repeated baseline orientations,
i.e. redundancy — which is measurable, and is measured in experiment 01.
[`tests/test_arrays.py::test_rational_angle_is_periodic_but_golden_angle_is_not`]

## What a fair comparison requires

Fixed across all models: $N$, $[r_{\min}, r_{\max}]$, $\lambda$, observing
geometry, UV grid and extent, and the metric definitions. The experiment
harness enforces this: one shared UV grid per experiment, sized from the most
extended model so that none is measured on a grid tailored to itself.

What is *not* controlled, and cannot be: these constructions differ in more
than one property at a time (radial law **and** angular law). A difference in
UV coverage between `golden` and `fibonacci_radial` cannot be attributed to
either property alone. Isolating them needs a factorial experiment — the same
radial law with two angular rules, and the same angular rule with two radial
laws — which is listed as open work in the README.

## What this repository does not claim

- Not that Fibonacci is better. The `radial` variant does markedly worse on
  every UV metric measured here, but from $N = 50$ it also violates the
  minimum-separation constraint, so it is not a valid array at all at that size
  — a disqualification, not a fair loss. Meanwhile the `golden_angle` variant,
  which is equally entitled to the name "Fibonacci", covers *more* unique UV
  cells than the golden spiral at every $N$ tested (9018 vs 6872 at $N=20$).
  "Fibonacci" is not one thing, and the choice of construction matters more
  than the Fibonacci-vs-golden framing.
- Not that the golden spiral is optimal. An area-uniform random layout with
  the same $N$ and extent covers *more* unique UV cells than the golden spiral
  in every experiment run here.
- Not that any of this settles the question. One radial-extent choice, one
  declination, one grid, one wavelength.
