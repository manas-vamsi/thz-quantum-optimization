# Mapping UV coverage onto a QUBO

This is the central research question of the repository:

> Can a global UV-cell occupancy objective be represented **exactly** as a
> quadratic function of antenna-selection variables?

The answer, derived below, is: **only in a special case**, and the general
case forces a choice between an approximation and a much larger problem.

---

## 1. The obstruction

Candidate pads $i = 1 \dots M$, binary selection $x_i$. A baseline exists iff
both ends are selected:

$$y_{ij} = x_i x_j$$

Let $C(i,j)$ be the set of UV cells that baseline $(i,j)$ touches — computed
offline, once, from geometry (snapshot or Earth-rotation track). The covered
set for a selection is

$$\mathrm{Cov}(x) = \bigcup_{i<j\,:\,y_{ij}=1} C(i,j)$$

and cell $c$ is covered iff

$$z_c = \bigvee_{(i,j)\,\ni\, c} y_{ij}
     = 1 - \prod_{(i,j)\,\ni\, c}\bigl(1 - y_{ij}\bigr)$$

Expanding that product gives monomials $y_{i_1j_1}y_{i_2j_2}\cdots$, i.e.
terms of degree up to $2m_c$ in $x$, where

$$m_c = \#\{\text{candidate pairs touching cell } c\}$$

So $\lvert\mathrm{Cov}(x)\rvert = \sum_c z_c$ is a **higher-order
pseudo-Boolean function**, not a quadratic form. The OR is the obstruction:
addition double-counts what a union does not.

### The one exact case

If $m_c \le 1$ for every cell — no two candidate pairs ever land in the same
cell — then each $z_c$ has a single term, the OR collapses, and

$$\lvert\mathrm{Cov}(x)\rvert = \sum_{i<j} \lvert C(i,j)\rvert\, x_i x_j$$

**exactly**. This is checked by `coverage_is_exactly_quadratic(cell_sets)` and
it is not a theoretical curiosity: it holds whenever the UV grid is fine
relative to the baseline spacing. In the toy problem (M=8, snapshot sampling)
it holds from $n_c = 32$ cells per axis upward.

But it comes with a catch, measured in experiment 05: at those resolutions
every feasible selection covers exactly the same number of cells
(each pair owns 2 cells, so 6 pairs always give 12). The objective becomes
*constant* and carries no information. **The regime where the quadratic form
is exact is exactly the regime where it is uninformative**, at least for
snapshot sampling of a small array. Under Earth-rotation sampling, where each
baseline sweeps a track through many cells, the toy problem never reaches
$m_c \le 1$ at any tested resolution — the quadratic form is never exact
there.

---

## 2. The four routes, with their costs

### A. Pairwise approximation (implemented, used by default)

Keep $H_{\rm uv} = -\sum_{i<j} w_{ij}x_ix_j$ and pick $w_{ij}$ to minimise the
damage. Three weightings are implemented:

| weighting | $w_{ij}$ | bias |
|---|---|---|
| `count` | $\lvert C(i,j)\rvert$ | upper bound: $\sum w_{ij}y_{ij}\ge\lvert\mathrm{Cov}\rvert$, over-counts every contested cell |
| `shared` | $\sum_{c\in C(i,j)} 1/m_c$ | credit for a contested cell split among all candidate claimants; exact when nothing is contested; neither bound in general |
| `radial` | `shared` with cell weight $1 + r_c/r_{\rm ref}$ | as `shared`, plus an explicit preference for the outer UV plane — a design choice, not a physical result |

**Measured behaviour** (toy problem, $\binom84 = 70$ feasible selections,
Spearman rank correlation against exact $\lvert\mathrm{Cov}(x)\rvert$, and the
coverage shortfall of the surrogate's own optimum):

Earth-rotation sampling, where the problem is never exactly quadratic:

| cells/axis | contested cells | `count` $\rho$ | `shared` $\rho$ | `radial` $\rho$ | best shortfall |
|---|---|---|---|---|---|
| 8 | 24 | 0.878 | 0.896 | 0.905 | 16.0 % |
| 16 | 48 | 0.952 | 0.964 | 0.944 | 7.7 % |
| 32 | 60 | 0.963 | 0.981 | 0.972 | 0 % (`count`) |
| 64 | 44 | 0.988 | 0.991 | 0.986 | 0 % |
| 256 | 20 | 0.999 | 0.999 | 0.989 | 0.5 % |

Snapshot sampling, coarse grids (the only informative snapshot regime):

| cells/axis | contested | `count` | `shared` | `radial` |
|---|---|---|---|---|
| 8 | 10 | $\rho=0.618$, 16.7 % short | $\rho=0.806$, 0 % short | $\rho=0.771$, 8.3 % short |
| 16 | 10 | surrogate is constant, 33.3 % short | $\rho=0.742$, 0 % short | $\rho=0.515$, 0 % short |

Readings that the data support:

1. The pairwise surrogate is a **good rank approximation** and a **poor value
   approximation**. Rank correlation is 0.88–0.999 under Earth-rotation
   sampling, but the surrogate's numerical value is not the cell count and
   must never be reported as one.
2. `shared` beats `count` almost everywhere, and by a lot where cells are
   heavily contested. Removing the systematic double-count is worth more than
   any clever weighting.
3. Even at $\rho \approx 0.96$ the surrogate's *argmax* can miss the true
   optimum by 5–8 % of coverage. Rank correlation does not imply optimum
   recovery.
4. `radial` is not an accuracy improvement — it is a different objective. It
   scores worse against unweighted cell count *by construction*, because it is
   optimising something else. That is a legitimate choice, but it must be
   argued from science requirements, not from these numbers.

### B. Auxiliary variables (implemented, exact, expensive)

Introduce $y_{ij}$ and $z_c$ as genuine variables:

- $y_{ij} = x_ix_j$ via the Rosenberg penalty
  $\lambda(x_ix_j - 2x_iy - 2x_jy + 3y)$, which is 0 iff the equality holds
  and positive otherwise;
- $z_c \le \sum_{(i,j)\ni c} y_{ij}$ via a binary-encoded slack $s_c$ and the
  equality penalty $\lambda\bigl(\sum y_{ij} - z_c - s_c\bigr)^2$;
- reward $-\sum_c z_c$.

The ground state then maximises $\lvert\mathrm{Cov}(x)\rvert$ **exactly**,
with no pairwise approximation.
[`tests/test_qubo.py::test_exact_auxiliary_qubo_reproduces_coverage`]

The price, for the toy problem alone ($M = 8$, $16^2$ grid):

| block | variables |
|---|---|
| $x$ (pads) | 8 |
| $y$ (pairs) | 28 |
| $z$ (cells) | 46 |
| slack | 56 |
| **total** | **138** |

138 variables to select 4 antennas from 8. The growth is
$O(M^2)$ for $y$ and $O(\lvert\text{cells}\rvert\log m_c)$ for $z$ and the
slack — and the cell count grows with both grid resolution and Earth-rotation
sampling. For $M = 100$ pads on a realistic grid this is far beyond current
hardware. The construction is nonetheless worth having: it is the ground truth
against which approximations are measured, and it is tractable for small
validation problems.

### C. Redundancy-aware penalty

Approximate uniqueness by penalising overlap directly:

$$H_{\rm uv} = -\sum_{i<j} \lvert C(i,j)\rvert x_ix_j
 + \mu \sum_{(i,j)<(k,l)} \lvert C(i,j)\cap C(k,l)\rvert\, y_{ij}y_{kl}$$

The second term is quartic in $x$, so it needs quadratisation too — it is
approach B by another route, with the same variable blow-up. The `shared`
weighting is the *mean-field* version of this idea: instead of penalising
overlap pair-by-pair, it discounts each cell by how contested it is on
average. That is why `shared` works as well as it does, and also why it has no
bound guarantee — the discount is computed over all *candidate* claimants, not
the selected ones.

### D. Alternative surrogate objectives

Not coverage at all, but a quantity that is genuinely pairwise: baseline
length diversity, angular diversity, or a target radial distribution matched
pair-by-pair. These are exactly quadratic by construction and need no
approximation — but they optimise a proxy chosen by the designer, and the link
from proxy to imaging performance would have to be established separately.
Not implemented; noted as an option.

---

## 3. Recommendation, and what would change it

For the *current* state of the project: **use approach A with the `shared`
weighting**, report the surrogate value as a surrogate, and validate the
selected configuration by computing its exact coverage and PSF afterwards. The
rank correlation is high enough to make the optimiser useful and the
verification step is cheap.

Use approach B for any problem small enough to afford it, precisely so that
the shortfall of approach A on that problem is known rather than assumed.

This recommendation should be revisited if:

- the science case turns out to need the *value* of coverage rather than a
  ranking (approach A cannot supply it);
- the grid resolution needed is fine enough that $m_c \le 1$ holds while the
  objective is still informative (then A is exact and the question dissolves);
- hardware makes a few thousand variables routine (then B is simply better).

## 4. Other terms

- **Selection** and **shadow** are natively quadratic; no approximation.
  Expansions in §12–13 of
  [mathematical_formulation.md](mathematical_formulation.md).
- **PWV** is linear and exact — but the cost vector is synthetic (see
  `thz_opt.constraints.pwv`), so the term is plumbing, not physics.
- **Phase** is the worst case: path existence from an outer station to the
  core is a *global* graph property, not quadratic and not even
  submodular-friendly. The implemented `pairwise_phase_penalty` penalises long
  pairs having no single available stepping stone — a strictly weaker,
  one-hop, selection-independent relaxation. It is a placeholder for a
  physically justified coherence model that the project does not yet have, and
  is labelled as such in code, docs and outputs.
