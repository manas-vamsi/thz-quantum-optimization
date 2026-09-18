# Reproducibility and provenance

Every number in the manuscript is written by a script in `code/experiments/`,
not entered by hand. This file maps each figure, table and equation to the code
behind it, so any claim in the paper can be traced to the line that produces it.

Paths are relative to this directory.

---

## Figures

| Figure | Shows | Produced by | Numbers written alongside |
|---|---|---|---|
| 1 | Formulation comparison: pad-space surrogate against baseline variables | `code/experiments/07_formulation_comparison.py` | `data/exp07_formulations.csv` |
| 2 | Classical optimiser comparison | `code/experiments/08_optimizer_benchmark.py` | `data/exp08_optimizer_benchmark.csv` |
| 3 | An optimised layout and its UV response | `code/experiments/08_optimizer_benchmark.py` | `data/exp08_summary.json` |
| 4 | Layouts selected by two science cases | `code/experiments/09_science_cases.py` | `data/exp09_science_cases.csv` |
| 5 | Multi-frequency UV filling | `code/experiments/09_science_cases.py` | `data/exp09_multifrequency.csv` |
| 6 | Multi-epoch coverage against reconfiguration cost | `code/experiments/10_multiepoch.py` | `data/exp10_multiepoch.csv` |

Each figure also carries its generating parameters printed in the corner of the
image, so a figure separated from this directory still states what produced it.

## Tables

| Section | Table | Produced by | Data file |
|---|---|---|---|
| 4.1 | Baseline-variable fidelity, M=8, N=4 | `07_formulation_comparison.py` | `data/exp07_formulations.csv` |
| 4.2 | Optimiser benchmark, M=120, N=20 | `08_optimizer_benchmark.py` | `data/exp08_optimizer_benchmark.csv` |
| 4.2.1 | MILP certification against problem size | `08_optimizer_benchmark.py` | `data/exp08_certified_instance.csv` |
| 4.2.1 | Heuristics against the proven optimum | `08_optimizer_benchmark.py` | `data/exp08_certified_instance.csv` |
| 4.3 | Science-case dependence | `09_science_cases.py` | `data/exp09_science_cases.csv` |
| 4.3 | Smooth preference against band restriction | `09_science_cases.py` | `data/exp09_soft_vs_hard.csv` |
| 4.4 | Multi-frequency synthesis | `09_science_cases.py` | `data/exp09_multifrequency.csv` |
| 4.5 | Multi-epoch reconfiguration | `10_multiepoch.py` | `data/exp10_multiepoch.csv` |
| 4.5 | Observing-window sweep | `10_multiepoch.py` | `data/exp10_window_sweep.csv` |
| 4.6 | Analytic layout families | `06_layout_shootout.py` | `data/exp06_shootout.csv` |

## Equations and claims

| Manuscript | Claim | Implementation |
|---|---|---|
| §2.1 | Baselines, UV coordinates, Earth-rotation tracks | `code/src/thz_opt/interferometry/baselines.py`, `uv.py`, `earth_rotation.py` |
| §2.2, Eq. (1) | Kolmogorov phase structure function, decorrelation factor | `code/src/thz_opt/constraints/coherence.py` |
| §2.2, Eq. (2) | Fixed antenna count; minimum separation | `code/src/thz_opt/qubo/objective.py`, `constraints/separation.py` |
| §2.3, Eq. (3) | The quartic UV-overlap obstruction | `code/src/thz_opt/qubo/coefficients.py` (derivation in the module docstring) |
| §2.4, Eq. (4) | Rosenberg quadratization | `code/src/thz_opt/qubo/baseline_qubo.py::build_baseline_qubo` |
| §2.4, Eq. (5) | Second-order Bonferroni coverage bound | `baseline_qubo.py::bonferroni_coverage_terms` |
| §2.4, Eq. (6) | Gridded sampling-density energy; Parseval relation | `baseline_qubo.py::sidelobe_energy_terms` |
| §2.4 | Coherence-weighted objective | `baseline_qubo.py::coherence_weighted_sidelobe_terms` |
| §3 | Mixed-integer certification and dual bounds | `code/src/thz_opt/optimize/exact.py` |
| §4.2 | Greedy, one-swap local search, simulated annealing | `code/src/thz_opt/optimize/heuristics.py` |
| §4.2 | Incremental UV occupancy used by every search | `code/src/thz_opt/optimize/state.py` |
| §4.3 | Science-case radial weightings | `code/src/thz_opt/optimize/science_cases.py` |
| §4.4 | Multi-frequency synthesis, `a_ck = sum_f a_ckf` | `code/src/thz_opt/interferometry/multifrequency.py` |
| §4.5, Eq. (8) | Multi-epoch selection and movement cost | `code/src/thz_opt/qubo/multiepoch.py` |
| §5 | Scenario-mean robustness; star cable cost | `code/src/thz_opt/qubo/robust.py`, `constraints/cable.py` |

---

## Verification suite

```bash
cd code
python -m pytest
```

144 checks. The ones that matter most for the manuscript's claims:

| Check | What it establishes |
|---|---|
| Matrix QUBO against a direct term-by-term objective over all 2^8 assignments | the assembled matrix is the objective, exactly |
| Upper-triangular against symmetric convention | the factor-of-two convention is consistent |
| Parseval identity on the gridded sampling function | Eq. (6) really is dirty-beam energy |
| Eq. (6) against direct gridding | the quadratic form is exact, to machine precision |
| Bonferroni expression never exceeds exact coverage | Eq. (5) is a genuine lower bound |
| Rosenberg penalty admits no cheaper auxiliary assignment | `y_k = x_i x_j` is enforced, not assumed |
| Earth rotation reduces to the zenith snapshot | the UV convention is correct |
| Dense and sparse QUBO assemblies agree | the large-instance path is the same model |
| Star cable cost never undercuts the minimum spanning tree | the linear cost term is an upper bound |

## Reproducing everything from scratch

```bash
cd code
pip install -r requirements.txt
python -m pytest

cd experiments
python 06_layout_shootout.py
python 07_formulation_comparison.py
python 08_optimizer_benchmark.py        # about ten minutes: two mixed-integer solves
python 09_science_cases.py
python 10_multiepoch.py
```

Outputs overwrite `figures/` and `data/results/` **inside the project
repository**, not this directory; run `python paper/sync_code.py` afterwards to
bring the refreshed figures and numbers back here, then `python
paper/build_pdf.py` to rebuild the document.

Determinism: all layouts are seeded from `code/configs/default.yaml`. The single
source of run-to-run variation is mixed-integer solver wall-clock time, which
can change whether a time-limited run closes its optimality gap. The manuscript
therefore reports the solver's `proven_optimal` flag rather than assuming
optimality, and distinguishes a proven optimum from an incumbent throughout.

## Scope of the data

All candidate pad coordinates and all atmospheric constants are synthetic,
chosen to make the comparisons controlled and fair. They are not a site survey,
and no result here should be read as a construction recommendation. Replacing
them with real pad positions and measured phase statistics is the stated next
step.
