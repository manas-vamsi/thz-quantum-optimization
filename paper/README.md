# Atmosphere-Aware QUBO Formulations for Discrete-Pad Terahertz Interferometric Array Design

Manuscript, source code, figures, data and verification suite for the study.
This directory is self-contained: everything needed to read the paper,
understand the approach, rerun every experiment and check every number is here.

```
paper/
├── manuscript.md              the paper
├── manuscript.pdf             rendered, six pages, two columns
├── README.md                  this file: approach, results, how to run
├── REPRODUCIBILITY.md         which script produced which figure, table and number
├── build_pdf.py               manuscript.md -> manuscript.pdf
├── sync_code.py               refresh code/ from the project repository
├── code/                      the implementation
│   ├── MANIFEST.md            checksum of every shipped file, and its source commit
│   ├── src/thz_opt/           library
│   ├── experiments/           the five scripts cited by the manuscript
│   ├── tests/                 144 checks on the mathematics
│   ├── configs/               parameters every experiment reads
│   └── requirements.txt
├── figures/                   the six figures used in the manuscript
└── data/                      the numbers behind every table
```

---

## 1. The problem

A terahertz interferometer is built on prepared concrete pads. There are `M`
candidate pads and budget for `N` antennas. Which `N`?

Every *pair* of occupied pads forms one baseline, and each baseline measures one
Fourier component of the sky. The quality of a layout therefore depends on
pairs, not on individual pads, and the number of candidate layouts is `C(M,N)` --
about 10^29 at `M=100, N=50`.

Terahertz adds three constraints that centimetre-wave arrays can largely ignore:
atmospheric water vapour, tropospheric phase decorrelation that worsens with
baseline length, and a minimum separation set by dish diameter.

## 2. The approach

**Pairwise structure makes this a natural QUBO.** With `x_i = 1` when pad `i` is
occupied, a baseline exists exactly when `x_i x_j = 1`, so a quadratic objective
is the physically correct shape:

```
H = H_select + H_shadow + H_phase + H_uv
```

**All physics is precomputed.** Baselines, Earth-rotation UV tracks and per-cell
occupancy are computed once, offline, and compiled into coefficients. The
optimiser never runs a simulation.

**The obstruction, and the resolution.** Objectives that compare one baseline
against another -- UV-track overlap, unique-cell coverage, sidelobe energy --
involve four pad indices and are therefore *quartic* in `x`, not quadratic. A
pad-space `Q_ij` cannot hold them.

Introducing one baseline-activation variable per candidate pair, `y_k = x_i x_j`,
and enforcing that identity with the Rosenberg penalty restores an exact
quadratic form at a cost of `M + M(M-1)/2` variables. In `y`, the gridded
sampling-density energy, Cornwell repulsion and target-density matching are
represented *exactly*; unique-cell coverage is represented by a second-order
Bonferroni *lower bound*, which is exact unless a cell is covered three or more
times.

**Atmospheric coherence is the terahertz-specific term** and, usefully, it is
already quadratic: the decorrelation factor `gamma = exp(-sigma_phi^2/2)`
depends only on a pad pair, so it enters `Q_ij` directly.

Full derivations: `manuscript.md` §2, and `code/src/thz_opt/qubo/baseline_qubo.py`.

## 3. Results

| Finding | Measurement |
|---|---|
| Baseline variables beat a pad-only surrogate | selects an exact-coverage optimum at every tested UV resolution; the surrogate misses it at four of eight |
| Search beats analytic layouts | 763 occupied UV cells against 417 for a golden spiral and 610 for the best analytic reference, same pads and constraints |
| Heuristics are near-optimal where that can be proven | on a certifiable instance the optimum is 154 cells; simulated annealing attains 154 in 6.3 s against 54.5 s for the certifying solver |
| Receiver bandwidth is free coverage | +43.2 % occupied cells at 30 % fractional bandwidth, no extra antennas and no extra variables |
| The science case determines the layout | compact-source and extended-emission weightings share only four of fifteen pads |
| A preference is not enough | smooth radial weighting barely moves the optimum; only restricting the UV band does |
| Reconfiguration substitutes for Earth rotation | 2.66x coverage gain for a 0.2 h window, 1.06x for 6 h, at fixed integration time |

The study is a validated classical and QUBO-compatible baseline. It is **not** a
claim of quantum advantage, and **not** a site-specific array recommendation:
all pad coordinates and atmospheric constants are synthetic controlled inputs.

## 4. Running it

```bash
cd paper/code
pip install -r requirements.txt

python -m pytest                 # 144 checks on the mathematics itself
```

Then reproduce the manuscript's figures and tables:

```bash
cd experiments
python 07_formulation_comparison.py     # Figure 1,   Table 4.1
python 08_optimizer_benchmark.py        # Figures 2-3, Tables 4.2 and 4.2.1
python 09_science_cases.py              # Figures 4-5, Tables 4.3 and 4.4
python 10_multiepoch.py                 # Figure 6,   Table 4.5
python 06_layout_shootout.py            # Table 4.6
```

Experiment 08 takes roughly ten minutes because it runs two mixed-integer
solves; the others are faster. Every run is deterministic apart from
mixed-integer wall-clock time, which is why the manuscript reports the solver's
`proven_optimal` flag rather than assuming optimality.

Rebuild the PDF after regenerating figures:

```bash
python paper/build_pdf.py
```

## 5. Verifying what is here

`code/MANIFEST.md` lists a checksum for every shipped file and the commit it was
taken from. To confirm this copy matches the project repository:

```bash
python paper/sync_code.py --check
```

`REPRODUCIBILITY.md` maps every figure, table and equation in the manuscript to
the file that produces or implements it.

## 6. What is still required before submission

Institutional affiliation, correspondence address, coauthor information,
funding statement, acknowledgements, and a repository DOI. Scientifically: real
candidate pad coordinates, measured site atmospheric phase statistics to replace
the assumed constants, and a defined science case from which the UV weighting
should be derived rather than chosen.
