# How to run any of it

Everything here has been run on a clean checkout. If a command below does not
work, that is a bug worth reporting rather than a local problem.

---

## Setup, once

```bash
git clone https://github.com/manas-vamsi/thz-quantum-optimization
cd thz-quantum-optimization
pip install -r requirements.txt
```

Python 3.10 or later. No commercial solver, no cloud service, no API key.

Check it works:

```bash
python -m pytest
```

Expect **272 passed**, in under half a minute.

---

## Design an array yourself

This is the quickest way to see what the work does.

```bash
python design_array.py --list-sites
```

**Pick 20 antennas from ALMA's real pads:**

```bash
python design_array.py --site alma --n 20 --m 80 --dish 12 -o runs/alma20
```

**Design a new array at Hanle, on real terrain:**

```bash
python design_array.py --site hanle --n 16 --dish 8 --max-baseline 3000 -o runs/hanle16
```

Add `--quick` to trade Earth rotation for a snapshot: seconds instead of
minutes, useful for exploring. A full run is 2 to 8 minutes, nearly all of it in
the UV-track precompute, which grows as the square of the candidate count.

### What you get

| File | Contents |
|---|---|
| `pads.csv` | Every pad: latitude, longitude, elevation, local east/north, radius |
| `baselines.csv` | Every pairwise distance — what a cable or transporter plan needs |
| `array.cfg` | CASA observatory configuration; loads directly into `simobserve` |
| `report.json` | Resolution, largest angular scale, primary beam, UV coverage, sidelobe |
| `run.yaml` | Every parameter used — feed it back with `--config` for an identical rerun |

### On penalty weights

There is deliberately no `--penalty` or `--lambda` option, and a test prevents
one being added later.

The classical search swaps one pad in and one out, so the antenna count never
changes and there is no cardinality constraint to penalise. Penalties exist only
on the QUBO route, and there they are **derived**:

```bash
python design_array.py --site alma --n 12 --m 40 --quick --show-penalties
```

prints the smallest weights that provably cannot be cheated. A hand-entered
value below them makes the formulation wrong **silently** — no error, just a bad
answer. Better to display it than to invite it.

---

## Reproduce the results

Each experiment is standalone and writes its figures to `figures/` and its
numbers to `data/results/`.

```bash
cd experiments

python 11_real_alma.py           # real ALMA pads + measured atmosphere
python 12_qubo_solver.py         # the QUBO benchmark against certified optima
python 13_science_case.py        # the Fisher-derived objective
python 14_imaging_validation.py  # end-to-end imaging: does the metric mean anything
python 15_instrumental_limits.py # primary beam, smearing, real sensitivity
python 16_factorial_layout.py    # radial law versus angular law
```

Experiments 01 to 10 are the earlier controlled studies on synthetic pad fields.

**Runtimes:** 11 and 12 take a few minutes each (12 includes MILP certification,
which is the slow part). 13 to 16 are one to three minutes. Fixed seeds
throughout, so cell counts and correlations reproduce exactly; only wall-clock
timings vary with machine load, which is why the manuscript quotes speed
*ratios* rather than absolute seconds.

---

## The Ladakh design

```bash
cd ladakh_array_design
python run_design.py
```

Takes about eight minutes and downloads terrain on the first run (cached
afterwards). Produces three configurations at two sites, with pad coordinates,
CASA files and figures in `outputs/`.

Read `ladakh_array_design/README.md` first — it states plainly what is measured,
what is modelled, and what is not modelled at all.

---

## The manuscript

```bash
python paper/build_pdf.py              # rebuild the PDF
python paper/verify_citations.py       # check every DOI and arXiv id resolves
python paper/sync_code.py --check      # confirm paper/ matches the repository
python paper/build_pdf.py --check-metadata   # what still blocks submission
```

`paper/` is self-contained: the code, figures and data behind every number are
copied into it with SHA-256 checksums, so the manuscript can be sent to a
reviewer as one folder without the repository.

---

## Where things live

```
docs/objective_function.md     the derivation you asked for
paper/manuscript.pdf           the write-up, 13 pages
paper/code/                    self-contained copy, checksummed
src/thz_opt/                   the library, 56 modules
  qubo/                        formulation, Rosenberg, penalty floors
  optimize/                    greedy, local search, annealing, MILP
  imaging/                     visibilities, CLEAN, primary beam, smearing
  science/                     the Fisher-derived objective
  sites/                       real site data and terrain
experiments/                   17 numbered, standalone
ladakh_array_design/           the Ladakh design and its outputs
tests/                         272 checks
handover/                      this folder
```

---

## If something disagrees with the paper

Re-run the experiment and trust the output. Every table in the manuscript is
generated, not typed, and `paper/sync_code.py --check` will tell you if the
paper's copy has drifted from the repository.
