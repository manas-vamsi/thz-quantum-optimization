# What remains, and how to continue

Three categories, kept separate because they need different things: some need a
decision, some need work, and some cannot be fixed by us at all.

---

## 1. Decisions, not work

### Which story does the paper tell?

The evidence supports: *"we formulated the problem properly, solved it, and on
this problem the quantum route does not currently pay."* That is honest,
defensible, and the kind of negative result reviewers trust.

The evidence does **not** support: *"quantum optimization of terahertz
interferometric architectures."* We ran the QUBO and it lost to a classical
heuristic by 35 to 72 times.

There is a third option worth considering: reframe the project around what it
actually produced — a validated formulation, certified optima, a science-derived
objective, and a terrain-constrained design for an Indian site. That is a
stronger paper than a quantum result we cannot support.

This is a strategy decision. It is yours and Manas's.

### Coauthorship

Your entry is prepared in `paper/metadata.yaml`, commented out. Nobody is added
to a paper without having read it and agreed. If there are others at SAC who
should be on it, they need naming and a chance to review.

### Submission metadata

Seven fields block submission. `python paper/build_pdf.py --check-metadata`
lists them: affiliation, correspondence email, funding, acknowledgements,
competing interests, target venue, repository DOI.

The DOI is the only one needing action rather than an answer — `.zenodo.json` is
prepared, so archiving the repository and pasting back the minted DOI is about
ten minutes.

### Target venue

Not chosen. The work fits an instrumentation venue (SPIE, JATIS) or an
optimisation-applications venue better than a pure astronomy journal, because
the contribution is the formulation and its evaluation rather than an
observation.

---

## 2. The one input I could not derive

**What should this array actually observe?**

Everything downstream depends on it and nothing in the code can supply it.

The project currently optimises occupied UV cells, which is a number invented
for convenience. Section 4.3.1 of the manuscript replaces it with Fisher
information on a protoplanetary disc — but that source was chosen by me as a
demonstration, not because anyone decided it is the science case.

Once a real target exists, the chain runs itself:

```
science target
  → the angular scales that matter
  → the baselines that deliver them
  → the UV weighting          (derived, not chosen)
  → the array
```

Concretely, the choice determines: which of the three Ladakh configurations to
build first, whether short spacings matter enough to justify a compact array at
all, and whether the flux-recovery failure in §5 of RESULTS.md is a problem or
an irrelevance.

**This is the highest-value thing you can supply.** A single sentence — "we want
to measure X in sources of type Y" — unblocks more than any amount of further
optimisation.

---

## 3. Technical work, in priority order

### 3.1 Optimise under the full model, not just measure under it

**The clearest remaining gap.** Imaging, primary beam and smearing are
implemented, but nothing is *optimised* under them. Arrays are still chosen by
coverage objectives and then measured under the fuller model.

Section 5 of RESULTS.md shows why this matters: the cell-maximising array
recovers 42% of a smooth source where a golden spiral recovers 91%. An objective
that carried flux recovery would not make that mistake.

**Difficulty:** the honest problem is that such an objective is **not
quadratic**, so it does not fit the QUBO. That is a real finding, not an
obstacle to route around, and it should be stated if the work continues in this
direction.

### 3.2 Real visibility data

No real observation has ever been processed. The imaging is simulated —
correctly, and validated to 2×10⁻¹⁶ on the transform pair, but simulated.

Running the pipeline on an archival ALMA or GMRT dataset would convert
"predicted image quality" into "measured image quality". If SAC has access to
archival visibilities, this is straightforward and would strengthen the paper
considerably.

### 3.3 Multi-scale CLEAN

Deconvolution uses Högbom, chosen for transparency. A multi-scale algorithm
recovers extended flux better and would narrow — though not close — the gap in
§5. Worth doing before anyone cites the 42% figure as a hard limit.

### 3.4 Quantum hardware, if the decision is to pursue it

The bar is set and documented: beat a 0.35-second classical heuristic already
within 2% of proven optimal. The embedding cost is the obstacle — a 40-pad
instance needs 820 variables and 33,689 couplings on a dense graph; the full
174-pad field would need roughly 15,000 variables.

If this is attempted, `dwave-system` and an API token are the only additions;
the QUBO object submits unchanged.

---

## 4. What cannot be fixed by us

### Phase stability in Ladakh does not exist as a measurement

ALMA's structure function comes from 17,000 observations by a working
interferometer over four years. Nothing comparable can exist at a site with no
interferometer.

Every coherence figure for Hanle or site A in this work is **transferred from
Chajnantor and labelled a model**. It must not be quoted as a site measurement.

The only way to fix this is to put a phase monitor on the site and run it for
years. If an Indian submillimetre array is genuinely under consideration, that
measurement is on the critical path and starting it early costs little.

### Site selection is partly political

The data says site A (34.25°N, 78.75°E) is four times drier than Hanle. Hanle
has the road, the power, the staff and the flatter ground. We can quantify the
trade; we cannot make the decision.

---

## 5. Deliberate limitations, not oversights

These are stated in the manuscript and should stay stated:

- **Sections 4.1 to 4.6 use synthetic pad fields.** This is by design: comparing
  layout strategies requires holding N, extent and UV grid fixed, which real pad
  lists do not allow. Sections 4.2.2, 4.3.1, 4.7 and 4.8 use real ALMA geometry.
- **No w-term during imaging**, no polarization, no mosaicking, no calibration or
  pointing error.
- **Noise is white and Gaussian on gridded visibilities** — the standard
  idealisation, and optimistic.
- **The Ladakh design models terrain only.** No land ownership, access roads,
  power routing, geotechnics, snow or wind. It is preliminary geometry and says
  so on every output file.

---

## 6. If you pick this up

Read in this order:

1. `handover/README.md` — the story
2. `docs/objective_function.md` — the derivation
3. `paper/manuscript.pdf` §2 — the formulation as written up
4. `experiments/12_qubo_solver.py` — the result that matters most

Then run `python design_array.py --site alma --n 20 --quick` to see the whole
thing work in a few seconds.

The code is documented at the level of *why*, not *what* — module docstrings
explain the reasoning, including the approaches that were tried and abandoned
and the mistakes that were caught. Those are the most useful parts if you intend
to extend it.
