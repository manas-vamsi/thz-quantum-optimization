"""Search, as opposed to scoring.

Before this package existed the repository could evaluate a configuration but
never look for a good one, which meant every "best layout" was the best of a
hand-made shortlist. These are the classical baselines that any quantum result
must be measured against:

* :func:`~thz_opt.optimize.heuristics.greedy_removal` -- drop pads one at a time
* :func:`~thz_opt.optimize.heuristics.greedy_addition` -- add pads one at a time
* :func:`~thz_opt.optimize.heuristics.local_search` -- 1-swap hill climbing
* :func:`~thz_opt.optimize.heuristics.simulated_annealing` -- the standard comparator
* :func:`~thz_opt.optimize.exact.max_coverage_milp` -- **provably optimal** at
  moderate M, via HiGHS

The heuristics share an incremental occupancy model
(:class:`~thz_opt.optimize.state.UVState`) so that moving one antenna costs
``O(N)`` rather than a full regrid, and they all preserve the cardinality and
shadowing constraints by construction rather than by penalty.
"""

from .exact import max_coverage_milp
from .heuristics import (
    OptimizerResult,
    greedy_addition,
    greedy_removal,
    local_search,
    simulated_annealing,
)
from .pareto import dominates, pareto_indices, sweep_weights
from .science_cases import (
    SCIENCE_CASES,
    compact_source_weights,
    extended_emission_weights,
    radial_weights,
    science_case_objective,
)
from .objectives import (
    coherence_weighted_cells,
    max_unique_cells,
    min_sidelobe_energy,
    weighted_objective,
)
from .state import UVState, build_pair_tables

__all__ = [
    "UVState",
    "build_pair_tables",
    "OptimizerResult",
    "greedy_removal",
    "greedy_addition",
    "local_search",
    "simulated_annealing",
    "max_coverage_milp",
    "max_unique_cells",
    "min_sidelobe_energy",
    "weighted_objective",
    "coherence_weighted_cells",
    "dominates",
    "pareto_indices",
    "sweep_weights",
    "SCIENCE_CASES",
    "compact_source_weights",
    "extended_emission_weights",
    "radial_weights",
    "science_case_objective",
]
