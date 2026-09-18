import numpy as np
import pytest

from thz_opt.optimize.heuristics import (
    greedy_addition,
    greedy_removal,
    local_search,
    simulated_annealing,
)
from thz_opt.optimize.objectives import max_unique_cells, min_sidelobe_energy
from thz_opt.optimize.state import UVState, build_pair_tables


def _state() -> UVState:
    # Pair order for four pads: 01, 02, 03, 12, 13, 23.
    multiplicities = [
        {0: 1, 1: 1}, {1: 1, 2: 1}, {2: 1, 3: 1},
        {0: 1, 3: 1}, {1: 1, 3: 1}, {0: 1, 2: 1},
    ]
    cells, counts = build_pair_tables(multiplicities)
    return UVState(cells, counts, n_pads=4, n_cells_total=4)


def _assert_consistent(state: UVState) -> None:
    unique, sum_sq = state.recompute()
    assert (state.unique_cells, state.sum_sq) == (unique, sum_sq)


def test_incremental_state_matches_recomputation_after_random_walk():
    state = _state()
    rng = np.random.default_rng(3)
    for _ in range(100):
        pad = int(rng.integers(state.n_pads))
        if state.selected[pad]:
            state.remove_pad(pad)
        else:
            state.add_pad(pad)
        _assert_consistent(state)


def test_classical_searches_keep_cardinality_and_improve_from_a_common_start():
    state = _state()
    objective = max_unique_cells()

    addition = greedy_addition(state, 3, objective)
    assert addition.selection.sum() == 3
    _assert_consistent(state)

    start_score = addition.score
    local = local_search(state, objective)
    assert local.selection.sum() == 3
    assert local.score <= start_score
    _assert_consistent(state)

    annealed = simulated_annealing(state, 3, objective, n_iterations=300, seed=7)
    assert annealed.selection.sum() == 3
    _assert_consistent(state)


def test_greedy_removal_respects_hard_forbidden_pairs():
    state = _state()
    forbidden = np.zeros((4, 4), dtype=bool)
    forbidden[0, 1] = forbidden[1, 0] = True
    forbidden[2, 3] = forbidden[3, 2] = True

    result = greedy_removal(state, 2, min_sidelobe_energy(), forbidden)
    assert result.selection.sum() == 2
    assert not np.any(np.outer(result.selection, result.selection) & forbidden)
    assert result.meta["start"] == "maximal feasible set"
    _assert_consistent(state)


def test_infeasible_cardinality_is_rejected_for_constrained_removal():
    state = _state()
    forbidden = np.ones((4, 4), dtype=bool)
    np.fill_diagonal(forbidden, False)
    with pytest.raises(ValueError, match="mutually compatible"):
        greedy_removal(state, 2, max_unique_cells(), forbidden)
