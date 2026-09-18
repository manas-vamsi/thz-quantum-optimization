"""Classical optimisers for pad selection.

These exist for one reason: **without a classical baseline, no quantum result
means anything.** A referee's first question about any quantum-annealing claim
is "did it beat simulated annealing?", and until now this repository could
score a configuration but never search for one.

All four optimisers preserve the cardinality constraint *by construction* --
they only ever consider selections with exactly ``n_select`` pads -- so no
selection penalty is needed and no infeasible state is ever evaluated. The
shadowing constraint is enforced the same way, by refusing moves that create a
violating pair, which is both faster and less arbitrary than a penalty weight.

Objectives are passed as a callable ``score(state) -> float`` that is
**minimised**. Ready-made ones are in :mod:`thz_opt.optimize.objectives`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .state import UVState

__all__ = [
    "OptimizerResult",
    "greedy_removal",
    "greedy_addition",
    "local_search",
    "simulated_annealing",
]


@dataclass
class OptimizerResult:
    """What an optimiser found, and what it cost to find it."""

    selection: np.ndarray            # boolean mask over pads
    score: float                     # objective value (lower is better)
    unique_cells: int
    sum_sq: int
    n_evaluations: int
    method: str
    history: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def indices(self) -> np.ndarray:
        return np.flatnonzero(self.selection)

    def as_dict(self) -> dict:
        return {
            "method": self.method,
            "score": float(self.score),
            "unique_cells": int(self.unique_cells),
            "sum_sq": int(self.sum_sq),
            "n_evaluations": int(self.n_evaluations),
            "n_selected": int(self.selection.sum()),
            "selection": self.indices().tolist(),
            **self.meta,
        }


def _forbidden_matrix(n_pads: int, forbidden_pairs) -> np.ndarray:
    """Boolean ``(n, n)`` matrix of pad pairs that may not both be selected."""
    f = np.zeros((n_pads, n_pads), dtype=bool)
    if forbidden_pairs is None:
        return f
    arr = np.asarray(forbidden_pairs, dtype=bool)
    if arr.shape == (n_pads, n_pads):
        f = arr.copy()
    else:
        for i, j in forbidden_pairs:
            f[int(i), int(j)] = f[int(j), int(i)] = True
    np.fill_diagonal(f, False)
    return f


def _compatible(p: int, selected: np.ndarray, forbidden: np.ndarray) -> bool:
    """True if pad ``p`` can join ``selected`` without a forbidden pair."""
    return not np.any(forbidden[p] & selected)


def _finalise(state: UVState, score, n_eval: int, method: str,
              history: list, meta: dict) -> OptimizerResult:
    return OptimizerResult(
        selection=state.copy_selection(),
        score=float(score(state)),
        unique_cells=state.unique_cells,
        sum_sq=state.sum_sq,
        n_evaluations=n_eval,
        method=method,
        history=history,
        meta=meta,
    )


# --------------------------------------------------------------------------
# greedy
# --------------------------------------------------------------------------

def greedy_removal(state: UVState, n_select: int, score, forbidden_pairs=None,
                   track_history: bool = False) -> OptimizerResult:
    """Start from all pads and drop them one at a time.

    At each step the pad whose removal leaves the best objective is removed,
    until ``n_select`` remain. This is the classical method of Panduranga Rao
    et al. (2009) applied to exactly this problem, and it is the natural
    baseline: it sees the full array first and decides what to give up.

    Cost: ``O((M - N) * M)`` objective evaluations, each ``O(N)`` with the
    incremental state.
    """
    n = state.n_pads
    if not 0 <= n_select <= n:
        raise ValueError("n_select must be between zero and the number of pads")
    forbidden = _forbidden_matrix(n, forbidden_pairs)
    start = "all pads"
    if forbidden.any():
        # Starting from every pad would violate a hard shadowing constraint.
        # A maximal independent set gives removal the largest feasible starting
        # point we can build cheaply while preserving feasibility throughout.
        _maximal_feasible(state, forbidden, n_select)
        start = "maximal feasible set"
    else:
        state.set_selection(np.arange(n))
    n_eval, history = 0, []

    while state.n_selected > n_select:
        best_p, best_val = -1, np.inf
        for p in state.selection_indices():
            state.remove_pad(int(p))
            val = score(state)
            n_eval += 1
            state.add_pad(int(p))
            if val < best_val:
                best_val, best_p = val, int(p)
        state.remove_pad(best_p)
        if track_history:
            history.append({"n_selected": state.n_selected, "score": float(best_val)})

    return _finalise(state, score, n_eval, "greedy_removal", history,
                     {"start": start})


def greedy_addition(state: UVState, n_select: int, score, forbidden_pairs=None,
                    seed_pad: int | None = None,
                    track_history: bool = False) -> OptimizerResult:
    """Start from one pad and add the best pad at each step.

    Cheaper than removal (``O(N * M)`` evaluations) and usually a little worse,
    because early choices are made with no view of the final array.
    """
    n = state.n_pads
    if not 0 <= n_select <= n:
        raise ValueError("n_select must be between zero and the number of pads")
    forbidden = _forbidden_matrix(n, forbidden_pairs)
    state.set_selection([])
    n_eval, history = 0, []

    if n_select == 0:
        return _finalise(state, score, n_eval, "greedy_addition", history,
                         {"seed_pad": None, "reached_target": True})

    first = 0 if seed_pad is None else int(seed_pad)
    state.add_pad(first)

    while state.n_selected < n_select:
        best_p, best_val = -1, np.inf
        for p in range(n):
            if state.selected[p] or not _compatible(p, state.selected, forbidden):
                continue
            state.add_pad(p)
            val = score(state)
            n_eval += 1
            state.remove_pad(p)
            if val < best_val:
                best_val, best_p = val, p
        if best_p < 0:
            break  # nothing compatible left
        state.add_pad(best_p)
        if track_history:
            history.append({"n_selected": state.n_selected, "score": float(best_val)})

    return _finalise(state, score, n_eval, "greedy_addition", history,
                     {"seed_pad": first, "reached_target": state.n_selected == n_select})


def _maximal_feasible(state: UVState, forbidden: np.ndarray, n_select: int) -> None:
    """Initialise a large feasible set for constrained greedy removal.

    The order of a greedy independent-set construction matters, so try a
    deterministic order plus fixed-seed shuffled orders and retain the largest
    result. This is a start-state heuristic, not a claim of maximum
    independent-set optimality.
    """
    rng = np.random.default_rng(0)
    best = np.empty(0, dtype=int)
    orders = [np.arange(state.n_pads)]
    orders.extend(rng.permutation(state.n_pads) for _ in range(63))
    for order in orders:
        chosen = np.zeros(state.n_pads, dtype=bool)
        for p in order:
            if _compatible(int(p), chosen, forbidden):
                chosen[int(p)] = True
        if chosen.sum() > best.size:
            best = np.flatnonzero(chosen)
    if best.size < n_select:
        raise ValueError(
            f"could not construct {n_select} mutually compatible pads; "
            f"largest sampled feasible set has {best.size}"
        )
    state.set_selection(best)


# --------------------------------------------------------------------------
# local search and annealing
# --------------------------------------------------------------------------

def local_search(state: UVState, score, forbidden_pairs=None,
                 max_passes: int = 50) -> OptimizerResult:
    """Best-improvement 1-swap hill climbing from the current selection.

    Repeatedly exchanges one selected pad for one unselected pad whenever that
    improves the objective, until no single swap helps. Terminates at a local
    optimum of the 1-swap neighbourhood.
    """
    forbidden = _forbidden_matrix(state.n_pads, forbidden_pairs)
    n_eval, passes = 0, 0
    current = score(state)

    improved = True
    while improved and passes < max_passes:
        improved = False
        passes += 1
        for out_p in state.selection_indices().copy():
            state.remove_pad(int(out_p))
            best_in, best_val = -1, current
            for in_p in range(state.n_pads):
                if state.selected[in_p] or in_p == out_p:
                    continue
                if not _compatible(in_p, state.selected, forbidden):
                    continue
                state.add_pad(in_p)
                val = score(state)
                n_eval += 1
                state.remove_pad(in_p)
                if val < best_val - 1e-12:
                    best_val, best_in = val, in_p
            if best_in >= 0:
                state.add_pad(best_in)
                current, improved = best_val, True
            else:
                state.add_pad(int(out_p))

    return _finalise(state, score, n_eval, "local_search", [],
                     {"passes": passes})


def simulated_annealing(state: UVState, n_select: int, score,
                        forbidden_pairs=None, n_iterations: int = 20000,
                        t_start: float | None = None, t_end: float | None = None,
                        seed: int = 0, track_history: bool = False) -> OptimizerResult:
    """Metropolis annealing over 1-swap moves at fixed cardinality.

    The move set (swap one pad in, one out) preserves ``n_select`` exactly, so
    the search never leaves the feasible set and no cardinality penalty is
    needed -- which also means the penalty-weight question of
    ``docs/objective_function.md`` §8 does not arise here at all. That is worth
    noting when comparing against a QUBO solver, which *does* have to pay it.

    Temperatures default to a geometric schedule calibrated from the spread of
    objective values seen in a short random probe, so the caller does not have
    to guess a scale.
    """
    rng = np.random.default_rng(seed)
    forbidden = _forbidden_matrix(state.n_pads, forbidden_pairs)

    if state.n_selected != n_select:
        _random_feasible(state, n_select, forbidden, rng)

    current = score(state)
    best_val, best_sel = current, state.copy_selection()
    n_eval = 0

    if t_start is None or t_end is None:
        probe = []
        for _ in range(min(200, 10 * state.n_pads)):
            mv = _propose(state, forbidden, rng)
            if mv is None:
                continue
            out_p, in_p = mv
            state.remove_pad(out_p)
            state.add_pad(in_p)
            probe.append(score(state))
            n_eval += 1
            state.remove_pad(in_p)
            state.add_pad(out_p)
        spread = float(np.std(probe)) if len(probe) > 2 else 1.0
        spread = max(spread, 1e-9)
        t_start = spread if t_start is None else t_start
        t_end = spread * 1e-3 if t_end is None else t_end

    history = []
    for it in range(n_iterations):
        t = t_start * (t_end / t_start) ** (it / max(n_iterations - 1, 1))
        mv = _propose(state, forbidden, rng)
        if mv is None:
            continue
        out_p, in_p = mv
        state.remove_pad(out_p)
        state.add_pad(in_p)
        candidate = score(state)
        n_eval += 1

        delta = candidate - current
        if delta <= 0 or rng.random() < np.exp(-delta / t):
            current = candidate
            if current < best_val:
                best_val, best_sel = current, state.copy_selection()
        else:
            state.remove_pad(in_p)
            state.add_pad(out_p)

        if track_history and it % max(n_iterations // 200, 1) == 0:
            history.append({"iteration": it, "temperature": float(t),
                            "current": float(current), "best": float(best_val)})

    state.set_selection(np.flatnonzero(best_sel))
    return _finalise(state, score, n_eval, "simulated_annealing", history,
                     {"n_iterations": n_iterations, "t_start": float(t_start),
                      "t_end": float(t_end), "seed": seed})


def _propose(state: UVState, forbidden: np.ndarray, rng):
    """Pick a legal (out, in) swap, or None if none found quickly."""
    sel = state.selection_indices()
    if sel.size == 0:
        return None
    for _ in range(20):
        out_p = int(rng.choice(sel))
        in_p = int(rng.integers(state.n_pads))
        if state.selected[in_p]:
            continue
        mask = state.selected.copy()
        mask[out_p] = False
        if np.any(forbidden[in_p] & mask):
            continue
        return out_p, in_p
    return None


def _random_feasible(state: UVState, n_select: int, forbidden: np.ndarray, rng) -> None:
    """Random shadow-free selection of the requested size."""
    for _ in range(200):
        state.set_selection([])
        order = rng.permutation(state.n_pads)
        for p in order:
            if state.n_selected >= n_select:
                break
            if _compatible(int(p), state.selected, forbidden):
                state.add_pad(int(p))
        if state.n_selected == n_select:
            return
    raise RuntimeError("could not build a feasible selection of the requested size")
