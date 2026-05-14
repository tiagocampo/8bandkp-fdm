---
date: 2026-05-13
topic: fit-tail-regression-direction
status: open
origin: code-review finding #1 from ce-code-review
---

# Fix: fit_tail_exponential Regression Direction for Right-Edge States

## Problem

`fit_tail_exponential` in `topological_analysis.f90` always runs its log-linear regression from `tail_start` forward to `n_points`, even when the backward search found the tail to the **left** of the density peak. For right-edge Majorana states (density peaked near the end of the wire), the backward search succeeds but the regression then sweeps over an **increasing** density region, producing a meaningless localization length.

### Walkthrough

For a right-edge Majorana with density `[0.01, 0.05, 0.2, 1.0]`:
1. `peak_idx = 4` (density 1.0)
2. Forward search from 5: fails (past end)
3. Backward search: `tail_start = 2` (density 0.05 < 0.1 * 1.0)
4. Regression from 2 to 4: density `[0.05, 0.2, 1.0]` — **increasing**, positive slope
5. `xi = |1/positive_slope|` — positive but meaningless

For a left-edge Majorana with density `[1.0, 0.2, 0.05, 0.01]` (current behavior):
1. `peak_idx = 1` (density 1.0)
2. Forward search: `tail_start = 2` (density 0.2 < 0.1 * 1.0)
3. Regression from 2 to 4: density `[0.2, 0.05, 0.01]` — **decreasing**, negative slope
4. `xi = |1/slope|` — correct

## Fix

Track which search direction found the tail. When backward search is used, run the regression from `tail_start` **down to 1** (the actual decay direction), measuring positions as distances from `positions(tail_start)` going backward.

### Files

- `src/physics/topological_analysis.f90` — `fit_tail_exponential` subroutine (lines ~467-555)
- `tests/unit/test_edge_states.pf` — add test for right-edge density profile

### Specific Changes

1. Add `logical :: forward_tail` variable
2. Set `forward_tail = .true.` on forward search success, `forward_tail = .false.` on backward search success
3. Guard `tail_start >= n_points` only applies to forward tails — backward tails go toward 1
4. Split regression loop:
   - Forward: `do i = tail_start, n_points` with `positions(i) - x_start` (current behavior)
   - Backward: `do i = tail_start, 1, -1` with `x_start - positions(i)` (distance going left)
5. Update `domain_extent` for R11 convergence check: use `positions(tail_start) - positions(1)` for backward tails

### Caller Impact

None. `fit_exponential_decay` and `compute_majorana_profile` call `fit_tail_exponential` and only consume `xi`, `n_fit_actual`, and `converged` — the fix is entirely internal.

## Test Scenarios

- Construct right-edge density: `exp(-(4-x)/5)` for x=1..4. Verify recovered xi is within 20% of 5.0.
- Construct left-edge density: `exp(-(x-1)/5)` for x=1..4 (existing path). Verify unchanged.
- Symmetric density centered at midpoint: both searches should work, verify consistent xi.

## Scope Boundaries

- No caller interface changes
- No changes to threshold, density floor, or convergence criteria
- No changes to `fit_exponential_decay` or `compute_majorana_profile`
