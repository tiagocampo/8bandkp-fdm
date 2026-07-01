# Issue 8: Minor carry-over fixes

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 23

## What to build

Two trivial fixes from the minor carry-over items in `docs/plans/BACKLOG.md`.

### Fix 1: Richardson absorption_edge JSON gap

The Richardson convergence helpers (`tests/integration/convergence_helpers.py`) are supposed to extract and store `absorption_edge` as an observable in the JSON results. The code to compute it exists, but the observable is absent from the stored JSON results.

Investigate why:
1. Read `convergence_helpers.py` and find the `extract_absorption_edge` function
2. Read the convergence test scripts that use it (e.g., `test_qw_grid_convergence.py`)
3. Check the JSON output format to see where `absorption_edge` should appear
4. Fix the gap — likely a missing key in the results dict or a conditional that skips the extraction

### Fix 2: Replace return None with RuntimeError

In `validation/bulk/test_bulk_zeeman.py` at line 44, there is a `return None` that should be `raise RuntimeError(...)` for consistency with the rest of the validation pipeline. If the function returns None instead of raising, downstream code may silently skip the test or produce confusing errors.

Find the `return None` at line 44 and replace with `raise RuntimeError("descriptive message")` following the pattern of other validation scripts.

## Acceptance criteria

- [ ] `absorption_edge` observable present in Richardson convergence JSON results after running convergence tests
- [ ] `test_bulk_zeeman.py:44` uses `raise RuntimeError` instead of `return None`
- [ ] All existing tests pass

## Blocked by

None — can start immediately.
