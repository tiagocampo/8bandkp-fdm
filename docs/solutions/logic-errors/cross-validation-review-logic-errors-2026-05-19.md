---
title: "Cross-validation code review: false-pass comparison, misleading test labels, fragile convergence"
date: "2026-05-19"
category: logic-errors
module: validation/shared
problem_type: logic_error
component: testing_framework
severity: high
symptoms:
  - "compare_eigenvalues returned passed=True when eigenvalue counts differed between Fortran and kdotpy"
  - "Four test files claimed cross-code validation but only compared against analytical formulas"
  - "SC convergence detection matched 'converged' in 'NOT converged'"
  - "param_mapper.py strain_Du docstring documented -1.5*b_dp while code used -0.75*b_dp"
root_cause: logic_error
resolution_type: code_fix
tags: [cross-validation, kdotpy, false-pass, code-review, eigenvalue-comparison, convergence]
related_components:
  - validation/shared/comparison.py
  - validation/shared/param_mapper.py
  - validation/landau/test_landau_bulk.py
  - validation/wire/test_wire_subbands.py
  - validation/gfactor/test_gfactor_qw.py
  - validation/selfconsistent/test_sc_qw.py
---

# Cross-validation code review: false-pass comparison, misleading test labels, fragile convergence

## Problem

The kdotpy cross-validation pipeline (12 Python tests comparing our Fortran 8-band k.p solver against kdotpy) contained five issues discovered during code review: a comparison function that silently passed on eigenvalue count mismatches, stale and misleading docstrings across 5 files, fragile convergence detection that could match failure messages, and a tolerance value that drifted from its documentation. Together these could cause the validation suite to report PASS when the codes actually disagreed.

## Symptoms

- `compare_eigenvalues()` in `comparison.py` compared only the overlapping subset when eigenvalue counts differed and returned `passed=True` if those matched
- `param_mapper.py:15` docstring claimed `strain_Du = -1.5 * b_dp * 1000` but code at line 248 used `-0.75 * b_dp * 1000`
- Four test files (Landau, wire, gfactor, SC) described themselves as "cross-code validation" but only compared against analytical formulas or hardcoded references
- `test_sc_qw.py:74` used `or "converged" in result.stdout.lower()` which would match "NOT converged"
- `test_sc_qw.py` docstring said "< 50 meV" but `TOL_CB1_MEV = 10.0`

## What Didn't Work

- The `strain_Du` formula was corrected from `-1.5` to `-0.75` in commit `b7d8841` (during the Bir-Pikus sign bug fix session) but the docstring was not updated, creating a code-vs-documentation inconsistency. This is a common pattern: focused bug fixes address runtime behavior but leave surrounding documentation stale (session history).
- The initial `fortran_runner.py` returned `None` on non-zero returncode, causing 5 test files to silently SKIP instead of FAIL when the Fortran executable crashed. Caught in an earlier code review round and fixed by raising `RuntimeError`.
- The `strain_Du` factor-of-2 error itself caused a 203 meV VB discrepancy in the strain bandedge test. Only CB matched before the fix; after fixing, all 8 eigenvalues matched (session history).

## Solution

### Fix 1: Eigenvalue count mismatch in comparison.py (lines 58-77)

**Before:**
```python
    n = min(len(ours_meV), len(theirs_meV))
    if n_ours != n_theirs:
        print(f"  Warning: eigenvalue count mismatch ...")
    ...
    passed = all(b["passed"] for b in per_band)
```

**After:**
```python
    count_mismatch = n_ours != n_theirs
    if count_mismatch:
        print(f"  Warning: eigenvalue count mismatch ...")
    ...
    passed = all(b["passed"] for b in per_band) and not count_mismatch
```

The function still compares the overlapping subset for diagnostic value but now correctly fails the test on count mismatch.

### Fix 2: Stale strain_Du docstring in param_mapper.py (line 15)

Updated from `strain_Du = -1.5 * b_dp * 1000` to `strain_Du = -0.75 * b_dp * 1000  (meV, convention: bs = b_dp/2, Du = -(3/2)*bs)` to match the implementation at line 248. The derivation chain is included so future readers can verify the convention without re-tracing kdotpy's `layerstack.py`.

### Fix 3: Misleading test docstrings (4 files)

Each module docstring was rewritten to accurately describe what the test validates:
- `test_landau_bulk.py`: "validation against analytical formulas; kdotpy LL comparison deferred"
- `test_wire_subbands.py`: "single-code consistency validation; kdotpy wire requires hzy() integration"
- `test_gfactor_qw.py`: "Bulk g-factor validation; QW comparison deferred"
- `test_sc_qw.py`: "single-code validation; kdotpy comparison deferred"

### Fix 4: Fragile SC convergence detection (test_sc_qw.py)

Removed `or "converged" in result.stdout.lower()` fallback, keeping only exact `"SC loop converged" in result.stdout`. The Fortran output is deterministic.

### Fix 5: Tolerance docstring (test_sc_qw.py)

Updated from "< 50 meV" to "< 10 meV" to match `TOL_CB1_MEV = 10.0`.

## Why This Works

1. **Count mismatch is a signal**, not noise. When two solvers return different numbers of eigenvalues, something fundamental differs. The diagnostic comparison of the overlapping subset is still useful for debugging, but the test must fail. This follows the project's established pattern: test functions must return `(bool, data)` tuples and assertions must check meaningful ranges.

2. **Docstring derivation chains prevent drift.** Including `bs = b_dp/2, Du = -(3/2)*bs` means the formula can be re-derived rather than copied. The derivation comes from kdotpy's `layerstack.py`: `bs = -(2/3)*strain_Du`, combined with `bs = b_dp/2` from the Q_eps/vv matching condition.

3. **Accurate test descriptions prevent false confidence.** A developer seeing "cross-validation" expects two independent codes agreeing. If the test only validates against analytical formulas, labeling it honestly prevents misinterpretation.

4. **Exact string matching for deterministic outputs** is both sufficient and safer than loose substring matching. The `or "converged" in stdout.lower()` fallback would match "NOT converged" or "did NOT converge" via Boolean short-circuit.

## Prevention

- **Grep for doc references when fixing formulas.** After correcting a formula in code, search the same file for any docstring or comment referencing the old value. The param_mapper docstring could have been caught by `grep -n '1.5' param_mapper.py`.
- **Comparison functions must not silently downgrade.** Structural mismatches (different array sizes, different counts) are findings, not noise. The diagnostic comparison of overlapping data is useful but should not influence the pass/fail verdict.
- **Update docstrings when test scope changes.** If practical constraints reduce a test from "cross-code" to "single-code against analytical formulas," update the docstring immediately. A test whose description does not match its behavior is worse than no test at all.
- **Exact string matching for deterministic executables.** Prefer exact substring matches over case-insensitive or partial matches. If the output format might change, define the expected format as a constant.
- **Docstring tolerance values should reference code constants** rather than duplicating numeric values. Say `Tolerance: TOL_CB1_MEV` instead of hardcoding "50 meV" in the docstring.

## Related Issues

- [[bir-pikus-qeps-sign-bug]] — The param_mapper docstring staleness is downstream of the Bir-Pikus correction session where param_mapper.py was modified but the docstring was not updated.
- [[lecture-test-pair-review-findings]] — Established project prevention rules for test quality (pass/fail tuples, meaningful assertions) that apply directly to the compare_eigenvalues count-mismatch bug.
- GitHub Issue #14 — Parent issue for the entire kdotpy cross-validation pipeline.
