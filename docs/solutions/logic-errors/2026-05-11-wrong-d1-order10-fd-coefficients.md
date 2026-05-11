---
title: Wrong 1st-derivative central coefficients at FDorder=10
date: "2026-05-11"
category: logic-errors
module: src/math/finitedifferences.f90
problem_type: logic_error
component: testing_framework
severity: high
symptoms:
  - "Order-10 1st-derivative convergence rate far below theoretical value 10"
  - "FDcentralCoeffs1st case(10) positions hf-1/hf+1 were 5/3 instead of 5/6 (2x error)"
  - "No runtime error -- silently wrong derivative operator"
root_cause: logic_error
resolution_type: code_fix
tags: [finite-differences, stencil-coefficients, convergence-testing, fortran]
---

# Wrong 1st-derivative central coefficients at FDorder=10

## Problem

Three coefficient pairs in `FDcentralCoeffs1st` case(10) (`src/math/finitedifferences.f90:380-391`) were wrong, causing `sum(c_k * k) = 163/63` instead of the analytically required value 1.0. Any simulation using FDorder=10 for the 1st-derivative operator (staggered D1, variable-coefficient D2, velocity matrices) produced silently incorrect results.

## Symptoms

- No crash or runtime error -- wrong coefficients produce a valid-looking but incorrect operator
- Richardson convergence fixture showed observed D1 rate far below 10
- The consistency condition `sum(c_k * k) = 1` fails with the wrong values

## What Didn't Work

- **Median-rate convergence**: At orders 8-10, the finest refinement levels show round-off contamination, so the median picks up degraded rates rather than the true asymptotic rate (session history)
- **pFUnit `@assertTrue` with `&` continuation**: pFUnit preprocessor doesn't handle Fortran continuation lines inside macro calls
- **Staggered D1 + half-point interpolation for variable-coeff tests**: Boundary artifacts from staggered operators masked the real convergence behavior

## Solution

Fixed six coefficient values in `FDcentralCoeffs1st`, case(10):

| Position | Wrong | Correct | Factor |
|----------|-------|---------|--------|
| hf-4 | 5/1008 | 5/504 | 2x |
| hf-3 | 5/126 | 5/84 | 1.5x |
| hf-1 | 5/3 | 5/6 | 2x |
| hf+1 | 5/3 | 5/6 | 2x |
| hf+3 | 5/126 | 5/84 | 1.5x |
| hf+4 | 5/1008 | 5/504 | 2x |

Added direct coefficient-value test in `tests/unit/test_finitedifferences.pf` (`test_d1_order10_coefficients`) that verifies all 11 values against analytical fractions.

## Why This Works

The 1st-derivative central stencil must satisfy the Vandermonde system `sum c_k * k^m = delta_{m,1}` for m=0..10. The m=1 constraint requires `sum c_k * k = 1`. The nearest-neighbor coefficient (offset +/-1) was 5/3 instead of 5/6, contributing 10/3 to the sum when it should contribute 5/3. With corrected values from the standard 11-point central stencil (Fornberg 1988), all moments vanish except m=1, giving true 10th-order accuracy.

## Prevention

- 27 Richardson convergence tests in `tests/unit/test_fd_convergence.pf` validate D1/D2 at all orders (2,4,6,8,10) -- any coefficient regression is caught by CI
- Direct coefficient-value test `test_d1_order10_coefficients` checks all 11 values against analytical fractions -- catches typos more precisely than convergence rate alone
- For future FD operators, verify hardcoded coefficients against the Vandermonde solver output before committing

## Related Issues

- `docs/plans/2026-03-29-pr-review-fixes-design.md` section "D1: Fix FDstencil order-10 bug" identified the same class of bug
- The bug was previously identified but the coefficient-value regression test was missing until the Richardson convergence fixture work
