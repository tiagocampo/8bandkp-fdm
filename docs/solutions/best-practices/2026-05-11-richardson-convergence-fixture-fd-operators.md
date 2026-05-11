---
title: Richardson convergence fixture methodology for FD operator validation
date: "2026-05-11"
category: best-practices
module: tests/unit/test_fd_convergence.pf
problem_type: best_practice
component: testing_framework
severity: medium
applies_when:
  - "Writing convergence tests for finite-difference stencil operators"
  - "Verifying FD coefficient correctness at orders 6-10 where round-off contaminates fine meshes"
  - "Extending interpolation or derivative stencils via Vandermonde system solving"
  - "Measuring convergence rates for operators with different boundary vs interior accuracy"
tags: [finite-differences, convergence-testing, vandermonde-solver, interpolation, round-off-handling]
---

# Richardson convergence fixture methodology for FD operator validation

## Context

Hand-computed FD coefficients have no runtime validation -- bugs go undetected until physics results look wrong. Naive convergence-rate measurement (take rate from finest two levels) fails at high orders because round-off dominates the error budget, producing rates far below the theoretical order. A robust methodology is needed to validate FD operators in isolation.

## Guidance

### Max-rate convergence strategy

Compute convergence rate for every consecutive refinement pair across 6 levels (`N = [21, 41, 81, 161, 321, 641]`), then take the **maximum** observed rate. This captures the true asymptotic rate from the "sweet spot" where truncation dominates but round-off has not yet contaminated. Round-off causes the observed rate to drop, never to increase, so the maximum is the best estimator.

Supplement with a **finest-pair rate check** for low orders (2, 4) where round-off is not an issue -- assert that the finest-pair rate exceeds a relaxed threshold (e.g., `> expected_order - 0.5`).

### Interior-only max-norm

Measure error only at interior points `i = half_bw+1 ... N-half_bw`, excluding boundary stencil distortion. This isolates central stencil accuracy from one-sided boundary effects. For half-point interpolation, further restrict to `i = half_bw+1 ... N-half_bw-1`.

### Order-adaptive tolerance

Widen tolerance at higher orders where round-off degrades convergence sooner:

```fortran
select case(order)
case(2:4);  tol = 0.05_dp
case(6:8);  tol = 0.2_dp
case(10);   tol = 0.5_dp
end select
```

### Variable-coefficient via product rule

Test `d/dz[a(z) * du/dz]` as `a * D2*u + D1*a * D1*u` rather than using staggered operators with half-point interpolation. This separates D1 and D2 accuracy and avoids boundary artifacts from the staggered approach.

### Legacy API sign convention

The legacy `FDstencil(order=1) / FDmatrixDense` Toeplitz construction produces `-f'(x)` due to antisymmetric stencil sign convention. Negate the result before comparing to the exact derivative.

### Interpolation convergence rate

Vandermonde-derived interpolation with FDorder+1 stencil points gives convergence at FDorder+1. Lagrange order-4 (4-point) gives rate 4, not 5.

## Why This Matters

Without convergence tests, coefficient bugs like the wrong D1 order-10 values (off by 2x at dominant positions) silently corrupt simulation results. The fixture catches these at the operator level before they propagate to physics. The methodology is robust against the main failure mode (round-off contamination at high orders) that defeats naive convergence measurement.

## When to Apply

- Adding new FD orders or stencil types
- Modifying `FDcentralCoeffs1st`, `FDcentralCoeffs2nd`, or any coefficient-generation routine
- Extending `interpolateToHalfPoints` or adding new interpolation schemes
- Verifying Vandermonde-derived coefficients against analytical values

## Examples

The 27-test fixture in `tests/unit/test_fd_convergence.pf` covers:

- Constant-coefficient D1 and D2 at orders 2, 4, 6, 8, 10
- Half-point interpolation at orders 2, 4, 6, 8, 10
- Variable-coefficient D2 via product rule at orders 2, 4, 6, 8, 10
- Legacy API (FDstencil/FDmatrixDense) D2 at orders 2-10 and D1 at orders 2, 8

Direct coefficient-value test in `tests/unit/test_finitedifferences.pf` (`test_d1_order10_coefficients`) validates all 11 D1 order-10 coefficients against analytical fractions.

## Related

- `src/math/finitedifferences.f90` -- `vandermonde_interp` for deriving interpolation stencils at arbitrary offsets
- `docs/lecture/11-convergence.md` -- Richardson extrapolation mathematical foundation (section 11.3.3)
- `docs/solutions/logic-errors/2026-05-11-wrong-d1-order10-fd-coefficients.md` -- the bug this methodology caught
