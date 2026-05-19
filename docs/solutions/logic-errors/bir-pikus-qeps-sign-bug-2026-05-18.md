---
title: "Bir-Pikus Q_eps sign error: wrong HH/LH ordering under compressive strain"
date: "2026-05-18"
category: logic-errors
module: src/physics/strain_solver.f90
problem_type: logic_error
component: service_object
severity: critical
symptoms:
  - "Under compressive strain (InAs/GaAs), HH shifted down instead of up; LH shifted up instead of down"
  - "Cross-code validation against kdotpy showed VB eigenvalues diverging by hundreds of meV while CB agreed within 0.001 meV"
  - "Strain verification test band ordering wrong: HH at eigenvalue index 2 instead of 4"
  - "No runtime error -- silently incorrect band structure in all strained simulations"
root_cause: logic_error
resolution_type: code_fix
tags: [bir-pikus, strain, sign-convention, cross-code-validation, kdotpy, deformation-potential]
related_components:
  - src/physics/hamiltonianConstructor.f90
  - tests/integration/star_helpers.py
  - tests/integration/verify_strain_rung5_bulk.py
---

# Bir-Pikus Q_eps sign error: wrong HH/LH ordering under compressive strain

## Problem

The `compute_bp_scalar` function in `strain_solver.f90` had a wrong sign in the Q_eps shear strain formula: `Q_eps = +(b_dp/2) * (eps_zz - 0.5*(eps_yy + eps_xx))` instead of `Q_eps = -(b_dp/2) * (...)`. This inverted the HH/LH band ordering under compressive strain -- HH shifted down when it should shift up. The error was silent (no crashes) and affected all strained simulations.

This was the third sign fix in the same Bir-Pikus formula chain across three separate commits:
1. Commit `dfed60d` (Apr 18): pseudomorphic strain convention `eps = (a_sub - a_well)/a_well`
2. Commit `fe1e243` (Apr 20): VB diagonal from `P_eps + Q_eps` to `-P_eps + Q_eps`
3. Current branch: Q_eps itself from `+b_dp/2` to `-b_dp/2`

## Symptoms

Under compressive strain (InAs pseudomorphic on GaAs, ~6.7% mismatch):

- HH shifted down by ~187 meV when it should shift up by ~65 meV
- LH shifted up when it should shift down
- CB shift was correct (~310 meV in both codes) -- unaffected because `delta_Ec = ac * Tr(eps)` has no P_eps or Q_eps
- HH-LH splitting was inverted: Fortran gave ~0.083 eV, kdotpy gave ~0.162 eV
- Strain verification test mapped eigenvalues to wrong bands (HH at index 2 instead of 4)

## What Didn't Work

- **Earlier sign fixes were incomplete.** The April 18 refactor corrected `P_eps` diagonal signs but left `Q_eps` with the wrong sign on `b_dp`. Each fix was necessary but insufficient alone. (session history)
- **Unit tests encoded the same wrong formula.** `test_strain_solver.pf` computed expected Q_eps using `+b_dp/2`, so tests passed despite the bug. Cross-code validation was essential to catch this. (session history)
- **Comments went stale.** `hamiltonianConstructor.f90:806` documented `Q_eps = b_dp/2 * (...)` (positive, no minus), misleading future readers even after the code delegated to `compute_bp_scalar`.

## Solution

### Core fix in strain_solver.f90:841

```fortran
! Before (WRONG):
Q_eps = params%b_dp * 0.5_dp * (eps_zz - 0.5_dp * (eps_yy + eps_xx))

! After (CORRECT):
Q_eps = -params%b_dp * 0.5_dp * (eps_zz - 0.5_dp * (eps_yy + eps_xx))
```

### Coordinated updates across 5 additional files

**Unit tests** (`tests/unit/test_strain_solver.pf`, lines ~452 and ~539):
```fortran
! Both test subroutines updated to match corrected formula:
Q_eps = -params(1)%b_dp * 0.5_dp * (eps_zz_val - 0.5_dp * 2.0_dp * eps_xx_val)
```

**Python analytical reference** (`tests/integration/star_helpers.py:315`):
```python
Q_eps = -b_dp / 2.0 * (eps_zz - 0.5 * (eps_xx + eps_yy))
```

**Eigenvalue index mapping** (`tests/integration/verify_strain_rung5_bulk.py:112-117`):

Under compressive strain with corrected Q_eps, ascending order becomes LHSO_low < LHSO_high < HH < CB (HH shifts above LH-SO mixed state):
```python
expected_map = {
    "LHSO_low":  (0, bp["E_LHSO_low"]),
    "LHSO_high": (2, bp["E_LHSO_high"]),
    "HH":        (4, bp["E_HH"]),
    "CB":        (6, bp["E_CB"]),
}
```

**Regression reference** (`tests/integration/verify_star_inas_gaas_qw.py:126`):
```python
HHLH_REF = 0.1623  # eV (updated from 0.0567 with corrected Q_eps)
```

**Stale comment** (`src/physics/hamiltonianConstructor.f90:806`):
```fortran
! Updated to match compute_bp_scalar:
!   Q_eps = -(b_dp/2) * (eps_zz - 0.5*(eps_xx+eps_yy)) (tetragonal shear)
```

## Why This Works

The standard Chuang/Winkler Bir-Pikus convention defines:

- `P_eps = -av * Tr(eps)` (hydrostatic, pushes all VB edges)
- `Q_eps = -(b/2) * (eps_zz - 0.5*(eps_xx + eps_yy))` (tetragonal shear)
- `delta_EHH = -P_eps + Q_eps`, `delta_ELH = -P_eps - Q_eps`, `delta_ESO = -P_eps`

Under compressive strain (InAs on GaAs: `Tr(eps) < 0`, `b_dp < 0`):
- `Q_eps > 0` (the minus on b_dp and negative b_dp give a positive result)
- HH shifts up (`+Q_eps`), LH shifts down (`-Q_eps`)
- Without the minus sign on b_dp, Q_eps was negative, inverting HH/LH ordering

The CB was unaffected because `delta_Ec = ac * Tr(eps)` has no Q_eps dependence, explaining why kdotpy cross-validation showed perfect CB agreement (< 0.001 meV) but VB diverged by hundreds of meV.

## Prevention

- **CLAUDE.md boundary rule**: The Bir-Pikus convention is documented as a NEVER-change boundary with the exact formula and expected physics (under compressive strain: HH up, LH down). Single source of truth: `compute_bp_scalar` in `strain_solver.f90`.
- **Cross-code validation pipeline**: `validation/strain/test_strain_bandedge.py` compares strained eigenvalues against kdotpy, catching sign errors that unit tests miss when both encode the same wrong formula.
- **Physical sign test**: Unit test `test_compute_bp_scalar_gaas_biaxial` verifies exact numerical values; `test_strain_sign_inas_on_gaas` verifies compressive eps_xx is negative.
- **Six-file coordination**: Any change to `compute_bp_scalar` requires coordinated updates in: strain_solver.f90, test_strain_solver.pf, star_helpers.py, verify_strain_rung5_bulk.py, verify_star_inas_gaas_qw.py, and hamiltonianConstructor.f90 comments.

## Related Issues

- [Archived strain sign fix plan](../../plans/archive/2026-04-20-strain-sign-and-buffer-fix.md) -- earlier P_eps diagonal sign fix (phase 2 of the three-stage correction)
- [8-band verification ladder](../best-practices/8band-verification-ladder-2026-05-09.md) -- strain rungs 5-6 extend this framework
- [Standard-star benchmark logic errors](standard-star-benchmark-suite-logic-errors-2026-05-09.md) -- testing framework patterns for silently wrong physics
- `docs/lecture/04-strain.md` -- updated to match corrected Q_eps formula and HH/LH ordering
- [Cross-validation code review fixes](cross-validation-review-logic-errors-2026-05-19.md) -- param_mapper docstring was left stale after the Bir-Pikus correction; caught in follow-up code review
