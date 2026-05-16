# FD Variable-Coefficient & Config Overlap Fix Design

**Date:** 2026-04-17
**Branch:** `feature/docs-overhaul`

## Problem

Two pre-existing bugs in the QW Hamiltonian construction:

1. **FDorder>2 variable coefficients:** `applyVariableCoeff` does naive row-scaling `kpterms(j,i) = -profile(j) * FD(j,i)`, which incorrectly discretizes `d/dz[g(z)·d/dz]` at heterointerfaces. Only FDorder=2 produces correct eigenvalues because `build_kpterm_block` implicitly averages profile values via `dgemv`.

2. **Config overlap (painter's algorithm):** `confinementInitialization_raw` assigns `profile(startPos:endPos)` sequentially — last layer wins at overlapping z-points. A symmetric `barrier(-200,200) / well(-50,50) / barrier(-200,200)` config erases the well entirely.

## Root Cause Analysis

### Bug 1: Variable Coefficients

The k.p Hamiltonian has terms like `d/dz[γ(z)·d/dz]` where γ differs between materials. The correct FD discretization requires **midpoint averaging**: `(γ_i + γ_j)/2` at the half-grid point between coupled cells i and j.

**FDorder=2 path (correct):**
- `dgemv('N', N, N, 1.0, forward, N, profile_vec, ...)` computes `offup(i) = profile(i) + profile(i+1)` — 2-point average
- `dgemv('N', N, N, 1.0, central, N, profile_vec, ...)` computes `diag(i) = profile(i-1) + profile(i) + profile(i+1)` — 3-point average
- After sign and scale (`1/(2·dz²)`), this gives the standard averaged-coefficient stencil

**FDorder>=4 path (wrong):**
- `kpterms(j,i) = -profile(j) * FD(j,i)` — every stencil weight scaled by `profile(j)` only
- At interfaces: row j uses `γ(j)` for all couplings, ignoring neighbors' different `γ`

### Bug 2: Config Overlap

Lines 148-161: sequential `profile(startPos:endPos) = ...` overwrites earlier layers. No overlap detection, no priority system, no masking.

## Fix Design

### Fix 1: Midpoint-Averaged Variable Coefficients

Replace `applyVariableCoeff` with element-wise multiply using a pre-averaged profile matrix:

```
G_avg(i,j) = profile(i)  if i == j  (diagonal: local value)
G_avg(i,j) = (profile(i) + profile(j)) / 2  if i ≠ j  (off-diagonal: midpoint average)
kpterms(:,:,term) = -G_avg .* FD
```

This matches the FDorder=2 path's implicit averaging and generalizes correctly to all orders.

**Why this is correct:** For the 2nd-derivative operator `d/dz[g·d/dz]`, the standard FD approximation at point i couples to neighbor j through the midpoint value `g_{i+1/2}`. The matrix `G_avg .* D2` captures exactly this: each off-diagonal entry `D2(i,j)` gets weighted by the average `(g_i + g_j)/2 ≈ g_{(i+j)/2}`.

For 1st-derivative terms `g(z)·d/dz`, the same averaging applies. The sign convention is preserved by the leading `-`.

### Fix 2: Mask-Based Layer Assignment

Add a boolean mask to `confinementInitialization_raw`:

```fortran
logical :: assigned(N)
assigned = .false.

do i = 1, nlayers
  do j = startPos(i), endPos(i)
    if (.not. assigned(j)) then
      profile(j,:) = params(i)%...
      kptermsProfile(j,:) = params(i)%...
      assigned(j) = .true.
    end if
  end do
end do

! Verify all points assigned
if (any(.not. assigned)) then
  print *, 'ERROR: grid points not covered by any material layer'
  stop 1
end if
```

First layer to claim a point wins. With natural config ordering (barrier, well, barrier), the well correctly occupies the interior.

### Fix 3: Config Validation

In `input_parser.f90`, after computing `startPos/endPos`, verify:
1. All grid points covered (no gaps)
2. Warn if layers overlap (now safe due to mask, but indicates sloppy config)

## Config Fixes

### qw_gaas_algaas_kpar.cfg (and optics variant)

**Before (broken):**
```
material1: Al30Ga70As -200 200 0   ! barrier covers everything
material2: GaAs -50 50 0           ! well overwritten by material3
material3: Al30Ga70As -200 200 0   ! barrier overwrites well
```

**After (correct):**
```
material1: Al30Ga70As -200 -50 0   ! left barrier
material2: GaAs -50 50 0           ! well (not overwritten)
material3: Al30Ga70As 50 200 0     ! right barrier
```

### qw_inas_gasb_broken_gap_kpar.cfg

Already non-overlapping. Will produce correct results once FDorder=4 variable-coefficient fix is applied.

## Physics Regeneration

After code and config fixes, re-run all three configs:
1. `qw_gaas_algaas_kpar.cfg` → new eigenvalues, E(k_parallel) dispersion
2. `qw_gaas_algaas_optics.cfg` → new optical matrix elements
3. `qw_inas_gasb_broken_gap_kpar.cfg` → new broken-gap eigenvalues, anticrossing data

Expected changes:
- GaAs/AlGaAs eigenvalues should shift significantly (FDorder=4 was giving ~1.8 eV too high)
- Optical matrix elements will use correct wavefunctions
- Broken-gap eigenvalues may shift slightly (InAsW/GaSbW interfaces also affected)

## Documentation Updates

Files requiring numerical value updates:
- `docs/lecture/02-quantum-well.md` — eigenvalue tables (A.3, B.3), optical table (A.6), config snippets, anticrossing prose (B.5)
- `docs/lecture/06-optical-properties.md` — Section 6.6.6 transition strengths, Section 6.8.1 FDorder recommendation
- `docs/lecture/11-convergence.md` — FDorder convergence data if affected

Figures requiring regeneration:
- `fig_qw_dispersion_gaas_algaas` — E(k_parallel) with new eigenvalues
- `fig_qw_optical_matrix_elements` — oscillator strengths with new data
- `fig_qw_potential_profile_gaas` — potential profile with new eigenvalues
- `fig_qw_dispersion_broken_gap` — broken-gap anticrossing with corrected data

## Files to Modify

### Code
- `src/physics/hamiltonianConstructor.f90` — `applyVariableCoeff`, `confinementInitialization_raw`
- `src/io/input_parser.f90` — overlap/gap validation

### Configs
- `tests/regression/configs/qw_gaas_algaas_kpar.cfg`
- `tests/regression/configs/qw_gaas_algaas_optics.cfg`

### Docs + Figures
- `docs/lecture/02-quantum-well.md`
- `docs/lecture/06-optical-properties.md`
- `scripts/plotting/generate_all_figures.py` (no code changes, but re-run for new data)
- `docs/figures/` — regenerated PNGs

### Tests
- `tests/unit/test_hamiltonian.pf` — add test for variable-coefficient FDorder>2
- Regression golden data may need updating (eigenvalues change)
