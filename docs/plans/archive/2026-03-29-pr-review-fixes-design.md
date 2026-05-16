# PR Review Fixes Design

Date: 2026-03-29
Branch: feature/codebase-improvement
PR: #2

## Overview

Fix all issues found in the comprehensive PR review of the codebase improvement branch.
Organized into 5 batches executed sequentially with a test-first approach for the critical
sign convention bug.

## Batch 0: FD Sign Verification Test

**File:** `tests/fd_sign_compare.f90` (NEW)

Standalone Fortran program that verifies the higher-order FD path produces identical
kpterms output to the order-2 path when called with FDorder=2.

### Procedure
1. Set up a 2-layer QW (GaAs/AlGaAs, N=41 points, z=[-10,10] nm)
2. Call `confinementInitialization` with `FDorder=2` → `kpterms_o2(:,:,:)
3. Call `confinementInitialization` with `FDorder=2` via higher-order path → `kpterms_ho(:,:,:)
4. Compare all 10 kpterms slices element-by-element
5. Print PASS/FAIL for each slice with max absolute difference
6. Also test FDorder=4 and FDorder=6, dumping kpterms for manual inspection

### Expected outcome
- All slices should match within `tolerance=1e-12` for order=2 comparison
- If mismatch found, the differing terms and their ratios guide the fix

### Build integration
Add to Makefile:
```makefile
test-fd-sign: tests/fd_sign_compare.f90
	$(FC) $(FFLAGS) -o $@ $< -Ibuild $(LDFLAGS)
```

## Batch A: Input Validation (input_parser.f90)

**File:** `src/io/input_parser.f90`

### A1: Add iostat to all read() calls
Add `iostat=status` parameter to every `read(data_unit, *)` call (lines 37-126).
On failure, print descriptive message with field name and `stop 1`.

Example transformation:
```fortran
! Before:
read(data_unit, *) label, cfg%waveVector
! After:
read(data_unit, *, iostat=status) label, cfg%waveVector
if (status /= 0) then
  print *, 'Error: Failed to read waveVector from input.cfg'
  stop 1
end if
```

Apply to all 12 read() calls.

### A2: fdStep >= FDorder + 1 guard
After FDorder validation (around line 63), add:
```fortran
if (cfg%confinement == 1 .and. cfg%fdStep < cfg%FDorder + 1) then
  print *, 'Error: fdStep must be >= FDorder + 1 for QW simulations'
  print *, '  fdStep=', cfg%fdStep, ' FDorder=', cfg%FDorder
  print *, '  (FDorder=', cfg%FDorder, ' requires at least ', cfg%FDorder+1, ' grid points)'
  stop 1
end if
```

### A3: fdStep >= 3 guard for QW mode
After fdStep is known and confinement==1, add:
```fortran
if (cfg%confinement == 1 .and. cfg%fdStep < 3) then
  print *, 'Error: QW mode requires fdStep >= 3, got:', cfg%fdStep
  stop 1
end if
```

### A4: z(1) /= 0 guard for external field
After z-grid construction and before external field setup:
```fortran
if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
  if (abs(cfg%z(1)) < tolerance) then
    print *, 'Error: Electric field requires z(1) /= 0 (grid origin must be nonzero)'
    print *, '  Adjust startPos/endPos so the grid does not start at z=0'
    stop 1
  end if
end if
```

## Batch B: Hamiltonian Sign Fix (hamiltonianConstructor.f90)

**File:** `src/physics/hamiltonianConstructor.f90`

### B1: Fix applyVariableCoeff sign/factor convention
Based on Batch 0 test results. Most likely fix:

Option A (if factor-of-2 is the only issue):
- Add `scale_factor` parameter to `applyVariableCoeff`
- For 2nd-derivative terms, pass `0.5_dp` to match order-2 path's `1/(2*delta^2)`
- For 1st-derivative terms, pass `0.5_dp` to match order-2 path's `1/(4*delta)`

Option B (if sign also differs):
- Change sign convention in `applyVariableCoeff` or add a `sign` parameter
- Adjust callers in `confinementInitialization` accordingly

The test from Batch 0 determines which option.

### B2: Guard z(1) in externalFieldSetup_electricField
```fortran
if (abs(z(1)) < tolerance) then
  print *, 'Error: electric field setup requires z(1) /= 0'
  stop 1
end if
```
Note: This is redundant with A4 but provides defense-in-depth if called directly.

### B3: Promote EV/EC=0 warning to error in QW mode
Change lines 112-116 from warning-only to:
```fortran
if (params(i)%EV == 0.0_dp .and. params(i)%EC == 0.0_dp) then
  print *, "ERROR: Material '", trim(material(i)), &
    & "' has EV=0 and EC=0. Band offsets are required."
  print *, "  Check that the material name is correct and EV/EC are in the database."
  stop 1
end if
```

## Batch C: gfactor Robustness (gfactor_functions.f90)

**File:** `src/physics/gfactor_functions.f90`

### C1: Fix nlayers == 2 → nlayers > 1
Replace all `nlayers == 2` conditions with `nlayers > 1` in:
- Line 205: `if (nlayers == 2 .and. .not. present(HT_csr)...)`
- Line 355: `if (nlayers == 2 .and. sparse)`
- Line 370: `if (nlayers == 2 .and. .not. sparse)`

### C2: Check mkl_sparse_z_mv return code
After line 219, add:
```fortran
if (info /= SPARSE_STATUS_SUCCESS) then
  print *, 'Error: mkl_sparse_z_mv failed, info =', info
  stop 1
end if
```

### C3: Validate direction in set_perturbation_direction
Add else clause:
```fortran
subroutine set_perturbation_direction(d, smallk)
  integer, intent(in) :: d
  type(wavevector), intent(inout) :: smallk

  smallk%kx = 0
  smallk%ky = 0
  smallk%kz = 0
  select case(d)
  case(1); smallk%kx = 1
  case(2); smallk%ky = 1
  case(3); smallk%kz = 1
  case default
    print *, 'Error: perturbation direction must be 1, 2, or 3, got:', d
    stop 1
  end select
end subroutine
```

### C4: Log near-zero denominator skips
In gfactorCalculation, add a skip counter and warning:
```fortran
integer :: skip_count
...
skip_count = 0
...
if (abs(denom) > tolerance) then
  tensor(ii, jj, d) = tensor(ii, jj, d) + (Pele1*Pele2 - Pele3*Pele4) / denom
else
  skip_count = skip_count + 1
  if (skip_count <= 5) then
    print *, 'WARNING: near-zero denominator at n=', n, 'm=', m, 'l=', l, 'denom=', denom
  end if
end if
...
if (skip_count > 0) print *, 'Total skipped contributions:', skip_count
```
Apply to both the valence-band and conduction-band loops.

### C5: Validate dz in sigmaElem and pMatrixEleCalc
Add at top of each function:
```fortran
if (abs(dz) < tolerance) then
  print *, 'Error: sigmaElem called with dz=0. Grid spacing is zero.'
  stop 1
end if
```

## Batch D: Minor Fixes

### D1: Fix FDstencil order-10 bug
**File:** `src/math/finitedifferences.f90` line 135
Change: `vector(hf-5) = 8./25200.` → `vector(hf+5) = 8.0_dp/25200.0_dp`

### D2: Remove dead deallocation
**File:** `src/physics/hamiltonianConstructor.f90` lines 220-221
Remove:
```fortran
if (allocated(ScnDer)) deallocate(ScnDer)
if (allocated(FstDer)) deallocate(FstDer)
```
(These arrays are never allocated in the current code.)

### D3: Document 2.00231 comment
**File:** `src/apps/main_gfactor.f90` lines 238, 247, 256
Change `!+ 2.00231` to `! free-electron g-factor already included via sigma tensor`

### D4: Improve diag error in main_gfactor.f90
**File:** `src/apps/main_gfactor.f90` line 174
Replace:
```fortran
if (info /= 0) stop "error diag"
```
With:
```fortran
if (info /= 0) then
  print *, "Diagonalization error in g-factor calculation, info = ", info
  if (info < 0) print *, "Parameter ", -info, " had illegal value"
  stop 1
end if
```

## Execution Order

1. **Batch 0** — Write and run test program → determine exact sign fix needed
2. **Batch B1** — Apply sign fix based on test results
3. **Batch A** — Input validation (independent of sign fix)
4. **Batch C** — gfactor robustness (independent of sign fix)
5. **Batch D** — Minor fixes
6. **Re-run test** — Verify Batch 0 test passes after Batch B1 fix
7. **Build & smoke test** — `make all` + run against known results

## Files Changed

| File | Batches | Lines Changed (est.) |
|------|---------|---------------------|
| tests/fd_sign_compare.f90 (NEW) | 0 | ~150 |
| src/io/input_parser.f90 | A | +50 |
| src/physics/hamiltonianConstructor.f90 | B | ~20 |
| src/physics/gfactor_functions.f90 | C | ~30 |
| src/math/finitedifferences.f90 | D | 1 |
| src/apps/main_gfactor.f90 | D | ~10 |
| Makefile | 0 | +3 |
