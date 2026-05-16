# Codebase Improvement Design

**Date:** 2026-03-28
**Status:** Draft
**Scope:** Bug fixes, code quality refactoring, physics enhancement

## Overview

Incremental improvement of the 8-band k.p FDM codebase across five phases.
Each phase is independently buildable and verifiable. Physics results must
match the previous phase's output (except Phase 5, which intentionally
changes accuracy).

---

## Phase 1: Bug Fixes

No structural changes. Pure correctness fixes.

### 1.1 Bulk g-factor fdStep guard

**File:** `src/apps/main_gfactor.f90`
**Problem:** When `confinement == 0` (bulk), the code allocates an
`8*fdStep x 8*fdStep` matrix but only fills an 8x8 sub-block. The
check at line 173 (`if (evnum /= N) stop`) crashes with a cryptic
message unless `fdStep == 1`.
**Fix:** After reading `fdStep` from input, add:

```fortran
if (confinement == 0 .and. fdStep /= 1) then
  print *, 'Warning: bulk mode requires fdStep=1. Forcing fdStep=1.'
  fdStep = 1
end if
```

### 1.2 Division-by-zero guard in g-factor tensor

**File:** `src/physics/gfactor_functions.f90`
**Problem:** Line 555 computes a tensor element with denominator
`(cb_value(n) - vb_value(l)) + (cb_value(m) - vb_value(l))`. When
conduction and valence band energies are degenerate (narrow-gap
materials, high-symmetry points), this is zero, producing NaN.
The conduction-band loop at line 578 has a guard but the
valence-band loop does not.
**Fix:** Wrap the tensor accumulation in a threshold check:

```fortran
denom = (cb_value(n) - vb_value(l)) + (cb_value(m) - vb_value(l))
if (abs(denom) > tolerance) then
  tensor(ii, jj, d) = tensor(ii, jj, d) + (Pele1*Pele2 - Pele3*Pele4) / denom
end if
```

### 1.3 Simpson's rule odd-point validation

**File:** `src/core/utils.f90`
**Problem:** `simpson` function assumes odd `size(f)` but does not
validate this. Even-size input silently produces incorrect results.
**Fix:** Add at subroutine entry:

```fortran
if (mod(size(f), 2) == 0) then
  print *, 'Error: Simpson integration requires odd number of points.'
  stop
end if
```

### 1.4 Missing EV/EC for alloy materials

**File:** `src/core/parameters.f90`
**Problem:** 12+ alloy materials (GaAsW, Al20Ga80As, Al30Ga70As,
Ga47In53AsW, Al47In53AsW, InAsW, AlAsW, AlSbW, InSbW, InPW, GaP, AlP)
have `EV == 0` and `EC == 0` by default because these fields are never
set. Used in quantum well mode, they produce incorrect band offsets.
**Fix:** Add `EV` and `EC` values for all affected materials using
Vegard's law (linear interpolation of binary endpoints) from published
band offset data. Add a runtime warning when any material with
`EV == 0 .and. EC == 0` is used in confinement mode.

Materials and interpolation strategy:
- Al20Ga80As: 80% GaAs + 20% AlAs band offsets
- Al30Ga70As: 70% GaAs + 30% AlAs band offsets
- Ga47In53AsW: 47% GaAs + 53% InAs band offsets
- Al47In53AsW: 47% AlAs + 53% InAs band offsets
- GaAsW, InAsW, AlAsW, AlSbW, InSbW, InPW: same as their non-W
  counterparts (W variants are "well" versions with same parameters)
- GaP, AlP, InP: add from published VBO references (e.g., Vurgaftman
  et al., J. Appl. Phys. 89, 5815 (2001))

### 1.5 Fix kronij implementation

**File:** `src/core/defs.f90`
**Problem:** Uses fragile floating-point arithmetic:
```fortran
kronij = int(float((i+j)-abs(i-j)) / float((i+j)+abs(i-j)))
```
**Fix:** Replace with:
```fortran
kronij = merge(1, 0, i == j)
```

### 1.6 Fix subroutine name typo

**File:** `src/physics/hamiltonianConstructor.f90`
**Problem:** Subroutine named `externalFiledSetup_electricField`
should be `externalFieldSetup_electricField`.
**Fix:** Rename subroutine and update call sites in `main.f90` and
`main_gfactor.f90`.

### Verification

Build both executables. Run with `input_bulk.cfg` and a QW config.
Record eigenvalues and g-factors as baseline for subsequent phases.

---

## Phase 2: Dead Code Cleanup

Pure removal. No physics changes.

### 2.1 Remove unused types and functions

| Item | File |
|------|------|
| `hilbertspace` type | `src/core/defs.f90` |
| `gcd_rec` function | `src/core/defs.f90` |
| `KPRODUCT` subroutine | `src/math/finitedifferences.f90` |
| `dnscsr` subroutine (real-valued) | `src/core/utils.f90` |
| `pm(128)` array | `src/physics/gfactor_functions.f90` |
| `type` unused parameter | `src/math/finitedifferences.f90` `FDstencil` |

### 2.2 Remove commented-out debug code

| File | Approximate lines removed |
|------|--------------------------|
| `src/physics/hamiltonianConstructor.f90` | ~40 lines of `write(102,...)` |
| `src/physics/gfactor_functions.f90` | ~50 lines of Hermiticity checks, derivative code |
| `src/apps/main_gfactor.f90` | ~30 lines of symmetry checks |

### 2.3 Remove unused variables

| Item | File |
|------|------|
| `bshift` array declaration and allocation | `src/apps/main.f90`, `main_gfactor.f90` |
| `bshift` argument | `src/physics/hamiltonianConstructor.f90` `confinementInitialization` |
| `ifail` deallocation | `src/apps/main_gfactor.f90` |
| Commented-out `derivative`/`der`/`der_csr` allocations | `src/physics/gfactor_functions.f90` |

### 2.4 Fix style inconsistencies in touched files

- Replace `dcmplx(a, b)` with `cmplx(a, b, kind=dp)` throughout
- Ensure numeric literals use `_dp` kind suffix where appropriate
- Use consistent `implicit none` at module level only (remove redundant
  subroutine-level declarations where module already has it)

### Verification

Build and run. Output must be byte-identical to Phase 1.

---

## Phase 3: Extract Shared Input/Setup Module

### 3.1 New derived type: `simulation_config`

**File:** `src/core/defs.f90`

```fortran
type :: simulation_config
  integer :: confinement = 0
  integer :: fdStep = 1
  integer :: numLayers = 1
  integer :: numcb = 2
  integer :: numvb = 6
  integer :: evnum = 8
  integer :: waveVectorStep = 100
  integer :: ExternalField = 0
  character(len=10) :: waveVector = 'k0'
  real(dp) :: waveVectorMax = 0.0_dp
  real(dp), allocatable :: startPos(:)
  real(dp), allocatable :: endPos(:)
  integer, allocatable :: intStartPos(:)
  integer, allocatable :: intEndPos(:)
  character(len=30), allocatable :: materialN(:)
  type(paramStruct), allocatable :: params(:)
  real(dp), allocatable :: EFParams(:)
end type
```

### 3.2 New module: `input_parser`

**File:** `src/io/input_parser.f90`

Single public subroutine:

```fortran
subroutine read_and_setup(cfg, fdStep_actual)
  type(simulation_config), intent(out) :: cfg
  integer, intent(out) :: fdStep_actual  ! for array sizing
```

Responsibilities:
- Open and parse `input.cfg`
- Populate `cfg%materialN` and call `paramDatabase` for each material
- Compute `cfg%intStartPos`/`cfg%intEndPos` from z-positions
- Call `confinementInitialization` and `externalFieldSetup`
- Handle bulk fdStep guard (from Phase 1.1)

### 3.3 Refactored mains

Both `main.f90` and `main_gfactor.f90` shrink to:

1. `call read_and_setup(cfg, N)` — shared setup
2. Allocate Hamiltonian `HT(N, N)` using `cfg%evnum`, `cfg%fdStep`
3. Build Hamiltonian (`ZB8bandBulk` or `ZB8bandQW`)
4. Diagonalize (program-specific: `zheevx` vs `zheev`/`zheevd`)
5. Post-process (band structure sweep vs g-factor tensor)
6. Write output
7. Deallocate

Expected line counts: ~100-150 lines each (down from 362/416).

### Module dependency update

```
defs.f90
  <- parameters.f90
  <- mkl_spblas.f90
       <- utils.f90
            <- finitedifferences.f90
                 <- hamiltonianConstructor.f90
                      <- gfactor_functions.f90
  <- outputFunctions.f90
  <- input_parser.f90          (NEW: uses defs, parameters, hamiltonianConstructor)
```

### Verification

Build and run. Output must match Phase 2 output exactly (same
eigenvalues, same g-factors to machine precision).

---

## Phase 4: Refactor Large Functions

### 4.1 `gfactorCalculation` decomposition

**File:** `src/physics/gfactor_functions.f90`

Current: single 403-line subroutine.
Target: main subroutine ~150 lines + 2 helpers.

**Helper 1: `compute_pele`**

```fortran
function compute_pele(band_idx, direction, HT, eigenvectors, &
                      N, nlayers, sparse, sparseHT) result(Pele)
  integer, intent(in) :: band_idx, direction, N, nlayers
  complex(dp), intent(in) :: HT(N, N), eigenvectors(N, N)
  logical, intent(in) :: sparse
  ! sparse handle type intent(in)
  real(dp) :: Pele
```

Encapsulates the `if (nlayers==1) ... else if (sparse) ... else ...`
branching pattern currently copy-pasted 4 times.

**Helper 2: `set_perturbation_direction`**

```fortran
subroutine set_perturbation_direction(d, smallk)
  integer, intent(in) :: d
  type(wavevector), intent(out) :: smallk
```

Sets `smallk` components based on direction index. Currently
copy-pasted 3+ times.

### 4.2 `confinementInitialization` decomposition

**File:** `src/physics/hamiltonianConstructor.f90`

Current: 206 lines with 5 identical FD matrix construction blocks.
Target: ~80 lines + 1 helper.

**Helper: `build_kpterm_block`**

```fortran
subroutine build_kpterm_block(kpterms, profile, central, forward, &
                               backward, fdStep, term_idx, scale_factor)
  real(dp), intent(inout) :: kpterms(:,:,:)
  real(dp), intent(in) :: profile(:,:)
  real(dp), intent(in) :: central(:,:), forward(:,:), backward(:,:)
  integer, intent(in) :: fdStep, term_idx
  real(dp), intent(in) :: scale_factor
```

Each of the 5 blocks calls this helper with different `term_idx` and
`scale_factor`.

### 4.3 `insertCOO_cmplx` optimization

**File:** `src/core/utils.f90`

Current: O(n^4) total (linear scan per insertion for duplicates).

**New approach:** Two-phase build:

1. `insertCOO_cmplx` appends entries without duplicate checking
   (O(1) per insertion)
2. New `finalizeCOO_cmplx` sorts by (row, col) and merges duplicates
   (O(n log n))

Update `dnscsr_z_mkl` to call `finalizeCOO_cmplx` after all insertions.

### Verification

Build and run. Output must match Phase 3 to machine precision.
Performance improvement measurable for large QW systems (fdStep > 100).

---

## Phase 5: Higher-Order Finite Differences

This is the only phase that intentionally changes physics results.
Higher-order stencils improve accuracy, so results will differ from
Phase 4 but should converge as order increases.

### 5.1 New input parameter

**File:** `input.cfg`

```
FDorder      4        ! Optional. Finite difference order: 2, 4, 6, 8, 10.
                       ! Default: 2 (backward compatible)
```

Parsed in `input_parser.f90` (Phase 3 module).

### 5.2 Extend forward/backward stencils

**File:** `src/math/finitedifferences.f90`

Currently `FDstencil` has forward/backward stencils only for 2nd-order.
Add forward and backward stencil coefficients for orders 4, 6, 8, 10.

Source: Fornberg, B. "Generation of Finite Difference Formulas on
Arbitrarily Spaced Grids", Math. Comp. 51, 699-706 (1988).

Standard asymmetric stencils (uniform grid, 2nd derivative):

**4th-order forward (3-point one-sided):**
```
stencil = [2, -5, 4, -1] / dx^2   (indices 0..3)
```

**4th-order backward (3-point one-sided):**
```
stencil = [-1, 4, -5, 2] / dx^2   (indices -3..0)
```

Higher orders follow the same pattern from Fornberg's algorithm.

Implementation: extend `FDstencil` select-case blocks with the new
coefficients for `stencil_type == 'forward'` and `stencil_type == 'backward'`.

### 5.3 Modify `build_kpterm_block` helper

**File:** `src/physics/hamiltonianConstructor.f90`

The helper (from Phase 4.2) currently receives pre-computed
forward/central/backward matrices. Modify `confinementInitialization`
to build these matrices using `FDstencil` at the user-specified order
instead of hardcoding 2nd-order.

Change:
```fortran
! Old: hardcoded 2nd-order
call FDmatrixDense('central', 2, fdStep, fdStep, fdStep, dz, cCentral)
call FDmatrixDense('forward',  2, fdStep, fdStep, fdStep, dz, cForward)
call FDmatrixDense('backward', 2, fdStep, fdStep, fdStep, dz, cBackward)
```
To:
```fortran
! New: user-specified order
call FDmatrixDense('central', FDorder, fdStep, fdStep, fdStep, dz, cCentral)
call FDmatrixDense('forward',  FDorder, fdStep, fdStep, fdStep, dz, cForward)
call FDmatrixDense('backward', FDorder, fdStep, fdStep, fdStep, dz, cBackward)
```

### 5.4 Optional: interface smoothing

New optional input parameter:
```
interfaceSmoothing  0    ! 0 = abrupt (default), 1 = linear average at boundaries
```

When enabled, at each interface between layer i and layer i+1,
replace the abrupt parameter jump with a linear average over the
adjacent grid points. This reduces spurious interface states.

Implementation in `confinementInitialization`: after computing the
step-function profile, apply a single-pass smoothing at interfaces
identified by `intStartPos`/`intEndPos`.

### 5.5 Verification

1. Run with `FDorder=2` — must match Phase 4 output exactly
2. Run with `FDorder=4`, `FDorder=6` — verify convergence:
   - Compare eigenvalue differences between orders
   - Higher order should give smoother, more converged results
   - Key test: GaAs/AlGaAs QW ground state energy should converge
     to within ~1 meV by order 4-6
3. Cross-reference against published k.p results:
   - GaAs/Al0.3Ga0.7As single QW (well-studied benchmark)
   - InAs/GaSb type-II QW (tests interband coupling)

---

## Risk Assessment

| Phase | Risk | Mitigation |
|-------|------|------------|
| 1 | Low — targeted fixes | Each fix is isolated and testable |
| 2 | Low — pure removal | If build breaks, the removed code was needed |
| 3 | Medium — structural refactor | Byte-comparison of output vs Phase 2 |
| 4 | Medium — function extraction | Byte-comparison of output vs Phase 3 |
| 5 | High — physics changes | Convergence testing, published benchmarks |

## File Change Summary

| Phase | New files | Modified files |
|-------|-----------|---------------|
| 1 | 0 | `defs.f90`, `parameters.f90`, `utils.f90`, `hamiltonianConstructor.f90`, `gfactor_functions.f90`, `main_gfactor.f90` |
| 2 | 0 | `defs.f90`, `utils.f90`, `finitedifferences.f90`, `hamiltonianConstructor.f90`, `gfactor_functions.f90`, `main.f90`, `main_gfactor.f90` |
| 3 | 1 (`input_parser.f90`) | `defs.f90`, `Makefile`, `main.f90`, `main_gfactor.f90` |
| 4 | 0 | `utils.f90`, `hamiltonianConstructor.f90`, `gfactor_functions.f90` |
| 5 | 0 | `finitedifferences.f90`, `hamiltonianConstructor.f90`, `input_parser.f90` |
