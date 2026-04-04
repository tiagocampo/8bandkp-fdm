# PR Review Fixes + Phase 1.5 LAPACK Interface Design

**Date**: 2026-04-04
**Author**: Tiago de Campos
**Status**: In Progress (3/4 commits done)
**Branch**: `feature/quantum-wire-2d`
**PR**: #8

## Motivation

Comprehensive PR review of the quantum wire branch identified 6 critical issues, 10 important issues, and 4 test coverage gaps. Combined with the outstanding Phase 1.5 (centralized LAPACK interface) from `docs/plans/2026-04-03-linalg-backend-design.md`, this plan addresses all findings in a structured sequence.

## Wire Geometry Convention

The quantum wire is confined in the **xy plane**, with **z as the free direction** (wire axis along z, k_z sweep). The 2D confinement grid uses:
- `grid%dx` — physical x-spacing of confinement plane
- `grid%dy` — physical y-spacing of confinement plane

## Sequencing

Phase 1.5 is done first because it touches the same files as several review fixes (`main.f90`, `eigensolver.f90`, `sc_loop.f90`). Doing it first means the review fixes modify clean `use linalg, only:` imports instead of `external ::` declarations.

## Commit Strategy

Four commits:

1. Phase 1.5: centralized LAPACK interface module
2. Critical fixes (C1-C6)
3. Important fixes (I1-I8, I10)
4. Test coverage additions (T1-T4)

---

## Commit 1: `refactor: add centralized LAPACK/FEAST interface module (Phase 1.5)`

### 1.1 Create `src/math/linalg.f90`

New module with explicit `interface` blocks for all LAPACK/FEAST routines currently declared as `external`:

```fortran
module linalg
  use definitions, only: dp
  implicit none
  private

  ! Standard LAPACK
  public :: zheevx
  public :: dgesv
  public :: ilaenv
  public :: dlamch

  ! MKL-specific
  public :: mkl_set_num_threads_local

  ! MKL FEAST (guarded)
#ifdef USE_MKL_FEAST
  public :: feastinit
  public :: zfeast_hcsrev
#endif

contains
  ! (empty — this module provides only interface blocks)
end module
```

Each routine needs an explicit interface block matching its LAPACK signature. The full signatures:

**`zheevx`** — 24 arguments (standard LAPACK hermitian eigensolver):
```fortran
interface
  subroutine zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
                    abstol, m, w, z, ldz, work, lwork, rwork, iwork, &
                    ifail, info)
    use definitions, only: dp
    character(len=1), intent(in) :: jobz, range, uplo
    integer, intent(in) :: n, lda, il, iu, ldz
    real(kind=dp), intent(in) :: vl, vu, abstol
    complex(kind=dp), intent(inout) :: a(lda, *)
    real(kind=dp), intent(out) :: w(*)
    complex(kind=dp), intent(out) :: z(ldz, *)
    complex(kind=dp), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    real(kind=dp), intent(out) :: rwork(*)
    integer, intent(out) :: iwork(*), ifail(*)
    integer, intent(out) :: m, info
  end subroutine
end interface
```

**`dgesv`** — 8 arguments (LAPACK general linear system):
```fortran
interface
  subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    use definitions, only: dp
    integer, intent(in) :: n, nrhs, lda, ldb
    real(kind=dp), intent(inout) :: a(lda, *), b(ldb, *)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
  end subroutine
end interface
```

**`ilaenv`** — 6 arguments, returns integer (LAPACK block size query):
```fortran
interface
  function ilaenv(ispec, name, opts, n1, n2, n3, n4) result(value)
    character(len=*), intent(in) :: name, opts
    integer, intent(in) :: ispec, n1, n2, n3, n4
    integer :: value
  end function
end interface
```

**`dlamch`** — 2 arguments, returns double (machine constants):
```fortran
interface
  function dlamch(cmach) result(value)
    use definitions, only: dp
    character(len=1), intent(in) :: cmach
    real(kind=dp) :: value
  end function
end interface
```

**`mkl_set_num_threads_local`** — returns previous thread count:
```fortran
interface
  function mkl_set_num_threads_local(nt) result(previous)
    integer, intent(in) :: nt
    integer :: previous
  end function
end interface
```

**`feastinit`** and **`zfeast_hcsrev`** — guarded with `#ifdef USE_MKL_FEAST`:
```fortran
#ifdef USE_MKL_FEAST
interface
  subroutine feastinit(fpm)
    integer, intent(inout) :: fpm(128)
  end subroutine
end interface

interface
  subroutine zfeast_hcsrev(uplo, n, a, ia, ja, fpm, epsout, loop, &
                           emin, emax, m0, e, x, m, res, info)
    use definitions, only: dp
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: n, m0
    complex(kind=dp), intent(in) :: a(*)
    integer, intent(in) :: ia(*), ja(*)
    integer, intent(inout) :: fpm(128)
    real(kind=dp), intent(out) :: epsout
    integer, intent(out) :: loop, m, info
    real(kind=dp), intent(in) :: emin, emax
    real(kind=dp), intent(inout) :: e(m0), res(m0)
    complex(kind=dp), intent(inout) :: x(n, m0)
  end subroutine
end interface
#endif
```

### 1.2 Delete `src/math/feast_call.f90`

The non-module wrappers existed to work around gfortran allocatable descriptor issues when passing arrays to `external` routines from inside module procedures. With explicit interface blocks in `linalg.f90`, gfortran gets the correct type signatures and generates proper calls without the wrapper indirection. The wrappers are no longer needed.

### 1.3 Update caller files

**`src/apps/main.f90`** (lines 36-37):
```fortran
! Before:
integer :: mkl_set_num_threads_local
external :: zheevx, mkl_set_num_threads_local

! After:
use linalg, only: zheevx, mkl_set_num_threads_local
```

**`src/physics/sc_loop.f90`** (lines 46-48):
```fortran
! Before:
integer :: ilaenv
real(kind=dp) :: dlamch
external :: ILAENV, DLAMCH, zheevx, dgesv

! After:
use linalg, only: zheevx, dgesv, ilaenv, dlamch
```

**`src/math/eigensolver.f90`** (lines 92-93, 204):
```fortran
! Before (solve_feast):
external :: feastinit_call
external :: feast_solve_hermitian_csr

! After:
use linalg, only: feastinit, zfeast_hcsrev
! (and inline the calls directly, removing the wrapper indirection)
```

```fortran
! Before (solve_dense_lapack):
external :: zheevx

! After:
use linalg, only: zheevx
```

### 1.4 Update CMake

**`src/math/CMakeLists.txt`**: Add `linalg.f90`, remove `feast_call.f90`.

### 1.5 Build and test

All 23+ existing tests must pass. The GCC 15 strict-aliasing note should be eliminated.

---

## Commit 2: `fix: address critical PR review issues (C1-C6)`

### C1: Replace `-ffast-math` with safe components

**File**: `CMakeLists.txt` line 25

**Problem**: `-ffast-math` enables `-ffinite-math-only`, which tells the compiler it can assume no NaN/Inf values exist. This allows it to optimize away guards like `if (abs(denom) > tolerance)` in `gfactor_functions.f90` (8+ locations). In a physics code with iterative solvers (FEAST, PARDISO, DIIS), NaN propagation is a real possibility these guards catch.

**Fix**: Replace `-ffast-math` with its individual safe components:
```cmake
# Before:
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=haswell ... -ffast-math ...")

# After:
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=haswell ... -fno-signed-zeros -freciprocal-math ...")
```

Kept:
- `-fno-signed-zeros`: safe, ~2% speedup
- `-freciprocal-math`: safe for physics, ~3% speedup
- `-fno-trapping-math`: safe, already implied by default in GCC 15

Dropped:
- `-ffinite-math-only`: **dangerous** — removes NaN/Inf guards
- `-funsafe-math-optimizations`: unnecessary risk

Expected performance loss: ~2-5% of the ~30% speedup from flags. Net still positive over the old `-O2`.

### C2: Add MKL threading guard to gfactor executable

**File**: `src/apps/main_gfactor.f90`

**Problem**: `main.f90` sets `mkl_set_num_threads_local(1)` inside OpenMP regions, but `main_gfactor.f90` does not. With `MKL_THREADING=intel_thread`, MKL routines in the gfactor executable can spawn threads that conflict with OpenMP or waste resources.

**Fix**: Add at program start (after `use linalg`):
```fortran
info = mkl_set_num_threads_local(1)
```

And before any OpenMP parallel region that calls MKL:
```fortran
!$omp parallel
info = mkl_set_num_threads_local(1)
!$omp do ...
```

### C3: PARDISO errors → `stop 1` instead of silent zero

**Files**: `src/physics/strain_solver.f90` lines 593-601, `src/physics/poisson.f90` lines 323-333

**Problem**: When PARDISO fails, both solvers print an error message and return with zero arrays. Downstream code proceeds with zero strain/potential, producing silently wrong physics. In batch jobs, the error message is often lost in output.

**Fix** (strain_solver.f90):
```fortran
! Before:
print *, '  Strain PDE solve failed. Setting strain to zero.'
return

! After:
print *, '  ERROR: Strain PDE solve failed (PARDISO error=', error, '). Aborting.'
stop 1
```

Same pattern for `solve_poisson_2d` in `poisson.f90`. Add PARDISO memory release (`phase=-1`) before the `stop 1` to avoid leaks.

### C4: FEAST non-convergence in wire k-sweep → `stop 1`

**File**: `src/apps/main.f90` lines 286-298 (serial k=1), lines 336-354 (OpenMP k=2..N)

**Problem**: When FEAST fails to converge, a warning is printed but wrong eigenvalues are stored in `eig_wire`. The output file has no indication which k-points failed.

**Fix** — for both serial and parallel paths:
```fortran
if (.not. eigen_res%converged) then
  print *, '  ERROR: FEAST did not converge at k-point', k
  stop 1
end if
```

Inside OpenMP parallel: use `!$omp critical` for the error print, then `stop 1`. This kills all threads, which is the correct behavior for a fatal error.

### C5: QW OpenMP diagonalization failure → `stop 1`

**File**: `src/apps/main.f90` lines 638-644

**Problem**: `zheevx` failure does `cycle`, leaving zero-initialized eigenvalues at that k-point. These appear as a band at E=0 in the output.

**Fix**:
```fortran
! Before:
print *, 'Error in diagonalization at k =', k, 'info =', info
cycle

! After:
print *, 'ERROR: diagonalization failed at k =', k, 'info =', info
stop 1
```

### C6: FEAST `info=3` → `converged=.false.`

**File**: `src/math/eigensolver.f90` line 160

**Problem**: `info=3` means the FEAST subspace was too small to contain all eigenvalues in the search window — eigenvalues are missing. Currently the code prints a warning but proceeds with partial results.

**Fix**:
```fortran
! Before:
result%converged = (info == 0) .or. (info == 2)

! After:
result%converged = (info == 0) .or. (info == 2)
! info == 3 is already converged=.false. since it's excluded above.
! But add explicit handling for clarity:
if (info == 3) then
  print *, '  WARNING: FEAST subspace too small (info=3). Missing eigenvalues.'
  print *, '  Increase feast_m0 or narrow the energy window.'
end if
```

No change to the `converged` logic — `info=3` already evaluates to `.false.` in the expression. But add the explicit warning with remediation advice.

---

## Commit 3: `fix: address important PR review issues (I1-I8, I10)`

### I1: Energy window zero-width edge case

**File**: `src/math/eigensolver.f90` lines 305-306

**Problem**: When Gershgorin bounds are near zero (e.g., bulk at k=0), the relative margin produces a zero-width window `[0, 0]`. FEAST finds nothing.

**Fix**: Add absolute margin floor:
```fortran
emin = rmin - max(0.1_dp * abs(rmin), 0.5_dp)
emax = rmax + max(0.1_dp * abs(rmax), 0.5_dp)
```

The 0.5 eV floor ensures a reasonable window even when bounds are tight.

### I2: Rename `dz_val → dx_val` in wire SC loop

**File**: `src/physics/sc_loop.f90` line 604

**Problem**: The variable `dz_val` is assigned `grid%dx`. In the wire convention, confinement is in the xy plane (grid x and y), and z is the free direction. The variable holds the x-spacing of the confinement grid, but the name suggests z-spacing.

**Fix**: Rename throughout `self_consistent_loop_wire`:
```fortran
! Before:
real(kind=dp) :: dy_val, dz_val
dz_val = grid%dx   ! AA (grid%dx for wire mode)

! After:
real(kind=dp) :: dy_val, dx_val
dx_val = grid%dx   ! AA — x-spacing of confinement plane
```

Update all references: `dy_val * dz_val` → `dy_val * dx_val`, `dz_val * ANGSTROM_TO_NM` → `dx_val * ANGSTROM_TO_NM`.

### I3: Cut-cell diagonal weighting — row-sum conservation

**File**: `src/math/sparse_matrices.f90` lines 886-891

**Problem**: For diagonal entries (col == row), the face-fraction weighting independently averages x and y face fractions. This does not guarantee the discrete Laplacian annihilates constants (row-sum = 0), which is required for box-integration accuracy at boundary cells.

**Fix**: Two-pass approach in `csr_apply_variable_coeff`:

Pass 1 — compute weighted off-diagonal entries and accumulate the negative sum:
```fortran
do k = rowptr(row), rowptr(row+1) - 1
  col = colind(k)
  if (col == row) then
    diag_pos = k
    cycle  ! skip diagonal for now
  end if
  ! ... compute face-fraction weighted value for off-diagonal ...
  values_out(k) = weighted_value
  diag_sum = diag_sum + weighted_value
end do
```

Pass 2 — set diagonal as negative sum (ensures row-sum = 0):
```fortran
if (diag_pos > 0) then
  values_out(diag_pos) = coeff(row) * cell_volume(row) * (-diag_sum)
end if
```

For interior cells (all faces = 1.0), the result is identical to the current formula. For boundary cells, the fix ensures exact conservation.

### I4: Deduplicate spin matrices

**Files**: `src/physics/gfactor_functions.f90`

**Problem**: The 8x8 spin matrices (SIGMA_X, SIGMA_Y, SIGMA_Z) are defined identically in both `sigmaElem` (QW, lines 83-110) and `sigmaElem_2d` (wire, lines 654-679). Any future correction to the phase convention must be applied consistently in both places.

**Fix**: Extract to module-level `save` parameters at the top of `gfactor_functions.f90`:
```fortran
! Module-level spin matrices (8-band basis)
complex(kind=dp), save :: SIGMA_X(8,8), SIGMA_Y(8,8), SIGMA_Z(8,8)
logical, save :: spin_matrices_initialized = .false.

private :: init_spin_matrices, spin_matrices_initialized, SIGMA_X, SIGMA_Y, SIGMA_Z

subroutine init_spin_matrices()
  if (spin_matrices_initialized) return
  SIGMA_X = cmplx(...)
  SIGMA_Y = cmplx(...)
  SIGMA_Z = cmplx(...)
  spin_matrices_initialized = .true.
end subroutine
```

Both `sigmaElem` and `sigmaElem_2d` call `init_spin_matrices()` and reference the shared arrays.

### I5: FEAST fallback warning

**File**: `src/math/eigensolver.f90` lines 53-58

**Problem**: When `USE_MKL_FEAST` is not defined, `solve_sparse_evp` silently falls through to `solve_dense_lapack` for FEAST/ARPACK requests. The dense solver ignores the energy window `[emin, emax]` and finds the smallest N eigenvalues instead — fundamentally different semantics.

**Fix**:
```fortran
#else
    case ('FEAST', 'ARPACK')
      print *, 'WARNING: FEAST/ARPACK not available. Falling back to dense LAPACK.'
      print *, '  Energy window [emin, emax] will be IGNORED.'
      print *, '  Results may differ from FEAST eigenvalue selection.'
      call solve_dense_lapack(H_csr, config, result)
#endif
```

### I6: `csr_set_values_from_coo` bounds check

**File**: `src/math/sparse_matrices.f90` lines 324-329

**Problem**: Out-of-range COO-to-CSR mappings are silently dropped. If the cache is corrupted (e.g., stale after grid resize), entries vanish from the Hamiltonian.

**Fix**:
```fortran
! Before:
if (csr_pos >= 1 .and. csr_pos <= mat%nnz) then
  mat%values(csr_pos) = mat%values(csr_pos) + vals(i)
end if

! After:
if (csr_pos < 1 .or. csr_pos > mat%nnz) then
  print *, 'ERROR: csr_set_values_from_coo: invalid mapping at COO index', i, &
    'csr_pos=', csr_pos, 'nnz=', mat%nnz
  stop 1
end if
mat%values(csr_pos) = mat%values(csr_pos) + vals(i)
```

### I7: DIIS fallback logging

**File**: `src/physics/sc_loop.f90` lines 357-363

**Problem**: When `dgesv` fails to solve the DIIS linear system, the code silently falls back to linear mixing. No diagnostic output.

**Fix**: Add one print statement:
```fortran
if (info /= 0) then
  print *, '  DIIS: linear system solve failed (info=', info, '). Falling back to linear mix.'
  deallocate(B, rhs, coeffs, ipiv)
  call linear_mix(phi_new, phi_old, phi_poisson, N, alpha)
  return
end if
```

### I8: Invalid `material_id` warning

**Files**: `src/physics/strain_solver.f90` lines 230-233, 130, 735

**Problem**: Grid points with `material_id = 0` (outside wire) or out-of-range values are silently skipped. If geometry initialization has a bug, strain has holes.

**Fix**: Add a warning on first occurrence:
```fortran
if (mid < 1 .or. mid > size(params)) then
  if (.not. warned_invalid_mat) then
    print *, 'WARNING: grid point', ij, 'has invalid material_id=', mid
    warned_invalid_mat = .true.
  end if
  cycle
end if
```

### I10: Fix Poisson mtype comment

**File**: `src/physics/poisson.f90` lines 296-305

**Problem**: Comment says "real symmetric positive definite, mtype=2" but code uses `mtype = 11` (real unsymmetric). The comment contradicts the code.

**Fix**:
```fortran
! Before:
! real symmetric positive definite, mtype=2

! After:
! Solve with MKL PARDISO (real unsymmetric, mtype=11)
! Box-integration with variable dielectric produces non-symmetric stencil
```

### I9: Deferred — `simulation_config` type bloat

The `simulation_config` god type (30+ fields) should be refactored by extracting a `wire_config` subtype, adding `grid_validate` and `config_validate` subroutines, and removing `fdStep = wire_ny` aliasing. This is tech debt, not a bug. Deferred to a dedicated refactoring PR after the wire feature is merged.

---

## Commit 4: `test: add wire strain PDE, wire SC loop, Fermi level, and cached CSR tests`

### T1: `test_strain_wire_biaxial`

**File**: `tests/unit/test_strain_solver.pf`

Tests `compute_strain` with `grid%ndim=2` (wire mode), dispatching to `compute_strain_wire` which solves the plane-strain PDE via PARDISO.

**Setup**: InAs core (r_inner=0, r_outer=25 A) in GaAs substrate (r_outer=50 A), 10x10 grid with circle geometry.

**Assertions**:
- `eps_xx` and `eps_yy` are non-zero at interior points (strain exists)
- `eps_xx` and `eps_yy` are continuous across the interface (displacement continuity)
- Strain is approximately zero far from the interface (stress-free BC)
- Tr(eps) is positive in the InAs core (hydrostatic tension from lattice mismatch, since InAs a0 > GaAs a0)

**Note**: This test requires PARDISO (MKL). It will only run when MKL is available. Guard with `#ifdef USE_MKL` if needed.

### T2: `test_sc_loop_wire_convergence`

**File**: `tests/unit/test_sc_loop.pf`

Tests `self_consistent_loop_wire` with a minimal wire configuration.

**Setup**: 3x3 grid, single material (GaAs), uniform ND doping = 1e18 cm^-3, 3 k-points, temperature = 300K.

**Assertions**:
- Returns `sc_converged = .true.` within max_iterations
- Eigenvalues at k=0 are non-zero and ordered (E1 < E2 < ...)
- `profile_2d` was modified from initial values (SC loop applied potential)

**Note**: This is an integration test. It exercises the full pipeline: grid → kpterms → Hamiltonian → FEAST → charge density → Poisson → DIIS → convergence.

### T3: `test_find_fermi_level_wire`

**File**: `tests/unit/test_sc_loop.pf`

Tests `find_fermi_level_wire` with known eigenvalues and kx grid.

**Setup**: 5 eigenvalues at known energies [-0.5, -0.3, 0.8, 1.0, 1.2] eV, 3 k-points (0.0, 0.05, 0.1) 1/A, 3x3 grid, temperature = 4K (step-like Fermi function).

**Assertions**:
- Fermi level is in the band gap (between highest VB and lowest CB)
- Charge neutrality holds within tolerance: `abs(total_charge - target_charge) < 1e-6`
- Fermi level shifts up with n-type doping

### T4: `test_csr_cached_values`

**File**: `tests/unit/test_csr_spmv.pf` or `tests/unit/test_hamiltonian_2d.pf`

Tests `csr_build_from_coo_cached` and `csr_set_values_from_coo` — the performance-critical path for wire k-point sweeps.

**Setup**: Build a CSR matrix from COO, then rebuild with different non-zero values using the cached path.

**Assertions**:
- Cached matrix has same sparsity pattern (same nnz, same rowptr, same colind)
- Cached matrix values match a freshly built matrix with the same COO values
- Building with the cache produces no out-of-range mappings

---

## Files Summary

### New files

| File | Purpose |
|---|---|
| `src/math/linalg.f90` | Centralized LAPACK/FEAST interface blocks |

### Deleted files

| File | Reason |
|---|---|
| `src/math/feast_call.f90` | Wrappers no longer needed with explicit interfaces |

### Modified files

| File | Changes |
|---|---|
| `CMakeLists.txt` | Replace `-ffast-math` with safe components (C1) |
| `src/math/CMakeLists.txt` | Add `linalg.f90`, remove `feast_call.f90` |
| `src/apps/main.f90` | `use linalg`, FEAST failure → `stop 1`, diag failure → `stop 1` |
| `src/apps/main_gfactor.f90` | `use linalg`, add MKL threading guard (C2) |
| `src/physics/sc_loop.f90` | `use linalg`, rename `dz_val→dx_val` (I2), DIIS logging (I7) |
| `src/math/eigensolver.f90` | `use linalg`, inline FEAST calls, `info=3` warning (C6), energy window floor (I1), fallback warning (I5) |
| `src/math/sparse_matrices.f90` | COO bounds check (I6), diagonal conservation (I3) |
| `src/physics/strain_solver.f90` | PARDISO → `stop 1` (C3), material_id warning (I8) |
| `src/physics/poisson.f90` | PARDISO → `stop 1` (C3), mtype comment fix (I10) |
| `src/physics/gfactor_functions.f90` | Deduplicate spin matrices (I4) |
| `tests/unit/test_strain_solver.pf` | Add wire strain PDE test (T1) |
| `tests/unit/test_sc_loop.pf` | Add wire SC loop and Fermi level tests (T2, T3) |
| `tests/unit/test_csr_spmv.pf` | Add cached CSR test (T4) |

## Progress

### Commit 1: `refactor: create centralized linalg module replacing external LAPACK/FEAST declarations` — DONE
- SHA: `8e8cc18`
- Created `src/math/linalg.f90` (117 lines) with 7 interface blocks
- Deleted `src/math/feast_call.f90`
- Updated callers: main.f90, main_gfactor.f90, eigensolver.f90, sc_loop.f90
- Bonus: fixed latent bug (vl/vu declared as integer, should be real(dp))
- 24/24 tests pass, zero warnings

### Commit 2: `fix: address critical PR review issues (C1-C6)` — DONE
- SHA: `9aad93e`
- C1: Replaced `-ffast-math` with safe components (kept -fno-signed-zeros, -freciprocal-math, -fno-trapping-math)
- C2: Added `mkl_set_num_threads_local(1)` in main_gfactor.f90
- C3: PARDISO errors → `stop 1` with cleanup in strain_solver.f90 and poisson.f90
- C4: FEAST non-convergence → `stop 1` when insufficient eigenvalues (nuanced: info=3 with enough eigenvalues continues with warning)
- C5: QW diag failure → `stop 1` instead of `cycle`
- C6: FEAST info=3 explicit warning with remediation advice
- 24/24 tests pass

### Commit 3: `fix: address important PR review issues (I1-I8, I10)` — DONE
- SHA: `9f5ff93`
- I1: Energy window floor `max(0.1*abs(), 0.5)` prevents zero-width windows
- I2: Renamed `dz_val→dx_val` in wire SC loop
- I3: Two-pass diagonal computation (row-sum=0) in `csr_apply_variable_coeff`
- I4: Deduplicated spin matrices to module-level save arrays
- I5: FEAST fallback warning in `#else` block
- I6: COO bounds check → `stop 1`
- I7: DIIS fallback logging on dgesv failure
- I8: Invalid material_id first-occurrence warning (3 locations)
- I10: Fixed Poisson mtype comment
- Updated test_hamiltonian_2d.pf for I3 diagonal convention change
- 24/24 tests pass

### Commit 4: `test: add cached CSR build and value mapping tests (T4)` — DONE (partial)
- SHA: `d193a2d`
- T4: `test_csr_cached_values` in test_csr_spmv.pf — PASS
- T1 (wire strain PDE), T2 (wire SC loop), T3 (Fermi level wire) — DEFERRED
  - T1: strain decay assertion failed on 10x10 grid; needs geometry/grid-size tuning
  - T2: required too many dependencies for integration test (eigensolver, geometry, kpterms); would need a mock/simplified approach
  - T3: depended on T2 infrastructure
  - Recommendation: implement T1-T3 as a separate follow-up PR with dedicated setup helpers

## Verification

1. All existing tests pass after each commit — DONE (24/24 after commits 1-3)
2. GCC 15 strict-aliasing note eliminated after Commit 1 — DONE
3. New tests (T1-T4) pass after Commit 4 — T4 PASS; T1-T3 deferred
4. Build with `cmake --build build` succeeds — DONE
4. Build with `cmake --build build` succeeds — DONE
5. `ctest --test-dir build` — all tests green — DONE (existing 24)
