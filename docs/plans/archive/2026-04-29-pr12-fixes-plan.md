# PR #12 Review Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix all 11 review issues found in PR #12 (modern Fortran migration) before merge.

**Architecture:** Mechanical fixes to 8 source files, 2 docs files. No algorithm changes. Each task is independent except Task 8 (PARDISO) which touches the same files as Task 7 (zdotc) in `linalg.f90`. Existing 39-test suite validates all changes.

**Tech Stack:** Fortran 2008, gfortran, Intel MKL (LP64, sequential), pFUnit for tests

---

## Task 1: Remove debug `print *` from eigensolver dispatch

Four `print *` statements fire per k-point (100-500 times per run). They are debug leftovers from the polymorphic dispatch implementation.

**Files:**
- Modify: `src/math/eigensolver.f90:146,149,866,877`

**Step 1: Remove the prints from old `solve_sparse_evp`**

At line 146, delete:
```fortran
      print *, '  Eigensolver: FEAST'
```

At line 149, delete:
```fortran
      print *, '  Eigensolver: dense LAPACK'
```

**Step 2: Remove the prints from polymorphic dispatch methods**

At line 866 (inside `feast_solve_dispatch`), delete:
```fortran
    print *, '  Eigensolver: FEAST'
```

At line 877 (inside `dense_lapack_solve_dispatch`), delete:
```fortran
    print *, '  Eigensolver: dense LAPACK'
```

**Step 3: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 4: Commit**

```bash
git add src/math/eigensolver.f90
git commit -m "fix: remove per-kpoint debug prints from eigensolver dispatch"
```

---

## Task 2: Remove `contiguous` from non-contiguous `eig` argument

`writeEigenvalues` declares `eig(:,:)` as `contiguous`, but `main.f90:448` passes `eig_wire(1:max_nev_found, :)` which is a non-contiguous row slice when `max_nev_found < size(eig_wire, 1)`. This violates the Fortran standard and creates hidden copy-in/copy-out temporaries.

**Files:**
- Modify: `src/io/outputFunctions.f90:149`

**Step 1: Remove the contiguous attribute**

Change line 149 from:
```fortran
      real(kind=dp), intent(in), contiguous :: eig(:,:)
```
to:
```fortran
      real(kind=dp), intent(in) :: eig(:,:)
```

**Step 2: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 3: Commit**

```bash
git add src/io/outputFunctions.f90
git commit -m "fix: remove contiguous from non-contiguous eig dummy argument"
```

---

## Task 3: Add `private` default to `gfactorFunctions`

This is the only module missing `private` default. Only `sigmaElem`, `sigmaElem_2d`, `gfactorCalculation`, and `gfactorCalculation_wire` are called externally (from `main_gfactor.f90`). Tests import the module only for transitive type access and don't call any procedures directly.

**Files:**
- Modify: `src/physics/gfactor_functions.f90:9-17`

**Step 1: Add private default and explicit public exports**

Change lines 9-17 from:
```fortran
  implicit none

  ! Module-level spin matrices (8-band basis).
  ! Shared between sigmaElem (QW) and sigmaElem_2d (wire).
  complex(kind=dp), save :: SIGMA_X(8,8), SIGMA_Y(8,8), SIGMA_Z(8,8)
  logical, save :: spin_matrices_initialized = .false.

  private :: init_spin_matrices, spin_matrices_initialized
  private :: SIGMA_X, SIGMA_Y, SIGMA_Z
```
to:
```fortran
  implicit none

  private
  public :: sigmaElem, sigmaElem_2d
  public :: gfactorCalculation, gfactorCalculation_wire

  ! Module-level spin matrices (8-band basis).
  ! Shared between sigmaElem (QW) and sigmaElem_2d (wire).
  complex(kind=dp), save :: SIGMA_X(8,8), SIGMA_Y(8,8), SIGMA_Z(8,8)
  logical, save :: spin_matrices_initialized = .false.
```

(The old `private ::` lines are removed because `private` default already covers everything; the `public ::` lines explicitly export the needed symbols.)

**Step 2: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass. If any test fails with "no implicit type" or similar, a symbol is missing from the `public` list — add it.

**Step 3: Commit**

```bash
git add src/physics/gfactor_functions.f90
git commit -m "refactor: add private default to gfactorFunctions module"
```

---

## Task 4: Fix README fpm section

The README says `fpm build --profile release` but `fpm.toml` doesn't encode MKL/FFTW3 link flags. The command won't actually build.

**Files:**
- Modify: `README.md:35-38`

**Step 1: Replace the fpm section**

Change lines 35-38 from:
```
**Alternative build with fpm** (optional, requires [fpm](https://fpm.fortran-lang.org/) v0.13.0+):
```bash
fpm build --profile release
```
See `fpm.toml` for the project manifest.
```
to:
```
**Experimental fpm manifest** — `fpm.toml` is provided as a project manifest but does not encode MKL/FFTW3 link flags. Use CMake for production builds.
```

**Step 2: Verify**

```bash
grep -A2 'fpm' README.md
```

Expected: the old `fpm build` command is gone, replaced by the experimental note.

**Step 3: Commit**

```bash
git add README.md
git commit -m "docs: mark fpm manifest as experimental in README"
```

---

## Task 5: Fix CLAUDE.md `value` contradiction

CLAUDE.md line 151 says "Scalars passed by reference (not `value`)" but `mkl_set_num_threads_local` correctly uses `value` because its C function takes `int` by value.

**Files:**
- Modify: `CLAUDE.md:151`

**Step 1: Update the convention**

Change line 151 from:
```
- `iso_c_binding` used for all MKL C APIs: PARDISO as `pardiso_c`, FEAST wrappers, `mkl_set_num_threads_local`. Scalars passed by reference (not `value`) since MKL C API passes everything as pointers.
```
to:
```
- `iso_c_binding` used for all MKL C APIs: PARDISO as `pardiso_c`, FEAST wrappers, `mkl_set_num_threads_local`. PARDISO/FEAST scalars passed by reference (MKL C API passes pointers). `mkl_set_num_threads_local` uses `value` since the C function takes `int` by value.
```

**Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: clarify iso_c_binding value/reference convention in CLAUDE.md"
```

---

## Task 6: Remove double finalization in `feast_solver_t`

When `feast_solver_t` is destroyed, `feast_solver_finalize` calls `feast_workspace_free(self%ws)`, then the `ws` component's own finalizer fires and calls `feast_workspace_free` again. Harmless but wasteful.

**Files:**
- Modify: `src/math/eigensolver.f90:881-884`

**Step 1: Remove explicit free from parent finalizer**

Change lines 881-884 from:
```fortran
  subroutine feast_solver_finalize(self)
    type(feast_solver_t), intent(inout) :: self
    call feast_workspace_free(self%ws)
  end subroutine feast_solver_finalize
```
to:
```fortran
  subroutine feast_solver_finalize(self)
    type(feast_solver_t), intent(inout) :: self
    ! ws component auto-finalizes via feast_workspace_finalize
  end subroutine feast_solver_finalize
```

**Step 2: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 3: Commit**

```bash
git add src/math/eigensolver.f90
git commit -m "refactor: remove redundant workspace free from feast_solver finalizer"
```

---

## Task 7: Fix `zdotc` interface — declare as subroutine

MKL's `zdotc` was compiled with ifort, which returns complex values via a hidden first-argument convention. Declaring it as a function in gfortran can cause ABI mismatch (garbage return values). The fix is to declare it as a subroutine with the result as the first output argument.

**Files:**
- Modify: `src/math/linalg.f90:104-112`
- Modify: `src/physics/optical_spectra.f90:162,266,782,976`

**Step 1: Change the zdotc interface in linalg.f90**

Replace lines 104-112:
```fortran
  ! zdotc - BLAS complex dot product (conjugated first vector)
  interface
    function zdotc(n, zx, incx, zy, incy) result(res)
      import :: dp
      integer, intent(in) :: n, incx, incy
      complex(kind=dp), intent(in) :: zx(*), zy(*)
      complex(kind=dp) :: res
    end function zdotc
  end interface
```
with:
```fortran
  ! zdotc - BLAS complex dot product (conjugated first vector)
  ! Declared as subroutine to match MKL's hidden-argument ABI for
  ! complex function returns (ifort convention on x86-64).
  interface
    subroutine zdotc(res, n, zx, incx, zy, incy)
      import :: dp
      complex(kind=dp), intent(out) :: res
      integer, intent(in) :: n, incx, incy
      complex(kind=dp), intent(in) :: zx(*), zy(*)
    end subroutine zdotc
  end interface
```

**Step 2: Update all call sites in optical_spectra.f90**

There are 4 call sites. Change each from function-call syntax to subroutine-call syntax.

Line 162, change:
```fortran
          Pele = zdotc(dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
```
to:
```fortran
          call zdotc(Pele, dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
```

Line 266, change:
```fortran
          Pele = zdotc(dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
```
to:
```fortran
          call zdotc(Pele, dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
```

Line 782, change:
```fortran
        pele_ij = zdotc(dim, eigvecs(1:dim,state_i), 1, Ytmp, 1)
```
to:
```fortran
        call zdotc(pele_ij, dim, eigvecs(1:dim,state_i), 1, Ytmp, 1)
```

Line 976, change:
```fortran
          Pele = zdotc(dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
```
to:
```fortran
          call zdotc(Pele, dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
```

**Step 3: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass. The optical regression tests compare numerical output against golden files — if the zdotc return was previously wrong, values may change. Check the regression test output carefully.

**Step 4: Commit**

```bash
git add src/math/linalg.f90 src/physics/optical_spectra.f90
git commit -m "fix: declare zdotc as subroutine to match MKL complex-return ABI"
```

---

## Task 8: Fix PARDISO interfaces — `c_intptr_t` for pt + real-valued interface for Poisson

PARDISO `pt` array should use `c_intptr_t` (pointer-sized integer) instead of `c_int64_t` for portability. Additionally, `poisson.f90` calls PARDISO through implicit interface with real-valued arrays — the migration is incomplete. Add a real-valued PARDISO interface and have `poisson.f90` use it.

**Files:**
- Modify: `src/math/linalg.f90:9,114-128`
- Modify: `src/math/eigensolver.f90:373`
- Modify: `src/physics/poisson.f90:1-10,203`

**Step 1: Update linalg.f90 imports**

Change line 9 from:
```fortran
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_double, c_double_complex, c_char
```
to:
```fortran
  use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_double_complex, c_char
```

**Step 2: Add `c_intptr_t` to the complex PARDISO interface import**

In the complex PARDISO interface (around line 119), change:
```fortran
      import :: c_int, c_int64_t, c_double_complex
```
to:
```fortran
      import :: c_int, c_intptr_t, c_double_complex
```

And change line 120:
```fortran
      integer(c_int64_t), intent(inout) :: pt(64)
```
to:
```fortran
      integer(c_intptr_t), intent(inout) :: pt(64)
```

**Step 3: Add real-valued PARDISO interface to linalg.f90**

After the complex PARDISO interface block (after line 128), add:

```fortran
  ! pardiso_real - MKL PARDISO for real-valued matrices (iso_c_binding)
  interface
    subroutine pardiso_real(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, &
                         nrhs, iparm, msglvl, b, x, error) bind(C, name="PARDISO")
      import :: c_int, c_intptr_t, c_double
      integer(c_intptr_t), intent(inout) :: pt(64)
      integer(c_int), intent(in) :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
      real(c_double), intent(in) :: a(*)
      integer(c_int), intent(in) :: ia(*), ja(*), perm(*)
      integer(c_int), intent(inout) :: iparm(64)
      real(c_double), intent(inout) :: b(*), x(*)
      integer(c_int), intent(out) :: error
    end subroutine pardiso_real
  end interface
```

**Step 4: Export `pardiso_real` from linalg.f90**

Add to the public declarations (around line 27), after the `pardiso_c` public:

```fortran
  public :: pardiso_real
```

And move the PARDISO publics outside the `#ifdef USE_ARPACK` guard since `pardiso_real` is needed by `poisson.f90` regardless of ARPACK:

Change lines 25-28 from:
```fortran
  ! MKL PARDISO (guarded)
#ifdef USE_ARPACK
  public :: pardiso_c
#endif
```
to:
```fortran
  ! MKL PARDISO
  public :: pardiso_real
#ifdef USE_ARPACK
  public :: pardiso_c
#endif
```

**Step 5: Update eigensolver.f90 pt declaration**

Change line 373 from:
```fortran
    integer(8) :: pt(64)
```
to:
```fortran
    integer(kind=c_intptr_t) :: pt(64)
```

And add the import at the top of `solve_arpack` (around line 355):
```fortran
    use, intrinsic :: iso_c_binding, only: c_intptr_t
```

**Step 6: Update poisson.f90 to use explicit PARDISO interface**

Add `use linalg, only: pardiso_real` to the module imports. The poisson module already has:
```fortran
  use definitions
```
Add after it:
```fortran
  use linalg, only: pardiso_real
```

Change line 203 from:
```fortran
    integer(8) :: pt(64)
```
to:
```fortran
    integer(kind=c_intptr_t) :: pt(64)
```

Add import at the top of `poisson_solve_2d` subroutine (or module level):
```fortran
  use, intrinsic :: iso_c_binding, only: c_intptr_t
```

Change all three `call pardiso(...)` to `call pardiso_real(...)` at lines 320, 327, 335.

**Step 7: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass. The Poisson regression tests verify the real-valued PARDISO path.

**Step 8: Commit**

```bash
git add src/math/linalg.f90 src/math/eigensolver.f90 src/physics/poisson.f90
git commit -m "fix: use c_intptr_t for PARDISO pt, add real-valued interface for Poisson"
```

---

## Task 9: Guard `feast_solver_t` in `#ifdef USE_MKL_FEAST`

The `feast_solver_t` type, its dispatch method, and its finalizer are defined unconditionally but reference FEAST-specific routines. If `USE_MKL_FEAST` is ever undefined, this causes compilation failure.

**Files:**
- Modify: `src/math/eigensolver.f90:83-91,859-868,881-884`

**Step 1: Wrap the feast_solver_t type definition**

Change lines 83-91 from:
```fortran
  ! ------------------------------------------------------------------
  ! Concrete solver types for polymorphic dispatch.
  ! ------------------------------------------------------------------
  type, extends(eigensolver_base) :: feast_solver_t
    type(feast_workspace) :: ws
  contains
    procedure :: solve => feast_solve_dispatch
    final :: feast_solver_finalize
  end type feast_solver_t
```
to:
```fortran
  ! ------------------------------------------------------------------
  ! Concrete solver types for polymorphic dispatch.
  ! ------------------------------------------------------------------
#ifdef USE_MKL_FEAST
  type, extends(eigensolver_base) :: feast_solver_t
    type(feast_workspace) :: ws
  contains
    procedure :: solve => feast_solve_dispatch
    final :: feast_solver_finalize
  end type feast_solver_t
#endif
```

**Step 2: Wrap the dispatch implementation and finalizer**

The `feast_solve_dispatch` subroutine (around line 859-868) and `feast_solver_finalize` (around line 881-884) are already inside `#ifdef USE_MKL_FEAST` — verify this. If `feast_solver_finalize` is NOT inside the `#ifdef`, wrap it:

```fortran
#ifdef USE_MKL_FEAST
  subroutine feast_solve_dispatch(self, H_csr, config, result)
    ...
  end subroutine feast_solve_dispatch
#endif

  subroutine feast_solver_finalize(self)
    type(feast_solver_t), intent(inout) :: self
    ! ws component auto-finalizes via feast_workspace_finalize
  end subroutine feast_solver_finalize
```

Wait — `feast_solver_finalize` references `feast_solver_t`, so it MUST also be inside `#ifdef USE_MKL_FEAST`. Wrap it too.

**Step 3: Update the public exports**

Change line 16 from:
```fortran
  public :: eigensolver_base, feast_solver_t, dense_lapack_solver_t
```
to:
```fortran
  public :: eigensolver_base, dense_lapack_solver_t
#ifdef USE_MKL_FEAST
  public :: feast_solver_t
#endif
```

**Step 4: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 5: Commit**

```bash
git add src/math/eigensolver.f90
git commit -m "fix: guard feast_solver_t in #ifdef USE_MKL_FEAST"
```

---

## Task 10: Fix `make_eigensolver` factory — restore ARPACK awareness and error on unknowns

The factory silently maps ARPACK and typos to dense LAPACK. The old `solve_sparse_evp` had explicit ARPACK handling and errored on unknown methods.

**Files:**
- Modify: `src/math/eigensolver.f90:886-903`

**Step 1: Replace the factory**

Change lines 886-903 from:
```fortran
  function make_eigensolver(config) result(solver)
    class(eigensolver_base), allocatable :: solver
    type(eigensolver_config), intent(in) :: config

    select case (trim(config%method))
#ifdef USE_MKL_FEAST
    case ('FEAST')
      allocate(feast_solver_t :: solver)
#else
    case ('FEAST')
      allocate(dense_lapack_solver_t :: solver)
#endif
    case ('DENSE')
      allocate(dense_lapack_solver_t :: solver)
    case default
      allocate(dense_lapack_solver_t :: solver)
    end select
  end function make_eigensolver
```
to:
```fortran
  function make_eigensolver(config) result(solver)
    class(eigensolver_base), allocatable :: solver
    type(eigensolver_config), intent(in) :: config

    select case (trim(config%method))
    case ('DENSE')
      allocate(dense_lapack_solver_t :: solver)
#ifdef USE_MKL_FEAST
    case ('FEAST')
      allocate(feast_solver_t :: solver)
#else
    case ('FEAST')
      print *, 'WARNING: FEAST requested but not compiled; using dense LAPACK'
      allocate(dense_lapack_solver_t :: solver)
#endif
    case ('ARPACK')
      print *, 'WARNING: ARPACK dispatch not yet polymorphic; using dense LAPACK'
      allocate(dense_lapack_solver_t :: solver)
    case default
      error stop 'Unknown eigensolver method: '//trim(config%method)
    end select
  end function make_eigensolver
```

**Step 2: Build and test**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass. All existing configs use 'FEAST' or 'DENSE', so the `error stop` won't fire.

**Step 3: Commit**

```bash
git add src/math/eigensolver.f90
git commit -m "fix: restore ARPACK awareness in make_eigensolver, error on unknown methods"
```

---

## Task 11: Document FEAST `bind(C)` x86-64 assumption

The FEAST interfaces use `bind(C)` which binds to MKL's C-exported symbols. This is safe on x86-64 Linux where C and Fortran pointer-passing ABIs are identical, but may not be portable.

**Files:**
- Modify: `src/math/linalg.f90:131-133`

**Step 1: Add comment before FEAST interfaces**

Change lines 131-133 from:
```fortran
#ifdef USE_MKL_FEAST
  ! feastinit - FEAST initialization (iso_c_binding)
  interface
```
to:
```fortran
#ifdef USE_MKL_FEAST
  ! FEAST interfaces use bind(C) binding to MKL's C-exported symbols.
  ! Safe on x86-64 Linux where C and Fortran pointer-passing ABIs are identical.
  ! May require adjustment for other platforms.

  ! feastinit - FEAST initialization (iso_c_binding)
  interface
```

**Step 2: Commit**

```bash
git add src/math/linalg.f90
git commit -m "docs: document x86-64 ABI assumption in FEAST bind(C) interfaces"
```

---

## Task 12: Final validation and CLAUDE.md update

Run the full test suite one final time and update CLAUDE.md to reflect the zdotc and PARDISO changes.

**Files:**
- Modify: `CLAUDE.md`

**Step 1: Run full test suite**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 2: Update CLAUDE.md conventions**

Add after the existing iso_c_binding convention (around line 151):

```
- `zdotc` declared as subroutine (not function) in `linalg.f90` to match MKL's hidden-argument ABI for complex function returns on gfortran. All call sites use `call zdotc(result, ...)` syntax.
- PARDISO has two `iso_c_binding` interfaces: `pardiso_c` (complex, used by ARPACK eigensolver) and `pardiso_real` (real-valued, used by Poisson solver). Both use `c_intptr_t` for the `pt` handle array.
```

**Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md with zdotc and PARDISO real interface conventions"
```
