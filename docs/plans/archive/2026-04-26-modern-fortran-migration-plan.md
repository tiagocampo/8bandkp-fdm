# Modern Fortran Migration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Migrate the current post-PR-11 codebase toward a cleaner Fortran 2008/2018 style while preserving the Hamiltonian performance refactor, numerical outputs, and cache-safety guarantees.

**Architecture:** Treat the merged Hamiltonian refactor as the baseline: `confinement_init.f90`, `hamiltonian_wire.f90`, and `hamiltonianConstructor.f90` remain separate modules, and `wire_workspace` plus `feast_workspace` remain explicit cache/workspace types. Start with low-risk F2008 mechanical cleanup, then encapsulate resource-management boundaries around CSR, FEAST, wire workspaces, parser state, and configuration validation without changing physics formulas, basis ordering, material parameters, FD coefficients, or Hamiltonian signs.

**Tech Stack:** Fortran 90/2003/2008 source, CMake/Ninja, Intel MKL FEAST/SpBLAS/PARDISO, FFTW3, pFUnit, shell/Python regression tests.

---

## Constraints From The Hamiltonian Performance Refactor

Do not undo these merged design choices:

- `src/physics/hamiltonian_wire.f90` owns 2D wire Hamiltonian assembly, `wire_workspace`, `wire_coo_cache`, commutator velocity matrices, and the table-driven strain COO insertion.
- `src/physics/confinement_init.f90` owns confinement initialization for 1D/2D grids.
- `src/physics/hamiltonianConstructor.f90` now focuses on bulk/QW Hamiltonians and shared dense strain insertion.
- `wire_workspace` caches normal wire-Hamiltonian CSR block structures and COO mappings. `g='g3'` derivative-Hamiltonian builds must not initialize or mutate it.
- `feast_workspace` caches FEAST upper-triangle CSR structure. Reuse must remain guarded by sparsity-pattern validation.
- Any OOP/finalizer work must preserve explicit `*_free` routines until all call sites and tests are migrated.

---

### Task 1: Make the Fortran 2008 Baseline Explicit

**Files:**
- Modify: `CMakeLists.txt`
- Test: full configure/build

**Step 1: Add target-level Fortran standard properties**

In `CMakeLists.txt`, after `project(8bandkp-fdm Fortran)`, add:

```cmake
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)
```

If `CMAKE_Fortran_EXTENSIONS OFF` breaks preprocessing or compiler-specific MKL flags, move the standard setting to targets in `src/CMakeLists.txt` instead:

```cmake
set_target_properties(8bandkp_common bandStructure gfactorCalculation opticalProperties
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    Fortran_STANDARD 2008
    Fortran_STANDARD_REQUIRED ON
)
```

**Step 2: Configure from scratch**

Run:

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
```

Expected: configure succeeds. If compiler-specific flags conflict with strict standard mode, keep `Fortran_STANDARD 2008` but allow extensions and document the reason in a code comment.

**Step 3: Build**

Run:

```bash
cmake --build build
```

Expected: all three executables and tests compile.

**Step 4: Run smoke tests**

Run:

```bash
ctest --test-dir build -R "test_defs|test_eigensolver|test_hamiltonian_2d" --output-on-failure
```

Expected: selected tests pass.

**Step 5: Commit**

```bash
git add CMakeLists.txt src/CMakeLists.txt
git commit -m "build: require Fortran 2008 standard"
```

---

### Task 2: Replace Legacy Intrinsics With Generic Intrinsics

**Files:**
- Modify: `src/core/defs.f90`
- Modify: `src/core/parameters.f90`
- Test: `tests/unit/test_defs.pf`, `tests/unit/test_parameters.pf`

**Step 1: Replace `dsqrt` in definitions constants**

In `src/core/defs.f90`, replace:

```fortran
real(kind=dp), parameter :: SQR3 = dsqrt(3.0_dp)
real(kind=dp), parameter :: SQR2 = dsqrt(2.0_dp)
real(kind=dp), parameter :: SQR2o3 = dsqrt(2.0_dp/3.0_dp)
real(kind=dp), parameter :: SQR3o2 = dsqrt(1.5_dp)
```

with:

```fortran
real(kind=dp), parameter :: SQR3 = sqrt(3.0_dp)
real(kind=dp), parameter :: SQR2 = sqrt(2.0_dp)
real(kind=dp), parameter :: SQR2o3 = sqrt(2.0_dp/3.0_dp)
real(kind=dp), parameter :: SQR3o2 = sqrt(1.5_dp)
```

**Step 2: Replace `dsqrt` in material-parameter calculations**

In `src/core/parameters.f90`, replace:

```fortran
params(i)%P = dsqrt(params(i)%EP*const)
P = dsqrt(EP*const)
```

with:

```fortran
params(i)%P = sqrt(params(i)%EP*const)
P = sqrt(EP*const)
```

Do not alter material parameter values.

**Step 3: Run tests**

Run:

```bash
ctest --test-dir build -R "test_defs|test_parameters" --output-on-failure
```

Expected: both tests pass.

**Step 4: Commit**

```bash
git add src/core/defs.f90 src/core/parameters.f90
git commit -m "refactor: replace legacy dsqrt intrinsics"
```

---

### Task 3: Replace `forall` Outside Hamiltonian Code First

**Files:**
- Modify: `src/core/utils.f90`
- Modify: `src/math/finitedifferences.f90`
- Test: `tests/unit/test_utils.pf`, `tests/unit/test_finitedifferences.pf`

**Step 1: Replace Simpson coefficient `forall` with loops**

In `src/core/utils.f90`, replace:

```fortran
forall(i=2:N-1:2) sc(i) = 4.0_dp
forall(i=3:N-1:2) sc(i) = 2.0_dp
```

with:

```fortran
do i = 2, N - 1, 2
  sc(i) = 4.0_dp
end do
do i = 3, N - 1, 2
  sc(i) = 2.0_dp
end do
```

Use `do concurrent` only if the loop body has no hidden dependencies and the compiler accepts it cleanly under the project flags.

**Step 2: Replace identity-matrix `forall` in finite differences**

In `src/math/finitedifferences.f90`, replace:

```fortran
forall(i = 1:n) matrix(i,i) = 1*factor
forall(i = 1:n) matrix(i,i) = 1
```

with:

```fortran
do i = 1, n
  matrix(i,i) = factor
end do
```

and:

```fortran
do i = 1, n
  matrix(i,i) = 1
end do
```

Keep coefficient values and FD stencil logic unchanged.

**Step 3: Run targeted tests**

Run:

```bash
ctest --test-dir build -R "test_utils|test_finitedifferences" --output-on-failure
```

Expected: both tests pass.

**Step 4: Commit**

```bash
git add src/core/utils.f90 src/math/finitedifferences.f90
git commit -m "refactor: replace simple forall constructs"
```

---

### Task 4: Replace `forall` In Refactored Confinement Initialization

**Files:**
- Modify: `src/physics/confinement_init.f90`
- Test: `tests/unit/test_hamiltonian.pf`, `tests/unit/test_hamiltonian_2d.pf`, wire regressions

**Step 1: Replace 1D `forall` blocks with explicit loops**

In `src/physics/confinement_init.f90`, replace each `forall` block with equivalent `do` loops. Preserve current indexing exactly.

Example pattern:

```fortran
forall(ii=2:N-1)
  kpterms(ii,ii+1,9) = ...
  kpterms(ii,ii-1,9) = ...
end forall
```

becomes:

```fortran
do ii = 2, N - 1
  kpterms(ii,ii+1,9) = ...
  kpterms(ii,ii-1,9) = ...
end do
```

**Step 2: Replace profile initialization `forall` blocks**

Keep the `profile(ii,*)` assignments in the same order. Do not change band-offset conventions.

**Step 3: Build**

Run:

```bash
cmake --build build
```

Expected: clean build.

**Step 4: Run targeted tests**

Run:

```bash
ctest --test-dir build -R "test_hamiltonian|test_hamiltonian_2d" --output-on-failure
```

Expected: both tests pass.

**Step 5: Run wire smoke regressions**

Run:

```bash
ctest --test-dir build -R "regression_wire_gaas_rectangle|regression_wire_dense_sparse_consistency" --output-on-failure
```

Expected: both pass.

**Step 6: Commit**

```bash
git add src/physics/confinement_init.f90
git commit -m "refactor: replace confinement forall constructs"
```

---

### Task 5: Replace Remaining `forall` In Bulk/QW Hamiltonian Construction

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90`
- Test: `tests/unit/test_hamiltonian.pf`, QW and g-factor regressions

**Step 1: Replace `forall(ii=1:N, jj=1:N)` blocks**

In `src/physics/hamiltonianConstructor.f90`, replace dense matrix assignment `forall` blocks with nested `do` loops:

```fortran
do jj = 1, N
  do ii = 1, N
    Q(ii,jj) = ...
    T(ii,jj) = ...
    S(ii,jj) = ...
    SC(jj,ii) = ...
    PZ(ii,jj) = ...
    A(ii,jj) = ...
  end do
end do
```

Keep the existing `SC(jj,ii)` transpose convention unchanged.

**Step 2: Replace diagonal-only `forall` blocks**

Replace:

```fortran
forall (ii=1:N)
  R(ii,ii) = ...
  RC(ii,ii) = ...
  PP(ii,ii) = ...
  PM(ii,ii) = ...
end forall
```

with:

```fortran
do ii = 1, N
  R(ii,ii) = ...
  RC(ii,ii) = ...
  PP(ii,ii) = ...
  PM(ii,ii) = ...
end do
```

**Step 3: Replace profile diagonal `forall`**

Preserve basis ordering:

- bands 1-4: valence
- bands 5-6: split-off
- bands 7-8: conduction

Do not change the array layout.

**Step 4: Run targeted tests**

Run:

```bash
ctest --test-dir build -R "test_hamiltonian|regression_qw_alsbw_gasbw_inasw|regression_gfactor_cb" --output-on-failure
```

Expected: all pass.

**Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "refactor: replace Hamiltonian forall constructs"
```

---

### Task 6: Remove `goto` From Optional Parser Blocks

**Files:**
- Modify: `src/io/input_parser.f90`
- Test: parser-sensitive regressions

**Step 1: Add helper for optional block status reset**

Add a small internal helper near the other parser helpers:

```fortran
subroutine reset_optional_block_status(status)
  integer, intent(inout) :: status
  if (status /= 0) status = 0
end subroutine reset_optional_block_status
```

This preserves the current behavior where an incomplete optional block falls through and parsing continues with defaults.

**Step 2: Refactor the SC block label `100`**

Replace each:

```fortran
if (status /= 0) goto 100
```

inside the SC block with a named block and `exit`:

```fortran
sc_optional: block
  read(data_unit, *, iostat=status) label, cfg%sc%max_iterations
  if (status /= 0) exit sc_optional
  print *, trim(label), cfg%sc%max_iterations

  ! repeat for the rest of the SC fields
end block sc_optional
call reset_optional_block_status(status)
```

Keep doping parsing behavior unchanged.

**Step 3: Refactor labels `101`, `102`, `103`, and `200`**

Use the same named-block pattern for:

- optics block
- exciton block
- scattering block
- strain block

Do not reorder optional block labels or field reads. The parser is sequence-sensitive even though block labels are name-aware.

**Step 4: Refactor trailing `goto` uses near FEAST/material options**

For any remaining `goto` instances in `input_parser.f90`, either use a named block or convert to local helper routines returning `status`.

**Step 5: Confirm no `goto` remains in parser**

Run:

```bash
rg -n "\bgoto\b" src/io/input_parser.f90
```

Expected: no matches.

**Step 6: Run parser-sensitive tests**

Run:

```bash
ctest --test-dir build -R "regression_wire_parser_optional_sections|regression_gfactor_no_optics|regression_qw_character_and_output_dir|regression_sc_gaas_alas_qw|regression_sc_qw_inas_alsb" --output-on-failure
```

Expected: all pass.

**Step 7: Commit**

```bash
git add src/io/input_parser.f90
git commit -m "refactor: replace parser goto fallthroughs"
```

---

### Task 7: Replace External BLAS/PARDISO Declarations With Explicit Interfaces

**Files:**
- Modify: `src/math/linalg.f90`
- Modify: `src/math/eigensolver.f90`
- Modify: `src/physics/optical_spectra.f90`
- Test: eigensolver and optics tests

**Step 1: Add `zdotc` interface to `linalg.f90`**

In `src/math/linalg.f90`, add a public interface for `zdotc` using the same ABI style as existing LAPACK/MKL interfaces. Export it from the module.

Expected interface shape:

```fortran
interface
  function zdotc(n, zx, incx, zy, incy) result(res)
    import :: dp
    integer, intent(in) :: n, incx, incy
    complex(kind=dp), intent(in) :: zx(*), zy(*)
    complex(kind=dp) :: res
  end function zdotc
end interface
```

If the compiler/MKL binding requires a subroutine form instead of a function form, match the currently working symbol usage before changing call sites.

**Step 2: Import `zdotc` in optical spectra**

In `src/physics/optical_spectra.f90`, replace every local:

```fortran
complex(kind=dp), external :: zdotc
```

with module import:

```fortran
use linalg, only: zdotc
```

Add it to the existing `use linalg` statement if one exists.

**Step 3: Add explicit PARDISO interface or central wrapper**

Move the `external :: pardiso` declaration out of `src/math/eigensolver.f90` and into `src/math/linalg.f90` as an explicit interface or wrapper subroutine. Keep `-fno-strict-aliasing` until the PARDISO binding is verified under GNU and ifx.

**Step 4: Build**

Run:

```bash
cmake --build build
```

Expected: no interface mismatch warnings or link errors.

**Step 5: Run tests**

Run:

```bash
ctest --test-dir build -R "test_eigensolver|test_optical|test_optical_qw|regression_wire_optical_selection|regression_qw_optics_commutator" --output-on-failure
```

Expected: all pass.

**Step 6: Commit**

```bash
git add src/math/linalg.f90 src/math/eigensolver.f90 src/physics/optical_spectra.f90
git commit -m "refactor: centralize external BLAS and PARDISO interfaces"
```

---

### Task 8: Add Finalizers While Preserving Explicit Free Routines

**Files:**
- Modify: `src/math/sparse_matrices.f90`
- Modify: `src/physics/hamiltonian_wire.f90`
- Modify: `src/math/eigensolver.f90`
- Test: workspace/unit tests

**Step 1: Add finalizer for `csr_matrix`**

In `src/math/sparse_matrices.f90`, update the type:

```fortran
type :: csr_matrix
  integer :: nrows = 0
  integer :: ncols = 0
  integer :: nnz   = 0
  complex(kind=dp), allocatable :: values(:)
  integer, allocatable          :: colind(:)
  integer, allocatable          :: rowptr(:)
contains
  final :: csr_finalize
end type csr_matrix
```

Add:

```fortran
subroutine csr_finalize(mat)
  type(csr_matrix), intent(inout) :: mat
  call csr_free(mat)
end subroutine csr_finalize
```

Keep `csr_free` public. Existing manual frees remain valid.

**Step 2: Add finalizer for `wire_workspace`**

In `src/physics/hamiltonian_wire.f90`, add:

```fortran
contains
  final :: wire_workspace_finalize
```

inside `type :: wire_workspace`, and implement:

```fortran
subroutine wire_workspace_finalize(ws)
  type(wire_workspace), intent(inout) :: ws
  call wire_workspace_free(ws)
end subroutine wire_workspace_finalize
```

Do not finalize `wire_coo_cache` until call-site behavior is audited; keep `wire_coo_cache_free`.

**Step 3: Add finalizer for `feast_workspace`**

In `src/math/eigensolver.f90`, add:

```fortran
contains
  final :: feast_workspace_finalize
```

inside `type :: feast_workspace`, and implement:

```fortran
subroutine feast_workspace_finalize(fw)
  type(feast_workspace), intent(inout) :: fw
  call feast_workspace_free(fw)
end subroutine feast_workspace_finalize
```

**Step 4: Run focused tests**

Run:

```bash
ctest --test-dir build -R "test_csr_spmv|test_hamiltonian_2d|test_eigensolver" --output-on-failure
```

Expected: all pass.

**Step 5: Run one wire regression**

Run:

```bash
ctest --test-dir build -R "regression_wire_gaas_rectangle" --output-on-failure
```

Expected: pass.

**Step 6: Commit**

```bash
git add src/math/sparse_matrices.f90 src/physics/hamiltonian_wire.f90 src/math/eigensolver.f90
git commit -m "refactor: add finalizers for core workspace types"
```

---

### Task 9: Introduce Type-Bound Convenience Methods Conservatively

**Files:**
- Modify: `src/math/sparse_matrices.f90`
- Modify call sites only where trivial
- Test: CSR/eigensolver/wire tests

**Step 1: Add type-bound wrappers without hiding internals**

In `type :: csr_matrix`, add:

```fortran
contains
  procedure :: free => csr_matrix_free_bound
  procedure :: clone_structure => csr_matrix_clone_structure_bound
  final :: csr_finalize
end type csr_matrix
```

Implement wrappers:

```fortran
subroutine csr_matrix_free_bound(this)
  class(csr_matrix), intent(inout) :: this
  call csr_free(this)
end subroutine csr_matrix_free_bound

subroutine csr_matrix_clone_structure_bound(this, dst)
  class(csr_matrix), intent(in) :: this
  type(csr_matrix), intent(out) :: dst
  call csr_clone_structure(this, dst)
end subroutine csr_matrix_clone_structure_bound
```

**Step 2: Do not make components private yet**

Do not change `values`, `colind`, `rowptr`, `nrows`, `ncols`, or `nnz` visibility in this task. The Hamiltonian fast paths use direct component access for performance and clarity. Hiding internals is a later migration after profiling.

**Step 3: Migrate one low-risk call site**

In one test or non-hot path, replace:

```fortran
call csr_free(H)
```

with:

```fortran
call H%free()
```

Keep broad call-site migration out of this task.

**Step 4: Run tests**

Run:

```bash
ctest --test-dir build -R "test_csr_spmv|test_hamiltonian_2d" --output-on-failure
```

Expected: pass.

**Step 5: Commit**

```bash
git add src/math/sparse_matrices.f90 tests/unit/test_csr_spmv.pf
git commit -m "refactor: add type-bound CSR convenience methods"
```

---

### Task 10: Add Validate-Once Entry Point For `simulation_config`

**Files:**
- Modify: `src/core/defs.f90`
- Modify: `src/io/input_parser.f90`
- Test: parser/unit/regression tests

**Step 1: Add validation result helper if not already present**

Keep this simple; do not introduce a full diagnostics framework.

```fortran
type :: config_validation_result
  logical :: ok = .true.
  character(len=256) :: message = ''
end type config_validation_result
```

**Step 2: Add `validate_simulation_config` procedure**

Prefer a module procedure in `input_parser.f90` rather than a type-bound method on `simulation_config` for this first step:

```fortran
subroutine validate_simulation_config(cfg)
  type(simulation_config), intent(in) :: cfg

  if (cfg%confDir /= 'z' .and. cfg%confDir /= 'x' .and. cfg%confDir /= 'y') then
    error stop 'validate_simulation_config: invalid confDir'
  end if

  if (cfg%fdStep <= 0) then
    error stop 'validate_simulation_config: FDstep must be positive'
  end if

  if (cfg%numcb < 0 .or. cfg%numvb < 0) then
    error stop 'validate_simulation_config: numcb/numvb must be non-negative'
  end if
end subroutine validate_simulation_config
```

Extend with only validations already assumed by downstream code. Do not change defaults.

**Step 3: Call validation once after parsing**

At the end of `read_input_config`, call:

```fortran
call validate_simulation_config(cfg)
```

**Step 4: Add/adjust tests**

Use existing parser regression coverage first. If adding a unit test is practical, add one malformed config test that expects failure. Do not commit `input.cfg`.

**Step 5: Run tests**

Run:

```bash
ctest --test-dir build -R "regression_bulk_gaas|regression_wire_parser_optional_sections|regression_gfactor_no_optics" --output-on-failure
```

Expected: pass.

**Step 6: Commit**

```bash
git add src/core/defs.f90 src/io/input_parser.f90 tests/unit/test_parameters.pf
git commit -m "refactor: validate simulation config after parsing"
```

---

### Task 11: Add `contiguous` Only To Proven Contiguous Hot-Path Arguments

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90`
- Modify: `src/physics/confinement_init.f90`
- Modify: `src/physics/optical_spectra.f90`
- Test: targeted wire/optics regressions

**Step 1: Identify safe candidates**

Only add `contiguous` to assumed-shape arrays that are always passed whole, allocated arrays in current call sites. Initial candidates:

- `profile_2d(:,:)` in `ZB8bandGeneralized` and helpers
- `coo_rows(:)`, `coo_cols(:)`, `coo_vals(:)` insertion helpers
- dense optical spectra arrays passed as whole arrays

Do not add `contiguous` to array slices, pointer targets, or generic interfaces without checking every call site.

**Step 2: Add attribute one procedure group at a time**

Example:

```fortran
real(kind=dp), intent(in), contiguous :: profile_2d(:,:)
```

For internal COO insertion helpers:

```fortran
integer, intent(inout), contiguous :: coo_r(:), coo_c(:)
complex(kind=dp), intent(inout), contiguous :: coo_v(:)
```

**Step 3: Build after each file**

Run:

```bash
cmake --build build
```

Expected: no interface errors.

**Step 4: Run targeted tests**

Run:

```bash
ctest --test-dir build -R "test_hamiltonian_2d|test_optical|regression_wire_optical_selection|regression_qw_absorption_polarization" --output-on-failure
```

Expected: all pass.

**Step 5: Commit**

```bash
git add src/physics/hamiltonian_wire.f90 src/physics/confinement_init.f90 src/physics/optical_spectra.f90
git commit -m "perf: mark proven contiguous hot-path arrays"
```

---

### Task 12: Document The New Modernization Baseline

**Files:**
- Modify: `README.md`
- Modify: `CLAUDE.md`
- Modify: `docs/plans/2026-04-26-modern-fortran-migration-design.md` if needed

**Step 1: Update build requirements**

In `README.md`, state that the code is built as Fortran 2008 and still requires:

- gfortran or ifx
- CMake >= 3.15
- Ninja optional
- Intel MKL
- FFTW3
- pFUnit for tests

Do not add stdlib or fpm as required dependencies unless those tasks have actually been implemented.

**Step 2: Update code conventions**

In `CLAUDE.md`, add:

- Prefer `do` / `do concurrent` over `forall`.
- Prefer generic intrinsics (`sqrt`) over legacy typed intrinsics (`dsqrt`).
- Keep direct workspace/free routines available even after adding finalizers.
- Do not hide CSR internals until wire fast paths are profiled and migrated.
- `g='g3'` derivative builds must stay isolated from `wire_workspace`.
- `feast_workspace` reuse must remain pattern-validated.

**Step 3: Update design doc status**

Change the design status from:

```markdown
**Status:** Approved design, pending implementation planning
```

to:

```markdown
**Status:** Implementation planned in `docs/plans/2026-04-26-modern-fortran-migration-plan.md`
```

**Step 4: Commit**

```bash
git add README.md CLAUDE.md docs/plans/2026-04-26-modern-fortran-migration-design.md docs/plans/2026-04-26-modern-fortran-migration-plan.md
git commit -m "docs: document modern Fortran migration baseline"
```

---

## Deferred Work

Do not include these in the first implementation pass:

- Making `csr_matrix` components private.
- Replacing `wire_workspace` with polymorphic/OOP cache types.
- Moving FEAST, ARPACK, and dense LAPACK into a full abstract solver hierarchy.
- Adding Fortran stdlib as a required dependency.
- Adding fpm as an alternate build system.
- GPU offloading or coarrays.
- Replacing `selected_real_kind` with `iso_fortran_env` kinds before confirming `real128` behavior on all supported compilers.

These remain valid design directions, but the merged Hamiltonian refactor made cache correctness and hot-path transparency more important than broad encapsulation.

---

## Final Verification

After all tasks in this plan are complete, run:

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected:

- Clean build.
- `39/39` tests pass.
- No `forall` or `goto` remains in production `src/` files:

```bash
rg -n "\bforall\b|\bgoto\b" src
```

- Remaining `external ::` declarations are either in tests/prototypes or justified by a central interface:

```bash
rg -n "external ::" src tests
```

If any numerical regression changes, stop and inspect the relevant Hamiltonian/optics/g-factor output before continuing. Do not update golden data for this migration unless the change is proven to be formatting-only or a test tolerance issue.
