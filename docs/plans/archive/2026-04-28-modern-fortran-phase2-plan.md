# Modern Fortran Phase 2: Encapsulation + `do concurrent`

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `private` defaults to all modules, upgrade scalar `pure` functions to `elemental`, complete finalizer coverage, and introduce `do concurrent` on proven-independent hot-path loops — all without changing physics formulas, basis ordering, material parameters, or numerical outputs.

**Architecture:** Incremental encapsulation on top of the Phase 1 baseline (`-std=f2008`, no `forall`/`goto`/`dsqrt`, centralized BLAS interfaces, finalizers on `csr_matrix`/`wire_workspace`/`feast_workspace`). Each task is independently buildable and testable. No new dependencies (stdlib/fpm deferred to Phase 3).

**Tech Stack:** Fortran 2008, CMake/Ninja, Intel MKL FEAST/SpBLAS/PARDISO, FFTW3, pFUnit, shell/Python regression tests.

---

## Constraints From Phase 1

Do not undo these:

- `-std=f2008` enforcement via `CMAKE_Fortran_FLAGS` in root `CMakeLists.txt`
- No `forall`, `goto`, or `dsqrt` in production `src/`
- External BLAS/LAPACK/PARDISO interfaces centralized in `linalg.f90`
- `contiguous` attribute on proven hot-path arrays in `hamiltonian_wire.f90`, `confinement_init.f90`
- `wire_workspace` cache semantics: `g='g3'` derivative builds must not initialize or mutate it
- `feast_workspace` reuse must remain pattern-validated
- Explicit `*_free` routines remain public alongside finalizers

---

### Task 1: Add `private` Default to Hamiltonian Modules

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90`
- Modify: `src/physics/confinement_init.f90`
- Modify: `src/physics/hamiltonian_wire.f90`
- Test: targeted Hamiltonian and wire tests

These three modules already have explicit `public ::` statements listing their exports, so adding `private` is safe — anything not explicitly listed is already intended to be internal.

**Step 1: Add `private` after `implicit none` in each module**

In `src/physics/hamiltonianConstructor.f90` (line 11, after `implicit none`):

```fortran
  implicit none
  private
```

The existing `public ::` declarations on lines 13-15 keep the intended exports visible.

In `src/physics/confinement_init.f90` (line 7, after `implicit none`):

```fortran
  implicit none
  private
```

The existing `public ::` on lines 14-15 keep `confinementInitialization` and `confinementInitialization_2d` visible. The generic interface `confinementInitialization` on lines 9-12 is already matched by the public statement.

In `src/physics/hamiltonian_wire.f90` (line 10, after `implicit none`):

```fortran
  implicit none
  private
```

The existing `public ::` on lines 89-93 keep the intended exports visible.

**Step 2: Build**

```bash
cmake --build build
```

Expected: clean build. If any "no specific procedure for the interface" errors appear, check that the generic interface `build_velocity_matrices` (lines 94-97) and `confinementInitialization` are still accessible.

**Step 3: Run targeted tests**

```bash
ctest --test-dir build -R "test_hamiltonian|test_hamiltonian_2d" --output-on-failure
```

**Step 4: Run wire regressions**

```bash
ctest --test-dir build -R "regression_wire_gaas_rectangle|regression_wire_dense_sparse_consistency" --output-on-failure
```

**Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 src/physics/confinement_init.f90 src/physics/hamiltonian_wire.f90
git commit -m "refactor: add private default to Hamiltonian modules"
```

---

### Task 2: Add `private` Default to Core Math Modules

**Files:**
- Modify: `src/core/parameters.f90`
- Modify: `src/core/utils.f90`
- Modify: `src/math/finitedifferences.f90`
- Test: targeted unit tests

These modules do NOT currently have explicit `public ::` lists, so every module procedure is implicitly public. Adding `private` requires enumerating all symbols that external code uses.

**Step 1: Enumerate and add public exports for `parameters.f90`**

Read `src/core/parameters.f90` and find all module procedures. Add after `implicit none`:

```fortran
  implicit none
  private
  public :: paramDatabase
```

`paramDatabase` is the only procedure in this module that external code calls.

**Step 2: Enumerate and add public exports for `utils.f90`**

Read `src/core/utils.f90`. The public routines are `dnscsr_z_mkl`, `simpson`, and `simpson_real`. Add:

```fortran
  implicit none
  private
  public :: dnscsr_z_mkl, simpson, simpson_real
```

**Step 3: Enumerate and add public exports for `finitedifferences.f90`**

Read `src/math/finitedifferences.f90`. Add `private` and enumerate all public routines. The key exports are: `FDmatrixDense`, `Identity`, `buildDenseFiniteDiff`, `buildVariableCoeff`, `buildVariableCoeffStaggered`, and any other routines called from outside the module. Read the file to identify all external-callable routines.

```fortran
  implicit none
  private
  public :: FDmatrixDense, Identity, buildDenseFiniteDiff
  public :: buildVariableCoeff, buildVariableCoeffStaggered
```

Verify this list by searching for `use finitedifferences` across the codebase.

**Step 4: Build**

```bash
cmake --build build
```

Expected: clean build. Any unresolved symbol errors indicate a missing `public ::` entry.

**Step 5: Run targeted tests**

```bash
ctest --test-dir build -R "test_finitedifferences|test_utils|test_parameters|test_defs" --output-on-failure
```

**Step 6: Commit**

```bash
git add src/core/parameters.f90 src/core/utils.f90 src/math/finitedifferences.f90
git commit -m "refactor: add private default to core math modules"
```

---

### Task 3: Add `private` Default to `definitions` Module

**Files:**
- Modify: `src/core/defs.f90`
- Test: full build + all tests

This is the highest-impact module — it defines all shared types and constants used throughout the codebase. Adding `private` requires enumerating every type, constant, and function that external code references.

**Step 1: Audit all imports from `definitions`**

Search the codebase for all `use definitions` and `use definitions, only:` statements to build a complete list of external references. Key exports include:

- Kinds: `sp`, `dp`, `qp`, `iknd`
- Constants: `pi_dp`, `a0`, `hbar`, `const`, `c`, `m0`, `e`, `e0`, `kB_eV`, `tolerance`, `SQR3`, `SQR2`, `SQR2o3`, `SQR3o2`, `hbar2O2m0`, `NUM_CB_STATES`, `NUM_VB_STATES`, `ge_free`
- Types: `paramStruct`, `simulation_config`, `optics_config`, `sc_config`, `strain_config`, `exciton_config`, `scattering_config`, `doping_spec`, `spatial_grid`, `wire_geometry`, `wire_region`, `group`, `bir_pikus_blocks`, `strain_entry`, `config_validation_result`, `eigensolver_result`
- Functions: `kronij`, `grid_ngrid`, `init_spatial_grid`

**Step 2: Add `private` and enumerate all `public` exports**

After `implicit none` on line 2:

```fortran
  implicit none
  private

  ! Kinds
  public :: sp, dp, qp, iknd

  ! Physical constants
  public :: pi_sp, pi_dp, pi_qp, a0, hbar, const, c, m0, e, e0
  public :: kB_eV, tolerance, renormalization, hbar2O2m0
  public :: NUM_CB_STATES, NUM_VB_STATES, ge_free
  public :: SQR3, SQR2, SQR2o3, SQR3o2

  ! Types
  public :: paramStruct, simulation_config, optics_config, sc_config
  public :: strain_config, exciton_config, scattering_config
  public :: doping_spec, spatial_grid, wire_geometry, wire_region
  public :: group, bir_pikus_blocks, strain_entry
  public :: config_validation_result, eigensolver_result

  ! Functions
  public :: kronij, grid_ngrid, init_spatial_grid
```

Adjust the list based on the actual audit from Step 1. If in doubt about a symbol, keep it public.

**Step 3: Build**

```bash
cmake --build build
```

Expected: clean build. Any "has no implicit type" errors mean a symbol was missed from the public list.

**Step 4: Run full test suite**

```bash
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 5: Commit**

```bash
git add src/core/defs.f90
git commit -m "refactor: add private default to definitions module"
```

---

### Task 4: Upgrade Scalar `pure` Functions to `elemental`

**Files:**
- Modify: `src/core/defs.f90`
- Modify: `src/physics/charge_density.f90`
- Modify: `src/math/geometry.f90`
- Modify: `src/physics/strain_solver.f90`
- Test: targeted unit tests

Six functions are scalar `pure` with scalar arguments — safe to upgrade to `elemental`.

**Step 1: Upgrade `kronij` in `src/core/defs.f90` (line 305)**

Replace:

```fortran
  pure function kronij(i,j)
    integer, intent(in) :: i,j
    integer :: kronij
    kronij = merge(1, 0, i == j)
  end function
```

with:

```fortran
  elemental pure function kronij(i,j)
    integer, intent(in) :: i,j
    integer :: kronij
    kronij = merge(1, 0, i == j)
  end function
```

Note: `elemental pure function` is the correct F2008 syntax. In F2018+, `elemental` implies `pure` by default, but for `-std=f2008` we keep both keywords.

**Step 2: Upgrade `fermi_dirac` in `src/physics/charge_density.f90` (line 29)**

Replace `pure function fermi_dirac(...)` with `elemental pure function fermi_dirac(...)`.

**Step 3: Upgrade `flat_idx` in `src/math/geometry.f90` (line 44)**

Replace `pure function flat_idx(...)` with `elemental pure function flat_idx(...)`.

**Step 4: Upgrade `segment_circle_fraction` in `src/math/geometry.f90` (line 191)**

Replace `pure function segment_circle_fraction(...)` with `elemental pure function segment_circle_fraction(...)`.

**Step 5: Upgrade `wire_flat_idx` in `src/physics/strain_solver.f90` (line 826)**

Replace `pure function wire_flat_idx(...)` with `elemental pure function wire_flat_idx(...)`.

**Step 6: Upgrade `compute_bp_scalar` in `src/physics/strain_solver.f90` (line 870)**

Replace `pure function compute_bp_scalar(...)` with `elemental pure function compute_bp_scalar(...)`.

**Step 7: Build**

```bash
cmake --build build
```

Expected: clean build.

**Step 8: Run targeted tests**

```bash
ctest --test-dir build -R "test_defs|test_hamiltonian|test_hamiltonian_2d" --output-on-failure
```

**Step 9: Commit**

```bash
git add src/core/defs.f90 src/physics/charge_density.f90 src/math/geometry.f90 src/physics/strain_solver.f90
git commit -m "refactor: upgrade scalar pure functions to elemental"
```

---

### Task 5: Add Finalizers for Types With Existing `_free` Routines

**Files:**
- Modify: `src/math/eigensolver.f90`
- Modify: `src/physics/strain_solver.f90`
- Modify: `src/physics/hamiltonian_wire.f90`
- Test: targeted tests

Four types already have manual `*_free` routines but lack `final ::` bindings. Adding finalizers is a one-line addition to each type definition plus a thin wrapper subroutine.

**Step 1: Add finalizer for `eigensolver_result` in `src/math/eigensolver.f90`**

The type `eigensolver_result` (line 38) has `eigenvalues(:)` and `eigenvectors(:,:)` allocatable components. Manual free: `eigensolver_result_free` (line 722).

Add to the type definition:

```fortran
contains
  final :: eigensolver_result_finalize
```

Add the wrapper:

```fortran
subroutine eigensolver_result_finalize(er)
  type(eigensolver_result), intent(inout) :: er
  call eigensolver_result_free(er)
end subroutine
```

Keep `eigensolver_result_free` public.

**Step 2: Add finalizer for `strain_result` in `src/physics/strain_solver.f90`**

The type `strain_result` (line 37) has 4 allocatable arrays. Manual free: `strain_result_free` (line 856).

Same pattern: add `final :: strain_result_finalize` and wrapper delegating to `strain_result_free`.

**Step 3: Add finalizer for `bir_pikus_blocks` in `src/physics/strain_solver.f90`**

The type `bir_pikus_blocks` has 7 allocatable arrays. Manual free: `bir_pikus_blocks_free` (line 808).

Same pattern: add `final :: bir_pikus_blocks_finalize` and wrapper.

**Step 4: Add finalizer for `wire_coo_cache` in `src/physics/hamiltonian_wire.f90`**

The type `wire_coo_cache` (line 27) has `coo_to_csr(:)` allocatable. Manual free: `wire_coo_cache_free` (line 1357).

Add to the type:

```fortran
contains
  final :: wire_coo_cache_finalize
```

And the wrapper.

**Step 5: Build and run targeted tests**

```bash
cmake --build build
ctest --test-dir build -R "test_eigensolver|test_hamiltonian_2d" --output-on-failure
```

**Step 6: Commit**

```bash
git add src/math/eigensolver.f90 src/physics/strain_solver.f90 src/physics/hamiltonian_wire.f90
git commit -m "refactor: add finalizers for eigensolver_result, strain_result, bir_pikus_blocks, wire_coo_cache"
```

---

### Task 6: Add Finalizers for Types Without `_free` Routines

**Files:**
- Modify: `src/core/defs.f90`
- Test: full build + targeted tests

Three types have allocatable components but no `_free` routine at all: `wire_geometry` (1 array), `spatial_grid` (8 arrays), and `simulation_config` (9 arrays). These are all in `module definitions`.

**Step 1: Add finalizer for `wire_geometry`**

The type has one allocatable component: `verts(:,:)`. Add:

```fortran
contains
  final :: wire_geometry_finalize
```

```fortran
subroutine wire_geometry_finalize(wg)
  type(wire_geometry), intent(inout) :: wg
  if (allocated(wg%verts)) deallocate(wg%verts)
end subroutine
```

**Step 2: Add finalizer for `spatial_grid`**

The type has 8 allocatable arrays: `x`, `z`, `coords`, `material_id`, `cell_volume`, `face_fraction_x`, `face_fraction_y`, `ghost_map`. Add:

```fortran
contains
  final :: spatial_grid_finalize
```

```fortran
subroutine spatial_grid_finalize(sg)
  type(spatial_grid), intent(inout) :: sg
  if (allocated(sg%x)) deallocate(sg%x)
  if (allocated(sg%z)) deallocate(sg%z)
  if (allocated(sg%coords)) deallocate(sg%coords)
  if (allocated(sg%material_id)) deallocate(sg%material_id)
  if (allocated(sg%cell_volume)) deallocate(sg%cell_volume)
  if (allocated(sg%face_fraction_x)) deallocate(sg%face_fraction_x)
  if (allocated(sg%face_fraction_y)) deallocate(sg%face_fraction_y)
  if (allocated(sg%ghost_map)) deallocate(sg%ghost_map)
end subroutine
```

**Step 3: Add finalizer for `simulation_config`**

The type has 9+ allocatable arrays. Add:

```fortran
contains
  final :: simulation_config_finalize
```

```fortran
subroutine simulation_config_finalize(cfg)
  type(simulation_config), intent(inout) :: cfg
  if (allocated(cfg%startPos)) deallocate(cfg%startPos)
  if (allocated(cfg%endPos)) deallocate(cfg%endPos)
  if (allocated(cfg%z)) deallocate(cfg%z)
  if (allocated(cfg%intStartPos)) deallocate(cfg%intStartPos)
  if (allocated(cfg%intEndPos)) deallocate(cfg%intEndPos)
  if (allocated(cfg%materialN)) deallocate(cfg%materialN)
  if (allocated(cfg%params)) deallocate(cfg%params)
  if (allocated(cfg%doping)) deallocate(cfg%doping)
  if (allocated(cfg%regions)) deallocate(cfg%regions)
end subroutine
```

**Important:** `simulation_config` lifetime typically spans the entire program. The finalizer is a safety net for error paths or future use in loops. The module already has `private` from Task 3, so the finalizer is automatically accessible.

**Step 4: Build and run full tests**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 5: Commit**

```bash
git add src/core/defs.f90
git commit -m "refactor: add finalizers for wire_geometry, spatial_grid, simulation_config"
```

---

### Task 7: Add `do concurrent` to Velocity Matrix Construction

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90`
- Test: wire regressions

The velocity matrix construction loops (`build_velocity_matrices_2d` lines 561-579, `build_velocity_matrices_1d` lines 615-628) iterate over all CSR nonzeros independently. Each `k` writes to a unique `vel_x%values(k)` and `vel_y%values(k)` — no loop-carried dependencies. These are O(NNZ) hot paths called per-kz in the sweep loop.

**Step 1: Replace nested `do` in `build_velocity_matrices_2d`**

Replace (lines 561-579):

```fortran
      do row = 1, H_csr%nrows
        do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
          col = H_csr%colind(k)
          sp_row = mod(row - 1, Ngrid) + 1
          sp_col = mod(col - 1, Ngrid) + 1
          dx_diff = grid%coords(1, sp_row) - grid%coords(1, sp_col)
          dy_diff = grid%coords(2, sp_row) - grid%coords(2, sp_col)
          vel_x%values(k) = cmplx(0.0_dp, -dx_diff, kind=dp) * H_csr%values(k)
          vel_y%values(k) = cmplx(0.0_dp, -dy_diff, kind=dp) * H_csr%values(k)
        end do
      end do
```

with a flattened `do concurrent` over all nonzeros:

```fortran
      do concurrent (k = 1:H_csr%nnz)
        col = H_csr%colind(k)
        row = csr_row_of(k, H_csr%rowptr, H_csr%nrows)
        sp_row = mod(row - 1, Ngrid) + 1
        sp_col = mod(col - 1, Ngrid) + 1
        dx_diff = grid%coords(1, sp_row) - grid%coords(1, sp_col)
        dy_diff = grid%coords(2, sp_row) - grid%coords(2, sp_col)
        vel_x%values(k) = cmplx(0.0_dp, -dx_diff, kind=dp) * H_csr%values(k)
        vel_y%values(k) = cmplx(0.0_dp, -dy_diff, kind=dp) * H_csr%values(k)
      end do
```

This requires a helper function to recover the row index from a flat `k` position. Add a `pure function` in the module:

```fortran
  pure function csr_row_of(k, rowptr, nrows) result(row)
    integer, intent(in) :: k, nrows
    integer, intent(in) :: rowptr(nrows + 1)
    integer :: row, lo, hi, mid
    lo = 1
    hi = nrows
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (rowptr(mid + 1) <= k) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    row = lo
  end function
```

This is a binary search on the CSR row pointer — O(log N) per call, called NNZ times total, but the compiler can vectorize the outer `do concurrent` across all `k` values.

**Alternative (simpler):** Keep the nested `do` but make the outer loop `do concurrent`:

```fortran
      do concurrent (row = 1:H_csr%nrows)
        do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
          col = H_csr%colind(k)
          sp_row = mod(row - 1, Ngrid) + 1
          sp_col = mod(col - 1, Ngrid) + 1
          dx_diff = grid%coords(1, sp_row) - grid%coords(1, sp_col)
          dy_diff = grid%coords(2, sp_row) - grid%coords(2, sp_col)
          vel_x%values(k) = cmplx(0.0_dp, -dx_diff, kind=dp) * H_csr%values(k)
          vel_y%values(k) = cmplx(0.0_dp, -dy_diff, kind=dp) * H_csr%values(k)
        end do
      end do
```

This is simpler and avoids the binary search. Each `row` writes to disjoint ranges of `vel_x%values` and `vel_y%values` (CSR structure guarantees no overlap). Use this approach.

**Step 2: Replace nested `do` in `build_velocity_matrices_1d`**

Same pattern for lines 615-628:

```fortran
      do concurrent (row = 1:H_csr%nrows)
        do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
          col = H_csr%colind(k)
          sp_row = mod(row - 1, Ngrid) + 1
          sp_col = mod(col - 1, Ngrid) + 1
          z_diff = grid%z(sp_row) - grid%z(sp_col)
          vel(3)%values(k) = cmplx(0.0_dp, -z_diff, kind=dp) * H_csr%values(k)
        end do
      end do
```

**Step 3: Build**

```bash
cmake --build build
```

Expected: clean build.

**Step 4: Run wire regressions to verify numerical equivalence**

```bash
ctest --test-dir build -R "regression_wire_gaas_rectangle|regression_wire_dense_sparse_consistency|regression_wire_optical_selection" --output-on-failure
```

Expected: all pass — `do concurrent` must produce identical numerical results to sequential `do`.

**Step 5: Commit**

```bash
git add src/physics/hamiltonian_wire.f90
git commit -m "perf: use do concurrent for velocity matrix construction"
```

---

### Task 8: Add `do concurrent` to Optics Finalization and Kpterms Init

**Files:**
- Modify: `src/physics/optical_spectra.f90`
- Modify: `src/physics/confinement_init.f90`
- Test: targeted regression tests

Two lower-impact but safe `do concurrent` targets: optics energy-grid prefactor application and kpterms diagonal initialization.

**Step 1: Replace optics finalization loop**

In `src/physics/optical_spectra.f90`, in `optics_finalize` (around line 357), replace:

```fortran
do ie = 1, nE
  if (E_grid(ie) > 0.0_dp) then
    prefactor_E = ...
  else
    prefactor_E = 0.0_dp
  end if
  alpha_te(ie) = prefactor_E * alpha_te(ie) * AA_TO_CM
  alpha_tm(ie) = prefactor_E * alpha_tm(ie) * AA_TO_CM
  alpha_isbt(ie) = prefactor_E * alpha_isbt(ie) * AA_TO_CM
end do
```

with:

```fortran
do concurrent (ie = 1:nE)
  if (E_grid(ie) > 0.0_dp) then
    prefactor_E = ...
  else
    prefactor_E = 0.0_dp
  end if
  alpha_te(ie) = prefactor_E * alpha_te(ie) * AA_TO_CM
  alpha_tm(ie) = prefactor_E * alpha_tm(ie) * AA_TO_CM
  alpha_isbt(ie) = prefactor_E * alpha_isbt(ie) * AA_TO_CM
end do
```

Each `ie` writes to unique array positions. No loop-carried dependencies.

**Step 2: Replace kpterms diagonal init loop**

In `src/physics/confinement_init.f90`, in `confinementInitialization_raw` (around line 135), replace:

```fortran
do ii = 1, N
  kpterms(ii,ii,1) = kptermsProfile(ii,1)
  kpterms(ii,ii,2) = kptermsProfile(ii,2)
  kpterms(ii,ii,3) = kptermsProfile(ii,3)
  kpterms(ii,ii,4) = kptermsProfile(ii,5)
  kpterms(ii,ii,10) = kptermsProfile(ii,4)
end do
```

with:

```fortran
do concurrent (ii = 1:N)
  kpterms(ii,ii,1) = kptermsProfile(ii,1)
  kpterms(ii,ii,2) = kptermsProfile(ii,2)
  kpterms(ii,ii,3) = kptermsProfile(ii,3)
  kpterms(ii,ii,4) = kptermsProfile(ii,5)
  kpterms(ii,ii,10) = kptermsProfile(ii,4)
end do
```

Each `ii` writes to unique diagonal positions `(ii,ii,*)`.

**Step 3: Build and run targeted tests**

```bash
cmake --build build
ctest --test-dir build -R "test_hamiltonian|regression_qw_absorption_polarization|regression_qw_optics_commutator" --output-on-failure
```

**Step 4: Commit**

```bash
git add src/physics/optical_spectra.f90 src/physics/confinement_init.f90
git commit -m "perf: use do concurrent for optics finalization and kpterms init"
```

---

### Task 9: Document Phase 2 Modernization

**Files:**
- Modify: `CLAUDE.md`
- Modify: `README.md`
- Modify: `docs/plans/2026-04-26-modern-fortran-migration-design.md`
- Test: full test suite

**Step 1: Update CLAUDE.md code conventions**

In the Code Conventions section, add/update:

- All modules now use `private` default with explicit `public` exports.
- Scalar `pure` functions upgraded to `elemental pure`: `kronij`, `fermi_dirac`, `flat_idx`, `wire_flat_idx`, `compute_bp_scalar`, `segment_circle_fraction`.
- All types with allocatable components now have finalizers (delegating to `*_free` routines where they exist).
- `do concurrent` used on proven-independent loops: velocity matrix construction, optics finalization, kpterms diagonal init.
- When adding new modules, use `private` default and enumerate `public ::` exports.
- When adding new scalar `pure` functions, use `elemental pure` by default.

**Step 2: Update README.md**

Update the project description to note modern Fortran 2008 features in use.

**Step 3: Update design doc status**

Update status to reflect Phase 2 completion.

**Step 4: Run full test suite**

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass.

**Step 5: Verify no regressions from Phase 2 patterns**

```bash
grep -rn "forall\|goto\|dsqrt\|external ::" src/ --include="*.f90"
```

Expected: no matches (Phase 1 guarantees preserved).

**Step 6: Commit**

```bash
git add CLAUDE.md README.md docs/plans/2026-04-26-modern-fortran-migration-design.md
git commit -m "docs: document Phase 2 modernization baseline"
```

---

## Deferred Work

These items from the design spec are intentionally excluded from Phase 2:

- **`csr_matrix` components made private** — 84 external access sites across eigensolver.f90 and hamiltonian_wire.f90. Requires adding accessor functions for all read/write patterns. High risk, high effort. Should be a dedicated task after profiling confirms no performance regression.
- **Polymorphic eigensolver dispatch** — Major refactor of the eigensolver interface. The current `select case` approach works. Deferred until GPU offloading needs it.
- **`spatial_grid` encapsulation with type-bound methods** — Would require changing Hamiltonian construction call sites. Too invasive for this round.
- **`iso_c_binding` for MKL interfaces** — Current hand-written interfaces in `linalg.f90` work. Deferred to Phase 3.
- **Fortran stdlib adoption** — Adds external dependency. Deferred to Phase 3.
- **fpm alternative build** — Deferred to Phase 3.
- **`iso_fortran_env` kinds** — Deferred to Phase 4 (needs `real128` compatibility testing).
- **GPU offloading, coarrays, F2023 features** — Deferred to Phase 4 (exploratory).

---

## Final Verification

After all tasks:

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass, no numerical regressions.

Verify Phase 2 patterns:

```bash
# All modules should have private default (except possibly geometry.f90 which already has it)
grep -L "private" src/**/*.f90

# All elemental functions confirmed
grep -n "elemental" src/**/*.f90

# do concurrent locations confirmed
grep -n "do concurrent" src/**/*.f90
```
