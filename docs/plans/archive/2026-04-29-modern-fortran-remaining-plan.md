# Modern Fortran Migration: Remaining Work

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the remaining Phase 2 items (spatial_grid accessors, polymorphic eigensolver, config validation) and deliver Phase 3 (contiguous attributes, iso_c_binding for MKL, stdlib, fpm).

**Architecture:** Incremental on top of the Phase 1+2 baseline (`-std=f2008`, `private` everywhere, `elemental` upgrades, finalizers, `do concurrent` in hot loops). Each task is independently buildable and testable. No physics formula, basis ordering, material parameter, or numerical output changes.

**Tech Stack:** Fortran 2008, CMake/Ninja, Intel MKL FEAST/SpBLAS/PARDISO, FFTW3, pFUnit, shell/Python regression tests. stdlib v0.8.0 and fpm v0.13.0 for ecosystem tasks.

---

## Constraints From Completed Phases

Do not undo these:

- `-std=f2008` enforcement via `CMAKE_Fortran_FLAGS` in root `CMakeLists.txt`
- No `forall`, `goto`, or `dsqrt` in production `src/`
- `private` default in all 22 modules with explicit `public` exports
- `elemental` on `kronij`, `fermi_dirac`, `flat_idx`, `wire_flat_idx`, `compute_bp_scalar`, `segment_circle_fraction`
- Finalizers on `csr_matrix`, `wire_workspace`, `feast_workspace`, `sc_workspace`, `bir_pikus_blocks`, `wire_geometry`, `spatial_grid`, `simulation_config`, `eigensolver_result`, `strain_result`, `wire_coo_cache`
- `contiguous` already on hot-path arrays in `confinement_init.f90`, `hamiltonian_wire.f90`
- `do concurrent` already in `build_velocity_matrices_1d`, `build_velocity_matrices_2d`, optics finalization, kpterms init
- `wire_workspace` cache semantics: `g='g3'` derivative builds must not initialize or mutate it
- `feast_workspace` reuse must remain pattern-validated
- Explicit `*_free` routines remain public alongside finalizers

---

### Task 1: Remove .bak Files and Resolve Dead Code

**Files:**
- Delete: `src/physics/poisson.f90.bak`
- Delete: `src/physics/strain_solver.f90.bak`
- Modify: `src/core/defs.f90:324-327` (remove `config_validation_result`)
- Modify: `src/core/defs.f90:24` (remove from public list)
- Test: `ctest --test-dir build`

- [ ] **Step 1: Delete .bak files**

```bash
rm src/physics/poisson.f90.bak src/physics/strain_solver.f90.bak
```

- [ ] **Step 2: Remove unused `config_validation_result` type**

In `src/core/defs.f90`, remove lines 324-327 (the entire type definition):

```fortran
! DELETE THIS:
type :: config_validation_result
  logical           :: ok      = .true.
  character(len=256) :: message = ''
end type config_validation_result
```

Also remove it from the public list (line 24):

```fortran
! CHANGE FROM:
public :: simulation_config, config_validation_result, group
! TO:
public :: simulation_config, group
```

- [ ] **Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "chore: remove .bak files and unused config_validation_result type"
```

---

### Task 2: Add `spatial_grid` Type-Bound Accessors

**Files:**
- Modify: `src/core/defs.f90` (type definition, `spatial_grid` type)
- Modify: `src/core/defs.f90` (move `grid_ngrid` to type-bound)
- Modify: all call sites replacing `grid%nx * grid%ny` with `grid%npoints()`
- Test: `ctest --test-dir build`

The `spatial_grid` type currently has only a finalizer. Add `npoints` as a type-bound procedure to eliminate the duplicated `grid%nx * grid%ny` pattern across 6+ files. Do NOT make components `private` — too many call sites access them directly and that refactoring is deferred.

- [ ] **Step 1: Add `npoints` type-bound procedure to `spatial_grid`**

In `src/core/defs.f90`, modify the `spatial_grid` type's `contains` block:

```fortran
type :: spatial_grid
  ! ... existing components unchanged ...
contains
  procedure :: npoints => spatial_grid_npoints
  final :: spatial_grid_finalize
end type spatial_grid
```

Add the implementation in the module's `contains` section (after `grid_ngrid`):

```fortran
elemental function spatial_grid_npoints(self) result(n)
  class(spatial_grid), intent(in) :: self
  integer :: n
  n = self%nx * self%ny
end function spatial_grid_npoints
```

Keep the standalone `grid_ngrid` function (it's a `pure function` taking `type(spatial_grid)`, not `class` — it's used in some contexts where the type-bound version can't be called).

- [ ] **Step 2: Replace `grid%nx * grid%ny` call sites**

Replace the pattern at these locations:

| File | Line | Change |
|---|---|---|
| `src/physics/gfactor_functions.f90` | 693 | `ngrid = grid%nx * grid%ny` → `ngrid = grid%npoints()` |
| `src/physics/gfactor_functions.f90` | 807 | `ngrid = cfg%grid%nx * cfg%grid%ny` → `ngrid = cfg%grid%npoints()` |
| `src/io/outputFunctions.f90` | 303 | `Ngrid = grid%nx * grid%ny` → `Ngrid = grid%npoints()` |
| `src/physics/hamiltonian_wire.f90` | 560 | `Ngrid = grid%nx * grid%ny` → `Ngrid = grid%npoints()` |
| `src/physics/hamiltonian_wire.f90` | 613 | `Ngrid = grid%ny` → leave as-is (this is just `grid%ny`, not the product) |

Also search `geometry.f90` for `ngrid = ... * ...` patterns and verify they use `grid%nx * grid%ny` before replacing.

- [ ] **Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass with no numerical changes.

- [ ] **Step 4: Commit**

```bash
git add src/core/defs.f90 src/physics/gfactor_functions.f90 src/io/outputFunctions.f90 src/physics/hamiltonian_wire.f90
git commit -m "refactor: add spatial_grid%npoints() type-bound accessor"
```

---

### Task 3: Add `contiguous` to Remaining Assumed-Shape Arrays

**Files:**
- Modify: `src/physics/optical_spectra.f90` (6 arrays)
- Modify: `src/physics/sc_loop.f90` (3 arrays)
- Modify: `src/physics/scattering.f90` (3 arrays)
- Modify: `src/physics/exciton.f90` (2 arrays)
- Modify: `src/physics/hamiltonianConstructor.f90` (1 array)
- Modify: `src/math/eigensolver.f90` (1 array)
- Modify: `src/io/outputFunctions.f90` (5 arrays)
- Modify: `src/apps/main.f90` (2 arrays)
- Test: `ctest --test-dir build`

This is a mechanical change — add `contiguous` to assumed-shape dummy arguments that always receive whole allocated arrays. Do NOT add it to `optional` arguments (they may be absent) or `allocatable` arguments (already guaranteed contiguous).

- [ ] **Step 1: Add `contiguous` to `optical_spectra.f90`**

In `src/physics/optical_spectra.f90`, add `contiguous` to these assumed-shape array arguments:

| Line | Subroutine | Change |
|---|---|---|
| 110 | `optics_accumulate` | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 218 | `optics_accumulate_spontaneous` | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 583 | `z_dipole` | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 667 | `compute_intersubband_transitions` | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 668 | `compute_intersubband_transitions` | `real(kind=dp), intent(in), contiguous :: z_grid(:)` |
| 731 | `compute_isbt_absorption` | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 902 | `compute_gain_qw` | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |

- [ ] **Step 2: Add `contiguous` to `sc_loop.f90`**

In `src/physics/sc_loop.f90`, lines 538-540 (`apply_potential_to_profile_2d`):

```fortran
real(kind=dp), intent(inout), contiguous :: profile_2d(:,:)
real(kind=dp), intent(in), contiguous :: profile_2d_base(:,:)
real(kind=dp), intent(in), contiguous :: phi(:,:)
```

- [ ] **Step 3: Add `contiguous` to `scattering.f90`**

In `src/physics/scattering.f90`:

| Line | Change |
|---|---|
| 52 | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 205 | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 404 | `real(kind=dp), intent(in), contiguous :: rate_em(:,:), rate_ab(:,:)` |

- [ ] **Step 4: Add `contiguous` to `exciton.f90`**

In `src/physics/exciton.f90`:

| Line | Change |
|---|---|
| 48 | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |
| 150 | `complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)` |

- [ ] **Step 5: Add `contiguous` to remaining files**

`src/physics/hamiltonianConstructor.f90` line 467:
```fortran
complex(kind=dp), intent(inout), contiguous :: HT(:,:)
```

`src/math/eigensolver.f90` line 791 (`sort_eigenpairs_ascending`):
```fortran
complex(kind=dp), intent(inout), contiguous :: evecs(:,:)
```

`src/io/outputFunctions.f90`:
| Line | Change |
|---|---|
| 33 | `complex(kind=dp), intent(in), contiguous :: A(:,:)` |
| 149 | `real(kind=dp), intent(in), contiguous :: eig(:,:)` |
| 258 | `complex(kind=dp), intent(in), contiguous :: eigenvectors(:,:)` |
| 289 | `complex(kind=dp), intent(in), contiguous :: eigenvectors(:,:)` |
| 357 | `complex(kind=dp), intent(in), contiguous :: eigenvectors(:,:)` |

`src/apps/main.f90` lines 774-776 (`reortho_eigenvectors`):
```fortran
complex(kind=dp), intent(in), contiguous :: prev_evec(:,:)
complex(kind=dp), intent(inout), contiguous :: curr_evec(:,:)
```

- [ ] **Step 6: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass. No numerical changes — `contiguous` is a compiler hint, not a behavior change.

- [ ] **Step 7: Commit**

```bash
git add src/physics/optical_spectra.f90 src/physics/sc_loop.f90 src/physics/scattering.f90 src/physics/exciton.f90 src/physics/hamiltonianConstructor.f90 src/math/eigensolver.f90 src/io/outputFunctions.f90 src/apps/main.f90
git commit -m "perf: add contiguous attribute to remaining assumed-shape hot-path arrays"
```

---

### Task 4: `iso_c_binding` for MKL C APIs

**Files:**
- Modify: `src/math/linalg.f90` (interface blocks for PARDISO, FEAST, mkl_set_num_threads_local)
- Test: `ctest --test-dir build`

The MKL routines `pardiso`, `zfeast_hcsrev`, `feastinit`, and `mkl_set_num_threads_local` are C APIs. Wrap them with `iso_c_binding` for type safety and portability. ARPACK routines are Fortran APIs — leave them as-is.

- [ ] **Step 1: Add `iso_c_binding` import to `linalg.f90`**

At the top of `linalg.f90`, add to the `use` statements:

```fortran
use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_double, c_double_complex, c_char
```

- [ ] **Step 2: Rewrite `mkl_set_num_threads_local` interface**

Replace the existing interface (lines ~95-100):

```fortran
interface
  function mkl_set_num_threads_local(nt) result(previous) bind(C, name="MKL_Set_Num_Threads_Local")
    import :: c_int
    integer(c_int), value :: nt
    integer(c_int) :: previous
  end function mkl_set_num_threads_local
end interface
```

- [ ] **Step 3: Rewrite PARDISO interface**

Replace the existing PARDISO interface (lines ~112-128). Use `c_int` for MKL_INT (LP64 interface), `c_int64_t` for handle arrays:

```fortran
interface
  subroutine pardiso_c(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, &
                       nrhs, iparm, msglvl, b, x, error) bind(C, name="PARDISO")
    import :: c_int, c_int64_t, c_double_complex
    integer(c_int64_t), intent(inout) :: pt(64)
    integer(c_int), value :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
    complex(c_double_complex), intent(in) :: a(*)
    integer(c_int), intent(in) :: ia(*), ja(*), perm(*)
    integer(c_int), intent(inout) :: iparm(64)
    complex(c_double_complex), intent(inout) :: b(*), x(*)
    integer(c_int), intent(out) :: error
  end subroutine pardiso_c
end interface
```

Then update the call site(s) in `eigensolver.f90` from `call pardiso(...)` to `call pardiso_c(...)`. The argument types remain the same because `dp` matches `c_double` and default `integer` matches `c_int` on LP64.

- [ ] **Step 4: Rewrite FEAST interface**

Replace the existing `zfeast_hcsrev` interface (lines ~131-155):

```fortran
interface
  subroutine zfeast_hcsrev(uplo, n, a, ia, ja, fpm, epsout, loop, &
                           emin, emax, m0, e, x, m, res, info) bind(C, name="zfeast_hcsrev")
    import :: c_int, c_double, c_double_complex, c_char
    character(c_char), intent(in) :: uplo
    integer(c_int), value :: n, m0
    complex(c_double_complex), intent(in) :: a(*)
    integer(c_int), intent(in) :: ia(*), ja(*)
    integer(c_int), intent(inout) :: fpm(128)
    real(c_double), intent(out) :: epsout
    integer(c_int), intent(out) :: loop, m, info
    real(c_double), intent(in) :: emin, emax
    real(c_double), intent(inout) :: e(m0), res(m0)
    complex(c_double_complex), intent(inout) :: x(n, m0)
  end subroutine zfeast_hcsrev
end interface
```

Similarly update `feastinit` interface with `bind(C)`.

- [ ] **Step 5: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass. If there are type mismatch errors, check that `dp` == `c_double` and default integer == `c_int` (true on LP64 x86_64).

- [ ] **Step 6: Commit**

```bash
git add src/math/linalg.f90 src/math/eigensolver.f90
git commit -m "refactor: use iso_c_binding for MKL C APIs (PARDISO, FEAST, mkl_set_num_threads_local)"
```

---

### Task 5: Polymorphic Eigensolver Dispatch

**Files:**
- Modify: `src/math/eigensolver.f90` (add abstract base + concrete types)
- Modify: call sites in `src/apps/main.f90`, `src/apps/main_gfactor.f90`, `src/apps/main_optics.f90`, `src/physics/sc_loop.f90`
- Test: `ctest --test-dir build`

This is the highest-risk task. Replace the `select case (trim(config%method))` dispatch with an abstract type hierarchy. Each solver type caches its own workspace.

- [ ] **Step 1: Define abstract base type**

In `src/math/eigensolver.f90`, add before the existing types:

```fortran
type, abstract :: eigensolver_base
contains
  procedure(solve_evp_interface), deferred :: solve
end type eigensolver_base

abstract interface
  subroutine solve_evp_interface(self, H_csr, config, result)
    import :: eigensolver_base, csr_matrix, eigensolver_config, eigensolver_result
    class(eigensolver_base), intent(inout) :: self
    type(csr_matrix), intent(in) :: H_csr
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result
  end subroutine solve_evp_interface
end interface
```

- [ ] **Step 2: Define concrete solver types**

```fortran
type, extends(eigensolver_base) :: feast_solver_t
  type(feast_workspace) :: ws
contains
  procedure :: solve => feast_solve
end type feast_solver_t

type, extends(eigensolver_base) :: dense_lapack_solver_t
contains
  procedure :: solve => dense_lapack_solve
end type dense_lapack_solver_t
```

Add the implementation subroutines that wrap the existing `solve_feast` and `solve_dense_lapack`:

```fortran
subroutine feast_solve(self, H_csr, config, result)
  class(feast_solver_t), intent(inout) :: self
  type(csr_matrix), intent(in) :: H_csr
  type(eigensolver_config), intent(in) :: config
  type(eigensolver_result), intent(out) :: result
  call solve_feast(H_csr, config, result, fw=self%ws)
end subroutine feast_solve

subroutine dense_lapack_solve(self, H_csr, config, result)
  class(dense_lapack_solver_t), intent(inout) :: self
  type(csr_matrix), intent(in) :: H_csr
  type(eigensolver_config), intent(in) :: config
  type(eigensolver_result), intent(out) :: result
  call solve_dense_lapack(H_csr, config, result)
end subroutine dense_lapack_solve
```

Export `eigensolver_base`, `feast_solver_t`, `dense_lapack_solver_t` in the public list.

- [ ] **Step 3: Add factory function**

```fortran
function make_eigensolver(config) result(solver)
  class(eigensolver_base), allocatable :: solver
  type(eigensolver_config), intent(in) :: config
  select case (trim(config%method))
  case ('FEAST')
    allocate(feast_solver_t :: solver)
  case ('DENSE')
    allocate(dense_lapack_solver_t :: solver)
  case default
    allocate(dense_lapack_solver_t :: solver)
  end select
end function make_eigensolver
```

Export `make_eigensolver` in the public list.

- [ ] **Step 4: Update call sites**

In `src/apps/main.f90`, `main_gfactor.f90`, `main_optics.f90`, and `sc_loop.f90` — wherever `solve_sparse_evp` or `solve_feast` is called directly, replace with:

```fortran
class(eigensolver_base), allocatable :: solver
solver = make_eigensolver(eigconfig)
call solver%solve(H_csr, eigconfig, eigresult)
```

This is the most invasive change. Do one call site at a time, building and testing between each.

- [ ] **Step 5: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass with no numerical changes.

- [ ] **Step 6: Commit**

```bash
git add src/math/eigensolver.f90 src/apps/main.f90 src/apps/main_gfactor.f90 src/apps/main_optics.f90 src/physics/sc_loop.f90
git commit -m "refactor: polymorphic eigensolver dispatch via abstract type hierarchy"
```

---

### Task 6: `simulation_config` Validate-Once

**Files:**
- Modify: `src/core/defs.f90` (add type-bound `validate` to `simulation_config`)
- Modify: `src/io/input_parser.f90` (move validation logic to type-bound method)
- Test: `ctest --test-dir build`

The validation logic currently lives in `input_parser.f90` as `validate_simulation_config`. Move it to a type-bound procedure so the config object can validate itself.

- [ ] **Step 1: Add `validate` type-bound procedure to `simulation_config`**

In `src/core/defs.f90`, add to the `simulation_config` type's `contains`:

```fortran
type :: simulation_config
  ! ... existing components ...
contains
  procedure :: validate => simulation_config_validate
  final :: simulation_config_finalize
end type simulation_config
```

- [ ] **Step 2: Move validation to type-bound implementation**

In `defs.f90`, add the implementation (move the logic from `input_parser.f90:37-95`):

```fortran
subroutine simulation_config_validate(self)
  class(simulation_config), intent(in) :: self
  associate(cfg => self)
    if (cfg%confinement > 0 .and. cfg%fdStep <= 0) error stop 'Error: fdStep must be > 0'
    if (cfg%numcb < 0) error stop 'Error: numcb must be >= 0'
    if (cfg%numvb < 0) error stop 'Error: numvb must be >= 0'
    if (cfg%evnum /= cfg%numcb + cfg%numvb) error stop 'Error: evnum must equal numcb+numvb'
    if (cfg%numLayers < 1) error stop 'Error: numLayers must be >= 1'
    if (cfg%confinement < 0 .or. cfg%confinement > 2) error stop 'Error: confinement must be 0, 1, or 2'
    ! ... remaining checks from input_parser ...
  end associate
end subroutine simulation_config_validate
```

- [ ] **Step 3: Update call site in `input_parser.f90`**

Replace `call validate_simulation_config(cfg)` (line ~863) with `call cfg%validate()`.

Remove the standalone `validate_simulation_config` subroutine from `input_parser.f90` and its `public` declaration.

- [ ] **Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90 src/io/input_parser.f90
git commit -m "refactor: move simulation_config validation to type-bound method"
```

---

### Task 7: stdlib Adoption (Selective)

**Files:**
- Modify: `CMakeLists.txt` (add `find_package(fortran_stdlib)`)
- Modify: `src/core/defs.f90` (use `stdlib_kinds` for `dp`, `sp`, `qp`)
- Modify: `src/core/utils.f90` (use `stdlib_quadrature` for Simpson integration if API compatible)
- Test: `ctest --test-dir build`

This task adds stdlib as a dependency and adopts it where it replaces hand-rolled code. MKL SpBLAS, FEAST, and PARDISO stay as direct calls.

- [ ] **Step 1: Install and configure stdlib**

Install stdlib v0.8.0:
```bash
pip install fypp
git clone https://github.com/fortran-lang/stdlib.git /tmp/stdlib
cd /tmp/stdlib && mkdir build && cd build
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=$HOME/.local \
      -DCMAKE_MAXIMUM_RANK=4 -DWITH_QP=ON ..
cmake --build . && cmake --install .
```

In root `CMakeLists.txt`, add after MKL find_package:
```cmake
find_package(fortran_stdlib CONFIG QUIET)
```

In `src/CMakeLists.txt`, link to stdlib:
```cmake
target_link_libraries(8bandkp_common PUBLIC fortran_stdlib::fortran_stdlib)
```

- [ ] **Step 2: Verify build**

Run: `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl && cmake --build build`
Expected: Clean build with stdlib linked.

- [ ] **Step 3: Commit**

```bash
git add CMakeLists.txt src/CMakeLists.txt
git commit -m "build: add fortran_stdlib as optional dependency"
```

---

### Task 8: fpm.toml (Alternative Build)

**Files:**
- Create: `fpm.toml`

Add `fpm.toml` alongside `CMakeLists.txt` for developer convenience. fpm v0.13.0 supports profiles and external module flags for MKL integration.

- [ ] **Step 1: Create `fpm.toml`**

```toml
name = "8bandkp-fdm"
version = "0.1.0"
license = "GPL-3.0"
author = "Tiago de Campos"

[dependencies]
stdlib = "*"

[build]
auto-executables = false
auto-tests = false

[[executable]]
name = "bandStructure"
source-dir = "src/apps"
main = "main.f90"

[[executable]]
name = "gfactorCalculation"
source-dir = "src/apps"
main = "main_gfactor.f90"

[[executable]]
name = "opticalProperties"
source-dir = "src/apps"
main = "main_optics.f90"

[library]
source-dir = "src"

[install]
library = false

[fortran]
implicit-typing = false
source-form = "free"
```

- [ ] **Step 2: Verify fpm can parse it**

Run: `fpm build --compiler gfortran --flag "-std=f2008 -I${MKLROOT}/include -L${MKLROOT}/lib -lmkl_rt -lfftw3"`

This will likely fail to fully link without the right MKL flags — the goal is just to verify fpm can parse the manifest and start compilation. CMake remains the primary build system.

- [ ] **Step 3: Commit**

```bash
git add fpm.toml
git commit -m "build: add fpm.toml as alternative build system"
```

---

### Task 9: Update Documentation

**Files:**
- Modify: `README.md`
- Modify: `CLAUDE.md`

Every phase requires updating documentation to reflect changes.

- [ ] **Step 1: Update CLAUDE.md**

Update the "Code Conventions" section to add:
- `contiguous` attribute on all assumed-shape hot-path arrays
- `spatial_grid%npoints()` accessor preferred over `grid%nx * grid%ny`
- `simulation_config%validate()` for validation instead of standalone subroutine
- `iso_c_binding` for MKL C APIs (PARDISO, FEAST, `mkl_set_num_threads_local`)
- Polymorphic eigensolver: `make_eigensolver(config)` factory, `solver%solve(...)` dispatch
- stdlib as optional dependency (`find_package(fortran_stdlib)`)

- [ ] **Step 2: Update README.md**

Add stdlib to prerequisites section. Add fpm as alternative build option. Update build commands section.

- [ ] **Step 3: Commit**

```bash
git add README.md CLAUDE.md
git commit -m "docs: update README and CLAUDE.md for Phase 2+3 modernization"
```

---

## Execution Order

Execute in task order. Tasks 1-3 are independent and could run in parallel. Task 4 (iso_c_binding) should precede Task 5 (polymorphic eigensolver) since the eigensolver touches the same interface blocks. Task 6 (config validation) is independent. Tasks 7-8 (stdlib, fpm) are independent of each other and of Tasks 4-6. Task 9 (docs) runs last.

```
Task 1 (cleanup) ──┐
Task 2 (grid)     ──┤
Task 3 (contiguous)┤── Task 4 (iso_c_binding) ── Task 5 (eigensolver)
Task 6 (config)   ──┤
Task 7 (stdlib)   ──┤
Task 8 (fpm)      ──┘
                                              ── Task 9 (docs)
```

## Risk Assessment

| Task | Risk | Mitigation |
|---|---|---|
| 1. Cleanup | Very low | Deleting unused files/types |
| 2. Grid accessors | Low | No behavior change, just redirect |
| 3. Contiguous | Very low | Compiler hint only, no behavior change |
| 4. iso_c_binding | Medium | Must verify `c_int` matches LP64 integer on this platform |
| 5. Eigensolver | High | Architecture change — do one call site at a time, test between each |
| 6. Config validation | Low | Moving existing logic to type-bound |
| 7. stdlib | Low | Optional dependency, existing code unchanged |
| 8. fpm | Low | Additive, CMake stays primary |
| 9. Docs | None | Documentation only |
