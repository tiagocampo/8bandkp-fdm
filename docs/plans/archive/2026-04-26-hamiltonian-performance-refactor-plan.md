# Hamiltonian Performance Refactor Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Eliminate per-k-point allocation overhead in wire kz-sweeps by pre-allocating CSR workspaces, caching FEAST upper-triangle extraction, replacing the unrolled strain COO loop with a data-driven table, and decomposing hamiltonianConstructor.f90 into focused modules.

**Architecture:** A `wire_workspace` type pre-allocates all intermediate CSR matrices and COO buffers. On the first kz-point, kp-terms are built normally via csr_add and the resulting CSR structures are stashed. On subsequent kz-points, values are computed directly from kpterms_2d using element-wise linear combinations with kz-dependent coefficients, skipping all csr_add/csr_scale allocation. A `feast_workspace` type caches the upper-triangle CSR structure across FEAST calls.

**Tech Stack:** Fortran 90, Intel MKL FEAST, pFUnit for unit tests, CMake/Ninja build system

---

## Phase 1: Wire Workspace + Compute-Values-Direct

### Task 1: Add `wire_workspace` type and init/free

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:31-44` (extend wire_coo_cache area)
- Test: `tests/unit/test_hamiltonian_2d.pf` (add workspace tests)

**Context:** The existing `wire_coo_cache` type at line 31 only stores the COO-to-CSR mapping. Extend it into a richer `wire_workspace` that holds pre-allocated CSR blocks and COO buffers.

**Step 1: Define the workspace type**

Add after the existing `wire_coo_cache` type (line 35) in `src/physics/hamiltonianConstructor.f90`:

```fortran
type :: wire_workspace
  ! Pre-allocated kp-term CSR blocks (structure from first call, values updated per kz)
  type(csr_matrix) :: blk_Q, blk_T, blk_S, blk_SC
  type(csr_matrix) :: blk_R, blk_RC, blk_PZ, blk_PP, blk_PM, blk_A
  type(csr_matrix) :: blk_diff, blk_temp

  ! Pre-allocated COO buffers
  integer, allocatable          :: coo_rows(:), coo_cols(:)
  complex(kind=dp), allocatable :: coo_vals(:)
  integer                       :: coo_capacity = 0

  ! COO-to-CSR mapping (replaces wire_coo_cache)
  integer, allocatable          :: coo_to_csr(:)
  integer                       :: coo_nnz_in = 0

  ! Diagonal position indices within operator-sparsity CSRs
  ! diag_pos(k) = CSR index of the k-th diagonal entry (row k, col k)
  integer, allocatable          :: diag_pos(:)

  logical :: initialized = .false.
end type wire_workspace
```

Add public declarations:
```fortran
public :: wire_workspace, wire_workspace_init, wire_workspace_free
```

**Step 2: Implement `wire_workspace_free`**

```fortran
subroutine wire_workspace_free(ws)
  type(wire_workspace), intent(inout) :: ws

  call csr_free(ws%blk_Q)
  call csr_free(ws%blk_T)
  call csr_free(ws%blk_S)
  call csr_free(ws%blk_SC)
  call csr_free(ws%blk_R)
  call csr_free(ws%blk_RC)
  call csr_free(ws%blk_PZ)
  call csr_free(ws%blk_PP)
  call csr_free(ws%blk_PM)
  call csr_free(ws%blk_A)
  call csr_free(ws%blk_diff)
  call csr_free(ws%blk_temp)
  if (allocated(ws%coo_rows)) deallocate(ws%coo_rows)
  if (allocated(ws%coo_cols)) deallocate(ws%coo_cols)
  if (allocated(ws%coo_vals)) deallocate(ws%coo_vals)
  if (allocated(ws%coo_to_csr)) deallocate(ws%coo_to_csr)
  if (allocated(ws%diag_pos)) deallocate(ws%diag_pos)
  ws%coo_capacity = 0
  ws%coo_nnz_in = 0
  ws%initialized = .false.
end subroutine wire_workspace_free
```

**Step 3: Write the failing test**

Add to `tests/unit/test_hamiltonian_2d.pf`:

```fortran
@test
subroutine test_wire_workspace_free()
  type(wire_workspace) :: ws
  call wire_workspace_free(ws)
  @assertFalse(ws%initialized)
end subroutine test_wire_workspace_free
```

**Step 4: Build and run the test**

Run: `cmake --build build && ctest --test-dir build -R test_hamiltonian_2d -V`
Expected: PASS

**Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 tests/unit/test_hamiltonian_2d.pf
git commit -m "refactor: add wire_workspace type definition and free routine"
```

---

### Task 2: Write regression test for workspace correctness

**Files:**
- Create: `tests/unit/test_wire_workspace.pf`
- Modify: `tests/unit/CMakeLists.txt` (add new test source)

**Context:** Before implementing the fast path, we need a test that builds a wire Hamiltonian at two kz-points using both the original (slow) and workspace (fast) paths, and verifies they produce identical CSR matrices.

**Step 1: Write the comparison test**

Create `tests/unit/test_wire_workspace.pf`:

```fortran
module test_wire_workspace
  use funit
  use definitions
  use parameters
  use sparse_matrices
  use geometry
  use hamiltonianConstructor
  implicit none

contains

  ! Helper: build kpterms_2d for a small GaAs wire grid
  subroutine build_test_kpterms(nx, ny, dx, dy, grid, params, regions, &
      profile_2d, kpterms_2d)
    integer, intent(in) :: nx, ny
    real(kind=dp), intent(in) :: dx, dy
    type(spatial_grid), intent(out) :: grid
    type(paramStruct), allocatable, intent(out) :: params(:)
    type(region_spec), allocatable, intent(out) :: regions(:)
    real(kind=dp), allocatable, intent(out) :: profile_2d(:,:)
    type(csr_matrix), allocatable, intent(out) :: kpterms_2d(:)
    character(len=255), allocatable :: material_names(:)
    integer :: mat2d(nx, ny)

    mat2d = 1
    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    allocate(regions(1))
    regions(1)%material = "GaAs"
    call confinementInitialization_2d(grid, params, regions, profile_2d, kpterms_2d, FDorder=2)
    deallocate(material_names)
  end subroutine build_test_kpterms

  ! Compare two CSR matrices element-by-element
  subroutine assert_csr_equal(A, B, label)
    type(csr_matrix), intent(in) :: A, B
    character(len=*), intent(in) :: label
    integer :: k

    @assertEqual(A%nrows, B%nrows, message=trim(label)//" nrows")
    @assertEqual(A%nnz, B%nnz, message=trim(label)//" nnz")
    do k = 1, A%nnz
      @assertEqual(A%colind(k), B%colind(k), message=trim(label)//" colind")
      @assertEqual(A%values(k), B%values(k), tolerance=1.0e-13_dp, message=trim(label)//" values")
    end do
  end subroutine assert_csr_equal

  @test
  subroutine test_workspace_matches_original_single_kz()
    ! Build wire Hamiltonian at one kz with and without workspace,
    ! verify identical CSR output.
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H_slow, H_fast
    type(wire_coo_cache) :: coo_cache
    real(kind=dp) :: kz
    type(simulation_config) :: cfg

    call build_test_kpterms(nx, ny, dx, dy, grid, params, regions, profile_2d, kpterms_2d)
    kz = 0.05_dp
    cfg%grid = grid

    ! Build without workspace (slow path)
    coo_cache%initialized = .false.
    call ZB8bandGeneralized(H_slow, kz, profile_2d, kpterms_2d, cfg, coo_cache)
    call wire_coo_cache_free(coo_cache)

    ! Build with workspace (first call = slow path, same result)
    coo_cache%initialized = .false.
    call ZB8bandGeneralized(H_fast, kz, profile_2d, kpterms_2d, cfg, coo_cache)
    call assert_csr_equal(H_slow, H_fast, "single kz slow vs first-call")

    call csr_free(H_slow)
    call csr_free(H_fast)
    call wire_coo_cache_free(coo_cache)
    do ij = 1, 17
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions)
    call grid_free(grid)
  end subroutine test_workspace_matches_original_single_kz
end module test_wire_workspace
```

Note: this test initially uses the existing slow path. Once the workspace fast path is implemented (Task 4), extend it to test kz-point-2 with workspace initialized from kz-point-1.

**Step 2: Register the test in CMakeLists.txt**

Add to `tests/unit/CMakeLists.txt`:

```cmake
create_test(test_wire_workspace test_wire_workspace.pf
    hamiltonianConstructor sparse_matrices geometry finitedifferences
    definitions parameters strain_solver input_parser utils)
```

Follow the pattern of existing test registrations in that file.

**Step 3: Build and run**

Run: `cmake --build build && ctest --test-dir build -R test_wire_workspace -V`
Expected: PASS (both paths are slow path currently)

**Step 4: Commit**

```bash
git add tests/unit/test_wire_workspace.pf tests/unit/CMakeLists.txt
git commit -m "test: add wire workspace comparison test skeleton"
```

---

### Task 3: Implement fast path for `build_kp_term_Q` (proof of concept)

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:1795-1812` (build_kp_term_Q)

**Context:** This is the proof of concept. `build_kp_term_Q` computes:
```
Q = -(0.25*kp7 + 0.75*kp8 + kz2*kp1 + (-2*kz2)*kp2)
```
kp7 and kp8 are Laplacian operator terms (same sparsity). kp1 and kp2 are diagonal (subset of Laplacian sparsity). The result has Laplacian sparsity.

On the first call (workspace not initialized), we build normally and stash the result structure. On subsequent calls, we compute values directly.

**Step 1: Add workspace argument and fast-path logic to `build_kp_term_Q`**

Change the subroutine signature to accept an optional workspace:

```fortran
subroutine build_kp_term_Q(kz2, kpterms_2d, blk, ws)
  real(kind=dp), intent(in) :: kz2
  type(csr_matrix), intent(in) :: kpterms_2d(:)
  type(csr_matrix), intent(inout) :: blk
  type(wire_workspace), intent(inout), optional :: ws

  if (present(ws) .and. ws%initialized) then
    ! Fast path: kp7 and kp8 have operator sparsity (= result sparsity)
    blk%values = cmplx(0.25_dp, 0.0_dp, kind=dp) * kpterms_2d(7)%values &
               + cmplx(0.75_dp, 0.0_dp, kind=dp) * kpterms_2d(8)%values
    ! Add diagonal contributions (kp1, kp2) at diagonal positions
    ! kp1 and kp2 are diagonal CSRs with NNZ = N, entries at (k,k)
    ! ws%diag_pos maps sequential index k to CSR position of (k,k) in blk
    block
      integer :: k
      do k = 1, kpterms_2d(1)%nnz
        blk%values(ws%diag_pos(k)) = blk%values(ws%diag_pos(k)) &
          + cmplx(kz2, 0.0_dp, kind=dp) * kpterms_2d(1)%values(k) &
          + cmplx(-2.0_dp * kz2, 0.0_dp, kind=dp) * kpterms_2d(2)%values(k)
      end do
    end block
    ! Apply negation
    blk%values = -blk%values
  else
    ! Original slow path
    type(csr_matrix) :: qxy, kz_diag
    call csr_add(kpterms_2d(7), kpterms_2d(8), qxy, &
      cmplx(0.25_dp, 0.0_dp, kind=dp), cmplx(0.75_dp, 0.0_dp, kind=dp))
    call csr_add(kpterms_2d(1), kpterms_2d(2), kz_diag, &
      cmplx(kz2, 0.0_dp, kind=dp), cmplx(-2.0_dp * kz2, 0.0_dp, kind=dp))
    call csr_add(qxy, kz_diag, blk, UM, UM)
    call csr_free(qxy)
    call csr_free(kz_diag)
    call negate_csr(blk)
  end if
end subroutine build_kp_term_Q
```

**Step 2: Build to verify compilation**

Run: `cmake --build build`
Expected: PASS (no functionality change yet, workspace never passed)

**Step 3: Run all tests**

Run: `ctest --test-dir build`
Expected: All 39 tests PASS (workspace optional arg is not used yet)

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "refactor: add optional workspace arg to build_kp_term_Q with fast path"
```

---

### Task 4: Extend fast path to all kp-term builders

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90` (all `build_kp_term_*` subroutines)

**Context:** Each kp-term builder needs the same treatment as Q. The coefficient tables are:

| Term | Operator kpterms (same sparsity) | Diagonal kpterms | kz power | Final sign |
|------|----------------------------------|-------------------|----------|------------|
| Q | 7*0.25 + 8*0.75 | 1*kz2 + 2*(-2*kz2) | - | negate |
| T | 7*0.75 + 8*0.25 | 1*kz2 + 2*(2*kz2) | - | negate |
| S | 14*(-2*sqrt3*kz) + 15*(i*2*sqrt3*kz) | none | kz^1 | none |
| SC | via S then conjugate_transpose | none | kz^1 | none |
| R | 16*(-sqrt3) + 11*(2i*sqrt3) | 2*(-sqrt3*kz2) | mixed | none |
| RC | via R then conjugate_transpose | none | mixed | none |
| PZ | none | 4*(kz) | kz^1 | none |
| PP | 12*(i/sqrt2) + 13*(-1/sqrt2) | none | kz^0 | none |
| PM | via PP then conjugate_transpose | none | kz^0 | none |
| A | 5*(1) | 10*(kz2) | mixed | none |

**Step 1: Add workspace argument to each `build_kp_term_*` subroutine**

For each subroutine, add `type(wire_workspace), intent(inout), optional :: ws` parameter and implement the fast path using element-wise operations.

Key patterns:
- **Operator-only terms** (S, PP, PM, SC, RC): All contributing kpterms have identical sparsity. Fast path is simple element-wise `blk%values = sum(coeff_i * kp_i%values)`.
- **Diagonal-only terms** (PZ): `blk%values = coeff * kp%values`.
- **Mixed terms** (Q, T, R, A): Add operator contributions element-wise, then add diagonal contributions at pre-computed `diag_pos` indices.

For `build_kp_term_S` (lines 1854-1878):
```fortran
if (present(ws) .and. ws%initialized) then
  blk%values = cmplx(-2.0_dp*SQR3*kz, 0.0_dp, kind=dp) * kpterms_2d(14)%values &
             + cmplx(0.0_dp, 2.0_dp*SQR3*kz, kind=dp) * kpterms_2d(15)%values
else
  ! ... existing code ...
end if
```

For `build_kp_term_R` (lines 1915-1942):
```fortran
if (present(ws) .and. ws%initialized) then
  blk%values = cmplx(-SQR3, 0.0_dp, kind=dp) * kpterms_2d(16)%values &
             + cmplx(0.0_dp, 2.0_dp*SQR3, kind=dp) * kpterms_2d(11)%values
  ! Add diagonal contribution from kp2 at diagonal positions
  block
    integer :: k
    do k = 1, kpterms_2d(2)%nnz
      blk%values(ws%diag_pos(k)) = blk%values(ws%diag_pos(k)) &
        + cmplx(-SQR3*kz2, 0.0_dp, kind=dp) * kpterms_2d(2)%values(k)
    end do
  end block
else
  ! ... existing code ...
end if
```

For `build_kp_term_A` (lines 2049-2056):
```fortran
if (present(ws) .and. ws%initialized) then
  blk%values = kpterms_2d(5)%values
  ! Add diagonal contribution from kp10 at diagonal positions
  block
    integer :: k
    do k = 1, kpterms_2d(10)%nnz
      blk%values(ws%diag_pos(k)) = blk%values(ws%diag_pos(k)) &
        + cmplx(kz2, 0.0_dp, kind=dp) * kpterms_2d(10)%values(k)
    end do
  end block
else
  ! ... existing code ...
end if
```

For `build_kp_term_PZ` (lines 1966-1972):
```fortran
if (present(ws) .and. ws%initialized) then
  blk%values = cmplx(kz, 0.0_dp, kind=dp) * kpterms_2d(4)%values
else
  call csr_scale(kpterms_2d(4), blk, cmplx(kz, 0.0_dp, kind=dp))
end if
```

For `build_kp_term_SC`, `build_kp_term_RC`, `build_kp_term_PM`: These are computed via `build_kp_term_X` then `csr_conjugate_transpose`. The fast path builds the non-conjugated version directly into a temporary workspace CSR, then copies the conjugate-transposed structure.

```fortran
! For build_kp_term_SC:
if (present(ws) .and. ws%initialized) then
  ! Build S values directly into a temp, then conjugate-transpose into blk
  ! ws%blk_temp holds S structure, blk holds SC structure
  call build_kp_term_S(kz, kpterms_2d, ws%blk_temp, ws)
  call csr_conjugate_transpose_inplace(ws%blk_temp, blk)
else
  ! ... existing code ...
end if
```

Note: `csr_conjugate_transpose_inplace` is a new helper that writes the conjugate transpose of A into a pre-allocated B (B already has the correct structure). This is a simple loop over A's nonzeros writing `conjg(A%values(k))` into B at the transposed position.

**Step 2: Build to verify compilation**

Run: `cmake --build build`
Expected: Clean build

**Step 3: Run all tests**

Run: `ctest --test-dir build`
Expected: All 39 tests PASS (workspace not yet wired in)

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "refactor: add optional workspace fast paths to all build_kp_term_* subroutines"
```

---

### Task 5: Implement workspace initialization and wire-up in ZB8bandGeneralized

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:1216-1665` (ZB8bandGeneralized)

**Context:** This is the integration step. `ZB8bandGeneralized` currently takes an optional `wire_coo_cache`. We extend it to also accept an optional `wire_workspace`. On the first call, it builds kp-terms via the slow path, stashes structures, and pre-computes `diag_pos`. On subsequent calls, it uses the fast path.

**Step 1: Change ZB8bandGeneralized signature**

```fortran
subroutine ZB8bandGeneralized(HT_csr, kz, profile_2d, kpterms_2d, cfg, coo_cache, ws, g)
  ! ... existing args ...
  type(wire_coo_cache), intent(inout), optional :: coo_cache
  type(wire_workspace), intent(inout), optional :: ws
  character(len=2), intent(in), optional  :: g
```

Keep `coo_cache` for backward compatibility. When both are present, prefer `ws`.

**Step 2: Add workspace initialization logic**

After building all kp-terms on the first call, stash their structures:

```fortran
if (present(ws)) then
  if (.not. ws%initialized) then
    ! First call: stash CSR structures from the slow-path build
    call csr_clone_structure(blk_Q, ws%blk_Q)
    call csr_clone_structure(blk_T, ws%blk_T)
    call csr_clone_structure(blk_S, ws%blk_S)
    call csr_clone_structure(blk_SC, ws%blk_SC)
    call csr_clone_structure(blk_R, ws%blk_R)
    call csr_clone_structure(blk_RC, ws%blk_RC)
    call csr_clone_structure(blk_PZ, ws%blk_PZ)
    call csr_clone_structure(blk_PP, ws%blk_PP)
    call csr_clone_structure(blk_PM, ws%blk_PM)
    call csr_clone_structure(blk_A, ws%blk_A)

    ! Pre-compute diagonal position mapping
    ! For operator-sparsity CSRs (e.g., blk_Q), find CSR index of each (k,k) entry
    allocate(ws%diag_pos(N))
    do i = 1, N
      do k = ws%blk_Q%rowptr(i), ws%blk_Q%rowptr(i+1) - 1
        if (ws%blk_Q%colind(k) == i) then
          ws%diag_pos(i) = k
          exit
        end if
      end do
    end do

    ! Stash COO buffer capacity for reuse
    ws%coo_capacity = coo_capacity

    ws%initialized = .true.
  else
    ! Fast path: reuse pre-allocated COO buffers
    ! (COO assembly still needed but into pre-allocated arrays)
    coo_capacity = ws%coo_capacity
    allocate(coo_rows(coo_capacity))
    allocate(coo_cols(coo_capacity))
    allocate(coo_vals(coo_capacity))
    coo_idx = 0
  end if
end if
```

**Step 3: Pass workspace to build_kp_term_* calls**

Change all calls like:
```fortran
call build_kp_term_Q(kz2, kpterms_2d, blk_Q)
```
to:
```fortran
call build_kp_term_Q(kz2, kpterms_2d, blk_Q, ws)
```

On the fast path, `blk_Q` points to `ws%blk_Q`. The build_kp_term_* subroutines detect `ws%initialized` and use the element-wise fast path.

**Step 4: On fast path, use workspace CSR blocks directly**

Replace local CSR declarations with references to workspace blocks:

```fortran
if (present(ws) .and. ws%initialized) then
  ! Use workspace CSR blocks directly (no local alloc)
  ! blk_Q etc. are aliases to ws%blk_Q etc.
  ! ... build kp-terms into workspace blocks ...
  ! ... insert into COO as before ...
  ! ... build CSR from COO ...
else
  ! Original path with local CSR variables
end if
```

In practice, since Fortran doesn't have pointer aliasing for derived type components, the cleanest approach is to use `associate`:
```fortran
associate(blk_Q => ws%blk_Q, blk_T => ws%blk_T, ...)
  ! ... existing block insertion code, unchanged ...
end associate
```

**Step 5: Build and run tests**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests PASS

**Step 6: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "feat: integrate wire_workspace into ZB8bandGeneralized"
```

---

### Task 6: Wire up workspace in main programs

**Files:**
- Modify: `src/apps/main.f90` (wire kz-sweep loop)
- Modify: `src/apps/main_optics.f90` (wire optics kz-sweep)

**Context:** The main programs currently use `wire_coo_cache` for the wire kz-sweep. Replace with `wire_workspace`.

**Step 1: In `main.f90`, replace `wire_coo_cache` with `wire_workspace`**

Find the wire kz-sweep loop (search for `coo_cache` usage). Replace:
```fortran
type(wire_coo_cache) :: coo_cache
```
with:
```fortran
type(wire_workspace) :: wire_ws
```

And the call to `ZB8bandGeneralized`:
```fortran
call ZB8bandGeneralized(HT_csr, kz, profile_2d, kpterms_2d, cfg, ws=wire_ws)
```

After the kz-sweep loop, clean up:
```fortran
call wire_workspace_free(wire_ws)
```

**Step 2: Same change in `main_optics.f90`**

**Step 3: Build and run full regression suite**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests PASS — results must be byte-identical to before

**Step 4: Commit**

```bash
git add src/apps/main.f90 src/apps/main_optics.f90
git commit -m "feat: wire up wire_workspace in main programs for kz-sweeps"
```

---

### Task 7: Pre-allocate COO buffers in workspace

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90` (ZB8bandGeneralized COO section)

**Context:** Currently `coo_rows`, `coo_cols`, `coo_vals` are allocated and freed every k-point. Pre-allocate them once in the workspace and reuse.

**Step 1: On first call, allocate COO buffers in workspace**

After computing `coo_capacity`, allocate into the workspace:
```fortran
if (present(ws)) then
  if (.not. ws%initialized) then
    ws%coo_capacity = coo_capacity
    allocate(ws%coo_rows(coo_capacity))
    allocate(ws%coo_cols(coo_capacity))
    allocate(ws%coo_vals(coo_capacity))
  end if
end if
```

**Step 2: On subsequent calls, reuse workspace COO buffers**

Instead of `allocate(coo_rows(coo_capacity))`, use the workspace arrays:
```fortran
if (present(ws) .and. ws%initialized) then
  coo_capacity = ws%coo_capacity
  ! Point local aliases to workspace arrays
  ! (Fortran: use associate or pass workspace arrays directly)
else
  allocate(coo_rows(coo_capacity), coo_cols(coo_capacity), coo_vals(coo_capacity))
end if
```

At cleanup, don't deallocate workspace buffers:
```fortran
if (.not. (present(ws) .and. ws%initialized)) then
  deallocate(coo_rows, coo_cols, coo_vals)
end if
```

**Step 3: Build and run tests**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests PASS

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "perf: pre-allocate COO buffers in wire_workspace across k-points"
```

---

## Phase 2: Data-Driven Strain COO Insertion

### Task 8: Define strain_entry type and build the table

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:2216-2424` (insert_strain_coo)

**Context:** The current `insert_strain_coo` is a 208-line manually unrolled loop with 32 COO entries per grid point. Replace it with a data-driven table of 32 entries.

**Step 1: Define the strain entry type**

Add before `insert_strain_coo`:

```fortran
type :: strain_entry
  integer :: row_band, col_band  ! 0-based band offsets (0-7)
  integer :: field_id            ! 1=EHH, 2=ELH, 3=ESO, 4=Ec, 5=S_eps, 6=R_eps, 7=QT2_eps, 8=P_eps_VBSO
  integer :: cflag               ! 0=plain, 1=conjg, 2=negate, 3=negate+conjg
end type strain_entry
```

**Step 2: Build the parameter table**

Construct the 32-entry table by reading the existing `insert_strain_coo` (lines 2228-2423). Each `coo_v(coo_idx) = <expression>` maps to a table entry with the appropriate field_id and cflag.

The table has 8 diagonal entries + 4 S_eps entries + 4 R_eps entries + 16 VB-SO coupling entries.

Enumerate all 32 entries as `type(strain_entry), parameter :: strain_table(32)`. The scale factors (like `-IU * RQS2`, `IU * SQR2`) are applied at runtime based on the field and flags.

Note: The VB-SO coupling entries use compound factors like `-IU * RQS2 * conjg(S_eps)` which combine cflag with runtime complex constants. Store the complex prefactor in the table as well:

```fortran
type :: strain_entry
  integer :: row_band, col_band
  integer :: field_id
  complex(kind=dp) :: prefactor  ! complex scale factor (IU*RQS2, etc.)
  logical :: use_conjg           ! whether to conjugate the field value
end type strain_entry
```

**Step 3: Build to verify compilation**

Run: `cmake --build build`
Expected: PASS

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "refactor: define strain_entry type and parameter table"
```

---

### Task 9: Implement table-driven insert_strain_coo and test

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:2216-2424` (replace insert_strain_coo)
- Modify: `tests/unit/test_strain_solver.pf` (add COO comparison test)

**Step 1: Write the test first**

Add to `tests/unit/test_strain_solver.pf`:

```fortran
@test
subroutine test_strain_coo_table_matches_original()
  ! Build a small bir_pikus_blocks, insert strain COO using both
  ! old (unrolled) and new (table) approaches, compare outputs.
  integer, parameter :: N = 5
  integer, parameter :: coo_cap = 32 * N + 100
  integer :: coo_r(coo_cap), coo_c(coo_cap), coo_r2(coo_cap), coo_c2(coo_cap)
  complex(kind=dp) :: coo_v(coo_cap), coo_v2(coo_cap)
  integer :: idx1, idx2, k
  type(bir_pikus_blocks) :: bp

  ! Build synthetic strain blocks
  allocate(bp%delta_EHH(N), bp%delta_ELH(N), bp%delta_ESO(N), bp%delta_Ec(N))
  allocate(bp%S_eps(N), bp%R_eps(N), bp%QT2_eps(N))
  bp%delta_EHH = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  bp%delta_ELH = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
  bp%delta_ESO = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp]
  bp%delta_Ec = [10.0_dp, 11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp]
  bp%S_eps = [(0.1_dp, 0.2_dp), (0.3_dp, 0.4_dp), (0.5_dp, 0.6_dp), (0.7_dp, 0.8_dp), (0.9_dp, 1.0_dp)]
  bp%R_eps = [(0.01_dp, 0.02_dp), (0.03_dp, 0.04_dp), (0.05_dp, 0.06_dp), (0.07_dp, 0.08_dp), (0.09_dp, 0.10_dp)]
  bp%QT2_eps = [(1.0_dp, 0.1_dp), (2.0_dp, 0.2_dp), (3.0_dp, 0.3_dp), (4.0_dp, 0.4_dp), (5.0_dp, 0.5_dp)]

  ! Insert with both approaches
  idx1 = 0
  call insert_strain_coo(coo_r, coo_c, coo_v, coo_cap, idx1, bp, N)
  idx2 = 0
  call insert_strain_coo_table(coo_r2, coo_c2, coo_v2, coo_cap, idx2, bp, N)

  ! Compare
  @assertEqual(idx1, idx2)
  do k = 1, idx1
    @assertEqual(coo_r(k), coo_r2(k))
    @assertEqual(coo_c(k), coo_c2(k))
    @assertEqual(coo_v(k), coo_v2(k), tolerance=1.0e-14_dp)
  end do

  call bir_pikus_blocks_free(bp)
end subroutine test_strain_coo_table_matches_original
```

**Step 2: Implement `insert_strain_coo_table`**

```fortran
subroutine insert_strain_coo_table(coo_r, coo_c, coo_v, coo_cap, coo_idx, bp, N)
  ! ... same interface as insert_strain_coo ...
  integer :: ii, e, g_row, g_col
  complex(kind=dp) :: field_val

  do ii = 1, N
    do e = 1, 32  ! size(strain_table)
      coo_idx = coo_idx + 1
      if (coo_idx > coo_cap) return

      g_row = strain_table(e)%row_band * N + ii
      g_col = strain_table(e)%col_band * N + ii

      ! Extract field value
      select case (strain_table(e)%field_id)
      case (1); field_val = cmplx(bp%delta_EHH(ii), 0.0_dp, kind=dp)
      case (2); field_val = cmplx(bp%delta_ELH(ii), 0.0_dp, kind=dp)
      case (3); field_val = cmplx(bp%delta_ESO(ii), 0.0_dp, kind=dp)
      case (4); field_val = cmplx(bp%delta_Ec(ii), 0.0_dp, kind=dp)
      case (5); field_val = bp%S_eps(ii)
      case (6); field_val = bp%R_eps(ii)
      case (7); field_val = bp%QT2_eps(ii)
      end select

      ! Apply conjugation
      if (strain_table(e)%use_conjg) field_val = conjg(field_val)

      coo_r(coo_idx) = g_row
      coo_c(coo_idx) = g_col
      coo_v(coo_idx) = strain_table(e)%prefactor * field_val
    end do
  end do
end subroutine insert_strain_coo_table
```

**Step 3: Run test to verify it fails (insert_strain_coo_table doesn't exist yet)**

Run: `cmake --build build && ctest --test-dir build -R test_strain -V`
Expected: FAIL (undefined symbol)

**Step 4: Implement and verify**

Run: `cmake --build build && ctest --test-dir build -R test_strain -V`
Expected: PASS

**Step 5: Replace `insert_strain_coo` body with table-driven version**

Once the test passes, replace the original 208-line `insert_strain_coo` with a call to `insert_strain_coo_table`. Or rename `insert_strain_coo_table` to `insert_strain_coo` and remove the old implementation.

**Step 6: Run full test suite**

Run: `ctest --test-dir build`
Expected: All 39 tests PASS

**Step 7: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 tests/unit/test_strain_solver.pf
git commit -m "refactor: replace 208-line unrolled strain COO with data-driven table"
```

---

## Phase 3: FEAST Workspace

### Task 10: Add `feast_workspace` type to eigensolver

**Files:**
- Modify: `src/math/eigensolver.f90:1-193` (solve_feast)

**Context:** `solve_feast` extracts the upper triangle of the CSR matrix on every call (lines 131-158). For wire kz-sweeps where the sparsity pattern doesn't change, this wastes one O(NNZ) scan + 3 allocations per k-point.

**Step 1: Define feast_workspace type**

Add to `src/math/eigensolver.f90` after the existing type definitions:

```fortran
type :: feast_workspace
  integer, allocatable          :: rowptr_loc(:)   ! (N+1) upper-triangle row pointers
  integer, allocatable          :: colind_loc(:)   ! (nnz_upper) column indices
  integer                       :: nnz_upper = 0
  ! Pre-allocated FEAST arrays
  real(kind=dp), allocatable    :: E(:)
  complex(kind=dp), allocatable :: X(:,:)
  real(kind=dp), allocatable    :: res(:)
  integer                       :: M0 = 0
  logical                       :: initialized = .false.
end type feast_workspace

public :: feast_workspace, feast_workspace_free
```

**Step 2: Implement feast_workspace_free**

```fortran
subroutine feast_workspace_free(fw)
  type(feast_workspace), intent(inout) :: fw
  if (allocated(fw%rowptr_loc)) deallocate(fw%rowptr_loc)
  if (allocated(fw%colind_loc)) deallocate(fw%colind_loc)
  if (allocated(fw%E)) deallocate(fw%E)
  if (allocated(fw%X)) deallocate(fw%X)
  if (allocated(fw%res)) deallocate(fw%res)
  fw%nnz_upper = 0
  fw%M0 = 0
  fw%initialized = .false.
end subroutine feast_workspace_free
```

**Step 3: Write the test**

Add to `tests/unit/test_eigensolver.pf`:

```fortran
@test
subroutine test_feast_workspace_cached_matches_uncached()
  ! Solve the same eigenvalue problem with and without workspace.
  ! Eigenvalues and eigenvectors must be identical.
  integer, parameter :: n = 30
  real(kind=dp) :: diag(n)
  type(csr_matrix) :: H
  type(eigensolver_config) :: cfg
  type(eigensolver_result) :: res1, res2
  type(feast_workspace) :: fw
  integer :: i

  diag = [(dble(i), i = 1, n)]
  call build_diagonal_csr(diag, n, H)

  cfg%method = 'FEAST'
  cfg%emin = 5.0_dp
  cfg%emax = 15.0_dp
  cfg%nev = 10

  ! Solve without workspace
  call solve_feast(H, cfg, res1)

  ! Solve with workspace (first call initializes cache)
  call solve_feast(H, cfg, res2, fw)
  @assertEqual(res1%nev_found, res2%nev_found)
  do i = 1, res1%nev_found
    @assertEqual(res1%eigenvalues(i), res2%eigenvalues(i), tolerance=1.0e-12_dp)
  end do
  call eigensolver_result_free(res1)
  call eigensolver_result_free(res2)

  ! Solve again with cached workspace (tests fast path)
  call solve_feast(H, cfg, res2, fw)
  @assertTrue(fw%initialized)
  ! Re-solve without workspace for comparison
  call solve_feast(H, cfg, res1)
  @assertEqual(res1%nev_found, res2%nev_found)
  do i = 1, res1%nev_found
    @assertEqual(res1%eigenvalues(i), res2%eigenvalues(i), tolerance=1.0e-12_dp)
  end do

  call eigensolver_result_free(res1)
  call eigensolver_result_free(res2)
  call feast_workspace_free(fw)
  call csr_free(H)
end subroutine test_feast_workspace_cached_matches_uncached
```

**Step 4: Implement caching in solve_feast**

Add optional `fw` argument to `solve_feast`:

```fortran
subroutine solve_feast(H_csr, config, result, fw)
  type(csr_matrix), intent(in)          :: H_csr
  type(eigensolver_config), intent(in)  :: config
  type(eigensolver_result), intent(out) :: result
  type(feast_workspace), intent(inout), optional :: fw
```

Replace the upper-triangle extraction block (lines 131-158) with:

```fortran
if (present(fw) .and. fw%initialized .and. fw%M0 == M0) then
  ! Fast path: reuse cached upper-triangle structure, just update values
  do i = 1, N
    k = fw%rowptr_loc(i)
    do j = H_csr%rowptr(i), H_csr%rowptr(i+1) - 1
      if (H_csr%colind(j) >= i) then
        val_loc(k) = H_csr%values(j)
        k = k + 1
      end if
    end do
  end do
else
  ! Original path: count, allocate, fill
  ! ... existing code from lines 131-158 ...

  ! Cache for future calls
  if (present(fw)) then
    call move_alloc(rowptr_loc, fw%rowptr_loc)
    call move_alloc(colind_loc, fw%colind_loc)
    fw%nnz_upper = nnz
    fw%M0 = M0
    fw%initialized = .true.
    ! Re-allocate local arrays that were moved
    allocate(val_loc(nnz))
    ! Refill val_loc (was already filled above)
    ! Actually: move rowptr/colind, keep val_loc local
    ! ... adjust logic to keep val_loc as local allocatable ...
  end if
end if
```

Also cache E, X, res arrays:
```fortran
if (present(fw) .and. fw%initialized) then
  E => fw%E
  X => fw%X
  res => fw%res
else
  allocate(E(M0), X(N, M0), res(M0))
  if (present(fw)) then
    call move_alloc(E, fw%E)
    call move_alloc(X, fw%X)
    call move_alloc(res, fw%res)
    fw%M0 = M0
  end if
end if
```

**Step 5: Build and run tests**

Run: `cmake --build build && ctest --test-dir build`
Expected: All tests PASS

**Step 6: Commit**

```bash
git add src/math/eigensolver.f90 tests/unit/test_eigensolver.pf
git commit -m "feat: add feast_workspace with upper-triangle caching"
```

---

### Task 11: Wire up feast_workspace in main programs

**Files:**
- Modify: `src/apps/main.f90` (wire kz-sweep)
- Modify: `src/apps/main_optics.f90` (wire optics kz-sweep)

**Context:** In the wire kz-sweep loops, declare a `feast_workspace` and pass it to `solve_sparse_evp` (or directly to `solve_feast`). The eigensolver dispatch needs to pass the workspace through.

**Step 1: Add feast_workspace parameter to solve_sparse_evp**

```fortran
subroutine solve_sparse_evp(H_csr, config, result, feast_ws)
  type(csr_matrix), intent(in)          :: H_csr
  type(eigensolver_config), intent(in)  :: config
  type(eigensolver_result), intent(out) :: result
  type(feast_workspace), intent(inout), optional :: feast_ws

  case ('FEAST')
    call solve_feast(H_csr, config, result, fw=feast_ws)
```

**Step 2: Declare and pass feast_workspace in main.f90**

In the wire kz-sweep:
```fortran
type(feast_workspace) :: feast_ws
! ... inside loop:
call solve_sparse_evp(HT_csr, eigen_cfg, eigen_res, feast_ws)
! ... after loop:
call feast_workspace_free(feast_ws)
```

**Step 3: Same in main_optics.f90**

**Step 4: Build and run full suite**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests PASS

**Step 5: Commit**

```bash
git add src/math/eigensolver.f90 src/apps/main.f90 src/apps/main_optics.f90
git commit -m "feat: wire up feast_workspace in wire kz-sweep loops"
```

---

## Phase 4: Three-Module Split

### Task 12: Extract `confinement_init` module

**Files:**
- Create: `src/physics/confinement_init.f90`
- Modify: `src/physics/hamiltonianConstructor.f90` (remove extracted code)
- Modify: `src/CMakeLists.txt` (add new source file)

**Context:** Extract the confinement initialization subroutines from `hamiltonianConstructor.f90` into a separate `confinement_init` module. This module has no circular dependency risk — it's pure infrastructure.

**What to extract:**
- `confinementInitialization_raw` (the raw-parameter overload)
- `confinementInitialization_cfg` (the config-based overload)
- `confinementInitialization_2d` (the 2D wire overload)
- `build_kpterm_block` (FD helper)
- `applyVariableCoeffStaggered` (variable coefficient helper)
- `build_diagonal_csr` (profile diagonal helper)
- Dirichlet boundary closure helpers (if any local subroutines)

**Estimated size:** ~650 lines

**Step 1: Create `src/physics/confinement_init.f90`**

```fortran
module confinement_init
  use definitions
  use finitedifferences
  use sparse_matrices
  use utils
  use strain_solver, only: compute_bp_scalar
  use input_parser, only: simulation_config
  implicit none

  interface confinementInitialization
    module procedure confinementInitialization_raw
    module procedure confinementInitialization_cfg
  end interface confinementInitialization

  public :: confinementInitialization, confinementInitialization_2d
  public :: build_kpterm_block, build_diagonal_csr

contains
  ! ... move extracted subroutines here ...
end module confinement_init
```

**Step 2: Update `hamiltonianConstructor.f90`**

- Add `use confinement_init` to the module
- Remove the extracted subroutines
- Remove the `confinementInitialization` interface block

**Step 3: Update `src/CMakeLists.txt`**

Add `physics/confinement_init.f90` to `COMMON_SOURCES` and the `set_source_files_properties` list.

**Step 4: Build and test**

Run: `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl && cmake --build build && ctest --test-dir build`
Expected: All 39 tests PASS (note: must reconfigure to pick up new source file)

**Step 5: Commit**

```bash
git add src/physics/confinement_init.f90 src/physics/hamiltonianConstructor.f90 src/CMakeLists.txt
git commit -m "refactor: extract confinement_init module from hamiltonianConstructor"
```

---

### Task 13: Extract `hamiltonian_wire` module

**Files:**
- Create: `src/physics/hamiltonian_wire.f90`
- Modify: `src/physics/hamiltonianConstructor.f90` (remove extracted wire code)
- Modify: `src/CMakeLists.txt` (add new source file)

**Context:** Extract all wire-specific code into a separate module. This module is self-contained — it never calls QW or bulk routines.

**What to extract:**
- `ZB8bandGeneralized` (the main wire builder)
- All `build_kp_term_*` subroutines (Q, T, S, SC, R, RC, PZ, PP, PM, A)
- `insert_csr_block`, `insert_csr_block_scaled`
- `insert_strain_coo` / `insert_strain_coo_table`
- `insert_profile_diagonal`
- `build_velocity_matrices_2d`
- `wire_workspace` type and `wire_workspace_free`
- `csr_conjugate_transpose`, `negate_csr`

**Estimated size:** ~1300 lines

**Step 1: Create `src/physics/hamiltonian_wire.f90`**

```fortran
module hamiltonian_wire
  use definitions
  use sparse_matrices
  use input_parser, only: simulation_config
  use strain_solver, only: bir_pikus_blocks, bir_pikus_blocks_free
  use confinement_init, only: confinementInitialization_2d
  implicit none

  public :: ZB8bandGeneralized
  public :: build_velocity_matrices_2d
  public :: wire_workspace, wire_workspace_free

  interface build_velocity_matrices
    module procedure build_velocity_matrices_2d
  end interface build_velocity_matrices

  ! ... type definitions and subroutines ...
end module hamiltonian_wire
```

**Step 2: Update `hamiltonianConstructor.f90`**

- Remove all extracted wire subroutines
- Remove `wire_coo_cache` type (replaced by `wire_workspace`)
- The module retains: `ZB8bandQW`, `ZB8bandBulk`, `add_bp_strain_dense`, `externalFieldSetup_electricField`, `build_velocity_matrices_1d`, profile helpers

**Step 3: Update app files**

In `main.f90`, `main_gfactor.f90`, `main_optics.f90`: add `use hamiltonian_wire` for wire mode paths.

In `optical_spectra.f90`: add `use hamiltonian_wire` if it calls `build_velocity_matrices_2d`.

**Step 4: Update CMakeLists.txt**

Add `physics/hamiltonian_wire.f90` to `COMMON_SOURCES` and `set_source_files_properties`.

**Step 5: Build and test**

Run: `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl && cmake --build build && ctest --test-dir build`
Expected: All 39 tests PASS

**Step 6: Commit**

```bash
git add src/physics/hamiltonian_wire.f90 src/physics/hamiltonianConstructor.f90 \
      src/CMakeLists.txt src/apps/main.f90 src/apps/main_gfactor.f90 \
      src/apps/main_optics.f90 src/physics/optical_spectra.f90
git commit -m "refactor: extract hamiltonian_wire module for wire-specific code"
```

---

### Task 14: Clean up and final verification

**Files:**
- All modified files
- Verify: `docs/plans/2026-04-26-hamiltonian-performance-refactor-design.md`

**Step 1: Remove stale `.mod` files and full rebuild**

```bash
rm -f *.mod
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl
cmake --build build
```

**Step 2: Run full test suite**

```bash
ctest --test-dir build -V
```
Expected: All tests PASS (39 total: 15 unit + 24 regression)

**Step 3: Verify line counts**

```bash
wc -l src/physics/confinement_init.f90 src/physics/hamiltonianConstructor.f90 src/physics/hamiltonian_wire.f90
```
Expected:
- confinement_init.f90: ~650 lines
- hamiltonianConstructor.f90: ~500 lines
- hamiltonian_wire.f90: ~1300 lines

**Step 4: Verify no circular dependencies**

```bash
grep -n "^  use " src/physics/hamiltonian_wire.f90
```
Expected: Uses `definitions`, `sparse_matrices`, `input_parser`, `strain_solver`, `confinement_init` — NOT `hamiltonianConstructor`.

**Step 5: Commit any final cleanup**

```bash
git add -A
git commit -m "refactor: final cleanup after 3-module split"
```

---

## Summary of Expected Impact

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Heap allocs per k-point (wire) | ~50 CSR + 3 COO | ~0 (workspace reuse) | ~98% reduction |
| FEAST allocs per k-point | 5 (rowptr, colind, val_loc, E, X) | 1 (val_loc) | ~80% reduction |
| insert_strain_coo lines | 208 (unrolled) | ~30 (table-driven) | ~85% reduction |
| hamiltonianConstructor.f90 | ~2450 lines | ~500 lines | 3-module split |
| Total allocs for 200 kz-points | ~10,600 | ~200 | ~98% reduction |

## Testing Strategy

Each task includes:
1. **Unit test** (pFUnit `.pf` file): verifies new code in isolation
2. **Regression suite** (39 tests): verifies no breakage across all configs
3. **Comparison test**: verifies byte-identical results between old and new paths

The comparison tests (Task 2, Task 9, Task 10) are the critical safety net — they prove the optimized code produces identical results to the original.
