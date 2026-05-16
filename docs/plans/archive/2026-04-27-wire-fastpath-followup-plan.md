# Wire Fast-Path Follow-Up Refactoring Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Eliminate per-kz COO sorting overhead and deduplicate the ~260-line block-insertion code in `hamiltonian_wire.f90`.

**Architecture:** Absorb `wire_coo_cache` into `wire_workspace` so callers pass a single object. Extract shared block-insertion subroutines that both fast and slow paths call with different COO backing stores. Move general CSR utilities to `sparse_matrices.f90`.

**Tech Stack:** Fortran 90, MKL SpBLAS, pFUnit tests, CMake/Ninja build.

---

### Task 1: Add `csr_find_diag_positions` to `sparse_matrices.f90`

Three locations in `hamiltonian_wire.f90` (lines 838-851, 854-867, 870-883) contain identical 13-line diagonal-position lookup loops. Extract to a reusable utility.

**Files:**
- Modify: `src/math/sparse_matrices.f90` (add new subroutine after `csr_scale` at ~line 806)
- Modify: `src/physics/hamiltonian_wire.f90:838-883` (replace 3 copies with calls)

**Step 1: Write the utility in `sparse_matrices.f90`**

Add after `csr_scale` (after line 806):

```fortran
  subroutine csr_find_diag_positions(mat, diag_pos)
    type(csr_matrix), intent(in)  :: mat
    integer, intent(out) :: diag_pos(mat%nrows)

    integer :: i, k

    diag_pos = 0
    do i = 1, mat%nrows
      do k = mat%rowptr(i), mat%rowptr(i + 1) - 1
        if (mat%colind(k) == i) then
          diag_pos(i) = k
          exit
        end if
      end do
    end do
    if (any(diag_pos == 0)) then
      print *, 'ERROR: csr_find_diag_positions: no diagonal entry for some rows'
      stop 1
    end if
  end subroutine csr_find_diag_positions
```

**Step 2: Add to public interface**

In `sparse_matrices.f90`, add `public :: csr_find_diag_positions`.

**Step 3: Replace 3 copies in `hamiltonian_wire.f90`**

Replace lines 838-851 (diag_pos from blk_Q):
```fortran
        allocate(ws%diag_pos(N))
        call csr_find_diag_positions(ws%blk_Q, ws%diag_pos)
```

Replace lines 854-867 (diag_pos_R from blk_R):
```fortran
        allocate(ws%diag_pos_R(N))
        call csr_find_diag_positions(ws%blk_R, ws%diag_pos_R)
```

Replace lines 870-883 (diag_pos_A from blk_A):
```fortran
        allocate(ws%diag_pos_A(N))
        call csr_find_diag_positions(ws%blk_A, ws%diag_pos_A)
```

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

**Step 5: Commit**

```bash
git add src/math/sparse_matrices.f90 src/physics/hamiltonian_wire.f90
git commit -m "refactor: extract csr_find_diag_positions to sparse_matrices"
```

---

### Task 2: Cache `build_strain_table()` with module-level `save`

`insert_strain_coo` (line 1537) calls `build_strain_table()` every kz-point. The table is a 32-entry compile-time constant.

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:1520-1570`

**Step 1: Add module-level cached table**

Near the top of the `contains` block (after line 92), add a cache variable and wrapper:

```fortran
  logical, save :: strain_table_cached = .false.
  type(strain_entry), save :: strain_table_cache(32)
```

**Step 2: Add cached accessor function**

```fortran
  function get_strain_table() result(table)
    type(strain_entry) :: table(32)
    if (.not. strain_table_cached) then
      strain_table_cache = build_strain_table()
      strain_table_cached = .true.
    end if
    table = strain_table_cache
  end function get_strain_table
```

**Step 3: Replace call in `insert_strain_coo`**

Line 1537: change `table = build_strain_table()` to `table = get_strain_table()`.

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

**Step 5: Commit**

```bash
git add src/physics/hamiltonian_wire.f90
git commit -m "perf: cache strain table across kz-points"
```

---

### Task 3: Move `csr_conjugate_transpose` utilities to `sparse_matrices.f90`

`csr_conjugate_transpose` (lines 1357-1381) and `csr_conjugate_transpose_to_preallocated` (lines 1383-1416) are general sparse linear algebra with no physics logic. They belong in `sparse_matrices.f90`.

**Files:**
- Modify: `src/math/sparse_matrices.f90` (add two subroutines)
- Modify: `src/physics/hamiltonian_wire.f90:1357-1416` (remove, use import)

**Step 1: Copy subroutines to `sparse_matrices.f90`**

Add `csr_conjugate_transpose` and `csr_conjugate_transpose_to_preallocated` after `csr_scale` in `sparse_matrices.f90`. Add `public :: csr_conjugate_transpose, csr_conjugate_transpose_to_preallocated`.

**Step 2: Remove from `hamiltonian_wire.f90`**

Delete lines 1357-1416. The `use sparse_matrices` at line 4 already provides these.

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

**Step 4: Commit**

```bash
git add src/math/sparse_matrices.f90 src/physics/hamiltonian_wire.f90
git commit -m "refactor: move csr_conjugate_transpose to sparse_matrices"
```

---

### Task 4: Absorb `wire_coo_cache` into `wire_workspace`

Neither `main.f90` nor `main_optics.f90` passes `coo_cache` alongside `wire_ws`. The O(NNZ log NNZ) sort runs every kz-point despite the caching infrastructure. Merge `coo_cache` into `wire_workspace` so the sort-avoidance is always active when the workspace is used.

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:27-64` (type definition)
- Modify: `src/physics/hamiltonian_wire.f90:190-228` (ZB8bandGeneralized signature)
- Modify: `src/physics/hamiltonian_wire.f90:489-507` (fast path finalization)
- Modify: `src/physics/hamiltonian_wire.f90:800-818` (slow path finalization)
- Modify: `src/physics/hamiltonian_wire.f90:893-898` (workspace init)
- Modify: `src/physics/hamiltonian_wire.f90:1587+` (wire_workspace_free)
- Modify: `src/apps/main.f90:46` (remove coo_cache declaration)
- Modify: `src/apps/main.f90:280-285` (remove coo_cache reset)
- Modify: `src/apps/main.f90:329-330` (remove coo_cache from call)
- Modify: `src/apps/main.f90:437` (remove coo_cache free)
- Modify: `src/apps/main_optics.f90` (remove coo_cache if present)
- Modify: `src/physics/sc_loop.f90:565` (update coo_cache usage)

**Step 1: Add `coo_cache` fields to `wire_workspace` type**

In the `wire_workspace` type (line 38), add before `logical :: initialized`:

```fortran
    ! COO-to-CSR sort cache (absorbed from wire_coo_cache)
    integer, allocatable          :: coo_to_csr(:)
    integer                       :: coo_nnz_in = 0
    logical                       :: coo_cache_valid = .false.
```

**Step 2: Update `wire_workspace_free`**

Add cleanup for the new fields:

```fortran
  if (allocated(ws%coo_to_csr)) deallocate(ws%coo_to_csr)
  ws%coo_cache_valid = .false.
```

**Step 3: Update slow-path workspace init (lines ~893)**

After `ws%initialized = .true.`, the COO cache is not yet valid — it will be populated on the first `csr_build_from_coo_cached` call in the finalization block.

**Step 4: Update fast-path COO finalization (lines 489-507)**

Replace the `if (present(coo_cache))` block with:

```fortran
if (ws%coo_cache_valid) then
  call csr_set_values_from_coo(HT_csr, coo_idx, &
    ws%coo_to_csr(1:coo_idx), coo_vals(1:coo_idx))
else
  call csr_build_from_coo_cached(HT_csr, Ntot, Ntot, coo_idx, &
    coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx), &
    ws%coo_to_csr)
  ws%coo_nnz_in = coo_idx
  ws%coo_cache_valid = .true.
end if
```

**Step 5: Update slow-path COO finalization (lines 800-818)**

Same pattern as Step 4 but using `ws` optional:

```fortran
if (present(ws)) then
  if (ws%coo_cache_valid) then
    call csr_set_values_from_coo(HT_csr, coo_idx, &
      ws%coo_to_csr(1:coo_idx), coo_vals(1:coo_idx))
  else
    call csr_build_from_coo_cached(HT_csr, Ntot, Ntot, coo_idx, &
      coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx), &
      ws%coo_to_csr)
    ws%coo_nnz_in = coo_idx
    ws%coo_cache_valid = .true.
  end if
else
  if (coo_idx > 0) then
    call csr_build_from_coo(HT_csr, Ntot, Ntot, coo_idx, &
      coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx))
  else
    call csr_init(HT_csr, Ntot, Ntot)
  end if
end if
```

**Step 6: Remove `coo_cache` parameter from `ZB8bandGeneralized` and internal calls**

- Remove `coo_cache` parameter from `ZB8bandGeneralized` (line 196), `zb8_generalized_fast` (line 239), and `zb8_generalized_slow` (line 522).
- Remove `type(wire_coo_cache)` from `main.f90` line 46 and all `coo_cache` references in `main.f90`.
- Update `sc_loop.f90` — the SC loop uses `coo_cache` without `ws`. Add a note that SC loop will need a separate follow-up to adopt workspace, or keep the `coo_cache` parameter for the SC path only.

**Step 7: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass. The `regression_wire_*` tests should show no numerical change.

**Step 8: Commit**

```bash
git add src/physics/hamiltonian_wire.f90 src/apps/main.f90 src/apps/main_optics.f90 src/physics/sc_loop.f90
git commit -m "perf: absorb coo_cache into wire_workspace for sort-free kz-sweeps"
```

---

### Task 5: Extract shared `insert_g3_blocks` subroutine

The g3-mode COO block insertion is duplicated between fast (lines 282-336) and slow (lines 574-628) paths. The insertion calls are identical; only the COO arrays differ.

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90` (add subroutine, replace both copies)

**Step 1: Write the shared subroutine**

Add as a private subroutine in `hamiltonian_wire`:

```fortran
  subroutine insert_g3_blocks(coo_r, coo_c, coo_v, coo_cap, coo_idx, &
      blk_S, blk_SC, blk_PZ, N)
    integer, intent(inout) :: coo_r(:), coo_c(:)
    complex(kind=dp), intent(inout) :: coo_v(:)
    integer, intent(in) :: coo_cap
    integer, intent(inout) :: coo_idx
    type(csr_matrix), intent(in) :: blk_S, blk_SC, blk_PZ
    integer, intent(in) :: N

    ! Row 1 (HH1)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 1, blk_SC, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 4, blk_SC, N, -IU*RQS2)
    ! Row 2 (HH2)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 0, blk_S, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 5, blk_SC, N, -IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 6, blk_PZ, N, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    ! Row 3 (LH1)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 3, blk_SC, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 4, blk_S, N, IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 7, blk_PZ, N, IU*SQR2*RQS3)
    ! Row 4 (LH2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 2, blk_S, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 5, blk_S, N, IU*RQS2)
    ! Row 5 (SO1)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 0, blk_S, N, IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 2, blk_SC, N, -IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 6, blk_PZ, N, IU*RQS3)
    ! Row 6 (SO2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 1, blk_S, N, IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 3, blk_SC, N, -IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 7, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))
    ! Row 7 (CB1)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 1, blk_PZ, N, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 4, blk_PZ, N, -IU*RQS3)
    ! Row 8 (CB2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 2, blk_PZ, N, -IU*SQR2*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 5, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))
  end subroutine insert_g3_blocks
```

**Step 2: Replace fast-path g3 insertion (lines 282-336)**

Replace with:
```fortran
        call insert_g3_blocks(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, blk_S, blk_SC, blk_PZ, N)
```

**Step 3: Replace slow-path g3 insertion (lines 574-628)**

Replace with:
```fortran
        call insert_g3_blocks(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, blk_S, blk_SC, blk_PZ, N)
```

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build -R wire`
Expected: All wire tests pass (regression_wire_insb_gfactor uses g3 mode).

**Step 5: Commit**

```bash
git add src/physics/hamiltonian_wire.f90
git commit -m "refactor: extract insert_g3_blocks shared between fast/slow paths"
```

---

### Task 6: Extract shared `insert_main_blocks` subroutine

The main 8x8 block insertion (~140 lines) is duplicated between fast (lines 345-478) and slow (lines 654-779). The insertion topology and prefactors are identical. The only difference is `blk_diff` and `blk_temp` are computed differently (element-wise vs csr_add), but by the time they reach the insertion calls, they hold the same values.

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90` (add subroutine, replace both copies)

**Step 1: Write the shared subroutine**

The subroutine takes all blocks + COO arrays as arguments. Both `blk_diff` and `blk_temp` are pre-computed before calling.

```fortran
  subroutine insert_main_blocks(coo_r, coo_c, coo_v, coo_cap, coo_idx, &
      blk_Q, blk_T, blk_S, blk_SC, blk_R, blk_RC, blk_PZ, blk_PP, &
      blk_PM, blk_A, blk_diff, blk_temp, N)
    integer, intent(inout) :: coo_r(:), coo_c(:)
    complex(kind=dp), intent(inout) :: coo_v(:)
    integer, intent(in) :: coo_cap
    integer, intent(inout) :: coo_idx
    type(csr_matrix), intent(in) :: blk_Q, blk_T, blk_S, blk_SC
    type(csr_matrix), intent(in) :: blk_R, blk_RC, blk_PZ, blk_PP
    type(csr_matrix), intent(in) :: blk_PM, blk_A, blk_diff, blk_temp
    integer, intent(in) :: N

    ! Row 1 (HH1)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 0, blk_Q, N)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 1, blk_SC, N)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 2, blk_RC, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 4, blk_SC, N, -IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 5, blk_RC, N, IU*SQR2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 0, 6, blk_PP, N, IU)
    ! Row 2 (HH2)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 0, blk_S, N)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 1, blk_T, N)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 3, blk_RC, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 4, blk_diff, N, IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 5, blk_SC, N, -IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 6, blk_PZ, N, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 1, 7, blk_PP, N, cmplx(-RQS3, 0.0_dp, kind=dp))
    ! Row 3 (LH1)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 0, blk_R, N)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 2, blk_T, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 3, blk_SC, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 4, blk_S, N, IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 5, blk_diff, N, IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 6, blk_PM, N, IU*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 2, 7, blk_PZ, N, IU*SQR2*RQS3)
    ! Row 4 (LH2)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 1, blk_R, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 2, blk_S, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 3, blk_Q, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 4, blk_R, N, IU*SQR2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 5, blk_S, N, IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 3, 7, blk_PM, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    ! Row 5 (SO1)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 0, blk_S, N, IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 1, blk_diff, N, -IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 2, blk_SC, N, -IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 3, blk_RC, N, -IU*SQR2)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 4, blk_temp, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 6, blk_PZ, N, IU*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 4, 7, blk_PP, N, IU*SQR2*RQS3)
    ! Row 6 (SO2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 0, blk_R, N, -IU*SQR2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 1, blk_S, N, IU*SQR3*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 2, blk_diff, N, -IU*RQS2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 3, blk_SC, N, -IU*RQS2)
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 5, blk_temp, N)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 6, blk_PM, N, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 5, 7, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))
    ! Row 7 (CB1)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 0, blk_PM, N, -IU)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 1, blk_PZ, N, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 2, blk_PP, N, -IU*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 4, blk_PZ, N, -IU*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 5, blk_PP, N, cmplx(SQR2*RQS3, 0.0_dp, kind=dp))
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 6, 6, blk_A, N)
    ! Row 8 (CB2)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 1, blk_PM, N, cmplx(-RQS3, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 2, blk_PZ, N, -IU*SQR2*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 3, blk_PP, N, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 4, blk_PM, N, -IU*SQR2*RQS3)
    call insert_csr_block_scaled(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 5, blk_PZ, N, cmplx(-RQS3, 0.0_dp, kind=dp))
    call insert_csr_block(coo_r, coo_c, coo_v, coo_cap, coo_idx, 7, 7, blk_A, N)
  end subroutine insert_main_blocks
```

**Step 2: Replace fast-path insertion (lines 345-478)**

Keep the `blk_diff%values = blk_Q%values - blk_T%values` and `blk_temp%values = ...` computation (lines 367, 417), then replace the insertion block with:

```fortran
        call insert_main_blocks(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, blk_Q, blk_T, blk_S, blk_SC, blk_R, blk_RC, blk_PZ, &
          blk_PP, blk_PM, blk_A, blk_diff, blk_temp, N)
```

**Step 3: Replace slow-path insertion (lines 654-779)**

Keep the `csr_add` calls that compute `blk_diff` and `blk_temp`, then replace with:

```fortran
        call insert_main_blocks(coo_rows, coo_cols, coo_vals, coo_capacity, &
          coo_idx, blk_Q, blk_T, blk_S, blk_SC, blk_R, blk_RC, blk_PZ, &
          blk_PP, blk_PM, blk_A, blk_diff, blk_temp, N)
```

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass. No numerical change.

**Step 5: Commit**

```bash
git add src/physics/hamiltonian_wire.f90
git commit -m "refactor: extract insert_main_blocks shared between fast/slow paths"
```

---

### Task 7: Pre-allocate `HT_csr_step` in `main.f90` kz-loop

`main.f90` lines 391-397 allocates and frees `HT_csr_step` every kz-point. Move outside the loop.

**Files:**
- Modify: `src/apps/main.f90:390-397`

**Step 1: Move `HT_csr_step` outside the loop**

Before the `do k = 2, cfg%waveVectorStep` loop, add:

```fortran
    type(csr_matrix) :: HT_csr_step
    call csr_clone_structure(HT_csr, HT_csr_step)
```

Inside the loop, replace the block:

```fortran
        call ZB8bandGeneralized(HT_csr_step, smallk(k)%kz, profile_2d, &
          & kpterms_2d, cfg, ws=wire_ws)
        call solve_sparse_evp(HT_csr_step, eigen_cfg, eigen_res, feast_ws)
```

After the loop, add:

```fortran
    call csr_free(HT_csr_step)
```

**Step 2: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

**Step 3: Commit**

```bash
git add src/apps/main.f90
git commit -m "perf: pre-allocate HT_csr_step across kz-points"
```

---

### Task 8: Extract COO finalization helper

The COO-to-CSR finalization block (20 lines) appears in both fast and slow paths. After Task 4, these will be slightly different (fast always has ws, slow may not), but the common pattern can still be extracted.

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90`

**Step 1: Write `finalize_coo_to_csr` helper**

```fortran
  subroutine finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
      coo_idx, ws)
    type(csr_matrix), intent(inout) :: HT_csr
    integer, intent(in) :: Ntot, coo_idx
    integer, intent(in) :: coo_rows(:), coo_cols(:)
    complex(kind=dp), intent(in) :: coo_vals(:)
    type(wire_workspace), intent(inout), optional :: ws

    if (present(ws)) then
      if (ws%coo_cache_valid) then
        call csr_set_values_from_coo(HT_csr, coo_idx, &
          ws%coo_to_csr(1:coo_idx), coo_vals(1:coo_idx))
      else
        call csr_build_from_coo_cached(HT_csr, Ntot, Ntot, coo_idx, &
          coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx), &
          ws%coo_to_csr)
        ws%coo_nnz_in = coo_idx
        ws%coo_cache_valid = .true.
      end if
    else
      if (coo_idx > 0) then
        call csr_build_from_coo(HT_csr, Ntot, Ntot, coo_idx, &
          coo_rows(1:coo_idx), coo_cols(1:coo_idx), coo_vals(1:coo_idx))
      else
        call csr_init(HT_csr, Ntot, Ntot)
      end if
    end if
  end subroutine finalize_coo_to_csr
```

**Step 2: Replace both finalization blocks**

In fast path and slow path, replace the 15-line if/else blocks with:

```fortran
      call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
        coo_idx, ws)
```

(For fast path, `ws` is always present. For slow path, it's optional.)

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 39 tests pass.

**Step 4: Commit**

```bash
git add src/physics/hamiltonian_wire.f90
git commit -m "refactor: extract finalize_coo_to_csr helper"
```

---

## Dependency Graph

```
Task 1 (csr_find_diag_positions)     -- independent
Task 2 (cache strain table)          -- independent
Task 3 (move CSR transpose utils)    -- independent
Task 4 (absorb coo_cache)            -- independent of 1-3
Task 5 (extract insert_g3_blocks)    -- independent
Task 6 (extract insert_main_blocks)  -- independent of 5 (but same file)
Task 7 (pre-allocate HT_csr_step)    -- independent
Task 8 (extract finalize_coo_to_csr) -- depends on Task 4
```

Recommended execution order: 1, 2, 3, 4, 8, 5, 6, 7.
Tasks 1-3 are small and safe. Task 4 is high-impact. Tasks 5-6 are large deduplications.
