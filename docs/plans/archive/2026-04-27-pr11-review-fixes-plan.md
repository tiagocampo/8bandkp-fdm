# PR #11 Review Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the sparsity-mismatch bug in wire_workspace fast paths (C1) plus defensive fixes (C2, C3, I2, I3, I4/I5) and add fast-path test coverage (I1).

**Architecture:** Add scatter maps (integer index arrays) to wire_workspace that map each kpterm's nonzero positions to the union-structure workspace block. Rewrite S, PP, R fast paths to scatter through these maps. Eliminate blk_temp dependency in SC/RC/PM by reading from the forward block directly.

**Tech Stack:** Fortran 90, pFUnit testing framework, MKL sparse BLAS.

---

### Task 1: Add scatter map fields to wire_workspace and build_scatter_map helper

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:38-58` (type definition)
- Modify: `src/physics/hamiltonian_wire.f90:86` (after contains, add helper)
- Modify: `src/physics/hamiltonian_wire.f90:1482-1505` (free routine)

**Step 1: Add 6 scatter map fields to wire_workspace type**

In `src/physics/hamiltonian_wire.f90`, after the `diag_pos` field (line 55), add:

```fortran
    ! Scatter maps: map each kpterm's nonzero index to the union-structure block position
    integer, allocatable :: scatter_S_14(:)   ! kp14 -> blk_S
    integer, allocatable :: scatter_S_15(:)   ! kp15 -> blk_S
    integer, allocatable :: scatter_PP_12(:)  ! kp12 -> blk_PP
    integer, allocatable :: scatter_PP_13(:)  ! kp13 -> blk_PP
    integer, allocatable :: scatter_R_16(:)   ! kp16 -> blk_R
    integer, allocatable :: scatter_R_11(:)   ! kp11 -> blk_R
```

Also **remove** the dead fields `coo_to_csr` and `coo_nnz_in` (lines 49-51).

**Step 2: Add `build_scatter_map` helper subroutine**

Add after line 86 (`contains`), before the first subroutine:

```fortran
    subroutine build_scatter_map(src, dst, scatter_map)
      type(csr_matrix), intent(in) :: src, dst
      integer, allocatable, intent(out) :: scatter_map(:)

      integer :: row, k_src, k_dst, lo, hi, mid, col_src

      allocate(scatter_map(src%nnz))
      do row = 1, src%nrows
        do k_src = src%rowptr(row), src%rowptr(row + 1) - 1
          col_src = src%colind(k_src)
          lo = dst%rowptr(row)
          hi = dst%rowptr(row + 1) - 1
          do while (lo <= hi)
            mid = (lo + hi) / 2
            if (dst%colind(mid) == col_src) then
              scatter_map(k_src) = mid
              exit
            else if (dst%colind(mid) < col_src) then
              lo = mid + 1
            else
              hi = mid - 1
            end if
          end do
          if (lo > hi) then
            print *, 'ERROR: build_scatter_map: entry not found for row', row, 'col', col_src
            stop 1
          end if
        end do
      end do
    end subroutine build_scatter_map
```

**Step 3: Add scatter map deallocation to wire_workspace_free**

In `wire_workspace_free` (line 1482), add before `ws%initialized = .false.`:

```fortran
    if (allocated(ws%scatter_S_14))  deallocate(ws%scatter_S_14)
    if (allocated(ws%scatter_S_15))  deallocate(ws%scatter_S_15)
    if (allocated(ws%scatter_PP_12)) deallocate(ws%scatter_PP_12)
    if (allocated(ws%scatter_PP_13)) deallocate(ws%scatter_PP_13)
    if (allocated(ws%scatter_R_16))  deallocate(ws%scatter_R_16)
    if (allocated(ws%scatter_R_11))  deallocate(ws%scatter_R_11)
```

Remove the dead field deallocation:
```fortran
    if (allocated(ws%coo_to_csr)) deallocate(ws%coo_to_csr)   ! REMOVE
    ws%coo_nnz_in = 0                                           ! REMOVE
```

**Step 4: Build scatter maps during workspace initialization**

In `zb8_generalized_slow`, after the `diag_pos` scan loop (line 806) and before `ws%coo_capacity = coo_capacity` (line 808), add:

```fortran
    call build_scatter_map(kpterms_2d(14), ws%blk_S, ws%scatter_S_14)
    call build_scatter_map(kpterms_2d(15), ws%blk_S, ws%scatter_S_15)
    call build_scatter_map(kpterms_2d(12), ws%blk_PP, ws%scatter_PP_12)
    call build_scatter_map(kpterms_2d(13), ws%blk_PP, ws%scatter_PP_13)
    call build_scatter_map(kpterms_2d(16), ws%blk_R, ws%scatter_R_16)
    call build_scatter_map(kpterms_2d(11), ws%blk_R, ws%scatter_R_11)
```

**Step 5: Build and verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: clean build, no errors.

**Step 6: Commit**

```
fix: add scatter map fields and builder to wire_workspace
```

---

### Task 2: Write failing fast-path vs slow-path comparison test

**Files:**
- Modify: `tests/unit/test_hamiltonian_2d.pf` (add test after `test_workspace_slow_path_reproducible`)

**Step 1: Add the test**

After `test_workspace_slow_path_reproducible` (before `end module`), add:

```fortran
  @test
  subroutine test_workspace_fast_path_matches_slow()
    ! Build wire Hamiltonian at kz=0.05 via slow path, then at kz=0.06
    ! via fast path (workspace initialized from first call).
    ! Verify the fast-path result matches a fresh slow-path build at kz=0.06.
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ngrid, ij, k
    integer :: mat2d(nx, ny)
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H_slow, H_fast
    type(wire_coo_cache) :: coo_cache
    type(wire_workspace) :: wire_ws
    real(kind=dp) :: kz1, kz2
    type(simulation_config) :: cfg

    ngrid = nx * ny
    mat2d = 1

    call grid_init_rect(grid, nx, ny, dx, dy, mat2d)
    allocate(material_names(1))
    material_names(1) = "GaAs"
    allocate(params(1))
    call paramDatabase(material_names, 1, params)
    allocate(regions(1))
    regions(1)%material = "GaAs"

    call confinementInitialization_2d(grid, params, regions, &
      profile_2d, kpterms_2d, FDorder=2)

    kz1 = 0.05_dp
    kz2 = 0.06_dp
    cfg%grid = grid

    ! First call at kz1: initializes workspace (slow path)
    call ZB8bandGeneralized(H_slow, kz1, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)
    call csr_free(H_slow)

    ! Second call at kz2: uses fast path
    call ZB8bandGeneralized(H_fast, kz2, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)

    ! Fresh slow-path at kz2 for comparison
    coo_cache%initialized = .false.
    call ZB8bandGeneralized(H_slow, kz2, profile_2d, kpterms_2d, cfg, &
      coo_cache)

    ! Compare: must match element-by-element
    @assertEqual(H_slow%nrows, H_fast%nrows, message="nrows match")
    @assertEqual(H_slow%nnz, H_fast%nnz, message="nnz match")
    do k = 1, H_slow%nnz
      @assertEqual(H_slow%colind(k), H_fast%colind(k), message="colind match")
    end do
    do k = 1, H_slow%nnz
      @assertEqual(H_slow%values(k), H_fast%values(k), tolerance=1.0e-13_dp, &
        message="values match at index " // trim(int_to_str(k)))
    end do

    call csr_free(H_slow)
    call csr_free(H_fast)
    call wire_coo_cache_free(coo_cache)
    call wire_workspace_free(wire_ws)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_workspace_fast_path_matches_slow
```

Note: `int_to_str` may not exist; if not, replace the message with a plain string like `"values match"`.

**Step 2: Build and run the test to verify it FAILS**

Run: `cmake --build build && ctest --test-dir build -R test_hamiltonian_2d -V 2>&1 | tail -20`
Expected: test_workspace_fast_path_matches_slow FAILS (values mismatch or crash due to corrupted CSR).

**Step 3: Commit (failing test)**

```
test: add fast-path vs slow-path comparison (currently fails)
```

---

### Task 3: Fix build_kp_term_S fast path with scatter maps

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:1028-1031` (S fast path)

**Step 1: Replace element-wise arithmetic with scatter**

Replace lines 1028-1031:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: kp14 and kp15 have gradient sparsity (same pattern)
        blk%values = cmplx(-2.0_dp*SQR3*kz, 0.0_dp, kind=dp) * kpterms_2d(14)%values &
                   + cmplx(0.0_dp, 2.0_dp*SQR3*kz, kind=dp) * kpterms_2d(15)%values
```

With:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: scatter kp14 and kp15 into union-structure block
        blk%values = cmplx(0.0_dp, 0.0_dp, kind=dp)
        blk%values(ws%scatter_S_14) = cmplx(-2.0_dp*SQR3*kz, 0.0_dp, kind=dp) &
          * kpterms_2d(14)%values
        blk%values(ws%scatter_S_15) = blk%values(ws%scatter_S_15) &
          + cmplx(0.0_dp, 2.0_dp*SQR3*kz, kind=dp) * kpterms_2d(15)%values
```

**Step 2: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_hamiltonian_2d -V 2>&1 | tail -20`
Expected: still fails (PP and R not yet fixed).

**Step 3: Commit**

```
fix: use scatter maps for build_kp_term_S fast path
```

---

### Task 4: Fix build_kp_term_PP fast path with scatter maps

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:1211-1214` (PP fast path)

**Step 1: Replace element-wise arithmetic with scatter**

Replace lines 1211-1214:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: kp12 and kp13 have gradient sparsity (same pattern)
        blk%values = cmplx(0.0_dp, RQS2, kind=dp) * kpterms_2d(12)%values &
                   + cmplx(-RQS2, 0.0_dp, kind=dp) * kpterms_2d(13)%values
```

With:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: scatter kp12 and kp13 into union-structure block
        blk%values = cmplx(0.0_dp, 0.0_dp, kind=dp)
        blk%values(ws%scatter_PP_12) = cmplx(0.0_dp, RQS2, kind=dp) &
          * kpterms_2d(12)%values
        blk%values(ws%scatter_PP_13) = blk%values(ws%scatter_PP_13) &
          + cmplx(-RQS2, 0.0_dp, kind=dp) * kpterms_2d(13)%values
```

**Step 2: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_hamiltonian_2d -V 2>&1 | tail -20`

**Step 3: Commit**

```
fix: use scatter maps for build_kp_term_PP fast path
```

---

### Task 5: Fix build_kp_term_R fast path with scatter maps

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:1109-1120` (R fast path)

**Step 1: Replace element-wise arithmetic with scatter**

Replace lines 1109-1120:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: operator part from kp16 and kp11
        blk%values = cmplx(-SQR3, 0.0_dp, kind=dp) * kpterms_2d(16)%values &
                   + cmplx(0.0_dp, 2.0_dp*SQR3, kind=dp) * kpterms_2d(11)%values
        ! Add diagonal contribution (kp2) at diagonal positions
        block
          integer :: k
          do k = 1, kpterms_2d(2)%nnz
            blk%values(ws%diag_pos(k)) = blk%values(ws%diag_pos(k)) &
              + cmplx(-SQR3*kz2, 0.0_dp, kind=dp) * kpterms_2d(2)%values(k)
          end do
        end block
```

With:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: scatter kp16 and kp11 into union-structure block
        blk%values = cmplx(0.0_dp, 0.0_dp, kind=dp)
        blk%values(ws%scatter_R_16) = cmplx(-SQR3, 0.0_dp, kind=dp) &
          * kpterms_2d(16)%values
        blk%values(ws%scatter_R_11) = blk%values(ws%scatter_R_11) &
          + cmplx(0.0_dp, 2.0_dp*SQR3, kind=dp) * kpterms_2d(11)%values
        ! Add diagonal contribution (kp2) at diagonal positions
        block
          integer :: k
          do k = 1, kpterms_2d(2)%nnz
            blk%values(ws%diag_pos(k)) = blk%values(ws%diag_pos(k)) &
              + cmplx(-SQR3*kz2, 0.0_dp, kind=dp) * kpterms_2d(2)%values(k)
          end do
        end block
```

**Step 2: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_hamiltonian_2d -V 2>&1 | tail -20`
Expected: test_workspace_fast_path_matches_slow now PASSES.

**Step 3: Commit**

```
fix: use scatter maps for build_kp_term_R fast path
```

---

### Task 6: Fix SC/RC/PM conjugate-transpose terms to read from forward blocks

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:1079-1082` (SC fast path)
- Modify: `src/physics/hamiltonian_wire.f90:1161-1164` (RC fast path)
- Modify: `src/physics/hamiltonian_wire.f90:1246-1249` (PM fast path)

**Step 1: Fix build_kp_term_SC fast path**

Replace lines 1079-1082:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: build S into temp, then conjugate transpose to blk
        call build_kp_term_S(kz, kpterms_2d, ws%blk_temp, ws)
        call csr_conjugate_transpose_to_preallocated(ws%blk_temp, blk)
```

With:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: conjugate transpose from already-built blk_S
        call csr_conjugate_transpose_to_preallocated(ws%blk_S, blk)
```

**Step 2: Fix build_kp_term_RC fast path**

Replace lines 1161-1164:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: build R into temp, then conjugate transpose to blk
        call build_kp_term_R(kz2, kpterms_2d, ws%blk_temp, ws)
        call csr_conjugate_transpose_to_preallocated(ws%blk_temp, blk)
```

With:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: conjugate transpose from already-built blk_R
        call csr_conjugate_transpose_to_preallocated(ws%blk_R, blk)
```

**Step 3: Fix build_kp_term_PM fast path**

Replace lines 1246-1249:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: build PP into temp, then conjugate transpose to blk
        call build_kp_term_PP(kz, kpterms_2d, ws%blk_temp, ws)
        call csr_conjugate_transpose_to_preallocated(ws%blk_temp, blk)
```

With:

```fortran
      if (present(ws) .and. ws%initialized) then
        ! Fast path: conjugate transpose from already-built blk_PP
        call csr_conjugate_transpose_to_preallocated(ws%blk_PP, blk)
```

**Step 4: Build and run all tests**

Run: `cmake --build build && ctest --test-dir build -L unit -V 2>&1 | tail -20`
Expected: all 15 unit tests pass.

**Step 5: Commit**

```
fix: eliminate blk_temp dependency in SC/RC/PM conjugate-transpose terms
```

---

### Task 7: Initialize diag_pos and harden COO overflow checks

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:797-806` (diag_pos init)
- Modify: `src/physics/hamiltonian_wire.f90:1330-1333` (insert_csr_block overflow)
- Modify: `src/physics/hamiltonian_wire.f90:1379-1383` (insert_csr_block_scaled overflow)
- Modify: `src/physics/hamiltonian_wire.f90:1417-1421` (insert_profile_diagonal overflow)
- Modify: `src/physics/hamiltonian_wire.f90:1447` (insert_strain_coo overflow)

**Step 1: Initialize diag_pos to 0 and add assertion**

Replace the diag_pos allocation and scan (lines 797-806) with:

```fortran
    allocate(ws%diag_pos(N))
    ws%diag_pos = 0
    do i = 1, N
      do k = ws%blk_Q%rowptr(i), ws%blk_Q%rowptr(i+1) - 1
        if (ws%blk_Q%colind(k) == i) then
          ws%diag_pos(i) = k
          exit
        end if
      end do
    end do
    if (any(ws%diag_pos == 0)) then
      print *, 'ERROR: diag_pos: no diagonal entry found for some rows'
      stop 1
    end if
```

**Step 2: Replace COO overflow WARNING with stop 1**

In `insert_csr_block` (around line 1330), replace:
```fortran
          print *, "WARNING: COO capacity exceeded in insert_csr_block, skipping entries"
          return
```
With:
```fortran
          print *, "ERROR: COO capacity exceeded in insert_csr_block"
          stop 1
```

Same replacement in `insert_csr_block_scaled` (line 1380-1382), `insert_profile_diagonal` (line 1418-1420), and `insert_strain_coo` (line 1447 — currently silent, add error message).

For `insert_strain_coo` line 1447, replace:
```fortran
          if (coo_idx > coo_cap) return
```
With:
```fortran
          if (coo_idx > coo_cap) then
            print *, "ERROR: COO capacity exceeded in insert_strain_coo"
            stop 1
          end if
```

**Step 3: Build and run tests**

Run: `cmake --build build && ctest --test-dir build -L unit 2>&1 | tail -5`
Expected: all pass.

**Step 4: Commit**

```
fix: initialize diag_pos, harden COO overflow to stop 1
```

---

### Task 8: Fix prev_wire_eval assignment guard (I2 / Codex P2)

**Files:**
- Modify: `src/apps/main.f90:369-376`

**Step 1: Move assignments inside the nev_found guard**

Replace lines 369-376:

```fortran
    if (allocated(prev_wire_eval)) deallocate(prev_wire_eval)
    if (allocated(prev_wire_evec)) deallocate(prev_wire_evec)
    if (eigen_res%nev_found > 0) then
      allocate(prev_wire_eval(eigen_res%nev_found))
      allocate(prev_wire_evec(Ntot, eigen_res%nev_found))
    end if
    prev_wire_eval = eigen_res%eigenvalues
    prev_wire_evec = eigen_res%eigenvectors
```

With:

```fortran
    if (allocated(prev_wire_eval)) deallocate(prev_wire_eval)
    if (allocated(prev_wire_evec)) deallocate(prev_wire_evec)
    if (eigen_res%nev_found > 0) then
      allocate(prev_wire_eval(eigen_res%nev_found))
      allocate(prev_wire_evec(Ntot, eigen_res%nev_found))
      prev_wire_eval = eigen_res%eigenvalues
      prev_wire_evec = eigen_res%eigenvectors
    end if
```

**Step 2: Build**

Run: `cmake --build build 2>&1 | tail -3`
Expected: clean build.

**Step 3: Commit**

```
fix: guard prev_wire_eval assignment against nev_found=0
```

---

### Task 9: Harden feast_workspace with N validation (I4/I5)

**Files:**
- Modify: `src/math/eigensolver.f90:50-56` (type definition)
- Modify: `src/math/eigensolver.f90:147` (fast-path check)
- Modify: `src/math/eigensolver.f90` around line 199 (cache population)

**Step 1: Add N field to feast_workspace**

After `M0 = 0` (line 54), add:

```fortran
    integer                       :: N = 0
```

**Step 2: Add N check to fast-path condition**

Replace line 147:

```fortran
    if (present(fw) .and. fw%initialized .and. fw%M0 == M0) then
```

With:

```fortran
    if (present(fw) .and. fw%initialized .and. fw%M0 == M0 .and. fw%N == N) then
```

**Step 3: Store N when caching**

In the slow path where `fw%nnz_upper` and `fw%M0` are set (around line 199), add:

```fortran
    fw%N = N
```

Also add to `feast_workspace_free`: `fw%N = 0`.

**Step 4: Build and run feast tests**

Run: `cmake --build build && ctest --test-dir build -R test_eigensolver -V 2>&1 | tail -10`
Expected: all eigensolver tests pass.

**Step 5: Commit**

```
fix: add N validation to feast_workspace cache
```

---

### Task 10: Small cleanups (suggestions from review)

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90:1452-1460` (add case default)
- Modify: `src/physics/hamiltonian_wire.f90:150` (remove dead nnz_est from outer scope — check if it exists there)
- Modify: `src/physics/hamiltonian_wire.f90:808` (guard double-init)

**Step 1: Add case default to insert_strain_coo**

After `case (7); field_val = bp%QT2_eps(ii)` (line 1459), add:

```fortran
    case default
      error stop "insert_strain_coo: invalid field_id"
```

**Step 2: Guard workspace double-init**

Before `ws%initialized = .true.` (line 809), add:

```fortran
    if (ws%initialized) then
      print *, 'ERROR: workspace already initialized in slow path'
      stop 1
    end if
```

**Step 3: Build and run all tests**

Run: `cmake --build build && ctest --test-dir build -L unit 2>&1 | tail -5`
Expected: all pass.

**Step 4: Commit**

```
refactor: add safety guards to strain COO and workspace init
```

---

### Task 11: Extend strain table test coverage (I6)

**Files:**
- Modify: `tests/unit/test_strain_solver.pf:626-639`

**Step 1: Add spot-checks for R_eps and VB-SO coupling entries**

After the existing spot-check for entry 9 (line 637), add:

```fortran
    ! Entry 13 (first R_eps entry): band 0->2, R_eps(1)
    @assertEqual(1, coo_r(13))
    @assertEqual(2*N + 1, coo_c(13))
    @assertEqual(bp%R_eps(1), coo_v(13), tolerance=1.0e-14_dp)

    ! Entry 17 (first VB-SO coupling): complex prefactor * conjg(S_eps(1))
    @assertEqual(1, coo_r(17))
    @assertEqual(4*N + 1, coo_c(17))
    ! VB-SO entry uses -IU*RQS2*conjg(S_eps) — verify it's nonzero and complex
    @assertTrue(abs(coo_v(17)) > 0.0_dp, message="VB-SO entry is nonzero")
```

**Step 2: Run test**

Run: `cmake --build build && ctest --test-dir build -R test_strain -V 2>&1 | tail -10`
Expected: pass.

**Step 3: Commit**

```
test: extend strain COO table coverage to R_eps and VB-SO entries
```

---

### Task 12: Run full regression suite and regenerate wire golden data

**Files:**
- Potentially: `tests/regression/data/wire_gaas_rectangle*`, `wire_inas_gaas_strain*`, etc.

**Step 1: Run all tests**

Run: `ctest --test-dir build -V 2>&1 | tail -30`

**Step 2: If wire regression tests fail, regenerate golden data**

The old fast-path produced wrong results due to C1. The golden data for wire configs needs regeneration. Run each wire config through the fixed code and update the reference data.

**Step 3: Verify all tests pass**

Run: `ctest --test-dir build 2>&1 | tail -10`
Expected: all tests pass (including the previously pre-existing failure if it was related).

**Step 4: Commit**

```
test: regenerate wire regression golden data after scatter-map fix
```
