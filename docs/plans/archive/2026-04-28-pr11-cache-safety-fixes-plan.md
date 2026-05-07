# PR #11 Cache Safety Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the remaining cache-safety issues from PR #11 review: prevent `wire_workspace` misuse in `g='g3'` velocity builds and make `feast_workspace` reuse validate the CSR sparsity pattern before taking the fast path.

**Architecture:** Treat `g='g3'` as a separate derivative-Hamiltonian build that must not initialize or mutate the normal wire Hamiltonian workspace. For FEAST, keep the cached upper-triangle arrays, but only reuse them when the current matrix has the exact same upper-triangle row counts and column indices as the cached structure; otherwise rebuild the cache through the existing slow path.

**Tech Stack:** Fortran 90, pFUnit, CMake/Ninja, Intel MKL FEAST.

---

### Task 1: Add a regression test for `g='g3'` with a workspace argument

**Files:**
- Modify: `tests/unit/test_hamiltonian_2d.pf`

**Step 1: Add the failing test**

Add a unit test near the existing wire workspace tests. It should verify that an accidental `ws=wire_ws` on the `g='g3'` path does not crash and does not corrupt an already initialized normal wire workspace.

```fortran
  @test
  subroutine test_workspace_g3_does_not_corrupt_normal_cache()
    integer, parameter :: nx = 4, ny = 4
    real(kind=dp), parameter :: dx = 5.0_dp, dy = 5.0_dp
    integer :: ij, k, normal_nnz, coo_nnz
    integer :: mat2d(nx, ny)
    type(spatial_grid) :: grid
    type(paramStruct), allocatable :: params(:)
    character(len=255), allocatable :: material_names(:)
    type(region_spec), allocatable :: regions(:)
    real(kind=dp), allocatable :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix) :: H_fast, H_slow, H_g3
    type(wire_workspace) :: wire_ws
    type(wire_coo_cache) :: coo_cache
    type(simulation_config) :: cfg

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
    cfg%grid = grid

    ! Initialize the normal workspace and its COO mapping.
    call ZB8bandGeneralized(H_fast, 0.05_dp, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)
    normal_nnz = H_fast%nnz
    coo_nnz = wire_ws%coo_nnz_in

    ! Accidental workspace use in g3 mode must not overwrite normal cache state.
    call ZB8bandGeneralized(H_g3, 1.0_dp, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws, g='g3')
    @assertTrue(H_g3%nnz > 0, message="g3 matrix built")
    @assertEqual(normal_nnz, H_fast%nnz, message="normal matrix retained")
    @assertEqual(coo_nnz, wire_ws%coo_nnz_in, message="normal COO mapping retained")

    ! Normal fast path must still match a fresh slow-path build after g3.
    call ZB8bandGeneralized(H_fast, 0.06_dp, profile_2d, kpterms_2d, cfg, &
      ws=wire_ws)
    call ZB8bandGeneralized(H_slow, 0.06_dp, profile_2d, kpterms_2d, cfg, &
      coo_cache)

    @assertEqual(H_slow%nrows, H_fast%nrows)
    @assertEqual(H_slow%nnz, H_fast%nnz)
    do k = 1, H_slow%nnz
      @assertEqual(H_slow%colind(k), H_fast%colind(k))
      @assertEqual(H_slow%values(k), H_fast%values(k), tolerance=1.0e-13_dp)
    end do

    call csr_free(H_fast)
    call csr_free(H_slow)
    call csr_free(H_g3)
    call wire_workspace_free(wire_ws)
    call wire_coo_cache_free(coo_cache)
    do ij = 1, size(kpterms_2d)
      call csr_free(kpterms_2d(ij))
    end do
    deallocate(kpterms_2d, profile_2d, params, regions, material_names)
    call grid_free(grid)
  end subroutine test_workspace_g3_does_not_corrupt_normal_cache
```

**Step 2: Run the test and verify current behavior**

Run:

```bash
ctest --test-dir build -R test_hamiltonian_2d --output-on-failure
```

Expected before the fix: failure or runtime error when the `g3` call tries to initialize or mutate the normal workspace path.

**Step 3: Commit only if the failing test is added in a separate TDD commit**

```bash
git add tests/unit/test_hamiltonian_2d.pf
git commit -m "test: cover g3 calls with wire workspace"
```

---

### Task 2: Make `g='g3'` ignore the normal wire workspace

**Files:**
- Modify: `src/physics/hamiltonian_wire.f90`

**Step 1: Change the `ZB8bandGeneralized` dispatch**

At the top-level dispatch in `ZB8bandGeneralized`, handle `gmode_dir == 3` before the normal workspace branch. Call the slow derivative build without passing `ws`, so it cannot clone unallocated full-Hamiltonian blocks or overwrite `ws%coo_to_csr`.

```fortran
      if (gmode_dir == 3) then
        ! g3 builds dH/dkz, whose sparsity and COO mapping differ from
        ! the normal Hamiltonian. Do not initialize or mutate wire_workspace.
        call zb8_generalized_slow(HT_csr, kz, kz2, profile_2d, kpterms_2d, &
          cfg, N=N, Ntot=Ntot, gmode_dir=gmode_dir)
      else if (present(ws) .and. ws%initialized) then
        call zb8_generalized_fast(HT_csr, kz, kz2, profile_2d, kpterms_2d, &
          cfg, ws, N, Ntot, gmode_dir)
      else
        call zb8_generalized_slow(HT_csr, kz, kz2, profile_2d, kpterms_2d, &
          cfg, coo_cache, ws, N, Ntot, gmode_dir)
      end if
```

If the compiler rejects the keyword-only call because required arguments follow optional ones, add a small dedicated helper instead:

```fortran
      call zb8_generalized_g3(HT_csr, profile_2d, kpterms_2d, cfg, N, Ntot)
```

and move the current `gmode_dir == 3` slow-path body into that helper.

**Step 2: Keep normal workspace initialization unchanged**

Do not add any `g3` structures to `wire_workspace`. The normal Hamiltonian workspace and the derivative-Hamiltonian build have different COO patterns, so sharing one mapping is the bug.

**Step 3: Run the targeted test**

Run:

```bash
cmake --build build
ctest --test-dir build -R test_hamiltonian_2d --output-on-failure
```

Expected: `test_workspace_g3_does_not_corrupt_normal_cache` passes.

**Step 4: Commit**

```bash
git add src/physics/hamiltonian_wire.f90 tests/unit/test_hamiltonian_2d.pf
git commit -m "fix: isolate g3 derivative builds from wire workspace"
```

---

### Task 3: Add a FEAST workspace sparsity-mismatch test

**Files:**
- Modify: `tests/unit/test_eigensolver.pf`

**Step 1: Add a helper to build a simple tridiagonal CSR if none exists**

Near the existing test helpers in `test_eigensolver.pf`, add:

```fortran
  subroutine build_tridiagonal_csr(n, H)
    integer, intent(in) :: n
    type(csr_matrix), intent(out) :: H
    integer, allocatable :: rows(:), cols(:)
    complex(kind=dp), allocatable :: vals(:)
    integer :: i, idx, nnz

    nnz = 3*n - 2
    allocate(rows(nnz), cols(nnz), vals(nnz))
    idx = 0
    do i = 1, n
      if (i > 1) then
        idx = idx + 1
        rows(idx) = i
        cols(idx) = i - 1
        vals(idx) = cmplx(-0.1_dp, 0.0_dp, kind=dp)
      end if
      idx = idx + 1
      rows(idx) = i
      cols(idx) = i
      vals(idx) = cmplx(real(i, kind=dp), 0.0_dp, kind=dp)
      if (i < n) then
        idx = idx + 1
        rows(idx) = i
        cols(idx) = i + 1
        vals(idx) = cmplx(-0.1_dp, 0.0_dp, kind=dp)
      end if
    end do

    call csr_build_from_coo(H, n, n, nnz, rows, cols, vals)
    deallocate(rows, cols, vals)
  end subroutine build_tridiagonal_csr
```

**Step 2: Add the failing test**

Add a test after `test_feast_workspace_cached_matches_uncached`:

```fortran
  @test
  subroutine test_feast_workspace_rebuilds_on_pattern_change()
    integer, parameter :: n = 20
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H_diag, H_tri
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: cached, uncached
    type(feast_workspace) :: fw
    integer :: i, min_nev

    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H_diag)
    call build_tridiagonal_csr(n, H_tri)

    cfg%method = 'FEAST'
    cfg%emin = 0.5_dp
    cfg%emax = 5.5_dp
    cfg%nev = 5

    ! Initialize cache with diagonal upper-triangle structure.
    call solve_feast(H_diag, cfg, cached, fw)
    call eigensolver_result_free(cached)
    @assertTrue(fw%initialized)

    ! Reuse same workspace on same-size but different sparsity matrix.
    call solve_feast(H_tri, cfg, cached, fw)
    call solve_feast(H_tri, cfg, uncached)

    @assertEqual(uncached%nev_found, cached%nev_found)
    min_nev = min(uncached%nev_found, cached%nev_found)
    do i = 1, min_nev
      @assertEqual(uncached%eigenvalues(i), cached%eigenvalues(i), tolerance=1.0e-10_dp)
    end do

    call eigensolver_result_free(cached)
    call eigensolver_result_free(uncached)
    call feast_workspace_free(fw)
    call csr_free(H_diag)
    call csr_free(H_tri)
  end subroutine test_feast_workspace_rebuilds_on_pattern_change
```

**Step 3: Run the test and verify current behavior**

Run:

```bash
ctest --test-dir build -R test_eigensolver --output-on-failure
```

Expected before the fix: failure, bounds error under checked builds, or eigenvalue mismatch.

**Step 4: Commit only if keeping TDD commits separate**

```bash
git add tests/unit/test_eigensolver.pf
git commit -m "test: cover FEAST workspace pattern changes"
```

---

### Task 4: Validate FEAST cached upper-triangle pattern before reuse

**Files:**
- Modify: `src/math/eigensolver.f90`

**Step 1: Add a private pattern-check helper**

Inside `module eigensolver`, before `solve_feast`, add:

```fortran
  logical function feast_workspace_matches_pattern(H_csr, fw, N, M0)
    type(csr_matrix), intent(in) :: H_csr
    type(feast_workspace), intent(in) :: fw
    integer, intent(in) :: N, M0

    integer :: i, j, k

    feast_workspace_matches_pattern = .false.
    if (.not. fw%initialized) return
    if (fw%N /= N .or. fw%M0 /= M0) return
    if (.not. allocated(fw%rowptr_loc)) return
    if (.not. allocated(fw%colind_loc)) return
    if (size(fw%rowptr_loc) /= N + 1) return

    do i = 1, N
      k = fw%rowptr_loc(i)
      do j = H_csr%rowptr(i), H_csr%rowptr(i+1) - 1
        if (H_csr%colind(j) >= i) then
          if (k >= fw%rowptr_loc(i + 1)) return
          if (k > fw%nnz_upper) return
          if (fw%colind_loc(k) /= H_csr%colind(j)) return
          k = k + 1
        end if
      end do
      if (k /= fw%rowptr_loc(i + 1)) return
    end do

    feast_workspace_matches_pattern = .true.
  end function feast_workspace_matches_pattern
```

Keep it private; do not add it to the module `public` list.

**Step 2: Use the helper in the fast-path condition**

Replace:

```fortran
    if (present(fw) .and. fw%initialized .and. fw%M0 == M0 .and. fw%N == N) then
```

with:

```fortran
    if (present(fw) .and. feast_workspace_matches_pattern(H_csr, fw, N, M0)) then
```

**Step 3: Rebuild stale/mismatched caches safely**

At the start of the slow-path `else` block, before allocating `rowptr_loc`, free stale cache arrays if a workspace is present but did not match:

```fortran
      if (present(fw)) then
        if (fw%initialized) call feast_workspace_free(fw)
      end if
```

This ensures the later `move_alloc(rowptr_loc, fw%rowptr_loc)` and `move_alloc(colind_loc, fw%colind_loc)` install a fresh pattern.

**Step 4: Run the targeted test**

Run:

```bash
cmake --build build
ctest --test-dir build -R test_eigensolver --output-on-failure
```

Expected: `test_feast_workspace_rebuilds_on_pattern_change` passes, along with existing FEAST workspace tests.

**Step 5: Commit**

```bash
git add src/math/eigensolver.f90 tests/unit/test_eigensolver.pf
git commit -m "fix: validate FEAST workspace sparsity before reuse"
```

---

### Task 5: Full verification and PR metadata cleanup

**Files:**
- No code changes expected.
- Optional: update PR #11 body on GitHub.

**Step 1: Run focused tests**

Run:

```bash
ctest --test-dir build -R "test_hamiltonian_2d|test_eigensolver" --output-on-failure
```

Expected: both tests pass.

**Step 2: Run full suite**

Run:

```bash
ctest --test-dir build --output-on-failure
```

Expected: 39/39 pass. If the old wire dense/sparse failure reappears, confirm it on `main` before calling it pre-existing.

**Step 3: Update PR body if needed**

The current PR body says `22/23 regression tests pass` with one pre-existing failure, but the review run on this branch produced `39/39` passing. Update the body to reflect the actual current result after these fixes.

**Step 4: Final commit if PR-body/docs changes are made locally**

```bash
git add docs/plans/2026-04-28-pr11-cache-safety-fixes-plan.md
git commit -m "docs: plan PR 11 cache safety fixes"
```

If the plan is not meant to be committed, leave it untracked with the rest of the local planning docs.
