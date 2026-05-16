# Topological & Magnetic Bug Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 5 open Codex review findings in topological analysis and magnetic field code paths, with new unit/regression tests for each fix.

**Architecture:** Five independent bug fixes applied sequentially. Each task follows TDD: write failing test first, fix, verify pass. No physics redesign — only correcting index calculations, loop structure, and coordinate extraction. Order by priority: P1 Peierls first, then BHZ hopping (most impactful P2), then remaining P2s.

**Tech Stack:** Fortran 2018, pFUnit, CMake/CTest, Python (regression verification).

---

## Task 1: Fix Peierls y-axis indexing for multi-column wires (F1, P1)

**Files:**
- Modify: `src/physics/magnetic_field.f90:86-91`
- Test: `tests/unit/test_magnetic_field.pf`

- [ ] **Step 1: Write the failing test**

In `tests/unit/test_magnetic_field.pf`, add a test that creates a 2D wire grid with `nx=3, ny=5`, builds a small off-diagonal COO entry between sites at different y-positions, and verifies that `add_peierls_coo` uses the correct y-coordinate (from `grid%coords(2,:)`) rather than a flat index into `grid%z`:

```fortran
@test
subroutine test_peierls_y_axis_2d_wire()
  type(spatial_grid) :: grid
  complex(kind=dp), allocatable :: coo_vals(:)
  integer, allocatable :: coo_row(:), coo_col(:)
  real(kind=dp) :: B_vec(3)
  integer :: nnz, i
  real(kind=dp) :: y_expected_i, y_expected_j, phase_expected
  complex(kind=dp) :: val_expected
  real(kind=dp) :: hbar_J

  ! 2D wire: 3 columns x 5 rows
  grid%ndim = 2
  grid%nx = 3
  grid%ny = 5
  grid%dx = 2.0_dp
  grid%dy = 2.0_dp
  allocate(grid%z(5))
  grid%z = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  allocate(grid%coords(2, 15))
  ! coords(1,:) = x, coords(2,:) = y for each of 15 sites
  do i = 1, 15
    grid%coords(1, i) = grid%dx * real(mod(i-1, 3), kind=dp)
    grid%coords(2, i) = grid%z(mod(i-1, 5) + 1)
  end do

  ! One off-diagonal entry: site 0 (x=0,y=1) to site 7 (x=1,y=3)
  ! In band-major with ngrid=15: flat site for row=1 is site 1, for row=8 is site 8
  nnz = 1
  allocate(coo_vals(1), coo_row(1), coo_col(1))
  coo_row(1) = 1  ! band 1, site 1 -> grid point 1
  coo_col(1) = 8  ! band 1, site 8 -> grid point 8
  coo_vals(1) = cmplx(1.0_dp, 0.0_dp, kind=dp)

  B_vec = [1.0_dp, 0.0_dp, 0.0_dp]  ! Bx = 1 T

  call add_peierls_coo(coo_vals, coo_row, coo_col, nnz, grid, B_vec)

  ! Verify the phase was computed using correct y-coordinates
  ! Site 1 -> y = grid%coords(2, 1) = 1.0 AA
  ! Site 8 -> y = grid%coords(2, 8) = 3.0 AA  (NOT grid%z(8) which is out of bounds)
  hbar_J = hbar * e
  phase_expected = e * 1.0_dp * (1.0e-10_dp - 3.0e-10_dp) * (grid%dy * 1.0e-10_dp) / hbar_J
  val_expected = cmplx(cos(phase_expected), -sin(phase_expected), kind=dp)

  @assertEqual(real(val_expected), real(coo_vals(1)), tolerance=1.0e-12_dp)
  @assertEqual(aimag(val_expected), aimag(coo_vals(1)), tolerance=1.0e-12_dp)

  deallocate(grid%z, grid%coords, coo_vals, coo_row, coo_col)
end subroutine test_peierls_y_axis_2d_wire
```

- [ ] **Step 2: Run test to verify it fails or exposes the bug**

Run: `OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_magnetic_field --output-on-failure`
Expected: The new test may crash with out-of-bounds access on `grid%z(8)` (since `grid%z` has length 5), or produce wrong phase values.

- [ ] **Step 3: Fix the Peierls y-axis indexing**

In `src/physics/magnetic_field.f90`, replace lines 86-91. Change the y-coordinate extraction to use `grid%coords` for 2D grids when available, and otherwise derive the y-axis index from the flat spatial site before indexing `grid%z`. This implements the spec's `iy = mod(flat_site - 1, ny) + 1` fallback and avoids indexing `grid%z` with a flat `nx*ny` site.

Replace the entire Peierls phase application block (lines 85-103) with:

```fortran
    ! Apply Peierls phase to all off-diagonal entries
    do idx = 1, nnz_offset
      if (coo_row(idx) /= coo_col(idx)) then
        block
          integer :: flat_i, flat_j, iy_i, iy_j

          flat_i = mod(coo_row(idx) - 1, ngrid) + 1
          flat_j = mod(coo_col(idx) - 1, ngrid) + 1

          ! Extract y-coordinate for each site. In 2D wire mode, coords(2,:)
          ! contains the y-coordinate for every flat spatial site. If coords
          ! is unavailable, grid%z has length ny, so convert flat site -> iy.
          if (grid%ndim == 2 .and. allocated(grid%coords)) then
            y_i = grid%coords(2, flat_i) * 1.0e-10_dp
            y_j = grid%coords(2, flat_j) * 1.0e-10_dp
          else
            iy_i = mod(flat_i - 1, grid%ny) + 1
            iy_j = mod(flat_j - 1, grid%ny) + 1
            y_i = grid%z(iy_i) * 1.0e-10_dp
            y_j = grid%z(iy_j) * 1.0e-10_dp
          end if
        end block

        ! Peierls phase: phi = e * Bx * (y_i - y_j) * dz / hbar
        phase = e * Bx * (y_i - y_j) * dy_m / hbar_J

        ! Phase factor: exp(-i * phase)
        exp_phase = cmplx(cos(phase), -sin(phase), kind=dp)

        ! Multiply the off-diagonal entry by the full complex phase factor
        coo_vals(idx) = coo_vals(idx) * exp_phase
      end if
    end do
```

- [ ] **Step 4: Build and run the test**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_magnetic_field --output-on-failure`
Expected: All magnetic field tests pass, including the new 2D wire test.

- [ ] **Step 5: Commit**

```bash
git add src/physics/magnetic_field.f90 tests/unit/test_magnetic_field.pf
git commit -m "fix: use grid coords for Peierls y-axis in multi-column wires"
```

---

## Task 2: Fix BHZ wire forward hopping (F2, P2)

**Files:**
- Modify: `src/physics/topological_analysis.f90:767-793`
- Test: `tests/unit/test_edge_states.pf`

- [ ] **Step 1: Write the failing test**

Add a test to `tests/unit/test_edge_states.pf` that verifies the BHZ wire Hamiltonian has correct off-diagonal coupling structure — specifically that forward A-term hopping connects site `i` to site `i+1`, not to a different orbital on site `i`:

```fortran
@test
subroutine test_bhz_forward_hopping_connects_neighbors()
  type(bhz_wire_params) :: params
  type(csr_matrix) :: H_csr
  real(kind=dp) :: val_fwd_A, val_fwd_BD
  integer :: i, j, col_found
  logical :: found

  params%A = 364.5_dp
  params%B = -686.0_dp
  params%D = -512.0_dp
  params%M = 10.0_dp
  params%d_wire = 58.0_dp
  params%N = 10
  params%dz = params%d_wire / real(params%N, kind=dp)

  call build_bhz_wire_hamiltonian(H_csr, params)

  ! Check site 1, band 1 forward coupling:
  ! Row (1-1)*4 + 1 = 1. Forward A-term should connect to site 2, not site 1.
  ! Site 2 starts at index (2-1)*4 + 1 = 5.
  ! So H(1, 5) should be nonzero (A-term forward), and H(1, 4) should be zero
  ! (old bug was connecting to same-site orbital 4 via 5-row=4).
  found = .false.
  do j = H_csr%rowptr(1), H_csr%rowptr(2) - 1
    if (H_csr%colind(j) == 5) then
      found = .true.
      exit
    end if
  end do
  @assertTrue(found, message="Forward A-term from site 1 should connect to site 2 (col=5)")

  ! Verify H(1, 4) is zero — old bug connected here
  found = .false.
  do j = H_csr%rowptr(1), H_csr%rowptr(2) - 1
    if (H_csr%colind(j) == 4) then
      found = .true.
      exit
    end if
  end do
  @assertFalse(found, message="Forward A-term should NOT connect site 1 band 1 to same-site band 4 (col=4)")

  call csr_free(H_csr)
end subroutine test_bhz_forward_hopping_connects_neighbors
```

- [ ] **Step 2: Run test to verify it fails**

Run: `OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_edge_states --output-on-failure`
Expected: Test fails because forward A-term connects to same-site `(i-1)*4 + 4` instead of neighbor `i*4 + 1`.

- [ ] **Step 3: Fix the BHZ forward hopping**

In `src/physics/topological_analysis.f90`, change lines 767-793.

**Forward A-term** (lines 767-776): Change column from `(i-1)*4 + 5 - row` to `i*4 + 5 - row`:

```fortran
        if (i < N) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = i*4 + 5 - row
          if (row == 1 .or. row == 4) then
            coo_vals(nnz_offset) = cmplx(A_over_2dz, 0.0_dp, kind=dp)
          else
            coo_vals(nnz_offset) = cmplx(-A_over_2dz, 0.0_dp, kind=dp)
          end if
        end if
```

**Forward B+D term** (lines 789-794): Change column from `(i-1)*4 + mod(row,4) + 1` to `i*4 + mod(row,4) + 1`:

```fortran
        if (i < N) then
          nnz_offset = nnz_offset + 1
          coo_row(nnz_offset) = (i-1)*4 + row
          coo_col(nnz_offset) = i*4 + mod(row, 4) + 1
          coo_vals(nnz_offset) = cmplx(B_plus_D_over_dz2, 0.0_dp, kind=dp)
        end if
```

- [ ] **Step 4: Build and run tests**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R "test_edge_states|test_z2_invariant|regression_topology_bhz" --output-on-failure`
Expected: New test passes. Existing BHZ Z2 tests still pass (trivial Z2=0, topological Z2=1). If Z2 results change, the fix is still correct — the old Hamiltonian was wrong and the new results are the physically correct ones.

- [ ] **Step 5: Run full suite to check for regressions**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure`
Expected: 66/66 pass (or existing BHZ tests may change values — verify physical correctness).

- [ ] **Step 6: Commit**

```bash
git add src/physics/topological_analysis.f90 tests/unit/test_edge_states.pf
git commit -m "fix: BHZ wire forward hopping must connect neighbor sites, not same-site orbitals"
```

---

## Task 3: Fix Majorana fit coordinates for 2D wires (F3, P2)

**Files:**
- Modify: `src/physics/topological_analysis.f90:655-662`
- Test: `tests/unit/test_edge_states.pf`

- [ ] **Step 1: Write the failing test**

Add a test that calls `compute_majorana_profile` with a 2D grid and verifies the coordinate array is sized correctly:

```fortran
@test
subroutine test_majorana_coords_2d_wire()
  type(spatial_grid) :: grid
  complex(kind=dp), allocatable :: evec(:)
  real(kind=dp) :: xi
  integer :: i, nspatial

  ! 2D wire: 3x5 = 15 spatial points
  grid%ndim = 2
  grid%nx = 3
  grid%ny = 5
  grid%dx = 2.0_dp
  grid%dy = 2.0_dp
  nspatial = 15
  allocate(grid%coords(2, nspatial))
  do i = 1, nspatial
    grid%coords(1, i) = real(mod(i-1, 3), kind=dp) * grid%dx
    grid%coords(2, i) = real((i-1)/3, kind=dp) * grid%dy
  end do

  ! Create a fake BdG eigenvector with a peak in the middle
  allocate(evec(2 * nspatial))
  evec = cmplx(0.0_dp, 0.0_dp, kind=dp)
  ! Peak at site 7 (middle of 3x5 grid)
  evec(7) = cmplx(1.0_dp, 0.0_dp, kind=dp)
  evec(nspatial + 7) = cmplx(1.0_dp, 0.0_dp, kind=dp)

  ! Should not crash with out-of-bounds on grid%z (which is not allocated)
  xi = compute_majorana_profile(evec, grid, 0.1_dp, nspatial)

  ! xi should be a valid number (not NaN or Inf)
  @assertTrue(xi >= 0.0_dp .or. xi == -1.0_dp, message="xi should be non-negative or -1 (fit failure)")

  deallocate(grid%coords, evec)
end subroutine test_majorana_coords_2d_wire
```

- [ ] **Step 2: Run test to verify it fails or crashes**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_edge_states --output-on-failure`
Expected: Test crashes or fails because `grid%z` is not allocated for 2D grids, or `xx = grid%z` produces wrong-size array.

- [ ] **Step 3: Fix the Majorana coordinate extraction**

In `src/physics/topological_analysis.f90`, replace lines 655-662 with:

```fortran
    if (grid%ndim == 2 .and. allocated(grid%coords)) then
      ! 2D wire: use y-coordinate at each spatial point
      allocate(xx(nspatial))
      xx(1:nspatial) = grid%coords(2, 1:nspatial)
    else if (allocated(grid%z)) then
      xx = grid%z
    else
      allocate(xx(nspatial))
      do i = 1, nspatial
        xx(i) = real(i - 1, kind=dp)
      end do
    end if
```

- [ ] **Step 4: Build and run tests**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_edge_states --output-on-failure`
Expected: All edge state tests pass, including new 2D wire test.

- [ ] **Step 5: Commit**

```bash
git add src/physics/topological_analysis.f90 tests/unit/test_edge_states.pf
git commit -m "fix: use grid coords for Majorana localization fit in 2D wires"
```

---

## Task 4: Fix QW Fu-Kane sweep to evaluate per grid point (F4, P2)

**Files:**
- Modify: `src/apps/main_topology.f90:1067-1105`
- Test: `tests/regression/configs/topology_qw_fukane.cfg` (verify sweep output changes)

- [ ] **Step 1: Understand the current broken behavior**

The current `compute_qw_fukane_gap_sweep` calls `compute_z2_fukane_qw_result` once at line 1087, then fills all grid points with the same values at lines 1095-1099. The wire BdG sweep at lines 1129-1135 shows the correct pattern: evaluate inside the loop with per-point parameters.

- [ ] **Step 2: Write the fix**

Replace `compute_qw_fukane_gap_sweep` (lines 1067-1105) with a per-point evaluation that follows the wire BdG sweep pattern. The Fu-Kane Z2 invariant depends on the QW Hamiltonian which in turn depends on Zeeman (B-field) and chemical potential (mu). Create a local config copy per point, set `cfg_local%bdg%enabled = .true.` so the Hamiltonian consumes the swept B-field, and set the swept chemical potential on `cfg_local%bdg%mu`.

```fortran
  subroutine compute_qw_fukane_gap_sweep(cfg_in, profile_in, kpterms_in, gap_threshold, &
      & z2_map, gap_map, transitions)
    type(simulation_config), intent(in) :: cfg_in
    real(kind=dp), contiguous, intent(in) :: profile_in(:,:)
    real(kind=dp), contiguous, intent(in) :: kpterms_in(:,:,:)
    real(kind=dp), intent(in) :: gap_threshold
    integer, allocatable, intent(out) :: z2_map(:,:)
    real(kind=dp), allocatable, intent(out) :: gap_map(:,:), transitions(:,:)

    integer :: iB, iMu, nB, nMu, z2_status, z2_val, n_occ
    real(kind=dp) :: min_gap, B_val, mu_val, dB, dmu
    type(simulation_config) :: cfg_local

    nB = cfg_in%topo%gap_sweep_nB
    nMu = cfg_in%topo%gap_sweep_nMu
    if (nB < 1 .or. nMu < 1) then
      allocate(z2_map(0,0), gap_map(0,0), transitions(0,2))
      return
    end if

    n_occ = 6 * size(profile_in, 1)
    allocate(z2_map(nMu, nB), gap_map(nMu, nB))

    dB = 0.0_dp
    if (nB > 1) dB = (cfg_in%topo%gap_sweep_B_max - cfg_in%topo%gap_sweep_B_min) / real(nB - 1, kind=dp)
    dmu = 0.0_dp
    if (nMu > 1) dmu = (cfg_in%topo%gap_sweep_mu_max - cfg_in%topo%gap_sweep_mu_min) / real(nMu - 1, kind=dp)

    do iB = 1, nB
      B_val = cfg_in%topo%gap_sweep_B_min + real(iB - 1, kind=dp) * dB
      do iMu = 1, nMu
        mu_val = cfg_in%topo%gap_sweep_mu_min + real(iMu - 1, kind=dp) * dmu

        cfg_local = cfg_in
        cfg_local%bdg%enabled = .true.
        cfg_local%bdg%B_vec = [0.0_dp, 0.0_dp, B_val]
        cfg_local%bdg%mu = mu_val

        call compute_z2_fukane_qw_result(cfg_local, profile_in, kpterms_in, n_occ, &
          & z2_val, min_gap, z2_status)
        if (z2_status /= 0) then
          print *, 'WARNING: QW Fu-Kane sweep failed at B=', B_val, ' mu=', mu_val, &
            & ' status=', z2_status
          z2_val = 0
          min_gap = huge(1.0_dp)
        end if

        z2_map(iMu, iB) = z2_val
        gap_map(iMu, iB) = min_gap
      end do
    end do

    call detect_z2_transitions(z2_map, gap_map, cfg_in%topo%gap_sweep_B_min, &
      & cfg_in%topo%gap_sweep_B_max, cfg_in%topo%gap_sweep_mu_min, &
      & cfg_in%topo%gap_sweep_mu_max, gap_threshold, transitions)
  end subroutine compute_qw_fukane_gap_sweep
```

- [ ] **Step 3: Build and verify**

Run: `cmake --build build`
Expected: Successful build.

- [ ] **Step 4: Run existing sweep regression tests**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 -R "regression_topology" --output-on-failure`
Expected: All topology regression tests pass. The Fu-Kane sweep test may produce different (correct) results — verify physical reasonableness.

- [ ] **Step 5: Manually verify the QW Fu-Kane sweep is no longer constant**

Run a QW Fu-Kane sweep config and inspect the output grid:

```bash
python3 - <<'PY'
from pathlib import Path
Path("input.cfg").write_text(Path("tests/regression/configs/topology_qw_fukane.cfg").read_text())
PY
OMP_NUM_THREADS=4 ./build/src/topologicalAnalysis
python3 - <<'PY'
from pathlib import Path
rows = []
for line in Path("output/z2_gap_sweep.dat").read_text().splitlines():
    if line.startswith("#") or not line.strip():
        continue
    b, mu, z2, gap = line.split()
    rows.append((float(b), float(mu), int(z2), float(gap)))
gaps = {round(r[3], 14) for r in rows}
z2s = {r[2] for r in rows}
print("points", len(rows), "unique_z2", sorted(z2s), "unique_gaps", len(gaps))
raise SystemExit(0 if len(gaps) > 1 or len(z2s) > 1 else 1)
PY
```

Expected: The script exits 0 and prints either more than one unique gap or more than one unique Z2 value. If the script exits 1, the sweep is still constant and the worker must trace whether `cfg_local%bdg%mu` is actually used by `compute_z2_fukane_qw_result`; do not commit a constant sweep.

- [ ] **Step 6: Run full suite**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure`
Expected: 66/66 pass.

- [ ] **Step 7: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "fix: evaluate QW Fu-Kane Z2 per sweep point instead of once"
```

---

## Task 5: Fix BdG wire gap to use min |E| instead of min adjacent spacing (F5, P2)

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `src/apps/main_topology.f90:534-544`
- Test: `tests/unit/test_bdg_hamiltonian.pf`

- [ ] **Step 1: Write the failing test**

Add a test that calls a production helper on eigenvalues where the minimum adjacent spacing gives the wrong gap but min |E| gives the correct one. This test must fail before the helper exists or before it is wired to the closest-to-zero definition; do not use a test-local implementation of the formula.

```fortran
@test
subroutine test_bdg_gap_uses_min_abs_energy()
  ! Simulate BdG eigenvalues: [-0.8, -0.3, 0.001, 0.3, 0.8] (meV)
  ! min adjacent spacing = 0.299 (between 0.001 and 0.3)
  ! min |E| = 0.001 (true superconducting gap)
  ! The gap should be 0.001, not 0.299
  real(kind=dp) :: evals(5)
  real(kind=dp) :: gap_min

  evals = [-0.8_dp, -0.3_dp, 0.001_dp, 0.3_dp, 0.8_dp]

  gap_min = bdg_zero_energy_gap(evals)

  @assertEqual(0.001_dp, gap_min, tolerance=1.0e-6_dp, &
    message="BdG gap should be min|E|, not min adjacent spacing")
end subroutine test_bdg_gap_uses_min_abs_energy
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_bdg --output-on-failure`
Expected: Fails to compile because `bdg_zero_energy_gap` is not exported yet, or fails behaviorally if the helper exists with the old adjacent-spacing definition.

- [ ] **Step 3: Add the production gap helper**

In `src/physics/topological_analysis.f90`, add `bdg_zero_energy_gap` to the module's public exports near the other gap helpers:

```fortran
  public :: bdg_zero_energy_gap
```

Then implement it in the module `contains` section:

```fortran
  pure function bdg_zero_energy_gap(eigenvalues) result(gap)
    real(kind=dp), intent(in), contiguous :: eigenvalues(:)
    real(kind=dp) :: gap

    if (size(eigenvalues) < 1) then
      gap = 0.0_dp
    else
      gap = minval(abs(eigenvalues))
    end if
  end function bdg_zero_energy_gap
```

- [ ] **Step 4: Use the helper in the BdG wire gap computation**

In `src/apps/main_topology.f90`, replace lines 533-544. Change from min adjacent spacing to min |E|:

```fortran
      ! Find minimum gap around zero energy
      result%min_gap = bdg_zero_energy_gap(eigvals_bdg)
```

- [ ] **Step 5: Build and run tests**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 -R test_bdg --output-on-failure`
Expected: All BdG tests pass.

- [ ] **Step 6: Run full suite**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure`
Expected: 66/66 pass.

- [ ] **Step 7: Commit**

```bash
git add src/physics/topological_analysis.f90 src/apps/main_topology.f90 tests/unit/test_bdg_hamiltonian.pf
git commit -m "fix: BdG wire gap uses min|E| (zero-energy distance) not min adjacent spacing"
```

---

## Task 6: Final verification and archive

**Files:**
- Archive: `docs/superpowers/specs/2026-05-07-topological-magnetic-bugfixes-design.md`
- Update: `docs/plans/REVIEW.md`
- Update: `docs/plans/BACKLOG.md`

- [ ] **Step 1: Run full test suite from clean build**

Run:
```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure
```
Expected: All tests pass.

- [ ] **Step 2: Update REVIEW.md**

Change group #38 status to reflect the Codex fixes. Mark the 5 open Codex findings as resolved.

- [ ] **Step 3: Update BACKLOG.md**

Add a new phase entry documenting the topological bug fixes.

- [ ] **Step 4: Archive the design spec**

Move the completed design spec out of active specs and into the archive:

```bash
mkdir -p docs/superpowers/specs/archive
git mv docs/superpowers/specs/2026-05-07-topological-magnetic-bugfixes-design.md \
  docs/superpowers/specs/archive/2026-05-07-topological-magnetic-bugfixes-design.md
```

- [ ] **Step 5: Commit**

```bash
git add docs/plans/REVIEW.md docs/plans/BACKLOG.md docs/superpowers/specs/archive/2026-05-07-topological-magnetic-bugfixes-design.md
git commit -m "docs: archive topological bug fixes, update review status"
```
