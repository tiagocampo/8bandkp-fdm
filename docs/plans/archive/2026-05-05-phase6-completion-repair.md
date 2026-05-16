# Phase 6 Completion Repair Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete and repair Phase 6 topological analysis so QW BdG, Fu-Kane Z2, Berry/Kubo conductance, spectral/LDOS, gap sweep, parser/app dispatch, docs, and tests are working and merge-ready.

**Architecture:** Keep existing module boundaries initially, but add small focused helpers inside the touched modules. Make every exposed executable mode compute real results or fail with a clear validation error. Use TDD at each checkpoint: write a failing unit/regression test for the reviewed defect, implement the smallest correct fix, then build and commit.

**Tech Stack:** Fortran 2018, CMake/Ninja, MKL LAPACK/PARDISO, pFUnit, shell/Python regression tests, existing `simulation_config`/`topological_result` types.

---

## File Structure

Modify:

- `src/physics/topological_analysis.f90`: band-major helpers, Majorana profile fix, Fu-Kane parity rewrite, Berry validation/performance, Kubo conductance, gap sweep evaluators.
- `src/physics/bdg_hamiltonian.f90`: import fix, QW BdG electron/hole assembly, Zeeman ownership cleanup, comments/stale variable cleanup.
- `src/physics/green_functions.f90`: spectral validation, LAPACK `info` checks, LDOS shifted matrix fix, bulk/wire spectral support, Landauer helper.
- `src/apps/main_topology.f90`: QW Fu-Kane dispatch, QW BdG dispatch, spectral/conductance/sweep dispatch, output fields, validation.
- `src/io/input_parser.f90`: name-aware Phase 6 topology parsing and config validation.
- `src/core/defs.f90`: result/config fields if status/error/method fields are missing or need clearer names.
- `docs/reference/input-reference.md`: document all Phase 6 modes, fields, outputs, units, and constraints.
- `tests/unit/CMakeLists.txt`: register any new pFUnit test files.
- `tests/unit/test_bdg_hamiltonian.pf`: QW BdG nonzero-k particle-hole, Zeeman, Majorana profile tests.
- `tests/unit/test_z2_invariant.pf`: Fu-Kane parity and product tests.
- `tests/unit/test_chern_number.pf`: Berry validation and Kubo integration tests.
- `tests/unit/test_green_functions.pf`: spectral/LDOS/Landauer tests.
- `tests/unit/test_phase_diagram.pf`: shared sweep/evaluator tests.
- `tests/unit/test_topology_parser.pf`: old/new parser field tests.
- `tests/regression/configs/*.cfg`: new Phase 6 example/regression configs.

Do not create new production modules during this plan. Keep helper routines inside the existing modules listed above. Use the existing regression harness for executable modes instead of adding `tests/unit/test_topology_app_modes.pf`.

## Common Commands

- Configure with tests:

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
```

- Build:

```bash
cmake --build build
```

- Run focused unit tests:

```bash
ctest --test-dir build -L unit --output-on-failure
```

- Run all tests:

```bash
OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure
```

If stale `.mod` errors appear, run:

```bash
rm -f *.mod
cmake --build build
```

---

### Task 1: Foundation Fixes And Guard Tests

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `src/physics/bdg_hamiltonian.f90`
- Modify: `src/physics/green_functions.f90`
- Modify: `src/apps/main_topology.f90`
- Test: `tests/unit/test_bdg_hamiltonian.pf`
- Test: `tests/unit/test_chern_number.pf`
- Test: `tests/unit/test_green_functions.pf`

- [ ] **Step 1: Write failing band-major Majorana profile test**

Add to `tests/unit/test_bdg_hamiltonian.pf` a test that constructs a band-major BdG vector with all density at site 2. Use `half_n = 8 * N`, electron rows `(ib-1)*N + isite`, and hole rows `half_n + (ib-1)*N + isite`.

```fortran
  @test
  subroutine test_majorana_profile_uses_band_major_order()
    type(spatial_grid) :: grid
    complex(kind=dp), allocatable :: psi(:)
    real(kind=dp), allocatable :: profile(:)
    real(kind=dp) :: xi
    integer :: N, half_n, ib, site

    N = 3
    half_n = 8 * N
    allocate(psi(2 * half_n), profile(N))
    psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
    site = 2
    do ib = 1, 8
      psi((ib - 1) * N + site) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      psi(half_n + (ib - 1) * N + site) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    end do

    grid%nx = 1
    grid%ny = N
    allocate(grid%z(N))
    grid%z = [0.0_dp, 1.0_dp, 2.0_dp]

    xi = compute_majorana_profile(psi, grid, 1.0e-8_dp, half_n, profile)

    @assertTrue(profile(2) > 0.99_dp, message="density peak must be at band-major site 2")
    @assertTrue(profile(1) < 1.0e-12_dp, message="site 1 should have no density")
    @assertTrue(profile(3) < 1.0e-12_dp, message="site 3 should have no density")

    deallocate(psi, profile, grid%z)
  end subroutine test_majorana_profile_uses_band_major_order
```

- [ ] **Step 2: Write failing Berry degenerate-grid test**

Add to `tests/unit/test_chern_number.pf`:

```fortran
  @test
  subroutine test_berry_curvature_rejects_degenerate_grid()
    complex(kind=dp), allocatable :: evecs(:,:,:,:)
    real(kind=dp), allocatable :: kx(:), ky(:), Omega(:,:)

    allocate(evecs(2, 1, 1, 2))
    allocate(kx(1), ky(2))
    evecs = cmplx(0.0_dp, 0.0_dp, kind=dp)
    evecs(1,1,:,:) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    kx = [0.0_dp]
    ky = [0.0_dp, 1.0_dp]

    Omega = compute_berry_curvature_lattice(evecs, kx, ky, 1)
    @assertEqual(0, size(Omega, 1), message="degenerate kx grid must return empty curvature")

    deallocate(evecs, kx, ky)
    if (allocated(Omega)) deallocate(Omega)
  end subroutine test_berry_curvature_rejects_degenerate_grid
```

- [ ] **Step 3: Write failing spectral eta validation test**

Add to `tests/unit/test_green_functions.pf` a small QW spectral call with `eta=0.0_dp`. If this codebase prefers `stop 1` for invalid input, put this validation in parser/app tests instead. Preferred module behavior is an empty output on invalid eta:

```fortran
  @test
  subroutine test_spectral_function_rejects_zero_eta()
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:), k_arr(:), E_arr(:), A(:,:)

    cfg%confinement = 1
    cfg%fdStep = 3
    allocate(profile(3,3), kpterms(3,3,10), k_arr(1), E_arr(3))
    profile = 0.0_dp
    kpterms = 0.0_dp
    k_arr = [0.0_dp]
    E_arr = [-1.0_dp, 0.0_dp, 1.0_dp]

    call compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, 0.0_dp, A)
    @assertTrue(allocated(A), message="invalid eta returns allocated empty output")
    @assertEqual(0, size(A, 1), message="eta=0 must not compute NaN spectrum")

    deallocate(profile, kpterms, k_arr, E_arr)
    if (allocated(A)) deallocate(A)
  end subroutine test_spectral_function_rejects_zero_eta
```

- [ ] **Step 4: Run focused tests and verify failures**

Run:

```bash
cmake --build build
ctest --test-dir build -R "bdg|chern|green" --output-on-failure
```

Expected: at least the new tests fail because the current code still uses the reviewed indexing and guard behavior.

- [ ] **Step 5: Implement band-major helper and Majorana profile fix**

In `src/physics/topological_analysis.f90`, add inside `contains` before consumers:

```fortran
  elemental integer function band_major_row(iband, isite, nsite) result(irow)
    implicit none
    integer, intent(in) :: iband, isite, nsite
    irow = (iband - 1) * nsite + isite
  end function band_major_row
```

Update `compute_majorana_profile` density loop:

```fortran
      do ib = 1, 8
        rho(i) = rho(i) + abs(evec_bdg(band_major_row(ib, i, nspatial)))**2
        rho(i) = rho(i) + abs(evec_bdg(half_n + band_major_row(ib, i, nspatial)))**2
      end do
```

When input is invalid or fit fails, return `xi = -1.0_dp` instead of `0.0_dp`.

- [ ] **Step 6: Restore imports and remove stale declarations**

In `src/physics/bdg_hamiltonian.f90`, change:

```fortran
  use magnetic_field, only: add_peierls_coo
```

to:

```fortran
  use magnetic_field, only: add_peierls_coo, compute_zeeman_vz
```

Remove `nnz_zeeman` from the declaration list. Fix pairing comments to bands 1-8.

In `src/physics/topological_analysis.f90`, remove unused `zheevd` from the `use linalg` line and replace touched `acos(-1.0_dp)` with `pi_dp`.

In `src/apps/main_topology.f90`, remove unused `confinementInitialization` if the build confirms it is unused.

- [ ] **Step 7: Implement Berry grid guard and reuse allocation**

At the start of `compute_berry_curvature_lattice`, before `dA`:

```fortran
    if (nkx < 2 .or. nky < 2 .or. n_occ < 1 .or. size(evecs_k, 2) < n_occ) then
      allocate(Omega(0, 0))
      return
    end if
```

Move `complex(kind=dp), allocatable :: overlap(:,:)` out of the inner `block`, allocate once before the loops, and deallocate after the loops.

- [ ] **Step 8: Implement spectral eta guard and LAPACK info check**

In `compute_spectral_function_qw`, before allocating the normal result:

```fortran
    if (eta <= 0.0_dp) then
      allocate(A_kE(0, 0))
      return
    end if
```

After each `zheev` call:

```fortran
      if (info /= 0) then
        A_kE = 0.0_dp
        return
      end if
```

Rename `lorntz` to `lorentz`.

- [ ] **Step 9: Run focused tests**

Run:

```bash
cmake --build build
ctest --test-dir build -R "bdg|chern|green" --output-on-failure
```

Expected: new guard/indexing tests pass.

- [ ] **Step 10: Commit**

```bash
git add src/physics/topological_analysis.f90 src/physics/bdg_hamiltonian.f90 src/physics/green_functions.f90 src/apps/main_topology.f90 tests/unit/test_bdg_hamiltonian.pf tests/unit/test_chern_number.pf tests/unit/test_green_functions.pf
git commit -m "fix(topo): repair phase 6 foundation defects"
```

---

### Task 2: Fu-Kane QW Parity Rewrite

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `src/apps/main_topology.f90`
- Test: `tests/unit/test_z2_invariant.pf`
- Test: `tests/regression/configs/topology_qw_fukane.cfg`

- [ ] **Step 1: Write parity operator unit tests**

In `tests/unit/test_z2_invariant.pf`, add tests for a synthetic band-major vector. Expected helper names are introduced in this task and should be `public` only if tests need direct access.

```fortran
  @test
  subroutine test_qw_parity_expectation_band_major_cb_odd()
    complex(kind=dp) :: psi(16)
    real(kind=dp) :: parity
    integer :: N

    N = 2
    psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
    psi((7 - 1) * N + 1) = cmplx(1.0_dp / sqrt(2.0_dp), 0.0_dp, kind=dp)
    psi((7 - 1) * N + 2) = cmplx(1.0_dp / sqrt(2.0_dp), 0.0_dp, kind=dp)

    parity = qw_inversion_expectation(psi, N)
    @assertEqual(-1.0_dp, parity, tolerance=1.0e-12_dp)
  end subroutine test_qw_parity_expectation_band_major_cb_odd

  @test
  subroutine test_qw_parity_expectation_envelope_odd()
    complex(kind=dp) :: psi(16)
    real(kind=dp) :: parity
    integer :: N

    N = 2
    psi = cmplx(0.0_dp, 0.0_dp, kind=dp)
    psi((1 - 1) * N + 1) = cmplx(1.0_dp / sqrt(2.0_dp), 0.0_dp, kind=dp)
    psi((1 - 1) * N + 2) = cmplx(-1.0_dp / sqrt(2.0_dp), 0.0_dp, kind=dp)

    parity = qw_inversion_expectation(psi, N)
    @assertEqual(-1.0_dp, parity, tolerance=1.0e-12_dp)
  end subroutine test_qw_parity_expectation_envelope_odd
```

- [ ] **Step 2: Write Fu-Kane Gamma inclusion product test**

Add a pure helper test:

```fortran
  @test
  subroutine test_fukane_product_includes_gamma()
    real(kind=dp) :: delta_trim(4)
    integer :: z2

    delta_trim = [-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    z2 = z2_from_trim_parities(delta_trim)
    @assertEqual(1, z2, message="Gamma parity inversion must affect Z2")
  end subroutine test_fukane_product_includes_gamma
```

- [ ] **Step 3: Run tests and verify failures**

Run:

```bash
cmake --build build
ctest --test-dir build -R z2 --output-on-failure
```

Expected: compile failure until helpers are implemented/exported for tests.

- [ ] **Step 4: Implement parity helpers**

In `src/physics/topological_analysis.f90`, add public testable helpers if the module test imports them:

```fortran
  public :: qw_inversion_expectation
  public :: z2_from_trim_parities
```

Add implementations:

```fortran
  elemental real(kind=dp) function band_inversion_parity(iband) result(p)
    implicit none
    integer, intent(in) :: iband
    if (iband <= 6) then
      p = 1.0_dp
    else
      p = -1.0_dp
    end if
  end function band_inversion_parity

  function qw_inversion_expectation(psi, nsite) result(parity)
    implicit none
    complex(kind=dp), intent(in) :: psi(:)
    integer, intent(in) :: nsite
    real(kind=dp) :: parity
    integer :: ib, isite, imirror

    parity = 0.0_dp
    do ib = 1, 8
      do isite = 1, nsite
        imirror = nsite + 1 - isite
        parity = parity + band_inversion_parity(ib) * real( &
          conjg(psi(band_major_row(ib, isite, nsite))) * &
          psi(band_major_row(ib, imirror, nsite)), kind=dp)
      end do
    end do
  end function qw_inversion_expectation

  integer function z2_from_trim_parities(delta_trim) result(z2)
    implicit none
    real(kind=dp), intent(in) :: delta_trim(4)
    real(kind=dp) :: product

    product = product(delta_trim)
    if (product < 0.0_dp) then
      z2 = 1
    else
      z2 = 0
    end if
  end function z2_from_trim_parities
```

Use a different local variable name than `product` if it shadows the intrinsic in your compiler:

```fortran
    trim_product = delta_trim(1) * delta_trim(2) * delta_trim(3) * delta_trim(4)
```

- [ ] **Step 5: Implement symmetry and status-aware Fu-Kane routine**

Add a status-aware routine:

```fortran
  subroutine compute_z2_fukane_qw_result(cfg, profile, kpterms, n_occ, z2, min_gap, status)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:)
    real(kind=dp), contiguous, intent(in) :: kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer, intent(out) :: z2, status
    real(kind=dp), intent(out) :: min_gap
```

Use status values:

```fortran
    integer, parameter :: topo_status_ok = 0
    integer, parameter :: topo_status_invalid = 1
    integer, parameter :: topo_status_asymmetric = 2
    integer, parameter :: topo_status_lapack = 3
```

Inside the routine:

- Validate `n_occ >= 2`, even `n_occ`, and `n_occ < dim_H`.
- Check `profile(i,:)` mirrors `profile(N+1-i,:)` within `1.0e-8_dp`.
- Use a lattice constant helper that reads `cfg%params(layer)%a0`. For QW, choose the central/well layer `layer = max(1, (cfg%numLayers + 1) / 2)` and fall back to the first layer with `a0 > 0.0_dp`. If no positive lattice constant exists, set `status = topo_status_invalid`.
- At each TRIM, diagonalize with `zheev`, check `info`.
- For occupied Kramers pairs `istate = 1, n_occ, 2`, compute parity expectation for one representative or the pair subspace. For the first implementation use the representative with larger norm stability:

```fortran
        parity_val = qw_inversion_expectation(HT(:, istate), N)
        if (parity_val >= 0.0_dp) then
          delta_trim(i_trim) = delta_trim(i_trim) * 1.0_dp
        else
          delta_trim(i_trim) = delta_trim(i_trim) * -1.0_dp
        end if
```

- Compute `min_gap = min(min_gap, evals(n_occ + 1) - evals(n_occ))`.
- Set `z2 = z2_from_trim_parities(delta_trim)`.

- [ ] **Step 6: Update compatibility wrapper**

Change `compute_z2_fukane_qw` to call the status-aware routine:

```fortran
    call compute_z2_fukane_qw_result(cfg, profile, kpterms, n_occ, z2, min_gap, status)
    if (status /= 0) z2 = 0
```

- [ ] **Step 7: Wire QSHE QW app dispatch**

In `src/apps/main_topology.f90`, replace the current QW branch in `case('qshe')` with a call that sets `topo_result%z2_invariant` and `topo_result%min_gap`. If status is not ok, print a clear error and `stop 1`.

- [ ] **Step 8: Add regression config**

Create `tests/regression/configs/topology_qw_fukane.cfg` based on a small symmetric QW config. Include:

```text
topology: T
mode: qshe
compute_z2: T
```

Use a small `FDstep` to keep runtime acceptable.

- [ ] **Step 9: Run tests**

Run:

```bash
cmake --build build
ctest --test-dir build -R z2 --output-on-failure
```

Expected: Fu-Kane helper tests pass.

- [ ] **Step 10: Commit**

```bash
git add src/physics/topological_analysis.f90 src/apps/main_topology.f90 tests/unit/test_z2_invariant.pf tests/regression/configs/topology_qw_fukane.cfg
git commit -m "feat(topo): implement QW Fu-Kane parity invariant"
```

---

### Task 3: Dense QW BdG Assembly And App Path

**Files:**
- Modify: `src/physics/bdg_hamiltonian.f90`
- Modify: `src/apps/main_topology.f90`
- Test: `tests/unit/test_bdg_hamiltonian.pf`
- Create: `tests/regression/configs/topology_qw_bdg.cfg`

- [ ] **Step 1: Write QW BdG nonzero-k particle-hole test**

Add to `tests/unit/test_bdg_hamiltonian.pf`:

```fortran
  @test
  subroutine test_bdg_qw_particle_hole_nonzero_k()
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:), evals(:), rwork(:)
    complex(kind=dp), allocatable :: H(:,:), work(:)
    integer :: dim, info, lwork, i

    call setup_minimal_qw_bdg_fixture(cfg, profile, kpterms)
    call build_bdg_hamiltonian_qw(H, cfg, profile, kpterms, 0.01_dp, 0.0_dp, 0.001_dp)
    dim = size(H, 1)
    allocate(evals(dim), rwork(max(1, 3*dim - 2)), work(1))
    call zheev('N', 'U', dim, H, dim, evals, work, -1, rwork, info)
    lwork = max(1, nint(real(work(1), kind=dp)))
    deallocate(work)
    allocate(work(lwork))
    call zheev('N', 'U', dim, H, dim, evals, work, lwork, rwork, info)

    @assertEqual(0, info, message="QW BdG diagonalization must succeed")
    do i = 1, dim / 2
      @assertEqual(-evals(i), evals(dim + 1 - i), tolerance=1.0e-7_dp)
    end do

    deallocate(profile, kpterms, H, evals, rwork, work)
  end subroutine test_bdg_qw_particle_hole_nonzero_k
```

Add this local helper in `tests/unit/test_bdg_hamiltonian.pf`:

```fortran
  subroutine setup_minimal_qw_bdg_fixture(cfg, profile, kpterms)
    type(simulation_config), intent(out) :: cfg
    real(kind=dp), allocatable, intent(out) :: profile(:,:), kpterms(:,:,:)
    integer :: N, i

    N = 3
    cfg%confinement = 1
    cfg%fdStep = N
    cfg%FDorder = 2
    cfg%numLayers = 1
    cfg%numcb = 2
    cfg%numvb = 6
    cfg%bdg%enabled = .true.
    cfg%bdg%delta_0 = 0.001_dp
    cfg%bdg%mu = 0.0_dp
    cfg%bdg%g_factor = 2.0_dp
    cfg%bdg%B_vec = [0.0_dp, 0.0_dp, 0.0_dp]
    allocate(profile(N,3), kpterms(N,N,10))
    profile = 0.0_dp
    kpterms = 0.0_dp
    do i = 1, N
      profile(i,1) = -0.1_dp
      profile(i,2) = -0.05_dp
      profile(i,3) = 0.5_dp
      kpterms(i,i,1) = 1.0_dp
    end do
  end subroutine setup_minimal_qw_bdg_fixture
```

- [ ] **Step 2: Write Zeeman no-double-counting test**

Add a test that builds QW BdG once with zero field and once with `cfg%bdg%B_field = [0,0,B]`, then checks the diagonal splitting equals one application of `compute_zeeman_vz`, not twice. Use the same row convention:

```fortran
row = (ib - 1) * cfg%fdStep + 1
```

- [ ] **Step 3: Run tests and verify failure**

Run:

```bash
cmake --build build
ctest --test-dir build -R bdg --output-on-failure
```

Expected: particle-hole or Zeeman test fails before assembly cleanup.

- [ ] **Step 4: Implement explicit electron/hole construction**

In `build_bdg_hamiltonian_qw`:

- Build `H_e` with `wv%kx = +k_par`.
- Build `H_h` with `wv%kx = -k_par`.
- Subtract/add `mu` only in final BdG blocks.
- Remove explicit post-`ZB8bandQW` Zeeman addition unless using a local config conversion path.

Core block assignment:

```fortran
    H_bdg(1:N8, 1:N8) = H_e
    do i = 1, N8
      H_bdg(i, i) = H_bdg(i, i) - cmplx(mu, 0.0_dp, kind=dp)
    end do

    do row = 1, N8
      do col = 1, N8
        H_bdg(N8 + row, N8 + col) = -conjg(H_h(row, col))
      end do
    end do
    do i = 1, N8
      H_bdg(N8 + i, N8 + i) = H_bdg(N8 + i, N8 + i) + cmplx(mu, 0.0_dp, kind=dp)
    end do
```

- [ ] **Step 5: Keep pairing blocks band-major**

Use:

```fortran
row = (ib - 1) * N + i
col = ((9 - ib) - 1) * N + i
```

Confirm `(1,2)` and `(2,1)` are Hermitian conjugates.

- [ ] **Step 6: Add QW BdG app dispatch**

In `case('bdg')`:

```fortran
      if (cfg%confinement == 2) then
        call run_bdg_wire(cfg, topo_result)
      else if (cfg%confinement == 1) then
        call run_bdg_qw(cfg, profile, kpterms, topo_result)
      else
        print *, 'Error: BdG mode requires confinement=1 (QW) or confinement=2 (wire)'
        stop 1
      end if
```

Implement `run_bdg_qw` using `build_bdg_hamiltonian_qw`, dense `zheev`, and Majorana/profile summary consistent with wire output.

- [ ] **Step 7: Add QW BdG regression config**

Create `tests/regression/configs/topology_qw_bdg.cfg` with small symmetric QW, `topology: T`, `mode: bdg`, and `bdg: T`.

- [ ] **Step 8: Run focused tests**

Run:

```bash
cmake --build build
ctest --test-dir build -R bdg --output-on-failure
```

Expected: QW BdG tests pass.

- [ ] **Step 9: Commit**

```bash
git add src/physics/bdg_hamiltonian.f90 src/apps/main_topology.f90 tests/unit/test_bdg_hamiltonian.pf tests/regression/configs/topology_qw_bdg.cfg
git commit -m "feat(topo): complete dense QW BdG path"
```

---

### Task 4: Berry, Kubo, Landauer, And Conductance Mode

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `src/physics/green_functions.f90`
- Modify: `src/apps/main_topology.f90`
- Modify: `src/core/defs.f90`
- Test: `tests/unit/test_chern_number.pf`
- Test: `tests/unit/test_green_functions.pf`
- Create: `tests/regression/configs/topology_conductance_qwz.cfg`

- [ ] **Step 1: Write Kubo integration test**

Add to `tests/unit/test_chern_number.pf`:

```fortran
  @test
  subroutine test_conductance_kubo_integrates_curvature()
    real(kind=dp), allocatable :: Omega(:,:), kx(:), ky(:)
    real(kind=dp) :: sigma

    allocate(Omega(2,2), kx(2), ky(2))
    kx = [0.0_dp, pi_dp]
    ky = [0.0_dp, pi_dp]
    Omega = 2.0_dp / pi_dp

    sigma = compute_conductance_kubo(Omega, kx, ky)
    @assertEqual(1.0_dp, sigma, tolerance=1.0e-12_dp)

    deallocate(Omega, kx, ky)
  end subroutine test_conductance_kubo_integrates_curvature
```

Update the function signature in tests to match the implementation after removing unused `n_occ`.

- [ ] **Step 2: Write Landauer analytic helper test**

In `tests/unit/test_green_functions.pf`, add:

```fortran
  @test
  subroutine test_landauer_single_perfect_channel()
    real(kind=dp) :: T, G

    T = compute_landauer_transmission_1d(E=0.0_dp, onsite=0.0_dp, hopping=-1.0_dp, eta=1.0e-6_dp)
    G = compute_conductance_landauer_from_transmission(T)

    @assertEqual(1.0_dp, T, tolerance=1.0e-5_dp)
    @assertEqual(1.0_dp, G, tolerance=1.0e-5_dp)
  end subroutine test_landauer_single_perfect_channel
```

- [ ] **Step 3: Run tests and verify failures**

Run:

```bash
cmake --build build
ctest --test-dir build -R "chern|green" --output-on-failure
```

Expected: missing/old conductance interfaces fail.

- [ ] **Step 4: Implement direct Kubo integration**

In `topological_analysis.f90`, replace or add:

```fortran
  function compute_conductance_kubo(berry_curvature, kx_arr, ky_arr) result(sigma_xy)
    real(kind=dp), contiguous, intent(in) :: berry_curvature(:,:)
    real(kind=dp), contiguous, intent(in) :: kx_arr(:), ky_arr(:)
    real(kind=dp) :: sigma_xy
    real(kind=dp) :: dkx, dky

    sigma_xy = 0.0_dp
    if (size(kx_arr) < 2 .or. size(ky_arr) < 2) return
    if (size(berry_curvature,1) /= size(kx_arr)) return
    if (size(berry_curvature,2) /= size(ky_arr)) return
    dkx = kx_arr(2) - kx_arr(1)
    dky = ky_arr(2) - ky_arr(1)
    sigma_xy = sum(berry_curvature) * dkx * dky / (2.0_dp * pi_dp)
  end function compute_conductance_kubo
```

Keep one Chern helper:

```fortran
  elemental real(kind=dp) function compute_hall_conductance_from_chern(C) result(sigma_xy)
    integer, intent(in) :: C
    sigma_xy = real(C, kind=dp)
  end function compute_hall_conductance_from_chern
```

- [ ] **Step 5: Implement Landauer effective-chain helpers**

In `green_functions.f90`, add public helpers:

```fortran
  public :: compute_landauer_transmission_1d
  public :: compute_conductance_landauer_from_transmission
```

Implement a single-site perfect-chain helper:

```fortran
  function compute_landauer_transmission_1d(E, onsite, hopping, eta) result(T)
    real(kind=dp), intent(in) :: E, onsite, hopping, eta
    real(kind=dp) :: T
    complex(kind=dp) :: z, root_term, sigma_l, sigma_r, G
    real(kind=dp) :: gamma_l, gamma_r

    if (eta <= 0.0_dp .or. abs(hopping) < 1.0e-30_dp) then
      T = 0.0_dp
      return
    end if
    z = cmplx(E - onsite, eta, kind=dp)
    root_term = sqrt(z*z - cmplx(4.0_dp * hopping*hopping, 0.0_dp, kind=dp))
    sigma_l = 0.5_dp * (z - root_term)
    sigma_r = sigma_l
    gamma_l = -2.0_dp * aimag(sigma_l)
    gamma_r = -2.0_dp * aimag(sigma_r)
    G = 1.0_dp / (z - sigma_l - sigma_r)
    T = max(0.0_dp, gamma_l * gamma_r * abs(G)**2)
  end function compute_landauer_transmission_1d

  elemental real(kind=dp) function compute_conductance_landauer_from_transmission(T) result(G)
    real(kind=dp), intent(in) :: T
    G = T
  end function compute_conductance_landauer_from_transmission
```

Document units in comments.

- [ ] **Step 6: Implement conductance mode prerequisites**

In `run_conductance`, branch on `cfg%topo%conductance_method`:

```fortran
select case(trim(cfg_in%topo%conductance_method))
case('kubo_chern')
  chern = compute_chern_qwz(cfg_in%topo%qwz_u, cfg_in%topo%berry_nk)
  result%chern_number = chern
  result%conductance_xy = compute_hall_conductance_from_chern(chern)
case('kubo_berry')
  ! build QWZ eigenvectors on configured grid, compute Omega, integrate
case('landauer')
  T = compute_landauer_transmission_1d(cfg_in%topo%landauer_energy, 0.0_dp, -1.0_dp, cfg_in%topo%spectral_eta)
  result%conductance_zz = compute_conductance_landauer_from_transmission(T)
case default
  print *, 'Error: unsupported conductance_method: ', trim(cfg_in%topo%conductance_method)
  stop 1
end select
```

If fields such as `berry_nk` or `landauer_energy` do not exist, add conservative defaults in `topology_config`.

- [ ] **Step 7: Add regression config**

Create `tests/regression/configs/topology_conductance_qwz.cfg`:

```text
topology: T
mode: conductance
compute_conductance: T
conductance_method: kubo_chern
qwz_u: -0.8
```

- [ ] **Step 8: Run tests**

Run:

```bash
cmake --build build
ctest --test-dir build -R "chern|green" --output-on-failure
```

Expected: Kubo and Landauer helper tests pass.

- [ ] **Step 9: Commit**

```bash
git add src/physics/topological_analysis.f90 src/physics/green_functions.f90 src/apps/main_topology.f90 src/core/defs.f90 tests/unit/test_chern_number.pf tests/unit/test_green_functions.pf tests/regression/configs/topology_conductance_qwz.cfg
git commit -m "feat(topo): complete conductance calculations"
```

---

### Task 5: Spectral, LDOS, Bulk/Wire Support

**Files:**
- Modify: `src/physics/green_functions.f90`
- Modify: `src/apps/main_topology.f90`
- Test: `tests/unit/test_green_functions.pf`
- Create: `tests/regression/configs/topology_spectral_qw.cfg`
- Create: `tests/regression/configs/topology_spectral_wire.cfg`
- Create: `tests/regression/configs/topology_spectral_bulk.cfg`

- [ ] **Step 1: Write LDOS off-diagonal shift regression test**

In `tests/unit/test_green_functions.pf`, add a test using a 2x2 CSR matrix with off-diagonal hopping. The expected LDOS should match dense inverse of `zI-H`, not `(z-H_ij)` for off-diagonals.

```fortran
  @test
  subroutine test_ldos_csr_shifts_diagonal_only()
    type(csr_matrix) :: H
    real(kind=dp) :: ldos(2)
    integer :: rows(4), cols(4)
    complex(kind=dp) :: vals(4)

    rows = [1, 1, 2, 2]
    cols = [1, 2, 1, 2]
    vals = [cmplx(0.0_dp,0.0_dp,kind=dp), cmplx(-0.5_dp,0.0_dp,kind=dp), &
            cmplx(-0.5_dp,0.0_dp,kind=dp), cmplx(0.0_dp,0.0_dp,kind=dp)]
    call csr_build_from_coo(H, 2, 2, 4, rows, cols, vals)

    call compute_ldos_csr(H, 0.0_dp, 0.1_dp, ldos)
    @assertTrue(all(ldos > 0.0_dp), message="coupled LDOS must be positive")
    @assertEqual(ldos(1), ldos(2), tolerance=1.0e-10_dp)

    call csr_free(H)
  end subroutine test_ldos_csr_shifts_diagonal_only
```

- [ ] **Step 2: Strengthen spectral sum-rule test**

Update existing QW spectral test to integrate:

```fortran
integral = sum(A_kE(1,:)) * (E_arr(2) - E_arr(1))
@assertTrue(abs(integral - real(8 * cfg%fdStep, kind=dp)) < 1.0_dp)
```

Use a broad energy window so the Lorentzian tails are captured.

- [ ] **Step 3: Run tests and verify failures**

Run:

```bash
cmake --build build
ctest --test-dir build -R green --output-on-failure
```

Expected: LDOS test fails before shifted-matrix fix if `USE_ARPACK` enables PARDISO LDOS tests.

- [ ] **Step 4: Fix CSR shifted matrix**

In `compute_ldos_csr`:

```fortran
    a_val = -H%values
    do r = 1, N
      found_diag = .false.
      do i = ia(r), ia(r + 1) - 1
        if (ja(i) == r) then
          a_val(i) = a_val(i) + shift
          found_diag = .true.
        end if
      end do
      if (.not. found_diag) then
        print *, 'ERROR: compute_ldos_csr: missing structural diagonal at row ', r
        stop 1
      end if
    end do
```

Declare `logical :: found_diag`.

- [ ] **Step 5: Add bulk spectral function**

Add `compute_spectral_function_bulk` in `green_functions.f90`, using an 8x8 dense Hamiltonian and the same Lorentzian broadening/`eta` validation pattern. Use the existing bulk Hamiltonian constructor call pattern from `src/apps/main.f90`: allocate an `8x8` complex matrix, set `wavevector` from `k_arr(ik)`, diagonalize with `zheev`, and accumulate the Lorentzian spectrum.

- [ ] **Step 6: Add wire spectral function**

Add `compute_spectral_function_wire` that builds CSR wire Hamiltonian for each `kz`, solves a configured small number of eigenvalues near the requested energy window with the existing `make_eigensolver` factory, and broadens eigenvalues into `A(k,E)`. Reuse the exact wire initialization sequence from `run_bdg_wire`: `confinementInitialization_2d`, `wire_workspace`, `ZB8bandGeneralized`, `auto_compute_energy_window`, `solver%solve`.

- [ ] **Step 7: Wire spectral app mode**

In `run_spectral`, branch:

```fortran
select case(cfg_in%confinement)
case(0)
  call compute_spectral_function_bulk(...)
case(1)
  call compute_spectral_function_qw(...)
case(2)
  call compute_spectral_function_wire(...)
case default
  print *, 'Error: unsupported confinement for spectral mode'
  stop 1
end select
```

Write `output/spectral_function.dat` with dimensions, units, `eta`, and k/E grid values.

- [ ] **Step 8: Add regression configs**

Create:

- `tests/regression/configs/topology_spectral_qw.cfg`
- `tests/regression/configs/topology_spectral_wire.cfg`
- `tests/regression/configs/topology_spectral_bulk.cfg`

Each should include:

```text
topology: T
mode: spectral
compute_spectral: T
spectral_k_grid: -0.01 0.01 3
spectral_E_grid: -0.1 0.1 51 0.002
```

Use these labels and implement them in Task 7: `compute_spectral`, `spectral_k_grid`, and `spectral_E_grid`.

- [ ] **Step 9: Run tests**

Run:

```bash
cmake --build build
ctest --test-dir build -R green --output-on-failure
```

Expected: green-function tests pass.

- [ ] **Step 10: Commit**

```bash
git add src/physics/green_functions.f90 src/apps/main_topology.f90 tests/unit/test_green_functions.pf tests/regression/configs/topology_spectral_qw.cfg tests/regression/configs/topology_spectral_wire.cfg tests/regression/configs/topology_spectral_bulk.cfg
git commit -m "feat(topo): complete spectral and LDOS paths"
```

---

### Task 6: Real Gap Sweep Evaluators

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `src/apps/main_topology.f90`
- Test: `tests/unit/test_phase_diagram.pf`
- Create: `tests/regression/configs/topology_sweep_bhz.cfg`
- Create: `tests/regression/configs/topology_sweep_qw.cfg`

- [ ] **Step 1: Replace heuristic tests with fake evaluator test**

In `tests/unit/test_phase_diagram.pf`, remove tests that assert `B > 0 => Z2=1`. Add a transition detector test:

```fortran
  @test
  subroutine test_transition_detection_scans_b_and_mu()
    integer :: z2_map(2,3)
    real(kind=dp) :: gap_map(2,3)
    real(kind=dp), allocatable :: transitions(:,:)

    z2_map = reshape([0, 0, 1, 1, 1, 0], [2,3])
    gap_map = 0.01_dp

    call detect_z2_transitions(z2_map, gap_map, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.001_dp, transitions)

    @assertTrue(size(transitions, 1) >= 2, message="transitions must be found along B and mu directions")
    deallocate(transitions)
  end subroutine test_transition_detection_scans_b_and_mu
```

- [ ] **Step 2: Add BHZ analytic evaluator test**

```fortran
  @test
  subroutine test_bhz_analytic_gap_sweep_has_nonzero_gaps()
    type(simulation_config) :: cfg
    integer, allocatable :: z2_map(:,:)
    real(kind=dp), allocatable :: gap_map(:,:), transitions(:,:)

    cfg%topo%sweep_model = 'bhz_analytic'
    cfg%topo%bhz_M = 0.0_dp
    call compute_z2_gap_sweep(cfg, -0.01_dp, 0.01_dp, 3, 0.0_dp, 0.0_dp, 1, &
      1.0e-6_dp, z2_map, gap_map, transitions)

    @assertTrue(any(gap_map > 0.0_dp), message="gap sweep must compute non-placeholder gaps")
    @assertTrue(any(z2_map == 0) .and. any(z2_map == 1), message="BHZ sweep must cross phases")

    deallocate(z2_map, gap_map, transitions)
  end subroutine test_bhz_analytic_gap_sweep_has_nonzero_gaps
```

Add `sweep_model` to `topology_config` if it does not exist.

- [ ] **Step 3: Run tests and verify failure**

Run:

```bash
cmake --build build
ctest --test-dir build -R phase --output-on-failure
```

Expected: new tests fail until evaluator-backed sweep is implemented.

- [ ] **Step 4: Implement transition detection helper**

In `topological_analysis.f90`, add public or test-visible:

```fortran
  subroutine detect_z2_transitions(z2_map, gap_map, B_min, B_max, mu_min, mu_max, gap_threshold, transitions)
```

Scan adjacent B and mu cells, append transition midpoint when Z2 changes or either endpoint gap is below threshold.

- [ ] **Step 5: Implement BHZ analytic evaluator**

Add:

```fortran
  subroutine eval_bhz_analytic(B_val, mu_val, cfg, z2, gap, status)
    real(kind=dp), intent(in) :: B_val, mu_val
    type(simulation_config), intent(in) :: cfg
    integer, intent(out) :: z2, status
    real(kind=dp), intent(out) :: gap
    real(kind=dp) :: M_eff

    M_eff = cfg%topo%bhz_M + B_val - mu_val
    gap = abs(2.0_dp * M_eff)
    if (M_eff < 0.0_dp) then
      z2 = 1
    else
      z2 = 0
    end if
    status = 0
  end subroutine eval_bhz_analytic
```

Confirm `bhz_M` units in `defs.f90` and normalize to eV consistently before using it.

- [ ] **Step 6: Implement wire and QW evaluators**

Add evaluator branches inside `compute_z2_gap_sweep`:

```fortran
select case(trim(cfg%topo%sweep_model))
case('bhz_analytic')
  call eval_bhz_analytic(...)
case('wire_bdg')
  call eval_wire_bdg_gap(...)
case('qw_fukane')
  call eval_qw_fukane_gap(...)
case default
  print *, 'ERROR: unknown sweep_model'
  stop 1
end select
```

`eval_wire_bdg_gap` uses existing wire BdG/eigensolver path and finite-wire Majorana classifier. `eval_qw_fukane_gap` calls `compute_z2_fukane_qw_result`.

- [ ] **Step 7: Update app output**

Ensure `run_gap_sweep` writes real `z2_phase_diagram.dat` columns:

```text
# B(T) mu(eV) z2 gap(eV)
```

Write transition records either to the same file header or `output/z2_transitions.dat`.

- [ ] **Step 8: Add regression configs**

Create:

- `tests/regression/configs/topology_sweep_bhz.cfg`
- `tests/regression/configs/topology_sweep_qw.cfg`

Both must use tiny grids such as `nB=3`, `nMu=2` for runtime.

- [ ] **Step 9: Run tests**

Run:

```bash
cmake --build build
ctest --test-dir build -R phase --output-on-failure
```

Expected: phase-diagram unit tests pass and no test asserts the old `B > 0` heuristic.

- [ ] **Step 10: Commit**

```bash
git add src/physics/topological_analysis.f90 src/apps/main_topology.f90 src/core/defs.f90 tests/unit/test_phase_diagram.pf tests/regression/configs/topology_sweep_bhz.cfg tests/regression/configs/topology_sweep_qw.cfg
git commit -m "feat(topo): replace gap sweep placeholders"
```

---

### Task 7: Parser, Output, Docs, Regression Wiring

**Files:**
- Modify: `src/io/input_parser.f90`
- Modify: `src/apps/main_topology.f90`
- Modify: `src/core/defs.f90`
- Modify: `docs/reference/input-reference.md`
- Modify: `tests/unit/test_topology_parser.pf`
- Modify: `tests/regression/configs/*.cfg`

- [ ] **Step 1: Write parser tests for old and new topology blocks**

In `tests/unit/test_topology_parser.pf`, add tests for:

```text
topology: T
mode: conductance
compute_conductance: T
conductance_method: kubo_chern
```

and old format:

```text
topology: T
mode: qhe
compute_chern: T
```

Assert parsed fields match defaults and requested methods.

- [ ] **Step 2: Run parser tests and verify failure**

Run:

```bash
cmake --build build
ctest --test-dir build -R topology_parser --output-on-failure
```

Expected: new method/range fields fail until parser is updated.

- [ ] **Step 3: Implement name-aware Phase 6 parsing**

In `input_parser.f90`, replace Phase 6 optional peek/backspace blocks with a loop that:

- Reads the next non-comment line.
- Splits label before `:`.
- Dispatches known labels:
  - `compute_gap_sweep`
  - `gap_sweep_B`
  - `gap_sweep_mu`
  - `gap_threshold`
  - `sweep_model`
  - `compute_conductance`
  - `conductance_method`
  - `compute_spectral`
  - `spectral_k_grid`
  - `spectral_E_grid`
- Backspaces only when the line belongs to a different top-level block.

Validation rules:

```fortran
if (cfg%topo%spectral_eta <= 0.0_dp) stop 1
if (cfg%topo%gap_sweep_nB < 1 .or. cfg%topo%gap_sweep_nMu < 1) stop 1
if (.not. is_supported_conductance_method(cfg%topo%conductance_method)) stop 1
```

- [ ] **Step 4: Update app output**

In `main_topology.f90`, ensure `topology_result.dat` writes:

```fortran
write(iounit, '(A,F12.6)') '# Conductance xy (e^2/h): ', topo_result%conductance_xy
write(iounit, '(A,F12.6)') '# Conductance zz: ', topo_result%conductance_zz
write(iounit, '(A,I0)') '# Majorana count: ', topo_result%n_majorana
write(iounit, '(A,I0)') '# Majorana fit failures: ', topo_result%n_majorana_fit_failed
```

Ensure `topological_result` contains these fields:

```fortran
integer :: n_majorana = 0
integer :: n_majorana_fit_failed = 0
```

This task adds no new allocatable result fields. Do not change `topological_result_finalize` in this step unless a compile error shows an existing allocatable field was renamed.

- [ ] **Step 5: Update input reference**

In `docs/reference/input-reference.md`, document:

- modes: `qhe`, `qshe`, `bdg`, `spectral`, `conductance`, `sweep`
- conductance methods: `kubo_chern`, `kubo_berry`, `landauer`, `longitudinal`
- spectral grids and `eta > 0`
- sweep grids, `sweep_model`, `gap_threshold`
- Fu-Kane QW symmetry requirement
- output files and units

- [ ] **Step 6: Add regression configs to test registry**

If the regression harness requires explicit registration, add new configs and expected behavior in the relevant CMake/test scripts. Ensure each new config is small enough for CI.

- [ ] **Step 7: Run parser and app build**

Run:

```bash
cmake --build build
ctest --test-dir build -R topology_parser --output-on-failure
```

Expected: parser tests pass.

- [ ] **Step 8: Run broader topology tests**

Run:

```bash
ctest --test-dir build -R "topology|bdg|chern|green|phase" --output-on-failure
```

Expected: focused Phase 6 tests pass.

- [ ] **Step 9: Commit**

```bash
git add src/io/input_parser.f90 src/apps/main_topology.f90 src/core/defs.f90 docs/reference/input-reference.md tests/unit/test_topology_parser.pf tests/regression/configs
git commit -m "feat(topo): wire phase 6 parser outputs and docs"
```

---

### Task 8: Full Verification And Cleanup

**Files:**
- Modify only files needed by failed verification.

- [ ] **Step 1: Clean stale root modules**

Run:

```bash
rm -f *.mod
```

- [ ] **Step 2: Fresh configure**

Run:

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
```

Expected: configure succeeds.

- [ ] **Step 3: Full build**

Run:

```bash
cmake --build build
```

Expected: build succeeds with all executables, including `build/src/topologicalAnalysis`.

- [ ] **Step 4: Full test suite**

Run:

```bash
OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure
```

Expected: all tests pass.

- [ ] **Step 5: Manual smoke test topologicalAnalysis modes**

For each new config, write it to `input.cfg` without using `cp -i`:

```bash
install -m 0644 tests/regression/configs/topology_conductance_qwz.cfg input.cfg
./build/src/topologicalAnalysis
```

Repeat for:

- `topology_qw_fukane.cfg`
- `topology_qw_bdg.cfg`
- `topology_spectral_qw.cfg`
- `topology_spectral_wire.cfg`
- `topology_spectral_bulk.cfg`
- `topology_sweep_bhz.cfg`

Expected: each mode exits successfully and writes its documented output files.

- [ ] **Step 6: Scan for placeholders and stale claims**

Run:

```bash
rg -n "placeholder|B > 0|Z2 = 0|not implemented|TODO|TBD|fallback" src/physics src/apps src/io docs/reference tests/unit tests/regression/configs
```

Expected: no Phase 6 executable path claims placeholder behavior. Legitimate historical docs must be reviewed and either updated or clearly labeled as archived.

- [ ] **Step 7: Commit verification fixes**

```bash
git add src/physics/topological_analysis.f90 src/physics/bdg_hamiltonian.f90 src/physics/green_functions.f90 src/apps/main_topology.f90 src/io/input_parser.f90 src/core/defs.f90 docs/reference/input-reference.md tests/unit tests/regression/configs
git commit -m "fix(topo): resolve phase 6 verification issues"
```

- [ ] **Step 8: Final status**

Run:

```bash
git status --short
```

Expected: only unrelated pre-existing user changes remain. Do not revert unrelated files.
