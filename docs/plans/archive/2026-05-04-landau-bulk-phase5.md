# Bulk Landau Levels (Confinement=3) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement confinement=3 mode for computing orbital Landau levels in bulk semiconductors with arbitrary B-field direction.

**Architecture:** New x-discretized confinement mode reuses QW machinery (kpterms, zheevx dense diagonalization, OpenMP k-sweep). A dedicated `ZB8bandLandau` assembly routine builds 8Nx8N Hamiltonians with position-dependent gauge-shifted wavevectors. Landau gauge A=(0, Bz·x, −By·x) captures orbital effects from B_y, B_z; B_x contributes Zeeman only.

**Tech Stack:** Fortran 2018, MKL LAPACK (zheevx), OpenMP, pFUnit tests, Python regression verification.

---

## File Structure

| File | Responsibility |
|------|---------------|
| `src/core/defs.f90` | Add `landau_nx`, `landau_width`, `landau_sweep` fields; update validator to allow confinement=3 |
| `src/physics/magnetic_field.f90` | New `compute_gauge_shifts` subroutine for Π_y(x), Π_z(x) |
| `src/physics/confinement_init.f90` | New `confinementInitialization_landau` for single-material kpterms setup |
| `src/physics/hamiltonianConstructor.f90` | New `ZB8bandLandau` assembly with position-dependent k-values |
| `src/io/input_parser.f90` | Parse confinement=3, landau_*, decouple b_field from bdg%enabled |
| `src/apps/main.f90` | confinement==3 dispatch branch |
| `tests/unit/test_landau.pf` | Unit tests for gauge shifts, hermiticity, bulk recovery |
| `tests/regression/configs/landau_InAs.cfg` | InAs at B=5T regression config |
| `tests/regression/configs/landau_GaAs.cfg` | GaAs at B=5T regression config |

---

### Task 1: Add Landau Fields to `simulation_config` and Update Validator

**Files:**
- Modify: `src/core/defs.f90:370-397` (simulation_config type)
- Modify: `src/core/defs.f90:512-515` (validator)
- Modify: `src/core/defs.f90:585-631` (init_grid_from_config)

- [ ] **Step 1: Add Landau fields to `simulation_config` type**

Add after the wire-specific fields block (after line 379, before the FEAST fields):

```fortran
    ! ---- Landau level fields (confinement=3) ----
    integer            :: landau_nx = 100          ! grid points in x-direction
    real(kind=dp)      :: landau_width = 2000.0_dp ! domain width in Angstrom
    character(len=4)   :: landau_sweep = 'ky'      ! sweep mode: ky, kz, or B
```

- [ ] **Step 2: Update the validator to allow confinement=3**

Change `src/core/defs.f90:513` from:

```fortran
      if (cfg%confinement < 0 .or. cfg%confinement > 2) then
        error stop 'validate_simulation_config: confinement must be 0, 1, or 2'
```

to:

```fortran
      if (cfg%confinement < 0 .or. cfg%confinement > 3) then
        error stop 'validate_simulation_config: confinement must be 0, 1, 2, or 3'
      end if
      if (cfg%confinement == 3) then
        if (cfg%landau_nx < 3) then
          error stop 'validate_simulation_config: landau_nx must be >= 3'
        end if
        if (cfg%landau_width <= 0.0_dp) then
          error stop 'validate_simulation_config: landau_width must be > 0'
        end if
```

- [ ] **Step 3: Add confinement=3 case to `init_grid_from_config`**

Add after the `case (2)` block (after line 630 in `init_grid_from_config`), before `end select`:

```fortran
    case (3)
      ! Landau: 1D x-discretization for orbital Landau quantization
      cfg%grid%ndim = 1
      cfg%grid%nx   = cfg%landau_nx
      cfg%grid%ny   = 1
      cfg%grid%dx   = cfg%landau_width / real(cfg%landau_nx - 1, kind=dp)
      cfg%grid%dy   = 0.0_dp

      ! Build x-coordinate array
      if (.not. allocated(cfg%grid%x)) then
        allocate(cfg%grid%x(cfg%landau_nx))
        do i = 1, cfg%landau_nx
          cfg%grid%x(i) = (i - 1) * cfg%grid%dx
        end do
      end if
```

- [ ] **Step 4: Build and verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build succeeds (new fields have defaults, no callers yet)

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90
git commit -m "feat: add landau_nx, landau_width, landau_sweep fields to simulation_config"
```

---

### Task 2: Implement `compute_gauge_shifts` in `magnetic_field.f90` (TDD)

**Files:**
- Modify: `src/physics/magnetic_field.f90` (add subroutine, update public)
- Create: `tests/unit/test_landau.pf`
- Modify: `tests/CMakeLists.txt` (register test)

- [ ] **Step 1: Write the failing unit test for gauge shifts**

Create `tests/unit/test_landau.pf`:

```fortran
module test_landau
  use funit
  use definitions
  use magnetic_field
  implicit none

contains

  @test
  subroutine test_gauge_shifts_zero_B()
    ! At B=0, gauge-shifted k-values should equal the input k-values
    real(kind=dp) :: x_grid(5)
    real(kind=dp) :: Pi_y(5), Pi_z(5)
    integer, parameter :: N = 5
    integer :: i

    x_grid = [0.0_dp, 10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]

    call compute_gauge_shifts(x_grid, [0.0_dp, 0.0_dp, 0.0_dp], 0.05_dp, 0.03_dp, Pi_y, Pi_z)

    do i = 1, N
      @assertEqual(0.05_dp, Pi_y(i), tolerance=1e-14_dp)
      @assertEqual(0.03_dp, Pi_z(i), tolerance=1e-14_dp)
    end do
  end subroutine

  @test
  subroutine test_gauge_shifts_nonzero_B()
    ! Verify Pi_y(i) = ky + e*Bz*x_i_m / hbar_J (converted to 1/AA)
    ! Gauge: A = (0, Bz*x, -By*x) => Pi_y = ky + Bz*x*1e-20/hbar
    real(kind=dp) :: x_grid(3)
    real(kind=dp) :: Pi_y(3), Pi_z(3)
    real(kind=dp) :: ky, kz, Bz, By
    real(kind=dp) :: expected_shift_y, expected_shift_z

    x_grid = [0.0_dp, 100.0_dp, 200.0_dp]  ! in Angstrom
    ky = 0.0_dp
    kz = 0.0_dp
    Bz = 5.0_dp   ! 5 Tesla
    By = 0.0_dp

    call compute_gauge_shifts(x_grid, [0.0_dp, By, Bz], ky, kz, Pi_y, Pi_z)

    ! At x=0: Pi_y = ky = 0
    @assertEqual(0.0_dp, Pi_y(1), tolerance=1e-14_dp)

    ! At x=100 AA: Pi_y = 0 + Bz * x * 1e-20 / hbar
    expected_shift_y = Bz * 100.0_dp * 1.0e-20_dp / hbar
    @assertEqual(expected_shift_y, Pi_y(2), tolerance=1e-14_dp)
    ! Expected: 5 * 100 * 1e-20 / 6.582e-16 = 7.596e-3 1/AA

    ! At x=200 AA: double the shift
    @assertEqual(2.0_dp * expected_shift_y, Pi_y(3), tolerance=1e-14_dp)
  end subroutine

  @test
  subroutine test_gauge_shifts_By_component()
    ! Verify Pi_z(i) = kz - By*x*1e-20/hbar
    real(kind=dp) :: x_grid(2)
    real(kind=dp) :: Pi_y(2), Pi_z(2)
    real(kind=dp) :: expected_shift_z

    x_grid = [0.0_dp, 100.0_dp]
    call compute_gauge_shifts(x_grid, [0.0_dp, 3.0_dp, 0.0_dp], 0.0_dp, 0.0_dp, Pi_y, Pi_z)

    ! Pi_y should be zero (Bz=0)
    @assertEqual(0.0_dp, Pi_y(2), tolerance=1e-14_dp)

    ! Pi_z = 0 - By * x * 1e-20 / hbar
    expected_shift_z = -3.0_dp * 100.0_dp * 1.0e-20_dp / hbar
    @assertEqual(expected_shift_z, Pi_z(2), tolerance=1e-14_dp)
  end subroutine

end module
```

- [ ] **Step 2: Register the test in CMake**

Add after the `test_green_functions` block (around line 144) in `tests/CMakeLists.txt`:

```cmake
    add_pfunit_ctest(test_landau
        TEST_SOURCES unit/test_landau.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )
```

- [ ] **Step 3: Build and verify test fails (subroutine not found)**

Run: `cmake --build build 2>&1 | grep -i error`
Expected: Compilation error — `compute_gauge_shifts` not found

- [ ] **Step 4: Implement `compute_gauge_shifts`**

Add to `src/physics/magnetic_field.f90` before `end module`, and add `compute_gauge_shifts` to the `public` list (line 8):

Change line 8:
```fortran
  public :: add_zeeman_coo, add_peierls_coo, compute_zeeman_vz, compute_gauge_shifts
```

Add the subroutine before `end module`:

```fortran
  subroutine compute_gauge_shifts(x_grid, B_vec, ky, kz, Pi_y, Pi_z)
    ! Computes gauge-shifted k-values for Landau mode.
    ! Landau gauge: A = (0, Bz*x, -By*x)
    !   Pi_y(i) = ky + Bz * x(i) * 1e-20 / hbar   [1/AA]
    !   Pi_z(i) = kz - By * x(i) * 1e-20 / hbar   [1/AA]
    !
    ! Unit derivation: Pi = k + eA/hbar
    !   e*B*x_SI / hbar_J = e*B*x*1e-10 / (hbar*e) = B*x*1e-10/hbar  [1/m]
    !   Convert to 1/AA: * 1e-10 => B*x*1e-20/hbar  [1/AA]
    real(kind=dp), intent(in) :: x_grid(:)
    real(kind=dp), intent(in) :: B_vec(3), ky, kz
    real(kind=dp), intent(out) :: Pi_y(:), Pi_z(:)

    real(kind=dp) :: By, Bz, inv_hbar
    integer :: i, N

    N = size(x_grid)
    By = B_vec(2)
    Bz = B_vec(3)
    inv_hbar = 1.0_dp / hbar

    do i = 1, N
      Pi_y(i) = ky + Bz * x_grid(i) * 1.0e-20_dp * inv_hbar
      Pi_z(i) = kz - By * x_grid(i) * 1.0e-20_dp * inv_hbar
    end do
  end subroutine compute_gauge_shifts
```

- [ ] **Step 5: Build and run the test**

Run: `cmake --build build && ctest --test-dir build -R test_landau --output-on-failure`
Expected: 3 tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/physics/magnetic_field.f90 tests/unit/test_landau.pf tests/CMakeLists.txt
git commit -m "feat: add compute_gauge_shifts for Landau gauge with unit tests"
```

---

### Task 3: Implement `confinementInitialization_landau`

**Files:**
- Modify: `src/physics/confinement_init.f90` (add subroutine, update public)

- [ ] **Step 1: Add `confinementInitialization_landau` to public exports**

Add to the `public` declarations after line 18:

```fortran
  public :: confinementInitialization_landau
```

- [ ] **Step 2: Implement the subroutine**

Add before `end module` in `confinement_init.f90`. This sets up kpterms for a single homogeneous material, reusing `build_kpterm_block`:

```fortran
  subroutine confinementInitialization_landau(x_grid, material, params, &
      & profile, kpterms, FDorder)
    real(kind=dp), intent(in), contiguous :: x_grid(:)
    character(len=255), intent(in) :: material(1)
    type(paramStruct), intent(in) :: params(1)
    real(kind=dp), intent(inout), allocatable :: profile(:,:)
    real(kind=dp), intent(inout), contiguous :: kpterms(:,:,:)
    integer, intent(in), optional :: FDorder

    real(kind=dp) :: delta
    integer :: N, order, ii

    N = size(x_grid)
    delta = abs(x_grid(2) - x_grid(1))

    if (present(FDorder)) then
      order = FDorder
    else
      order = 2
    end if

    ! Constant profile (single material)
    if (allocated(profile)) deallocate(profile)
    allocate(profile(N, 3))
    profile(:, 1) = params(1)%EV
    profile(:, 2) = params(1)%EV - params(1)%DeltaSO
    profile(:, 3) = params(1)%EC

    ! Diagonal material parameters (constant across all grid points)
    do concurrent (ii = 1:N)
      kpterms(ii, ii, 1) = params(1)%gamma1
      kpterms(ii, ii, 2) = params(1)%gamma2
      kpterms(ii, ii, 3) = params(1)%gamma3
      kpterms(ii, ii, 4) = params(1)%P
      kpterms(ii, ii, 10) = params(1)%A
    end do

    ! Build FD stencil terms using the same machinery as QW
    ! For a single material, the profile vector is constant
    block
      real(kind=dp), allocatable :: kptermsProfile(:,:)
      real(kind=dp), allocatable :: forward(:,:), backward(:,:), central(:,:)
      real(kind=dp), allocatable :: diag(:), offup(:), offdown(:)
      real(kind=dp), allocatable :: D_inner(:,:), D_outer(:,:), g_half(:)

      allocate(kptermsProfile(N, 5))
      kptermsProfile(:, 1) = params(1)%gamma1
      kptermsProfile(:, 2) = params(1)%gamma2
      kptermsProfile(:, 3) = params(1)%gamma3
      kptermsProfile(:, 4) = params(1)%A
      kptermsProfile(:, 5) = params(1)%P

      allocate(forward(N,N), backward(N,N), central(N,N))
      forward = 0.0_dp; backward = 0.0_dp; central = 0.0_dp
      do ii = 1, N - 1
        forward(ii, ii) = 1; forward(ii, ii+1) = 1
      end do
      forward(N, N) = 1
      backward = transpose(forward)
      central = backward + forward

      if (order == 2) then
        allocate(diag(N), offup(N), offdown(N))

        ! A*kx**2 (term 5)
        call build_kpterm_block(kpterms, kptermsProfile(:,4), central, forward, &
          & backward, diag, offup, offdown, N, 5, 1.0_dp/(2.0_dp*delta**2), .True.)
        ! Q (term 7): (gamma1-2*gamma2)
        call build_kpterm_block(kpterms, kptermsProfile(:,1) - 2.0_dp*kptermsProfile(:,2), &
          & central, forward, backward, diag, offup, offdown, N, 7, 1.0_dp/(2.0_dp*delta**2), .True.)
        ! T (term 8): (gamma1+2*gamma2)
        call build_kpterm_block(kpterms, kptermsProfile(:,1) + 2.0_dp*kptermsProfile(:,2), &
          & central, forward, backward, diag, offup, offdown, N, 8, 1.0_dp/(2.0_dp*delta**2), .True.)
        ! P*kx (term 6)
        call build_kpterm_block(kpterms, kptermsProfile(:,5), central, forward, &
          & backward, diag, offup, offdown, N, 6, 1.0_dp/(4.0_dp*delta), .False.)
        ! S -> gamma3*kx (term 9)
        call build_kpterm_block(kpterms, kptermsProfile(:,3), central, forward, &
          & backward, diag, offup, offdown, N, 9, 1.0_dp/(4.0_dp*delta), .False.)
      else
        ! Higher-order FD (same as QW path)
        call buildStaggeredD1Inner(N, delta, 2, D_inner)
        call buildStaggeredD1Outer(N, delta, 2, D_outer)
        allocate(g_half(N - 1))
        allocate(diag(N), offup(N), offdown(N))

        call interpolateToHalfPoints(kptermsProfile(:,4), N, order, g_half)
        call applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, 5)

        call interpolateToHalfPoints(kptermsProfile(:,1) - 2.0_dp*kptermsProfile(:,2), N, order, g_half)
        call applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, 7)

        call interpolateToHalfPoints(kptermsProfile(:,1) + 2.0_dp*kptermsProfile(:,2), N, order, g_half)
        call applyVariableCoeffStaggered(kpterms, g_half, D_inner, D_outer, N, 8)

        call build_kpterm_block(kpterms, kptermsProfile(:,5), central, forward, &
          & backward, diag, offup, offdown, N, 6, 1.0_dp/(4.0_dp*delta), .False.)
        call build_kpterm_block(kpterms, kptermsProfile(:,3), central, forward, &
          & backward, diag, offup, offdown, N, 9, 1.0_dp/(4.0_dp*delta), .False.)

        deallocate(D_inner, D_outer, g_half)
      end if

      deallocate(kptermsProfile, forward, backward, central, diag, offup, offdown)
    end block
  end subroutine confinementInitialization_landau
```

- [ ] **Step 3: Build and verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build succeeds

- [ ] **Step 4: Commit**

```bash
git add src/physics/confinement_init.f90
git commit -m "feat: add confinementInitialization_landau for single-material kpterms"
```

---

### Task 4: Implement `ZB8bandLandau` Assembly (TDD)

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90` (add subroutine)
- Modify: `tests/unit/test_landau.pf` (add hermiticity and bulk-recovery tests)

- [ ] **Step 1: Write failing tests for hermiticity and bulk recovery**

Append to `tests/unit/test_landau.pf`, before `end module`:

```fortran
  @test
  subroutine test_landau_hermitian()
    ! ZB8bandLandau should produce a Hermitian matrix at nonzero B
    use parameters
    use confinement_init, only: confinementInitialization_landau
    use magnetic_field, only: compute_zeeman_vz

    type(paramStruct) :: params(1)
    character(len=255) :: material(1)
    type(simulation_config) :: cfg
    type(wavevector) :: wv

    integer, parameter :: N = 21
    real(kind=dp) :: x_grid(N), delta_x
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    complex(kind=dp), allocatable :: HT(:,:)
    real(kind=dp) :: max_err
    integer :: i, j

    material(1) = "InAs"
    call paramDatabase(material, 1, params)

    delta_x = 10.0_dp  ! 10 AA spacing, 200 AA total
    do i = 1, N
      x_grid(i) = (i - 1) * delta_x
    end do

    allocate(profile(N, 3), kpterms(N, N, 10))
    kpterms = 0.0_dp
    call confinementInitialization_landau(x_grid, material, params, profile, kpterms)

    cfg%confinement = 3
    cfg%bdg%B_vec = [0.0_dp, 0.0_dp, 5.0_dp]
    cfg%bdg%g_factor = 2.0_dp
    cfg%bdg%enabled = .false.

    wv%kx = 0.0_dp; wv%ky = 0.0_dp; wv%kz = 0.0_dp

    allocate(HT(8*N, 8*N))
    HT = (0.0_dp, 0.0_dp)
    call ZB8bandLandau(HT, wv, profile, kpterms, x_grid, cfg)

    max_err = 0.0_dp
    do j = 1, 8*N
      do i = 1, 8*N
        max_err = max(max_err, abs(HT(i,j) - conjg(HT(j,i))))
      end do
    end do
    @assertTrue(max_err < 1.0e-12_dp, message="Landau Hamiltonian is Hermitian")

    deallocate(HT, profile, kpterms)
  end subroutine

  @test
  subroutine test_landau_recovers_bulk_at_B0()
    ! At B=0, the Landau Hamiltonian should give the same eigenvalues as bulk
    use parameters
    use confinement_init, only: confinementInitialization_landau
    use linalg, only: zheevx

    type(paramStruct) :: params(1)
    character(len=255) :: material(1)
    type(simulation_config) :: cfg
    type(wavevector) :: wv

    integer, parameter :: N = 21
    real(kind=dp) :: x_grid(N), delta_x
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    complex(kind=dp), allocatable :: HT(:,:), HT_bulk(:,:)

    real(kind=dp) :: eig_landau(8), eig_bulk(8)
    real(kind=dp), allocatable :: work(:), rwork(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info, M, lwork, i

    material(1) = "InAs"
    call paramDatabase(material, 1, params)

    delta_x = 10.0_dp
    do i = 1, N
      x_grid(i) = (i - 1) * delta_x
    end do

    allocate(profile(N, 3), kpterms(N, N, 10))
    kpterms = 0.0_dp
    call confinementInitialization_landau(x_grid, material, params, profile, kpterms)

    cfg%confinement = 3
    cfg%bdg%B_vec = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%bdg%g_factor = 2.0_dp
    cfg%bdg%enabled = .false.

    ! k at Gamma (k_y=0, k_z=0)
    wv%kx = 0.0_dp; wv%ky = 0.0_dp; wv%kz = 0.0_dp

    ! Solve Landau Hamiltonian
    allocate(HT(8*N, 8*N))
    HT = (0.0_dp, 0.0_dp)
    call ZB8bandLandau(HT, wv, profile, kpterms, x_grid, cfg)

    ! Extract the central 8x8 block (grid point at center)
    ! At B=0, all blocks are identical => central block has bulk eigenvalues
    block
      integer :: ic
      real(kind=dp) :: eig_block(8)
      complex(kind=dp) :: Hblock(8,8)
      ic = N / 2  ! center grid point
      Hblock = HT((ic-1)*8+1:ic*8, (ic-1)*8+1:ic*8)
      ! Diagonalize
      allocate(work(1), rwork(7*8), iwork(5*8), ifail(8))
      lwork = -1
      call zheevx('N', 'A', 'U', 8, Hblock, 8, 0.0_dp, 0.0_dp, 0, 0, &
        & 0.0_dp, M, eig_block, Hblock, 8, work, lwork, rwork, iwork, ifail, info)
      lwork = int(real(work(1)))
      deallocate(work); allocate(work(lwork))
      call zheevx('N', 'A', 'U', 8, Hblock, 8, 0.0_dp, 0.0_dp, 0, 0, &
        & 0.0_dp, M, eig_block, Hblock, 8, work, lwork, rwork, iwork, ifail, info)
      eig_landau = eig_block
      deallocate(work, rwork, iwork, ifail)
    end block

    ! Solve bulk Hamiltonian for comparison
    allocate(HT_bulk(8, 8))
    HT_bulk = (0.0_dp, 0.0_dp)
    call ZB8bandBulk(HT_bulk, wv, params)
    allocate(work(1), rwork(7*8), iwork(5*8), ifail(8))
    lwork = -1
    call zheevx('N', 'A', 'U', 8, HT_bulk, 8, 0.0_dp, 0.0_dp, 0, 0, &
      & 0.0_dp, M, eig_bulk, HT_bulk, 8, work, lwork, rwork, iwork, ifail, info)
    lwork = int(real(work(1)))
    deallocate(work); allocate(work(lwork))
    call zheevx('N', 'A', 'U', 8, HT_bulk, 8, 0.0_dp, 0.0_dp, 0, 0, &
      & 0.0_dp, M, eig_bulk, HT_bulk, 8, work, lwork, rwork, iwork, ifail, info)
    deallocate(work, rwork, iwork, ifail)

    ! The diagonal block should match bulk eigenvalues
    ! (within FD discretization error)
    do i = 1, 8
      @assertTrue(abs(eig_landau(i) - eig_bulk(i)) < 0.01_dp, &
        & message="Landau B=0 recovers bulk eigenvalues")
    end do

    deallocate(HT, HT_bulk, profile, kpterms)
  end subroutine
```

- [ ] **Step 2: Build and verify tests fail (ZB8bandLandau not found)**

Run: `cmake --build build 2>&1 | grep -i error`
Expected: Compilation error — `ZB8bandLandau` not found

- [ ] **Step 3: Implement `ZB8bandLandau`**

Add to `src/physics/hamiltonianConstructor.f90`, after the `ZB8bandQW` subroutine (after its `end subroutine`), before the `ZB8bandBulk` subroutine. The subroutine follows the same block structure as `ZB8bandQW` but with position-dependent k-values:

```fortran
  subroutine ZB8bandLandau(HT, wv, profile, kpterms, x_grid, cfg)
    ! Landau Hamiltonian for bulk in magnetic field.
    ! Discretizes x-direction with FD; k_y, k_z are parameters with gauge shifts.
    ! Landau gauge: A = (0, Bz*x, -By*x)
    !   Pi_y(i) = ky + Bz*x(i)*1e-20/hbar
    !   Pi_z(i) = kz - By*x(i)*1e-20/hbar
    ! B_x contributes Zeeman only (parallel to discretization direction).
    use magnetic_field, only: compute_zeeman_vz

    complex(kind=dp), intent(inout) :: HT(:,:)
    type(wavevector), intent(in) :: wv
    real(kind=dp), intent(in) :: profile(:,:), kpterms(:,:,:), x_grid(:)
    type(simulation_config), intent(in), optional :: cfg

    integer :: i, j, N
    real(kind=dp) :: By, Bz, inv_hbar
    real(kind=dp), allocatable :: k2(:), kx2(:), ky2(:), kxky(:)
    complex(kind=dp), allocatable :: kplus(:), kminus(:)
    complex(kind=dp), allocatable :: kplus_md(:), kminus_md(:)
    real(kind=dp), allocatable :: Pi_y_md(:), Pi_z_md(:)

    complex(kind=dp), allocatable :: Q(:,:), T(:,:), S(:,:), SC(:,:), R(:,:), RC(:,:)
    complex(kind=dp), allocatable :: PZ(:,:), PP(:,:), PM(:,:), A(:,:)

    N = size(HT, dim=1) / 8

    By = 0.0_dp; Bz = 0.0_dp
    if (present(cfg)) then
      By = cfg%bdg%B_vec(2)
      Bz = cfg%bdg%B_vec(3)
    end if
    inv_hbar = 1.0_dp / hbar

    ! --- Per-grid-point k-values (diagonal blocks) ---
    allocate(k2(N), kx2(N), ky2(N), kxky(N), kplus(N), kminus(N))
    do i = 1, N
      kx2(i) = (wv%ky + Bz * x_grid(i) * 1.0e-20_dp * inv_hbar)**2
      ky2(i) = (wv%kz - By * x_grid(i) * 1.0e-20_dp * inv_hbar)**2
      k2(i) = kx2(i) + ky2(i)
      kxky(i) = (wv%ky + Bz * x_grid(i) * 1.0e-20_dp * inv_hbar) &
            * (wv%kz - By * x_grid(i) * 1.0e-20_dp * inv_hbar)
      kplus(i) = (wv%ky + Bz * x_grid(i) * 1.0e-20_dp * inv_hbar) &
            + IU * (wv%kz - By * x_grid(i) * 1.0e-20_dp * inv_hbar)
      kminus(i) = (wv%ky + Bz * x_grid(i) * 1.0e-20_dp * inv_hbar) &
            - IU * (wv%kz - By * x_grid(i) * 1.0e-20_dp * inv_hbar)
    end do

    ! --- Midpoint k-values (off-diagonal FD hopping) ---
    allocate(kplus_md(N-1), kminus_md(N-1), Pi_y_md(N-1), Pi_z_md(N-1))
    do i = 1, N - 1
      block
        real(kind=dp) :: x_mid
        x_mid = 0.5_dp * (x_grid(i) + x_grid(i+1))
        Pi_y_md(i) = wv%ky + Bz * x_mid * 1.0e-20_dp * inv_hbar
        Pi_z_md(i) = wv%kz - By * x_mid * 1.0e-20_dp * inv_hbar
        kplus_md(i) = Pi_y_md(i) + IU * Pi_z_md(i)
        kminus_md(i) = Pi_y_md(i) - IU * Pi_z_md(i)
      end block
    end do

    ! --- Build k.p matrices ---
    allocate(Q(N,N), T(N,N), S(N,N), SC(N,N), R(N,N), RC(N,N))
    allocate(PZ(N,N), PP(N,N), PM(N,N), A(N,N))
    Q = ZERO; T = ZERO; S = ZERO; SC = ZERO
    R = ZERO; RC = ZERO; PZ = ZERO; PP = ZERO; PM = ZERO; A = ZERO

    ! Diagonal terms (position-dependent k)
    do i = 1, N
      Q(i,i) = -((kpterms(i,i,1) + kpterms(i,i,2))*k2(i) + kpterms(i,i,7))
      T(i,i) = -((kpterms(i,i,1) - kpterms(i,i,2))*k2(i) + kpterms(i,i,8))
      A(i,i) = kpterms(i,i,5) + k2(i)*kpterms(i,i,10)
      PZ(i,i) = kpterms(i,i,6) * (-IU)
      R(i,i) = -SQR3*(kpterms(i,i,2)*(kx2(i) - ky2(i)) - 2.0_dp*IU*kpterms(i,i,3)*kxky(i))
      RC(i,i) = -SQR3*(kpterms(i,i,2)*(kx2(i) - ky2(i)) + 2.0_dp*IU*kpterms(i,i,3)*kxky(i))
      PP(i,i) = kpterms(i,i,4) * kplus(i) * RQS2
      PM(i,i) = kpterms(i,i,4) * kminus(i) * RQS2
    end do

    ! Off-diagonal FD stencil terms
    ! kpterms(ii,jj,1/2/3/10) = 0 for ii/=jj, so k2 factor doesn't matter
    ! S/SC use midpoint k-values
    do i = 1, N - 1
      ! Upper: (i, i+1)
      j = i + 1
      Q(i,j) = -kpterms(i,j,7)
      T(i,j) = -kpterms(i,j,8)
      A(i,j) = kpterms(i,j,5)
      PZ(i,j) = kpterms(i,j,6) * (-IU)
      S(i,j) = 2.0_dp * SQR3 * kminus_md(i) * kpterms(i,j,9)
      SC(j,i) = 2.0_dp * SQR3 * kplus_md(i) * kpterms(i,j,9)

      ! Lower: (i+1, i)
      Q(j,i) = -kpterms(j,i,7)
      T(j,i) = -kpterms(j,i,8)
      A(j,i) = kpterms(j,i,5)
      PZ(j,i) = kpterms(j,i,6) * (-IU)
      S(j,i) = 2.0_dp * SQR3 * kminus_md(i) * kpterms(j,i,9)
      SC(i,j) = 2.0_dp * SQR3 * kplus_md(i) * kpterms(j,i,9)
    end do

    ! --- Assemble 8Nx8N Hamiltonian (identical block structure to QW) ---
    HT = ZERO

    ! col 1
    HT(1+0*N:1*N, 1+0*N:1*N) = Q
    HT(1+0*N:1*N, 1+1*N:2*N) = SC
    HT(1+0*N:1*N, 1+2*N:3*N) = RC
    HT(1+0*N:1*N, 1+4*N:5*N) = -IU*RQS2*SC
    HT(1+0*N:1*N, 1+5*N:6*N) = IU*SQR2*RC
    HT(1+0*N:1*N, 1+6*N:7*N) = IU*PP

    ! col 2
    HT(1+1*N:2*N, 1+0*N:1*N) = S
    HT(1+1*N:2*N, 1+1*N:2*N) = T
    HT(1+1*N:2*N, 1+3*N:4*N) = RC
    HT(1+1*N:2*N, 1+4*N:5*N) = IU*RQS2*(Q - T)
    HT(1+1*N:2*N, 1+5*N:6*N) = -IU*SQR3*RQS2*SC
    HT(1+1*N:2*N, 1+6*N:7*N) = SQR2*RQS3*PZ
    HT(1+1*N:2*N, 1+7*N:8*N) = -RQS3*PP

    ! col 3
    HT(1+2*N:3*N, 1+0*N:1*N) = R
    HT(1+2*N:3*N, 1+2*N:3*N) = T
    HT(1+2*N:3*N, 1+3*N:4*N) = -SC
    HT(1+2*N:3*N, 1+4*N:5*N) = IU*SQR3*RQS2*S
    HT(1+2*N:3*N, 1+5*N:6*N) = IU*RQS2*(Q - T)
    HT(1+2*N:3*N, 1+6*N:7*N) = IU*RQS3*PM
    HT(1+2*N:3*N, 1+7*N:8*N) = IU*SQR2*RQS3*PZ

    ! col 4
    HT(1+3*N:4*N, 1+1*N:2*N) = R
    HT(1+3*N:4*N, 1+2*N:3*N) = -S
    HT(1+3*N:4*N, 1+3*N:4*N) = Q
    HT(1+3*N:4*N, 1+4*N:5*N) = IU*SQR2*R
    HT(1+3*N:4*N, 1+5*N:6*N) = IU*RQS2*S
    HT(1+3*N:4*N, 1+7*N:8*N) = -PM

    ! col 5
    HT(1+4*N:5*N, 1+0*N:1*N) = IU*RQS2*S
    HT(1+4*N:5*N, 1+1*N:2*N) = -IU*RQS2*(Q - T)
    HT(1+4*N:5*N, 1+2*N:3*N) = -IU*SQR3*RQS2*SC
    HT(1+4*N:5*N, 1+3*N:4*N) = -IU*SQR2*RC
    HT(1+4*N:5*N, 1+4*N:5*N) = 0.5_dp*(Q + T)
    HT(1+4*N:5*N, 1+6*N:7*N) = IU*RQS3*PZ
    HT(1+4*N:5*N, 1+7*N:8*N) = IU*SQR2*RQS3*PP

    ! col 6
    HT(1+5*N:6*N, 1+0*N:1*N) = -IU*SQR2*R
    HT(1+5*N:6*N, 1+1*N:2*N) = IU*SQR3*RQS2*S
    HT(1+5*N:6*N, 1+2*N:3*N) = -IU*RQS2*(Q - T)
    HT(1+5*N:6*N, 1+3*N:4*N) = -IU*RQS2*SC
    HT(1+5*N:6*N, 1+5*N:6*N) = 0.5_dp*(Q + T)
    HT(1+5*N:6*N, 1+6*N:7*N) = SQR2*RQS3*PM
    HT(1+5*N:6*N, 1+7*N:8*N) = -RQS3*PZ

    ! col 7
    HT(1+6*N:7*N, 1+0*N:1*N) = -IU*PM
    HT(1+6*N:7*N, 1+1*N:2*N) = SQR2*RQS3*PZ
    HT(1+6*N:7*N, 1+2*N:3*N) = -IU*RQS3*PP
    HT(1+6*N:7*N, 1+4*N:5*N) = -IU*RQS3*PZ
    HT(1+6*N:7*N, 1+5*N:6*N) = SQR2*RQS3*PP
    HT(1+6*N:7*N, 1+6*N:7*N) = A

    ! col 8
    HT(1+7*N:8*N, 1+1*N:2*N) = -RQS3*PM
    HT(1+7*N:8*N, 1+2*N:3*N) = -IU*SQR2*RQS3*PZ
    HT(1+7*N:8*N, 1+3*N:4*N) = -PP
    HT(1+7*N:8*N, 1+4*N:5*N) = -IU*SQR2*RQS3*PM
    HT(1+7*N:8*N, 1+5*N:6*N) = -RQS3*PZ
    HT(1+7*N:8*N, 1+7*N:8*N) = A

    ! Profile (band edges)
    do i = 1, N
      HT(      i,      i) = HT(      i,      i) + profile(i,1)
      HT(  N + i,  N + i) = HT(  N + i,  N + i) + profile(i,1)
      HT(2*N + i,2*N + i) = HT(2*N + i,2*N + i) + profile(i,1)
      HT(3*N + i,3*N + i) = HT(3*N + i,3*N + i) + profile(i,1)
      HT(4*N + i,4*N + i) = HT(4*N + i,4*N + i) + profile(i,2)
      HT(5*N + i,5*N + i) = HT(5*N + i,5*N + i) + profile(i,2)
      HT(6*N + i,6*N + i) = HT(6*N + i,6*N + i) + profile(i,3)
      HT(7*N + i,7*N + i) = HT(7*N + i,7*N + i) + profile(i,3)
    end do

    ! Zeeman splitting (full |B| magnitude)
    if (present(cfg)) then
      block
        real(kind=dp) :: B_mag, Vz(8)
        B_mag = sqrt(sum(cfg%bdg%B_vec**2))
        if (B_mag > 1.0e-12_dp) then
          call compute_zeeman_vz(cfg%bdg%g_factor, mu_B, B_mag, Vz)
          do i = 1, N
            HT(      i,      i) = HT(      i,      i) + Vz(1)
            HT(  N + i,  N + i) = HT(  N + i,  N + i) + Vz(2)
            HT(2*N + i,2*N + i) = HT(2*N + i,2*N + i) + Vz(3)
            HT(3*N + i,3*N + i) = HT(3*N + i,3*N + i) + Vz(4)
            HT(4*N + i,4*N + i) = HT(4*N + i,4*N + i) + Vz(5)
            HT(5*N + i,5*N + i) = HT(5*N + i,5*N + i) + Vz(6)
            HT(6*N + i,6*N + i) = HT(6*N + i,6*N + i) + Vz(7)
            HT(7*N + i,7*N + i) = HT(7*N + i,7*N + i) + Vz(8)
          end do
        end if
      end block
    end if

    deallocate(Q, T, S, SC, R, RC, PZ, PP, PM, A)
    deallocate(k2, kx2, ky2, kxky, kplus, kminus)
    deallocate(kplus_md, kminus_md, Pi_y_md, Pi_z_md)
  end subroutine ZB8bandLandau
```

- [ ] **Step 4: Build and run unit tests**

Run: `cmake --build build && ctest --test-dir build -R test_landau --output-on-failure`
Expected: All 5 tests PASS (3 gauge shift + hermiticity + bulk recovery)

- [ ] **Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 tests/unit/test_landau.pf
git commit -m "feat: add ZB8bandLandau assembly with hermiticity and bulk-recovery tests"
```

---

### Task 5: Update Input Parser for Confinement=3

**Files:**
- Modify: `src/io/input_parser.f90`

This task adds parsing for `landau_nx`, `landau_width`, `landau_sweep`, accepts `confinement: 3`, and decouples `b_field` from `bdg%enabled` when confinement=3.

- [ ] **Step 1: Accept confinement=3**

In the validation block after reading confinement (around line 211), change:

```fortran
    if (cfg%confinement < 0 .or. cfg%confinement > 2) then
      print *, 'Error: confinement must be 0, 1, or 2, got:', cfg%confinement
      stop 1
    end if
```

to:

```fortran
    if (cfg%confinement < 0 .or. cfg%confinement > 3) then
      print *, 'Error: confinement must be 0, 1, 2, or 3, got:', cfg%confinement
      stop 1
    end if
```

- [ ] **Step 2: Add Landau-specific parsing after wire block**

After the `confinement == 2` block ends (after line 466 where `end if` closes the wire block), add:

```fortran
    else if (cfg%confinement == 3) then
      ! ---- Landau mode (1D x-discretization for orbital Landau quantization) ----
      cfg%confDir = 'x'

      ! Landau grid parameters
      read(data_unit, *, iostat=status) label, cfg%landau_nx
      if (status /= 0) then
        print *, 'Error: Failed to read landau_nx from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%landau_nx

      read(data_unit, *, iostat=status) label, cfg%landau_width
      if (status /= 0) then
        print *, 'Error: Failed to read landau_width from input.cfg'
        stop 1
      end if
      print *, trim(label), cfg%landau_width

      ! Landau sweep mode (optional, default 'ky')
      call read_next_data_line(data_unit, line, status)
      if (status == 0) then
        colon_pos = index(line, ':')
        if (colon_pos > 0) then
          label = adjustl(line(:colon_pos-1))
          if (trim(to_lower_ascii(label)) == 'landau_sweep') then
            read(line(colon_pos+1:), *, iostat=status) cfg%landau_sweep
            if (status /= 0) then
              print *, 'Error: Failed to parse landau_sweep value'
              stop 1
            end if
            print *, 'landau_sweep:', trim(cfg%landau_sweep)
          else
            backspace(data_unit)
          end if
        else
          backspace(data_unit)
        end if
      end if
```

- [ ] **Step 3: Decouple b_field from bdg%enabled for Landau mode**

In the b_field parsing path (around lines 531-557), change the line that sets `cfg%bdg%enabled = .true.` to be conditional:

For Path B (standalone b_field parsing, around line 543), change:
```fortran
              cfg%bdg%B_vec = cfg%b_field
              cfg%bdg%enabled = .true.
```
to:
```fortran
              cfg%bdg%B_vec = cfg%b_field
              if (cfg%confinement /= 3) cfg%bdg%enabled = .true.
```

Apply the same change to Path A (Evalue collision, around line 516):
```fortran
              cfg%bdg%B_vec = cfg%b_field
              if (cfg%confinement /= 3) cfg%bdg%enabled = .true.
```

- [ ] **Step 4: Set up the spatial grid for confinement=3**

In `read_and_setup` where `init_grid_from_config` is called (around line 1118), the grid initialization already handles case(3) from Task 1. The parser also needs to set up the material layers. After reading `numLayers` and material data, ensure the Landau path (confinement=3) sets up `cfg%z` (reuse for x-grid) and `cfg%fdStep`:

After the Landau parsing block added in Step 2, add:

```fortran
      ! Set up grid and fdStep for compatibility with existing infrastructure
      cfg%fdStep = cfg%landau_nx
      allocate(cfg%z(cfg%landau_nx))
      do i = 1, cfg%landau_nx
        cfg%z(i) = (i - 1) * cfg%landau_width / real(cfg%landau_nx - 1, kind=dp)
      end do
      cfg%dz = cfg%landau_width / real(cfg%landau_nx - 1, kind=dp)
      cfg%totalSize = cfg%landau_width
      cfg%intStartPos = [1]
      cfg%intEndPos = [cfg%landau_nx]
```

- [ ] **Step 5: Build and verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build succeeds

- [ ] **Step 6: Commit**

```bash
git add src/io/input_parser.f90
git commit -m "feat: parse confinement=3, landau_nx, landau_width, landau_sweep"
```

---

### Task 6: Add confinement=3 Dispatch in `main.f90`

**Files:**
- Modify: `src/apps/main.f90`

- [ ] **Step 1: Add confinement=3 initialization after wire block**

After the `if (cfg%confinement == 2) then ... end if` block in main.f90 (after line 474), add the Landau initialization. The initialization needs to call `confinementInitialization_landau`, which is done in `read_and_setup` for QW but must be done in main for Landau (since `confinementInitialization_cfg` only handles `confDir=='z'`):

```fortran
  ! ---- Landau mode initialization ----
  if (cfg%confinement == 3) then
    block
      character(len=255) :: lmat(1)
      type(paramStruct) :: lparams(1)

      lmat = cfg%materialN(1:1)
      call paramDatabase(lmat, 1, lparams)
      cfg%params(1) = lparams(1)

      N = cfg%landau_nx * 8
      allocate(kpterms(cfg%landau_nx, cfg%landau_nx, 10))
      kpterms = 0.0_dp
      call confinementInitialization_landau(cfg%grid%x, lmat, lparams, &
        & profile, kpterms, cfg%FDorder)
    end block

    ! Eigenvalue index range
    il = NUM_VB_STATES*cfg%landau_nx - cfg%numvb + 1
    iuu = NUM_VB_STATES*cfg%landau_nx + cfg%numcb
  end if
```

Add `use confinement_init, only: confinementInitialization_2d, confinementInitialization_landau` to the module imports (update the existing `use confinement_init` line).

- [ ] **Step 2: Add confinement=3 k-sweep dispatch**

In the dense k-sweep section, after the QW OpenMP block and before `call writeEigenvalues`, add a Landau branch. The Landau path parallels the QW path but calls `ZB8bandLandau` instead of `ZB8bandQW`:

After the QW `else if (cfg%confDir == 'z') then` block ends, before `end if` of the outer confDir check, add:

```fortran
    else if (cfg%confinement == 3) then
    ! --- LANDAU (8NxN, OpenMP parallel) ---
    info = mkl_set_num_threads_local(1)

    print '(A,I0,A,I0,A)', ' Landau: N=', cfg%landau_nx, ', ', cfg%waveVectorStep, ' k-points'

    ! Workspace query
    call ZB8bandLandau(HT, smallk(1), profile, kpterms, cfg%grid%x, cfg=cfg)
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,1), &
               HT, N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print *, 'Error: zheevx Landau workspace query failed, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))

    ! OpenMP k-sweep
    block
      complex(kind=dp), allocatable :: HT_loc(:,:), work_loc(:)
      real(kind=dp), allocatable    :: rwork_loc(:)
      integer, allocatable          :: iwork_loc(:), ifail_loc(:)
      integer :: info_loc, M_loc

      !$omp parallel private(k, HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc, info_loc, M_loc)
      allocate(HT_loc(N, N), work_loc(lwork), rwork_loc(7*N), iwork_loc(5*N), ifail_loc(N))
      HT_loc = ZERO

      !$omp do schedule(static)
      do k = 1, cfg%waveVectorStep
        call ZB8bandLandau(HT_loc, smallk(k), profile, kpterms, cfg%grid%x, cfg=cfg)
        call zheevx('V', 'I', 'U', N, HT_loc, N, vl, vu, il, iuu, abstol, M_loc, &
                   eig(:,k), HT_loc, N, work_loc, lwork, rwork_loc, iwork_loc, &
                   ifail_loc, info_loc)
        if (info_loc /= 0) then
          !$omp critical
          print *, "ERROR: Landau diagonalization failed at k=", k, "info=", info_loc
          !$omp end critical
          stop 1
        end if
        eigv(:,:,k) = HT_loc(:, 1:iuu-il+1)
      end do
      !$omp end do

      deallocate(HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc)
    !$omp end parallel
    end block

    ! Write eigenfunctions at start, middle, end
    do k = 1, cfg%waveVectorStep
      if (k == 1 .or. k == int(cfg%waveVectorStep/2) .or. k == cfg%waveVectorStep) then
        call writeEigenfunctions(N, iuu-il+1, eigv(:,1:iuu-il+1,k), &
          & k, cfg%landau_nx, cfg%grid%x, .false.)
      end if
    end do
  end if
```

- [ ] **Step 3: Build and verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build succeeds

- [ ] **Step 4: Run existing tests to verify no regressions**

Run: `OMP_NUM_THREADS=4 ctest --test-dir build -j4 --output-on-failure`
Expected: All 53 tests PASS (no existing tests affected)

- [ ] **Step 5: Commit**

```bash
git add src/apps/main.f90
git commit -m "feat: add confinement=3 Landau dispatch with OpenMP k-sweep"
```

---

### Task 7: Create Regression Configs and Run End-to-End Tests

**Files:**
- Create: `tests/regression/configs/landau_bulk_InAs.cfg`
- Create: `tests/regression/configs/landau_bulk_GaAs.cfg`
- Modify: `input.cfg` (temporary, for testing — never committed)

- [ ] **Step 1: Create InAs Landau regression config**

Create `tests/regression/configs/landau_bulk_InAs.cfg`:

```
waveVector: ky
waveVectorMax: 0.05
waveVectorStep: 50
confinement: 3
landau_nx: 100
landau_width: 2000.0
numLayers: 1
material1: InAs
numcb: 4
numvb: 4
FDorder: 2
b_field: 0.0 0.0 5.0
g_factor: 2.0
```

- [ ] **Step 2: Create GaAs Landau regression config**

Create `tests/regression/configs/landau_bulk_GaAs.cfg`:

```
waveVector: ky
waveVectorMax: 0.05
waveVectorStep: 50
confinement: 3
landau_nx: 100
landau_width: 2000.0
numLayers: 1
material1: GaAs
numcb: 4
numvb: 4
FDorder: 2
b_field: 0.0 0.0 5.0
g_factor: 2.0
```

- [ ] **Step 3: Run InAs Landau test manually**

```bash
cp tests/regression/configs/landau_bulk_InAs.cfg input.cfg
mkdir -p output
OMP_NUM_THREADS=4 ./build/src/bandStructure
```

Expected: Completes without error, produces `output/eigenvalues.dat`.

- [ ] **Step 4: Verify Landau level spacing**

Run a quick check of the eigenvalues. At B=5T, InAs m*=0.026m₀:
- ℏω_c = ℏ·e·B/m* = 6.582e-16 * 1.602e-19 * 5 / (0.026 * 9.109e-31) / 1.602e-19 = 0.02226 eV = 22.26 meV
- The CB eigenvalue spacing should be approximately 22 meV

Examine `output/eigenvalues.dat`:
```bash
head -3 output/eigenvalues.dat
```

Expected: The first few eigenvalues show CB states separated by ~22 meV.

- [ ] **Step 5: Verify k_y degeneracy (eigenvalues should be flat vs ky)**

The eigenvalues at different ky values should be nearly identical (Landau level degeneracy). Check:

```bash
# Compare eigenvalues at k=1 and k=25
awk 'NR==2 {print "k=1:", $0}' output/eigenvalues.dat
awk 'NR==26 {print "k=25:", $0}' output/eigenvalues.dat
```

Expected: Eigenvalue differences < 0.1 meV (small FD discretization error only).

- [ ] **Step 6: Run GaAs Landau test**

```bash
cp tests/regression/configs/landau_bulk_GaAs.cfg input.cfg
OMP_NUM_THREADS=4 ./build/src/bandStructure
```

Expected: Completes without error. CB spacing should be ~8.67 meV (heavier mass).

- [ ] **Step 7: Run full test suite to confirm no regressions**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure`
Expected: All 53 tests PASS + any new Landau tests

- [ ] **Step 8: Commit**

```bash
git add tests/regression/configs/landau_bulk_InAs.cfg tests/regression/configs/landau_bulk_GaAs.cfg
git commit -m "test: add InAs and GaAs Landau regression configs for confinement=3"
```

---

### Task 8: Add B-Sweep Support and Fan Diagram Output

**Files:**
- Modify: `src/apps/main.f90` (add B-sweep loop)
- Create: `tests/regression/configs/landau_bulk_InAs_Bsweep.cfg`

- [ ] **Step 1: Create B-sweep config**

Create `tests/regression/configs/landau_bulk_InAs_Bsweep.cfg`:

```
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 3
landau_nx: 100
landau_width: 2000.0
landau_sweep: B
numLayers: 1
material1: InAs
numcb: 4
numvb: 4
FDorder: 2
b_field: 0.0 0.0 0.5
g_factor: 2.0
b_sweep: 0.5 10.0 0.5
```

- [ ] **Step 2: Add B-sweep logic in main.f90**

Inside the `if (cfg%confinement == 3)` block in main.f90, after the standard k-sweep, add B-sweep support. When `cfg%landau_sweep == 'B'`, the code loops over B values instead of k-points:

After the Landau k-sweep block, before the `end if` of the confinement==3 check:

```fortran
    ! --- B-sweep for fan diagram ---
    if (trim(cfg%landau_sweep) == 'B') then
      block
        real(kind=dp) :: B_val, B_min, B_max, B_step
        integer :: nB, iB, nL
        real(kind=dp), allocatable :: fan_evals(:,:)
        type(wavevector) :: wv0
        integer :: iounit_fan

        B_min = cfg%bdg%B_sweep(1)
        B_max = cfg%bdg%B_sweep(2)
        B_step = cfg%bdg%B_sweep(3)
        nB = int((B_max - B_min) / B_step) + 1
        nL = iuu - il + 1

        wv0%kx = 0.0_dp; wv0%ky = 0.0_dp; wv0%kz = 0.0_dp
        allocate(fan_evals(nL, nB))

        print '(A,I0,A,g8.2,A,g8.2)', ' B-sweep: ', nB, ' steps from ', B_min, ' to ', B_max

        do iB = 1, nB
          B_val = B_min + (iB - 1) * B_step
          cfg%bdg%B_vec(3) = B_val  ! sweep Bz

          call ZB8bandLandau(HT, wv0, profile, kpterms, cfg%grid%x, cfg=cfg)
          call zheevx('N', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, &
                     fan_evals(:, iB), HT, N, work, lwork, rwork, iwork, ifail, info)
          if (info /= 0) then
            print *, 'ERROR: Landau B-sweep failed at B=', B_val, 'info=', info
            stop 1
          end if
        end do

        ! Write fan diagram
        open(newunit=iounit_fan, file='output/landau_fan.dat', status='replace', action='write')
        write(iounit_fan, '(A)', advance='no') '# B[T]'
        do i = 1, nL
          write(iounit_fan, '(A,I0)', advance='no') '  E_', i - 1
        end do
        write(iounit_fan, *)
        do iB = 1, nB
          B_val = B_min + (iB - 1) * B_step
          write(iounit_fan, '(*(1x,g14.6))') B_val, fan_evals(:, iB)
        end do
        close(iounit_fan)

        print '(A)', ' Fan diagram written to output/landau_fan.dat'
        deallocate(fan_evals)
      end block
    end if
```

- [ ] **Step 3: Run B-sweep test**

```bash
cp tests/regression/configs/landau_bulk_InAs_Bsweep.cfg input.cfg
OMP_NUM_THREADS=4 ./build/src/bandStructure
```

Expected: Produces `output/landau_fan.dat` with B values and eigenvalues. Check that eigenvalues scale linearly with B.

- [ ] **Step 4: Verify fan diagram linearity**

```bash
head -5 output/landau_fan.dat
```

Expected: CB eigenvalues increase linearly with B.

- [ ] **Step 5: Run full test suite**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure`
Expected: All tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/apps/main.f90 tests/regression/configs/landau_bulk_InAs_Bsweep.cfg
git commit -m "feat: add B-sweep fan diagram output for Landau mode"
```

---

## Self-Review Checklist

**1. Spec coverage:**
- [x] Gauge A=(0, Bz*x, -By*x) → Task 2 (compute_gauge_shifts), Task 4 (assembly)
- [x] Position-dependent Π_y, Π_z → Task 2, Task 4
- [x] Midpoint evaluation for FD hopping → Task 4 (kplus_md, kminus_md)
- [x] Zeeman splitting (full |B|) → Task 4 (copied from ZB8bandQW pattern)
- [x] Open boundaries → Reuses QW machinery (no change needed)
- [x] confinement=3 type/validator → Task 1
- [x] confinementInitialization_landau → Task 3
- [x] ZB8bandLandau assembly → Task 4
- [x] main.f90 dispatch → Task 6
- [x] Input parser → Task 5
- [x] landau_nx, landau_width, landau_sweep fields → Task 1, Task 5
- [x] B-field decoupled from bdg%enabled → Task 5
- [x] k_y degeneracy sweep → Task 7 (InAs config with ky sweep)
- [x] k_z dispersion sweep → Supported via waveVector: kz in config
- [x] B-sweep fan diagram → Task 8
- [x] Unit tests (4) → Task 2 (3 gauge) + Task 4 (2 assembly)
- [x] Regression configs (3) → Task 7 (2 materials) + Task 8 (B-sweep)

**2. Placeholder scan:** No TBDs, TODOs, or "implement later" patterns found.

**3. Type consistency:**
- `compute_gauge_shifts(x_grid, B_vec, ky, kz, Pi_y, Pi_z)` — all `real(dp)`, consistent across test and implementation
- `ZB8bandLandau(HT, wv, profile, kpterms, x_grid, cfg)` — `HT` complex(dp), `wv` wavevector, `cfg` optional simulation_config
- `confinementInitialization_landau(x_grid, material, params, profile, kpterms, FDorder)` — matches build_kpterm_block signatures
- Grid fields: `cfg%grid%x(landau_nx)`, `cfg%landau_nx`, `cfg%landau_width` — consistent across defs, parser, main
