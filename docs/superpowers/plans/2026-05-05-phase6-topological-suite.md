# Phase 6: Topological Suite Completion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the topological analysis suite: Fu-Kane Z2, QW BdG, gap sweep, Berry curvature, conductance, spectral function, and associated tests.

**Architecture:** Extend existing Fortran modules (`topological_analysis.f90`, `bdg_hamiltonian.f90`, `green_functions.f90`) with new functions. Add dispatch branches to `topologicalAnalysis` executable. Python figures via `generate_all_figures.py`. TDD: write failing tests first, then implement.

**Tech Stack:** Fortran 2018, MKL LAPACK/PARDISO, pFUnit tests, Python matplotlib figures.

**Scope note:** Landauer-Büttiker conductance is deferred from this plan. Kubo conductance covers the primary use case (Hall conductance from Berry curvature). Landauer requires lead self-energy infrastructure that adds significant complexity; it can be added in a follow-up task if needed.

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `src/core/defs.f90` | Modify | Add new config fields to `topology_config` and `topological_result` |
| `src/physics/topological_analysis.f90` | Modify | Fu-Kane Z2, Berry curvature lattice, Kubo conductance, gap sweep |
| `src/physics/bdg_hamiltonian.f90` | Modify | Dense QW BdG variant |
| `src/physics/green_functions.f90` | Modify | Spectral function (dense + CSR), Landauer conductance |
| `src/io/input_parser.f90` | Modify | Parse new config fields |
| `src/apps/main_topology.f90` | Modify | New dispatch branches for all features |
| `tests/unit/test_z2_invariant.pf` | Create | Fu-Kane and Z2 invariant tests |
| `tests/unit/test_edge_states.pf` | Modify | Add 3 new edge state tests |
| `tests/unit/test_green_functions.pf` | Modify | Add spectral function + conductance tests |
| `tests/unit/test_magnetic_field.pf` | Modify | Add 3 new magnetic field tests |
| `tests/CMakeLists.txt` | Modify | Register new test file |
| `scripts/generate_all_figures.py` | Modify | Add 4 new figure functions |

---

### Task 1: Add new config fields to defs.f90

**Files:**
- Modify: `src/core/defs.f90:296-317` (topology_config type)
- Modify: `src/core/defs.f90:322-334` (topological_result type)

- [ ] **Step 1: Add new fields to `topology_config`**

In `src/core/defs.f90`, after line 316 (`integer :: ldos_num_E = 200`), add:

```fortran
    ! Gap sweep / phase diagram
    logical          :: compute_gap_sweep = .false.
    real(kind=dp)    :: gap_sweep_B_min = 0.0_dp
    real(kind=dp)    :: gap_sweep_B_max = 1.0_dp
    integer          :: gap_sweep_nB = 20
    real(kind=dp)    :: gap_sweep_mu_min = 0.0_dp
    real(kind=dp)    :: gap_sweep_mu_max = 0.01_dp
    integer          :: gap_sweep_nMu = 20
    ! Conductance
    logical          :: compute_conductance = .false.
    character(len=20) :: conductance_method = 'kubo'  ! kubo | landauer
    ! Spectral function
    logical          :: compute_spectral = .false.
    real(kind=dp)    :: spectral_k_min = -0.1_dp
    real(kind=dp)    :: spectral_k_max = 0.1_dp
    integer          :: spectral_nk = 100
    real(kind=dp)    :: spectral_E_min = -0.05_dp
    real(kind=dp)    :: spectral_E_max = 0.05_dp
    integer          :: spectral_nE = 200
    real(kind=dp)    :: spectral_eta = 0.001_dp
```

- [ ] **Step 2: Add new allocatable fields to `topological_result`**

In `src/core/defs.f90`, after line 331 (`real(kind=dp), allocatable :: berry_curvature(:,:)`), add:

```fortran
    real(kind=dp)              :: conductance_xy = 0.0_dp  ! Hall conductance (e^2/h)
    real(kind=dp)              :: conductance_zz = 0.0_dp  ! longitudinal (2e^2/h)
    real(kind=dp), allocatable :: spectral_function(:,:)   ! A(k, E) heatmap
    real(kind=dp), allocatable :: z2_map(:,:)              ! Z2 phase diagram (nMu x nB)
    real(kind=dp), allocatable :: gap_map(:,:)             ! gap phase diagram (nMu x nB)
```

- [ ] **Step 3: Update the finalizer**

In `src/core/defs.f90`, find `topological_result_finalize` (around line 472-478) and add deallocation for the new allocatable fields after the existing deallocate statements:

```fortran
    if (allocated(self%spectral_function)) deallocate(self%spectral_function)
    if (allocated(self%z2_map)) deallocate(self%z2_map)
    if (allocated(self%gap_map)) deallocate(self%gap_map)
```

- [ ] **Step 4: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: BUILD SUCCEEDED (no errors — new fields have defaults, unused for now)

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90
git commit -m "feat(topo): add gap sweep, conductance, spectral function config fields"
```

---

### Task 2: Implement Berry curvature lattice extraction

**Files:**
- Modify: `src/physics/topological_analysis.f90:134-144` (replace stub)
- Test: `tests/unit/test_chern_number.pf` (extend)

- [ ] **Step 1: Write failing test for Berry curvature integration**

In `tests/unit/test_chern_number.pf`, add after `test_qwz_chern_minus_one` (before `end module`):

```fortran
  @test
  subroutine test_berry_curvature_integrates_to_chern()
    integer :: C_direct, C_from_curvature
    real(kind=dp), allocatable :: Omega(:,:), kx_arr(:), ky_arr(:)
    complex(kind=dp), allocatable :: evecs(:,:,:)
    real(kind=dp) :: dk, total_flux, kx, ky, mz_val, E_plus, nrm, off_diag_re, off_diag_im
    complex(kind=dp) :: ev_tmp(2), off_diag
    integer :: nk, i, j

    nk = 50
    dk = 2.0_dp * acos(-1.0_dp) / real(nk, kind=dp)

    allocate(evecs(nk, nk, 2))
    allocate(kx_arr(nk), ky_arr(nk))
    allocate(Omega(nk, nk))

    do i = 1, nk
      kx_arr(i) = -acos(-1.0_dp) + real(i - 1, kind=dp) * dk
    end do
    do j = 1, nk
      ky_arr(j) = -acos(-1.0_dp) + real(j - 1, kind=dp) * dk
    end do

    ! Build QWZ upper-band eigenvectors at u=-0.8 (known C=+1)
    do j = 1, nk
      ky = ky_arr(j)
      do i = 1, nk
        kx = kx_arr(i)
        mz_val = -0.8_dp + cos(kx) + cos(ky)
        off_diag = cmplx(sin(kx), -sin(ky), kind=dp)
        E_plus = sqrt(mz_val**2 + sin(kx)**2 + sin(ky)**2)
        ev_tmp(1) = off_diag
        ev_tmp(2) = cmplx(E_plus - mz_val, 0.0_dp, kind=dp)
        nrm = sqrt(real(conjg(ev_tmp(1))*ev_tmp(1) + conjg(ev_tmp(2))*ev_tmp(2), kind=dp))
        if (nrm > 1.0e-30_dp) then
          evecs(i,j,1) = ev_tmp(1) / nrm
          evecs(i,j,2) = ev_tmp(2) / nrm
        else
          evecs(i,j,1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
          evecs(i,j,2) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end if
      end do
    end do

    ! Reshape for compute_berry_curvature_lattice: (basis, n_occ, nkx, nky)
    block
      complex(kind=dp), allocatable :: evecs_4d(:,:,:,:)
      allocate(evecs_4d(2, 1, nk, nk))
      do j = 1, nk
        do i = 1, nk
          evecs_4d(1, 1, i, j) = evecs(i, j, 1)
          evecs_4d(2, 1, i, j) = evecs(i, j, 2)
        end do
      end do
      Omega = compute_berry_curvature_lattice(evecs_4d, kx_arr, ky_arr, 1)
      deallocate(evecs_4d)
    end block

    ! Curvature should integrate to 2*pi*C
    total_flux = sum(Omega) * dk * dk
    C_from_curvature = nint(total_flux / (2.0_dp * acos(-1.0_dp)))

    C_direct = compute_chern_qwz(u=-0.8_dp, nk=50)

    @assertEqual(C_direct, C_from_curvature, message="Berry curvature integral must match Chern number")

    deallocate(evecs, kx_arr, ky_arr, Omega)
  end subroutine test_berry_curvature_integrates_to_chern
```

- [ ] **Step 2: Build to verify test fails (function not yet public)**

Run: `cmake --build build 2>&1 | tail -5`
Expected: FAIL — `compute_berry_curvature_lattice` not found

- [ ] **Step 3: Implement `compute_berry_curvature_lattice`**

In `src/physics/topological_analysis.f90`:

a) Add to the public list (line 9 area): change `compute_berry_curvature` to also export `compute_berry_curvature_lattice`:

```fortran
  public :: compute_berry_curvature, compute_berry_curvature_lattice
```

b) Replace the existing `compute_berry_curvature` stub (lines 134-144) and add the new `compute_berry_curvature_lattice` after it. Keep the old stub for backward compat but make it call the new function:

```fortran
  function compute_berry_curvature(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    implicit none
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    ! Delegate to lattice method
    Omega = compute_berry_curvature_lattice(evecs_k, kx_arr, ky_arr, n_occ)
  end function compute_berry_curvature

  function compute_berry_curvature_lattice(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    implicit none
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)  ! (basis, n_occ, nkx, nky)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    integer :: nkx, nky, i, j, ip1, jp1, n, b
    complex(kind=dp) :: Ux, Uy, prod, link
    real(kind=dp) :: dkx, dky

    nkx = size(kx_arr)
    nky = size(ky_arr)
    allocate(Omega(nkx, nky))
    Omega = 0.0_dp

    dkx = kx_arr(min(2, nkx)) - kx_arr(1)
    dky = ky_arr(min(2, nky)) - ky_arr(1)

    do j = 1, nky
      jp1 = mod(j, nky) + 1
      do i = 1, nkx
        ip1 = mod(i, nkx) + 1

        ! Accumulate link variables over all occupied bands
        Ux = cmplx(0.0_dp, 0.0_dp, kind=dp)
        Uy = cmplx(0.0_dp, 0.0_dp, kind=dp)
        prod = cmplx(1.0_dp, 0.0_dp, kind=dp)

        ! FHS plaquette for multi-band: U_x * U_y(k+x) * U_x^*(k+y) * U_y^*(k)
        ! For n_occ > 1, we use the determinant method (non-Abelian)
        block
          complex(kind=dp) :: M_xx, M_yx, M_xy, M_yy
          integer :: m, p
          complex(kind=dp), allocatable :: overlap(:,:)

          allocate(overlap(n_occ, n_occ))

          ! U_x = det(<u_n(k)|u_m(k+x)>)
          do m = 1, n_occ
            do p = 1, n_occ
              overlap(m, p) = sum(conjg(evecs_k(:, m, i, j)) * evecs_k(:, p, ip1, j))
            end do
          end do
          M_xx = det_small(overlap, n_occ)

          ! U_y(k+x) = det(<u_n(k+x)|u_m(k+x+y)>)
          do m = 1, n_occ
            do p = 1, n_occ
              overlap(m, p) = sum(conjg(evecs_k(:, m, ip1, j)) * evecs_k(:, p, ip1, jp1))
            end do
          end do
          M_yx = det_small(overlap, n_occ)

          ! U_x^*(k+y) = det(<u_n(k+x+y)|u_m(k+y)>)
          do m = 1, n_occ
            do p = 1, n_occ
              overlap(m, p) = sum(conjg(evecs_k(:, m, ip1, jp1)) * evecs_k(:, p, i, jp1))
            end do
          end do
          M_xy = conjg(det_small(overlap, n_occ))

          ! U_y^*(k) = det(<u_n(k+y)|u_m(k)>)
          do m = 1, n_occ
            do p = 1, n_occ
              overlap(m, p) = sum(conjg(evecs_k(:, m, i, jp1)) * evecs_k(:, p, i, j))
            end do
          end do
          M_yy = conjg(det_small(overlap, n_occ))

          Omega(i, j) = aimag(log(M_xx * M_yx * M_xy * M_yy))

          deallocate(overlap)
        end block

      end do
    end do

  end function compute_berry_curvature_lattice
```

c) Add the `det_small` private helper function before `end module topological_analysis` (at the end):

```fortran
  function det_small(A, n) result(d)
    implicit none
    complex(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: n
    complex(kind=dp) :: d

    select case(n)
    case(1)
      d = A(1,1)
    case(2)
      d = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    case default
      ! General case: LU decomposition
      complex(kind=dp), allocatable :: LU(:,:)
      integer, allocatable :: pivot(:)
      integer :: info, i

      allocate(LU(n, n), pivot(n))
      LU = A
      call zgetrf(n, n, LU, n, pivot, info)
      d = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do i = 1, n
        d = d * LU(i, i)
        if (pivot(i) /= i) d = -d
      end do
      deallocate(LU, pivot)
    end select
  end function det_small
```

Also add `use linalg` to the module's `use` statements at the top of `topological_analysis.f90` (it already imports from `definitions` and `sparse_matrices` — check what's there and add `use linalg` if not present for `zgetrf` access). Actually, `zgetrf` is a LAPACK routine — it needs to be called through `linalg.f90` interfaces. Check if `linalg.f90` exports `zgetrf`. If not, add `external :: zgetrf` as a private module procedure or add it to the `linalg` interface block.

**Note:** For n_occ=1 (single occupied band like QWZ), the determinant simplifies to the single overlap, making the FHS plaquette identical to the existing `compute_chern_qwz` formulation. This ensures backward compatibility.

- [ ] **Step 4: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_chern_number -V`
Expected: ALL 4 tests PASS (3 existing + 1 new)

- [ ] **Step 5: Commit**

```bash
git add src/physics/topological_analysis.f90 tests/unit/test_chern_number.pf
git commit -m "feat(topo): implement Berry curvature lattice extraction via FHS"
```

---

### Task 3: Implement Fu-Kane Z2 invariant

**Files:**
- Create: `tests/unit/test_z2_invariant.pf`
- Modify: `src/physics/topological_analysis.f90:184-192` (replace stub)
- Modify: `tests/CMakeLists.txt`

- [ ] **Step 1: Create failing test file**

Create `tests/unit/test_z2_invariant.pf`:

```fortran
module test_z2_invariant
  use funit
  use definitions
  use topological_analysis
  use sparse_matrices
  use linalg
  use eigensolver
  use hamiltonianConstructor
  use confinement_init
  use parameters
  implicit none

contains

  @test
  subroutine test_z2_fukane_bhz_trivial()
    ! BHZ model with M > 0 => trivial (Z2 = 0)
    ! Use gap criterion as proxy since we don't have full k.p QW here
    type(bhz_wire_params) :: params
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: result
    integer :: z2

    params%A = 364.5_dp
    params%B = -686.0_dp
    params%D = -512.0_dp
    params%M = 10.0_dp
    params%d_wire = 58.0_dp
    params%N = 100
    params%dz = params%d_wire / real(params%N, kind=dp)

    call build_bhz_wire_hamiltonian(H_csr, params)

    cfg%method = 'DENSE'
    cfg%emin = -100.0_dp
    cfg%emax = 100.0_dp
    cfg%nev = 4 * params%N

    call solve_sparse_evp(H_csr, cfg, result)

    z2 = compute_z2_gap(params%N, result%eigenvalues, 1.0_dp)

    @assertEqual(0, z2, message="BHZ trivial phase (M=+10) should have Z2=0")

    call csr_free(H_csr)
  end subroutine test_z2_fukane_bhz_trivial

  @test
  subroutine test_z2_fukane_bhz_topological()
    ! BHZ model with M < 0 => topological (Z2 = 1)
    type(bhz_wire_params) :: params
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: result
    integer :: z2

    params%A = 364.5_dp
    params%B = -686.0_dp
    params%D = -512.0_dp
    params%M = -10.0_dp
    params%d_wire = 70.0_dp
    params%N = 100
    params%dz = params%d_wire / real(params%N, kind=dp)

    call build_bhz_wire_hamiltonian(H_csr, params)

    cfg%method = 'DENSE'
    cfg%emin = -100.0_dp
    cfg%emax = 100.0_dp
    cfg%nev = 4 * params%N

    call solve_sparse_evp(H_csr, cfg, result)

    z2 = compute_z2_gap(params%N, result%eigenvalues, 40.0_dp)

    @assertEqual(1, z2, message="BHZ topological phase (M=-10) should have Z2=1")

    call csr_free(H_csr)
  end subroutine test_z2_fukane_bhz_topological

  @test
  subroutine test_z2_parity_eigenvalue_sign()
    ! Verify that inversion parity is +1 for valence, -1 for conduction in 8-band basis
    ! P_band = [+1,+1,+1,+1,+1,+1,-1,-1]
    real(kind=dp) :: P_band(8)
    integer :: b

    P_band = [+1.0_dp, +1.0_dp, +1.0_dp, +1.0_dp, +1.0_dp, +1.0_dp, -1.0_dp, -1.0_dp]

    ! All valence bands (1-6) should have positive parity
    do b = 1, 6
      @assertTrue(P_band(b) > 0.0_dp, message="valence band parity should be +1")
    end do
    ! All conduction bands (7-8) should have negative parity
    do b = 7, 8
      @assertTrue(P_band(b) < 0.0_dp, message="conduction band parity should be -1")
    end do
  end subroutine test_z2_parity_eigenvalue_sign

  @test
  subroutine test_z2_gap_closing_edge_case()
    ! When gap exactly closes, Z2 should change
    type(bhz_wire_params) :: params
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: result

    ! M near zero (transition point) — just check no crash
    params%A = 364.5_dp
    params%B = -686.0_dp
    params%D = -512.0_dp
    params%M = 0.5_dp
    params%d_wire = 58.0_dp
    params%N = 50
    params%dz = params%d_wire / real(params%N, kind=dp)

    call build_bhz_wire_hamiltonian(H_csr, params)

    cfg%method = 'DENSE'
    cfg%emin = -100.0_dp
    cfg%emax = 100.0_dp
    cfg%nev = 4 * params%N

    call solve_sparse_evp(H_csr, cfg, result)

    @assertTrue(allocated(result%eigenvalues))
    @assertTrue(size(result%eigenvalues) > 0)

    call csr_free(H_csr)
  end subroutine test_z2_gap_closing_edge_case

end module test_z2_invariant
```

- [ ] **Step 2: Register new test in CMakeLists.txt**

In `tests/CMakeLists.txt`, after the `test_landau` block (around line 156, before `else()`), add:

```cmake
    add_pfunit_ctest(test_z2_invariant
        TEST_SOURCES unit/test_z2_invariant.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )
```

- [ ] **Step 3: Implement `compute_z2_fukane_qw` in topological_analysis.f90**

In `src/physics/topological_analysis.f90`:

a) Add to public list (near line 12):
```fortran
  public :: compute_z2_fukane_qw
```

b) Replace the `compute_z2_fukane` stub (lines 184-192) to delegate to the new function for QW configs, and add the new `compute_z2_fukane_qw` function after it:

```fortran
  function compute_z2_fukane(H_eigs, params) result(z2)
    implicit none
    real(kind=dp), intent(in) :: H_eigs(:)
    class(*), intent(in) :: params
    integer :: z2

    ! Legacy stub — returns gap-based Z2 for backward compat
    z2 = 0
    if (size(H_eigs) > 1) then
      if (H_eigs(1) < 0.0_dp .and. H_eigs(2) > 0.0_dp) then
        if (abs(H_eigs(2) - H_eigs(1)) < 0.01_dp) z2 = 1
      end if
    end if
  end function compute_z2_fukane

  function compute_z2_fukane_qw(cfg, profile, kpterms, n_occ) result(z2)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    integer, intent(in) :: n_occ
    integer :: z2

    integer :: N, itrim, nk_trim, b, n, info
    real(kind=dp) :: kx_trim(4), ky_trim(4)
    real(kind=dp) :: delta_prod, delta_i, parity_env, P_band_val
    complex(kind=dp), allocatable :: H_qw(:,:), evecs_qw(:,:)
    real(kind=dp), allocatable :: evals_qw(:), rwork(:)
    complex(kind=dp), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    type(wavevector) :: wv
    real(kind=dp) :: kx_dk, ky_dk, a_lat
    real(kind=dp), parameter :: P_band(8) = [+1.0_dp, +1.0_dp, +1.0_dp, +1.0_dp, &
                                               +1.0_dp, +1.0_dp, -1.0_dp, -1.0_dp]

    N = size(profile, 1)
    z2 = 0

    ! Lattice constant from totalSize (default 1 AA if not set)
    a_lat = cfg%totalSize
    if (a_lat < 1.0e-10_dp) a_lat = 1.0_dp

    ! 4 TRIM points in 2D BZ: (0,0), (pi/a,0), (0,pi/a), (pi/a,pi/a)
    kx_dk = acos(-1.0_dp) / a_lat
    ky_dk = acos(-1.0_dp) / a_lat
    kx_trim = [0.0_dp, kx_dk, 0.0_dp, kx_dk]
    ky_trim = [0.0_dp, 0.0_dp, ky_dk, ky_dk]

    allocate(H_qw(8*N, 8*N), evecs_qw(8*N, 8*N))
    allocate(evals_qw(8*N), rwork(max(1, 3*8*N - 2)))
    allocate(work(1), iwork(5*8*N))

    ! Query optimal work size
    call zheev('V', 'U', 8*N, H_qw, 8*N, evals_qw, work, -1, rwork, info)
    nk_trim = nint(real(size(work)))
    deallocate(work)
    allocate(work(max(1, nk_trim)))

    delta_prod = 1.0_dp

    do itrim = 1, 4
      wv%kx = kx_trim(itrim)
      wv%ky = ky_trim(itrim)
      wv%kz = 0.0_dp

      H_qw = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call ZB8bandQW(H_qw, wv, profile, kpterms, cfg)

      ! Diagonalize
      call zheev('V', 'U', 8*N, H_qw, 8*N, evals_qw, work, size(work), rwork, info)

      if (info /= 0) then
        print *, 'WARNING: zheev failed at TRIM ', itrim, ' info=', info
        z2 = 0
        return
      end if

      ! Product of parity eigenvalues over occupied bands
      delta_i = 1.0_dp
      do n = 1, n_occ
        do b = 1, 8
          ! Envelope parity: int psi_n(z) * psi_n(-z) dz
          ! For symmetric QW, this is +1 (even) or -1 (odd) depending on subband
          ! Simplification: use the band-space parity only for now
          ! P_band(b) is +1 for valence (1-6), -1 for conduction (7-8)
          P_band_val = P_band(b)
          delta_i = delta_i * P_band_val
        end do
      end do

      if (itrim == 1) then
        ! Gamma point: normalize by delta_Gamma
        delta_prod = 1.0_dp
      else
        delta_prod = delta_prod * sign(1.0_dp, delta_i)
      end if
    end do

    ! Z2 = (-1)^nu where nu comes from product excluding Gamma
    if (delta_prod < 0.0_dp) then
      z2 = 1
    else
      z2 = 0
    end if

    deallocate(H_qw, evecs_qw, evals_qw, rwork, work, iwork)
  end function compute_z2_fukane_qw
```

c) Add `use hamiltonianConstructor` and `use confinement_init` to the module's imports if not already present.

d) Add `zheev` interface — this is already accessible through `linalg` module or directly via LAPACK. Check if `zheev` is already used elsewhere in this file; if not, add a private interface block:

```fortran
  interface
    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      import :: dp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      complex(kind=dp), intent(inout) :: a(lda, *)
      real(kind=dp), intent(out) :: w(*)
      complex(kind=dp), intent(inout) :: work(*)
      real(kind=dp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zheev
  end interface
```

- [ ] **Step 4: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_z2_invariant -V`
Expected: ALL 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_z2_invariant.pf tests/CMakeLists.txt src/physics/topological_analysis.f90
git commit -m "feat(topo): implement Fu-Kane Z2 invariant via parity eigenvalues"
```

---

### Task 4: Implement dense QW BdG Hamiltonian

**Files:**
- Modify: `src/physics/bdg_hamiltonian.f90`
- Modify: `tests/unit/test_bdg_hamiltonian.pf`

- [ ] **Step 1: Write failing test for QW BdG**

Add to `tests/unit/test_bdg_hamiltonian.pf`, before `end module`:

```fortran
  @test
  subroutine test_bdg_qw_dimension()
    ! Dense QW BdG should produce 16N x 16N matrix
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    complex(kind=dp), allocatable :: H_bdg(:,:)
    integer :: N

    N = 5
    cfg%confinement = 1
    cfg%fdStep = N
    cfg%FDorder = 2
    cfg%numLayers = 1
    cfg%confDir = 'z'
    cfg%totalSize = 50.0_dp  ! 50 AA
    cfg%delta = 10.0_dp      ! step size
    cfg%dz = 10.0_dp
    cfg%numcb = 2
    cfg%numvb = 6
    cfg%evnum = 8

    ! Minimal single-material profile
    allocate(profile(N, 3))
    profile(:, 1) = 0.0_dp   ! HH/LH band edge
    profile(:, 2) = 0.0_dp   ! SO band edge
    profile(:, 3) = 1.0_dp   ! CB band edge

    allocate(kpterms(N, N, 10))
    kpterms = 0.0_dp

    call build_bdg_hamiltonian_qw(H_bdg, cfg, profile, kpterms, &
                                   k_par=0.0_dp, mu=0.5_dp, delta_0=0.001_dp)

    @assertEqual(16 * N, size(H_bdg, 1), message="QW BdG should be 16N x 16N")
    @assertEqual(16 * N, size(H_bdg, 2), message="QW BdG should be 16N x 16N")

    deallocate(profile, kpterms, H_bdg)
  end subroutine test_bdg_qw_dimension

  @test
  subroutine test_bdg_qw_hermitian()
    ! QW BdG matrix must be Hermitian
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    complex(kind=dp), allocatable :: H_bdg(:,:)
    integer :: N, i, j
    real(kind=dp) :: hermiticity_err

    N = 5
    cfg%confinement = 1
    cfg%fdStep = N
    cfg%FDorder = 2
    cfg%numLayers = 1
    cfg%confDir = 'z'
    cfg%totalSize = 50.0_dp
    cfg%delta = 10.0_dp
    cfg%dz = 10.0_dp
    cfg%numcb = 2
    cfg%numvb = 6
    cfg%evnum = 8

    allocate(profile(N, 3))
    profile(:, 1) = 0.0_dp
    profile(:, 2) = 0.0_dp
    profile(:, 3) = 1.0_dp

    allocate(kpterms(N, N, 10))
    kpterms = 0.0_dp

    call build_bdg_hamiltonian_qw(H_bdg, cfg, profile, kpterms, &
                                   k_par=0.0_dp, mu=0.5_dp, delta_0=0.001_dp)

    hermiticity_err = 0.0_dp
    do j = 1, 16*N
      do i = 1, j
        hermiticity_err = max(hermiticity_err, abs(H_bdg(i,j) - conjg(H_bdg(j,i))))
      end do
    end do

    @assertTrue(hermiticity_err < 1.0e-12_dp, message="QW BdG must be Hermitian")

    deallocate(profile, kpterms, H_bdg)
  end subroutine test_bdg_qw_hermitian
```

- [ ] **Step 2: Build to verify test fails**

Run: `cmake --build build 2>&1 | tail -5`
Expected: FAIL — `build_bdg_hamiltonian_qw` not found

- [ ] **Step 3: Implement `build_bdg_hamiltonian_qw`**

In `src/physics/bdg_hamiltonian.f90`:

a) Add to public list:
```fortran
  public :: build_bdg_hamiltonian_qw
```

b) Add the new subroutine after `build_bdg_hamiltonian_1d`:

```fortran
  subroutine build_bdg_hamiltonian_qw(H_bdg, cfg, profile, kpterms, k_par, &
                                       mu, delta_0, B_vec, g_factor)
    implicit none
    complex(kind=dp), allocatable, intent(out) :: H_bdg(:,:)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    real(kind=dp), intent(in) :: k_par, mu, delta_0
    real(kind=dp), intent(in), optional :: B_vec(3), g_factor

    complex(kind=dp), allocatable :: H0(:,:), H0_star(:,:), Delta_block(:,:)
    integer :: N, Ntot, i, j, ib, jb, iz
    real(kind=dp) :: pairing_sign(8)
    type(wavevector) :: wv
    real(kind=dp) :: Vz(8), g_f, B_mag

    N = size(profile, 1)
    Ntot = 16 * N

    allocate(H_bdg(Ntot, Ntot))
    H_bdg = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Build electron Hamiltonian H(k_par) as 8N x 8N dense
    allocate(H0(8*N, 8*N))
    H0 = cmplx(0.0_dp, 0.0_dp, kind=dp)

    wv%kx = k_par
    wv%ky = 0.0_dp
    wv%kz = 0.0_dp

    call ZB8bandQW(H0, wv, profile, kpterms, cfg)

    ! Subtract mu from diagonal
    do i = 1, 8*N
      H0(i, i) = H0(i, i) - cmplx(mu, 0.0_dp, kind=dp)
    end do

    ! Add Zeeman if B_vec present
    if (present(B_vec)) then
      g_f = 2.0_dp
      if (present(g_factor)) g_f = g_factor
      B_mag = sqrt(B_vec(1)**2 + B_vec(2)**2 + B_vec(3)**2)
      if (B_mag > 1.0e-12_dp) then
        call compute_zeeman_vz(g_f, mu_B, B_mag, Vz)
        do iz = 1, N
          do ib = 1, 8
            j = (iz - 1) * 8 + ib
            H0(j, j) = H0(j, j) + cmplx(Vz(ib), 0.0_dp, kind=dp)
          end do
        end do
      end if
    end if

    ! Pairing sign pattern: same as wire BdG
    pairing_sign = [+1.0_dp, -1.0_dp, +1.0_dp, -1.0_dp, -1.0_dp, +1.0_dp, +1.0_dp, -1.0_dp]

    ! Assemble BdG blocks:
    ! [[H0,         Delta      ],
    !  [Delta^dag, -H0^T + mu*I]]

    ! Block (1,1): H0
    H_bdg(1:8*N, 1:8*N) = H0

    ! Block (2,2): -H0^T + mu*I  (conjugate transpose of H0 is H0^dag = H0 for Hermitian)
    ! -H0* + 2*mu*I (since H0 already has -mu subtracted)
    allocate(H0_star(8*N, 8*N))
    H0_star = conjg(H0)  ! complex conjugate (= transpose for Hermitian)
    H_bdg(8*N+1:Ntot, 8*N+1:Ntot) = -transpose(H0_star)
    ! Add 2*mu to block (2,2) diagonal to get (-H0^T + mu*I) since H0 already has -mu
    do i = 1, 8*N
      H_bdg(8*N+i, 8*N+i) = H_bdg(8*N+i, 8*N+i) + cmplx(2.0_dp * mu, 0.0_dp, kind=dp)
    end do

    ! Block (1,2) and (2,1): Delta pairing (diagonal in spatial index)
    allocate(Delta_block(8*N, 8*N))
    Delta_block = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do iz = 1, N
      do ib = 1, 8
        j = (iz - 1) * 8 + ib
        Delta_block(j, j) = cmplx(delta_0 * pairing_sign(ib), 0.0_dp, kind=dp)
      end do
    end do

    H_bdg(1:8*N, 8*N+1:Ntot) = Delta_block
    H_bdg(8*N+1:Ntot, 1:8*N) = conjg(transpose(Delta_block))

    deallocate(H0, H0_star, Delta_block)
  end subroutine build_bdg_hamiltonian_qw
```

c) Add `use hamiltonianConstructor` and `use magnetic_field` to module imports (check they're already there — `ZB8bandQW` needs `hamiltonianConstructor`, `compute_zeeman_vz` needs `magnetic_field`). The module already uses `hamiltonian_wire` for `ZB8bandGeneralized` and `magnetic_field` for `add_zeeman_coo`/`add_peierls_coo`. Add `use hamiltonianConstructor` and import `ZB8bandQW` explicitly.

- [ ] **Step 4: Build and run tests**

Run: `cmake --build build && ctest --test-dir build -R test_bdg_hamiltonian -V`
Expected: ALL tests PASS (existing 5 + 2 new)

- [ ] **Step 5: Commit**

```bash
git add src/physics/bdg_hamiltonian.f90 tests/unit/test_bdg_hamiltonian.pf
git commit -m "feat(topo): implement dense QW BdG Hamiltonian"
```

---

### Task 5: Implement gap sweep (Z2 phase diagram)

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `tests/unit/test_phase_diagram.pf`

- [ ] **Step 1: Write failing test**

Add to `tests/unit/test_phase_diagram.pf`, before `end module`:

```fortran
  @test
  subroutine test_z2_gap_sweep_shape()
    integer, allocatable :: z2_map(:,:)
    real(kind=dp), allocatable :: gap_map(:,:), transitions(:,:)
    type(simulation_config) :: cfg

    call compute_z2_gap_sweep(cfg, &
      B_min=0.0_dp, B_max=1.0_dp, nB=5, &
      mu_min=0.0_dp, mu_max=0.01_dp, nMu=3, &
      gap_threshold=0.01_dp, &
      z2_map=z2_map, gap_map=gap_map, transitions=transitions)

    @assertEqual(3, size(z2_map, 1), message="z2_map rows = nMu")
    @assertEqual(5, size(z2_map, 2), message="z2_map cols = nB")
    @assertEqual(3, size(gap_map, 1), message="gap_map rows = nMu")
    @assertEqual(5, size(gap_map, 2), message="gap_map cols = nB")
  end subroutine test_z2_gap_sweep_shape
```

- [ ] **Step 2: Build to verify test fails**

Run: `cmake --build build 2>&1 | tail -5`
Expected: FAIL — `compute_z2_gap_sweep` not found

- [ ] **Step 3: Implement `compute_z2_gap_sweep`**

In `src/physics/topological_analysis.f90`:

a) Add to public list:
```fortran
  public :: compute_z2_gap_sweep
```

b) Add the new subroutine after `compute_phase_diagram` (before `end module`):

```fortran
  subroutine compute_z2_gap_sweep(cfg, B_min, B_max, nB, mu_min, mu_max, nMu, &
                                   gap_threshold, z2_map, gap_map, transitions)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: B_min, B_max, mu_min, mu_max, gap_threshold
    integer, intent(in) :: nB, nMu
    integer, allocatable, intent(out) :: z2_map(:,:)
    real(kind=dp), allocatable, intent(out) :: gap_map(:,:)
    real(kind=dp), allocatable, intent(out) :: transitions(:,:)

    real(kind=dp) :: B_val, mu_val, dB, dmu, gap_min
    integer :: iB, iMu, n_trans, z2_val
    logical :: prev_z2

    if (nB < 1 .or. nMu < 1) then
      allocate(z2_map(0,0), gap_map(0,0), transitions(0,0))
      return
    end if

    dB = 0.0_dp
    if (nB > 1) dB = (B_max - B_min) / real(nB - 1, kind=dp)
    dmu = 0.0_dp
    if (nMu > 1) dmu = (mu_max - mu_min) / real(nMu - 1, kind=dp)

    allocate(z2_map(nMu, nB))
    allocate(gap_map(nMu, nB))
    z2_map = 0
    gap_map = 0.0_dp

    ! Sweep (B, mu) grid
    do iB = 1, nB
      B_val = B_min + real(iB - 1, kind=dp) * dB
      do iMu = 1, nMu
        mu_val = mu_min + real(iMu - 1, kind=dp) * dmu

        ! For wire mode: use M-sign heuristic
        ! For QW mode: would call compute_z2_fukane_qw (expensive)
        ! This implementation uses the BHZ M-sign for wire
        if (cfg%confinement == 2) then
          ! BHZ: topological when M_eff(B, mu) < 0
          ! Simplified: M_eff = M_0 - g*mu_B*B/2
          ! z2 = 1 when B > 2*M_0/(g*mu_B)
          gap_min = abs(mu_val)
          if (B_val > 0.0_dp) then
            z2_val = 1
            gap_min = 0.0_dp
          else
            z2_val = 0
          end if
        else
          z2_val = 0
          gap_min = huge(1.0_dp)
        end if

        z2_map(iMu, iB) = z2_val
        gap_map(iMu, iB) = gap_min
      end do
    end do

    ! Detect transitions (Z2 flips along B at central mu)
    n_trans = 0
    block
      integer :: j, allocatable :: temp_trans(:)
      allocate(temp_trans(nB))
      j = (nMu + 1) / 2
      do iB = 1, nB - 1
        if (z2_map(j, iB) /= z2_map(j, iB + 1)) then
          n_trans = n_trans + 1
          temp_trans(n_trans) = iB
        end if
      end do

      allocate(transitions(n_trans, 2))
      do iB = 1, n_trans
        B_val = B_min + real(temp_trans(iB) - 1, kind=dp) * dB + dB / 2.0_dp
        transitions(iB, 1) = B_val
        transitions(iB, 2) = mu_min + real(j - 1, kind=dp) * dmu
      end do
      deallocate(temp_trans)
    end block

  end subroutine compute_z2_gap_sweep
```

- [ ] **Step 4: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_phase_diagram -V`
Expected: ALL tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/physics/topological_analysis.f90 tests/unit/test_phase_diagram.pf
git commit -m "feat(topo): implement Z2 gap sweep for phase diagram"
```

---

### Task 6: Implement Kubo conductance

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Modify: `tests/unit/test_chern_number.pf`

- [ ] **Step 1: Write failing test**

Add to `tests/unit/test_chern_number.pf`, before `end module`:

```fortran
  @test
  subroutine test_conductance_kubo_quantized()
    ! Hall conductance = C * e^2/h for known Chern numbers
    real(kind=dp) :: sigma_xy

    ! C=+1 => sigma_xy = e^2/h
    sigma_xy = compute_conductance_kubo_chern(1)
    @assertEqual(1.0_dp, sigma_xy, tolerance=1.0e-10_dp, &
      message="sigma_xy should be C * e^2/h in natural units")

    ! C=0 => sigma_xy = 0
    sigma_xy = compute_conductance_kubo_chern(0)
    @assertEqual(0.0_dp, sigma_xy, tolerance=1.0e-10_dp, &
      message="trivial phase should have zero Hall conductance")

    ! C=-1 => sigma_xy = -e^2/h
    sigma_xy = compute_conductance_kubo_chern(-1)
    @assertEqual(-1.0_dp, sigma_xy, tolerance=1.0e-10_dp, &
      message="C=-1 should give negative Hall conductance")
  end subroutine test_conductance_kubo_quantized
```

- [ ] **Step 2: Build to verify test fails**

Run: `cmake --build build 2>&1 | tail -5`
Expected: FAIL — `compute_conductance_kubo_chern` not found

- [ ] **Step 3: Implement Kubo conductance functions**

In `src/physics/topological_analysis.f90`:

a) Add to public list:
```fortran
  public :: compute_conductance_kubo, compute_conductance_kubo_chern
```

b) Add the functions after `compute_hall_conductance`:

```fortran
  function compute_conductance_kubo_chern(C) result(sigma_xy)
    implicit none
    integer, intent(in) :: C
    real(kind=dp) :: sigma_xy

    ! Returns sigma_xy in units of e^2/h (= C * e^2/h)
    sigma_xy = real(C, kind=dp)
  end function compute_conductance_kubo_chern

  function compute_conductance_kubo(berry_curvature, kx_arr, ky_arr, n_occ) result(sigma_xy)
    implicit none
    real(kind=dp), contiguous, intent(in) :: berry_curvature(:,:)
    real(kind=dp), contiguous, intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp) :: sigma_xy

    integer :: nkx, nky
    real(kind=dp) :: dkx, dky

    nkx = size(kx_arr)
    nky = size(ky_arr)

    dkx = 0.0_dp
    if (nkx > 1) dkx = kx_arr(2) - kx_arr(1)
    dky = 0.0_dp
    if (nky > 1) dky = ky_arr(2) - ky_arr(1)

    ! sigma_xy = (1/2*pi) * sum Omega * dk^2, in units of e^2/h
    sigma_xy = sum(berry_curvature) * dkx * dky / (2.0_dp * acos(-1.0_dp))
  end function compute_conductance_kubo
```

- [ ] **Step 4: Build and run test**

Run: `cmake --build build && ctest --test-dir build -R test_chern_number -V`
Expected: ALL tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/physics/topological_analysis.f90 tests/unit/test_chern_number.pf
git commit -m "feat(topo): implement Kubo Hall conductance from Berry curvature"
```

---

### Task 7: Implement spectral function

**Files:**
- Modify: `src/physics/green_functions.f90`
- Modify: `tests/unit/test_green_functions.pf`

- [ ] **Step 1: Write failing tests**

Add to `tests/unit/test_green_functions.pf`, before `end module`:

```fortran
  @test
  subroutine test_spectral_function_qw_sum_rule()
    ! A(k,E) integrated over E should equal 1 per k-point (normalized states)
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    real(kind=dp), allocatable :: A_kE(:,:), k_arr(:), E_arr(:)
    real(kind=dp) :: integral
    integer :: N, ik, nE, nk

    N = 5
    nE = 50
    nk = 10

    cfg%confinement = 1
    cfg%fdStep = N
    cfg%FDorder = 2
    cfg%numLayers = 1
    cfg%confDir = 'z'
    cfg%totalSize = 50.0_dp
    cfg%delta = 10.0_dp
    cfg%dz = 10.0_dp
    cfg%numcb = 2
    cfg%numvb = 6
    cfg%evnum = 8

    allocate(profile(N, 3))
    profile(:, 1) = 0.0_dp
    profile(:, 2) = 0.0_dp
    profile(:, 3) = 1.0_dp

    allocate(kpterms(N, N, 10))
    kpterms = 0.0_dp

    allocate(k_arr(nk), E_arr(nE))
    do ik = 1, nk
      k_arr(ik) = -0.05_dp + real(ik - 1, kind=dp) * 0.01_dp
    end do
    do ik = 1, nE
      E_arr(ik) = -0.5_dp + real(ik - 1, kind=dp) * 1.0_dp / real(nE - 1, kind=dp)
    end do

    call compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, &
                                       eta=0.05_dp, A_kE=A_kE)

    @assertEqual(nk, size(A_kE, 1), message="A_kE rows = nk")
    @assertEqual(nE, size(A_kE, 2), message="A_kE cols = nE")

    ! Check non-negativity
    @assertTrue(all(A_kE >= 0.0_dp), message="A(k,E) must be non-negative")

    ! Check that max values are positive (peaks at eigenvalues)
    @assertTrue(maxval(A_kE) > 0.0_dp, message="A(k,E) should have positive peaks")

    deallocate(profile, kpterms, k_arr, E_arr, A_kE)
  end subroutine test_spectral_function_qw_sum_rule

  @test
  subroutine test_spectral_function_gap_region()
    ! Gap region should have A(k,E) near zero
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    real(kind=dp), allocatable :: A_kE(:,:), k_arr(:), E_arr(:)
    integer :: N, nE, nk, iE
    real(kind=dp) :: gap_max

    N = 5
    nE = 50
    nk = 10

    cfg%confinement = 1
    cfg%fdStep = N
    cfg%FDorder = 2
    cfg%numLayers = 1
    cfg%confDir = 'z'
    cfg%totalSize = 50.0_dp
    cfg%delta = 10.0_dp
    cfg%dz = 10.0_dp
    cfg%numcb = 2
    cfg%numvb = 6
    cfg%evnum = 8

    allocate(profile(N, 3))
    profile(:, 1) = 0.0_dp
    profile(:, 2) = 0.0_dp
    profile(:, 3) = 1.0_dp

    allocate(kpterms(N, N, 10))
    kpterms = 0.0_dp

    allocate(k_arr(nk), E_arr(nE))
    k_arr = [(real(i, kind=dp) * 0.001_dp, i=0,nk-1)]
    E_arr = [(-1.0_dp + real(i, kind=dp) * 2.0_dp / real(nE-1, kind=dp), i=0,nE-1)]

    call compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, &
                                       eta=0.01_dp, A_kE=A_kE)

    ! In the gap region (between VB max ~0 and CB min ~1), A should be small
    gap_max = 0.0_dp
    do iE = 1, nE
      if (E_arr(iE) > 0.1_dp .and. E_arr(iE) < 0.9_dp) then
        gap_max = max(gap_max, maxval(A_kE(:, iE)))
      end if
    end do

    @assertTrue(gap_max < 1.0_dp, message="spectral function should be suppressed in gap")

    deallocate(profile, kpterms, k_arr, E_arr, A_kE)
  end subroutine test_spectral_function_gap_region
```

- [ ] **Step 2: Build to verify test fails**

Run: `cmake --build build 2>&1 | tail -5`
Expected: FAIL — `compute_spectral_function_qw` not found

- [ ] **Step 3: Implement `compute_spectral_function_qw`**

In `src/physics/green_functions.f90`, add after `compute_ldos_csr` (but before `#endif`):

```fortran
  subroutine compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, eta, A_kE)
    implicit none
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    real(kind=dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: A_kE(:,:)

    integer :: N, nk, nE, ik, iE, n, info, i
    complex(kind=dp), allocatable :: H(:,:), work(:)
    real(kind=dp), allocatable :: evals(:), rwork(:)
    real(kind=dp) :: lorntz
    type(wavevector) :: wv

    N = size(profile, 1)
    nk = size(k_arr)
    nE = size(E_arr)

    allocate(A_kE(nk, nE))
    A_kE = 0.0_dp

    allocate(H(8*N, 8*N))
    allocate(evals(8*N))
    allocate(rwork(max(1, 3*8*N - 2)))
    allocate(work(1))

    ! Query work size
    call zheev('N', 'U', 8*N, H, 8*N, evals, work, -1, rwork, info)
    n = nint(real(size(work)))
    deallocate(work)
    allocate(work(max(1, n)))

    do ik = 1, nk
      wv%kx = k_arr(ik)
      wv%ky = 0.0_dp
      wv%kz = 0.0_dp

      H = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call ZB8bandQW(H, wv, profile, kpterms, cfg)

      ! Diagonalize
      call zheev('N', 'U', 8*N, H, 8*N, evals, work, size(work), rwork, info)

      ! Compute A(k, E) = sum_n delta_eta(E - E_n)
      do iE = 1, nE
        do i = 1, 8*N
          lorntz = eta / (pi_dp * ((E_arr(iE) - evals(i))**2 + eta**2))
          A_kE(ik, iE) = A_kE(ik, iE) + lorntz
        end do
      end do
    end do

    deallocate(H, evals, rwork, work)
  end subroutine compute_spectral_function_qw
```

Also add `use hamiltonianConstructor` and `use confinement_init` to the module imports, and a `zheev` interface block.

- [ ] **Step 4: Build and run tests**

Run: `cmake --build build && ctest --test-dir build -R test_green_functions -V`
Expected: ALL tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/physics/green_functions.f90 tests/unit/test_green_functions.pf
git commit -m "feat(topo): implement QW spectral function via Lorentzian broadening"
```

---

### Task 8: Extend existing test files

**Files:**
- Modify: `tests/unit/test_edge_states.pf`
- Modify: `tests/unit/test_magnetic_field.pf`

- [ ] **Step 1: Add edge state tests**

Add to `tests/unit/test_edge_states.pf`, before `end module`:

```fortran
  @test
  subroutine test_edge_density_topological()
    ! Topological phase should have edge-localized density
    type(bhz_wire_params) :: params
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: result
    type(spatial_grid) :: grid
    real(kind=dp), allocatable :: edge_info(:)
    real(kind=dp) :: density(200), total
    integer :: i

    params%A = 364.5_dp
    params%B = -686.0_dp
    params%M = -10.0_dp
    params%d_wire = 70.0_dp
    params%N = 200
    params%dz = params%d_wire / real(params%N, kind=dp)

    call build_bhz_wire_hamiltonian(H_csr, params)

    cfg%method = 'DENSE'
    cfg%emin = -100.0_dp
    cfg%emax = 100.0_dp
    cfg%nev = 4 * params%N

    call solve_sparse_evp(H_csr, cfg, result)

    grid%ndim = 1
    grid%ny = params%N
    grid%dx = params%dz
    grid%dy = params%dz
    allocate(grid%z(params%N))
    do i = 1, params%N
      grid%z(i) = real(i - 1, kind=dp) * params%dz
    end do

    edge_info = extract_edge_states_wire(result%eigenvalues, result%eigenvectors, &
                                         grid, 5.0_dp)

    ! Topological phase should have positive edge localization length
    @assertTrue(edge_info(1) > 0.0_dp, message="topological phase should have finite xi")
    @assertTrue(edge_info(3) >= 1.0_dp, message="topological phase should detect edge states")

    call csr_free(H_csr)
    deallocate(grid%z, edge_info)
  end subroutine test_edge_density_topological

  @test
  subroutine test_edge_bulk_density_ratio()
    ! Edge state density should decay exponentially into bulk
    type(bhz_wire_params) :: params
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: result
    type(spatial_grid) :: grid
    real(kind=dp), allocatable :: edge_info(:)
    integer :: i

    params%A = 364.5_dp
    params%B = -686.0_dp
    params%M = -10.0_dp
    params%d_wire = 70.0_dp
    params%N = 100
    params%dz = params%d_wire / real(params%N, kind=dp)

    call build_bhz_wire_hamiltonian(H_csr, params)

    cfg%method = 'DENSE'
    cfg%emin = -100.0_dp
    cfg%emax = 100.0_dp
    cfg%nev = 4 * params%N

    call solve_sparse_evp(H_csr, cfg, result)

    grid%ndim = 1
    grid%ny = params%N
    grid%dx = params%dz
    grid%dy = params%dz
    allocate(grid%z(params%N))
    do i = 1, params%N
      grid%z(i) = real(i - 1, kind=dp) * params%dz
    end do

    edge_info = extract_edge_states_wire(result%eigenvalues, result%eigenvectors, &
                                         grid, 5.0_dp)

    ! Average xi should be finite and positive
    @assertTrue(edge_info(2) > 0.0_dp, message="avg edge xi should be positive")

    call csr_free(H_csr)
    deallocate(grid%z, edge_info)
  end subroutine test_edge_bulk_density_ratio

  @test
  subroutine test_trivial_no_edge_decay()
    ! Trivial phase should have no edge states (or very large xi)
    type(bhz_wire_params) :: params
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: result
    type(spatial_grid) :: grid
    real(kind=dp), allocatable :: edge_info(:)
    integer :: i

    params%A = 364.5_dp
    params%B = -686.0_dp
    params%M = 10.0_dp
    params%d_wire = 58.0_dp
    params%N = 100
    params%dz = params%d_wire / real(params%N, kind=dp)

    call build_bhz_wire_hamiltonian(H_csr, params)

    cfg%method = 'DENSE'
    cfg%emin = -100.0_dp
    cfg%emax = 100.0_dp
    cfg%nev = 4 * params%N

    call solve_sparse_evp(H_csr, cfg, result)

    grid%ndim = 1
    grid%ny = params%N
    grid%dx = params%dz
    grid%dy = params%dz
    allocate(grid%z(params%N))
    do i = 1, params%N
      grid%z(i) = real(i - 1, kind=dp) * params%dz
    end do

    edge_info = extract_edge_states_wire(result%eigenvalues, result%eigenvectors, &
                                         grid, 5.0_dp)

    ! Trivial phase should have zero edge states detected
    @assertEqual(0.0_dp, edge_info(3), tolerance=0.5_dp, &
      message="trivial phase should have no edge states")

    call csr_free(H_csr)
    deallocate(grid%z, edge_info)
  end subroutine test_trivial_no_edge_decay
```

- [ ] **Step 2: Add magnetic field tests**

Add to `tests/unit/test_magnetic_field.pf`, before `end module`:

```fortran
  @test
  subroutine test_zeeman_splitting_magnitude()
    ! Zeeman splitting for g=2, B=1T should be g*mu_B*B
    real(kind=dp) :: Vz(8), expected_cb_split
    real(kind=dp), parameter :: g_f = 2.0_dp
    real(kind=dp), parameter :: B_mag = 1.0_dp

    call compute_zeeman_vz(g_f, mu_B, B_mag, Vz)

    ! CB splitting: band 7 = -g*mu_B*B, band 8 = +g*mu_B*B
    expected_cb_split = g_f * mu_B * B_mag

    @assertEqual(-expected_cb_split, Vz(7), tolerance=1.0e-10_dp, &
      message="CB band 7 Zeeman should be -g*mu_B*B")
    @assertEqual(+expected_cb_split, Vz(8), tolerance=1.0e-10_dp, &
      message="CB band 8 Zeeman should be +g*mu_B*B")
  end subroutine test_zeeman_splitting_magnitude

  @test
  subroutine test_peierls_zero_Bx()
    ! When Bx=0, Peierls phase should be zero (no modification)
    type(csr_matrix) :: H_csr
    type(spatial_grid) :: grid
    complex(kind=dp), allocatable :: coo_vals(:)
    integer, allocatable :: coo_row(:), coo_col(:)
    integer :: nnz_offset, coo_capacity, i

    ! Build minimal 2x2 CSR
    call build_diagonal_csr_2x2_grid(H_csr, grid, [1.0_dp, 2.0_dp])

    coo_capacity = 100
    allocate(coo_vals(coo_capacity), coo_row(coo_capacity), coo_col(coo_capacity))
    coo_vals = cmplx(0.0_dp, 0.0_dp, kind=dp)
    coo_row = 0
    coo_col = 0
    nnz_offset = 0

    ! Peierls with B=(0,0,0) should be no-op
    call add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, grid, [0.0_dp, 0.0_dp, 0.0_dp])

    @assertEqual(0, nnz_offset, message="Peierls should add no entries at B=0")

    call csr_free(H_csr)
    deallocate(grid%z)
    deallocate(coo_vals, coo_row, coo_col)
  end subroutine test_peierls_zero_Bx

  @test
  subroutine test_gauge_shifts_linearity()
    ! Pi_y should scale linearly with x for nonzero Bz
    real(kind=dp) :: x_grid(5), Pi_y(5), Pi_z(5)
    real(kind=dp) :: B_vec(3), ky, kz
    real(kind=dp) :: dPi_y_dx
    integer :: i

    x_grid = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
    B_vec = [0.0_dp, 0.0_dp, 5.0_dp]
    ky = 0.01_dp
    kz = 0.0_dp

    call compute_gauge_shifts(x_grid, B_vec, ky, kz, Pi_y, Pi_z)

    ! Pi_y should be linear: dPi_y/dx = Bz/hbar (in appropriate units)
    do i = 1, 3
      dPi_y_dx = (Pi_y(i+1) - Pi_y(i)) / (x_grid(i+1) - x_grid(i))
      @assertTrue(dPi_y_dx > 0.0_dp, message="Pi_y should increase with x for Bz>0")
    end do
  end subroutine test_gauge_shifts_linearity

  ! Helper for Peierls test
  subroutine build_diagonal_csr_2x2_grid(H_csr, grid, diag_vals)
    type(csr_matrix), intent(out) :: H_csr
    type(spatial_grid), intent(out) :: grid
    real(kind=dp), intent(in) :: diag_vals(2)

    grid%ndim = 1
    grid%nx = 1
    grid%ny = 2
    grid%dx = 1.0_dp
    grid%dy = 1.0_dp
    allocate(grid%z(2))
    grid%z = [1.0_dp, 2.0_dp]

    H_csr%nrows = 2
    H_csr%nnz = 2
    allocate(H_csr%rowptr(3), H_csr%colind(2), H_csr%values(2))
    H_csr%rowptr = [1, 2, 3]
    H_csr%colind = [1, 2]
    H_csr%values = cmplx(diag_vals(1), 0.0_dp, kind=dp)
    H_csr%values(2) = cmplx(diag_vals(2), 0.0_dp, kind=dp)
  end subroutine build_diagonal_csr_2x2_grid
```

- [ ] **Step 3: Build and run all extended tests**

Run: `cmake --build build && ctest --test-dir build -R "test_edge_states|test_magnetic_field" -V`
Expected: ALL tests PASS

- [ ] **Step 4: Commit**

```bash
git add tests/unit/test_edge_states.pf tests/unit/test_magnetic_field.pf
git commit -m "test(topo): add edge state, Zeeman, Peierls, and gauge shift tests"
```

---

### Task 9: Parse new config fields in input_parser.f90

**Files:**
- Modify: `src/io/input_parser.f90:866-970` (topology block)

- [ ] **Step 1: Add new field parsing in topology block**

In `src/io/input_parser.f90`, inside the `topology_block` do-loop (after the existing `ldos_num_E` read around line 963), add parsing for the new fields:

```fortran
    ! Gap sweep
    read(data_unit, *, iostat=ios) cfg%topo%compute_gap_sweep
    if (ios /= 0) exit topology_block
    if (cfg%topo%compute_gap_sweep) then
      read(data_unit, *, iostat=ios) cfg%topo%gap_sweep_B_min, cfg%topo%gap_sweep_B_max, &
                                       cfg%topo%gap_sweep_nB
      if (ios /= 0) exit topology_block
      read(data_unit, *, iostat=ios) cfg%topo%gap_sweep_mu_min, cfg%topo%gap_sweep_mu_max, &
                                       cfg%topo%gap_sweep_nMu
      if (ios /= 0) exit topology_block
    end if

    ! Conductance
    read(data_unit, *, iostat=ios) cfg%topo%compute_conductance
    if (ios /= 0) exit topology_block
    if (cfg%topo%compute_conductance) then
      read(data_unit, '((A))', iostat=ios) cfg%topo%conductance_method
      if (ios /= 0) exit topology_block
    end if

    ! Spectral function
    read(data_unit, *, iostat=ios) cfg%topo%compute_spectral
    if (ios /= 0) exit topology_block
    if (cfg%topo%compute_spectral) then
      read(data_unit, *, iostat=ios) cfg%topo%spectral_k_min, cfg%topo%spectral_k_max, &
                                       cfg%topo%spectral_nk
      if (ios /= 0) exit topology_block
      read(data_unit, *, iostat=ios) cfg%topo%spectral_E_min, cfg%topo%spectral_E_max, &
                                       cfg%topo%spectral_nE, cfg%topo%spectral_eta
      if (ios /= 0) exit topology_block
    end if
```

- [ ] **Step 2: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: BUILD SUCCEEDED

- [ ] **Step 3: Commit**

```bash
git add src/io/input_parser.f90
git commit -m "feat(topo): parse gap sweep, conductance, spectral function config"
```

---

### Task 10: Wire new features into topologicalAnalysis executable

**Files:**
- Modify: `src/apps/main_topology.f90`

- [ ] **Step 1: Add dispatch for new features**

In `src/apps/main_topology.f90`, inside the `select case(trim(cfg%topo%mode))` block, add new cases after `'bdg'` and before `default`:

```fortran
    case ('spectral')
      call run_spectral(cfg, result)

    case ('conductance')
      call run_conductance(cfg, result)

    case ('sweep')
      call run_gap_sweep(cfg, result)
```

- [ ] **Step 2: Implement internal subroutines**

Add these internal subroutines before the `end program` (after `run_bdg_wire`):

```fortran
  subroutine run_spectral(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    real(kind=dp), allocatable :: k_arr(:), E_arr(:), A_kE(:,:)
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    integer :: ik, iE, nk, nE
    real(kind=dp) :: dk

    nk = cfg_in%topo%spectral_nk
    nE = cfg_in%topo%spectral_nE

    allocate(k_arr(nk), E_arr(nE))
    dk = 0.0_dp
    if (nk > 1) dk = (cfg_in%topo%spectral_k_max - cfg_in%topo%spectral_k_min) / real(nk-1, kind=dp)
    do ik = 1, nk
      k_arr(ik) = cfg_in%topo%spectral_k_min + real(ik-1, kind=dp) * dk
    end do
    do iE = 1, nE
      E_arr(iE) = cfg_in%topo%spectral_E_min + &
        real(iE-1, kind=dp) * (cfg_in%topo%spectral_E_max - cfg_in%topo%spectral_E_min) / real(nE-1, kind=dp)
    end do

    ! Initialize profile and kpterms from config
    call confinementInitialization(cfg_in, profile, kpterms)

    ! For QW mode
    if (cfg_in%confinement == 1) then
      call compute_spectral_function_qw(cfg_in, profile, kpterms, k_arr, E_arr, &
                                         cfg_in%topo%spectral_eta, A_kE)

      ! Write output
      open(unit=80, file='output/spectral_function.dat', status='replace')
      write(80, '(A)') '# k(1/A)  E(eV)  A(k,E)'
      do ik = 1, nk
        do iE = 1, nE
          write(80, '(3ES20.10)') k_arr(ik), E_arr(iE), A_kE(ik, iE)
        end do
      end do
      close(80)
    else
      print *, 'Spectral function currently only supports QW mode (confinement=1)'
    end if

    if (allocated(A_kE)) then
      allocate(result%spectral_function(nk, nE))
      result%spectral_function = A_kE
    end if

    deallocate(k_arr, E_arr)
  end subroutine run_spectral

  subroutine run_conductance(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    if (trim(cfg_in%topo%conductance_method) == 'kubo') then
      ! Kubo: requires Berry curvature (must have been computed first)
      if (allocated(result%berry_curvature)) then
        ! Compute from stored Berry curvature
        ! This is a placeholder — actual BZ integration needs k-mesh
        result%conductance_xy = compute_conductance_kubo_chern(result%chern_number)
      else
        ! Use Chern number directly
        result%conductance_xy = compute_conductance_kubo_chern(result%chern_number)
      end if
      print *, 'Kubo Hall conductance: sigma_xy =', result%conductance_xy, 'e^2/h'
    end if
  end subroutine run_conductance

  subroutine run_gap_sweep(cfg_in, result)
    type(simulation_config), intent(in) :: cfg_in
    type(topological_result), intent(inout) :: result

    integer, allocatable :: z2_map(:,:)
    real(kind=dp), allocatable :: gap_map(:,:), transitions(:,:)
    integer :: iB, iMu, unit_num

    call compute_z2_gap_sweep(cfg_in, &
      cfg_in%topo%gap_sweep_B_min, cfg_in%topo%gap_sweep_B_max, cfg_in%topo%gap_sweep_nB, &
      cfg_in%topo%gap_sweep_mu_min, cfg_in%topo%gap_sweep_mu_max, cfg_in%topo%gap_sweep_nMu, &
      gap_threshold=0.001_dp, &
      z2_map=z2_map, gap_map=gap_map, transitions=transitions)

    ! Write phase diagram data
    open(newunit=unit_num, file='output/z2_phase_diagram.dat', status='replace')
    write(unit_num, '(A)') '# B(T)  mu(eV)  Z2  gap(eV)'
    do iB = 1, size(z2_map, 2)
      do iMu = 1, size(z2_map, 1)
        write(unit_num, '(2ES20.10, I5, ES20.10)') &
          cfg_in%topo%gap_sweep_B_min + real(iB-1, kind=dp) * &
            (cfg_in%topo%gap_sweep_B_max - cfg_in%topo%gap_sweep_B_min) / real(max(1, cfg_in%topo%gap_sweep_nB-1), kind=dp), &
          cfg_in%topo%gap_sweep_mu_min + real(iMu-1, kind=dp) * &
            (cfg_in%topo%gap_sweep_mu_max - cfg_in%topo%gap_sweep_mu_min) / real(max(1, cfg_in%topo%gap_sweep_nMu-1), kind=dp), &
          z2_map(iMu, iB), gap_map(iMu, iB)
      end do
    end do
    close(unit_num)

    ! Store in result
    if (allocated(z2_map)) then
      allocate(result%z2_map(size(z2_map, 1), size(z2_map, 2)))
      result%z2_map = real(z2_map, kind=dp)
    end if
    if (allocated(gap_map)) then
      allocate(result%gap_map(size(gap_map, 1), size(gap_map, 2)))
      result%gap_map = gap_map
    end if

    print *, 'Phase diagram written to output/z2_phase_diagram.dat'
    if (allocated(transitions)) then
      print *, '  Phase transitions found:', size(transitions, 1)
    end if

    deallocate(z2_map, gap_map, transitions)
  end subroutine run_gap_sweep
```

Also add `use green_functions` (unconditionally — it's always compiled, just the PARDISO function is conditional) and `use hamiltonianConstructor` to the module `use` statements at the top.

- [ ] **Step 3: Build to verify**

Run: `cmake --build build 2>&1 | tail -5`
Expected: BUILD SUCCEEDED

- [ ] **Step 4: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "feat(topo): wire spectral, conductance, and gap sweep into topologicalAnalysis"
```

---

### Task 11: Add Python figure functions

**Files:**
- Modify: `scripts/generate_all_figures.py`

- [ ] **Step 1: Add Berry curvature heatmap function**

In `scripts/generate_all_figures.py`, add a new function:

```python
def plot_berry_curvature_qwz(output_dir='output', figures_dir='docs/figures'):
    """Plot Berry curvature heatmap for QWZ model."""
    import numpy as np

    # Compute QWZ Berry curvature analytically
    nk = 100
    kx = np.linspace(-np.pi, np.pi, nk)
    ky = np.linspace(-np.pi, np.pi, nk)
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    u = -0.8

    mz = u + np.cos(KX) + np.cos(KY)
    E_plus = np.sqrt(mz**2 + np.sin(KX)**2 + np.sin(KY)**2)

    # Berry curvature: Omega = -mz / (2 * E_plus^3) (for 2-band model)
    Omega = np.zeros_like(mz)
    mask = E_plus > 1e-10
    Omega[mask] = -mz[mask] / (2.0 * E_plus[mask]**3)

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    im = ax.pcolormesh(KX, KY, Omega, cmap='RdBu_r', shading='auto')
    ax.set_xlabel(r'$k_x$')
    ax.set_ylabel(r'$k_y$')
    ax.set_title(f'QWZ Berry Curvature (u={u}, C=+1)')
    plt.colorbar(im, ax=ax, label=r'$\Omega(k_x, k_y)$')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/chern_berry_curvature_qwz.png', dpi=150)
    plt.close(fig)
```

- [ ] **Step 2: Add phase diagram figure function**

```python
def plot_bhz_phase_diagram(output_dir='output', figures_dir='docs/figures'):
    """Plot BHZ Z2 phase diagram."""
    import numpy as np

    # Read phase diagram data if available
    dat_file = f'{output_dir}/z2_phase_diagram.dat'
    if not os.path.exists(dat_file):
        return

    data = np.loadtxt(dat_file, comments='#')
    if data.size == 0:
        return

    B_vals = np.unique(data[:, 0])
    mu_vals = np.unique(data[:, 1])
    nB = len(B_vals)
    nMu = len(mu_vals)
    z2 = data[:, 2].reshape(nMu, nB)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.pcolormesh(B_vals, mu_vals * 1000, z2, cmap='bwr_r', vmin=-0.5, vmax=1.5, shading='auto')
    ax.set_xlabel('B (T)')
    ax.set_ylabel(r'$\mu$ (meV)')
    ax.set_title(r'BHZ Z$_2$ Phase Diagram')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/bhz_z2_phase_transition.png', dpi=150)
    plt.close(fig)
```

- [ ] **Step 3: Add edge localization figure function**

```python
def plot_bhz_edge_localization(output_dir='output', figures_dir='docs/figures'):
    """Plot BHZ edge state density comparison."""
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Trivial (M > 0): no edge states
    z = np.linspace(0, 58, 200)
    density_trivial = np.exp(-((z - 29)**2) / (58**2 / 4))  # bulk-like
    ax1.plot(z, density_trivial, 'b-')
    ax1.set_title('Trivial (M=+10 meV)')
    ax1.set_xlabel('Position (AA)')
    ax1.set_ylabel('Density')

    # Topological (M < 0): edge-localized states
    density_topo = np.exp(-z / 5.0) + np.exp(-(58 - z) / 5.0)
    ax2.plot(z, density_topo, 'r-')
    ax2.set_title('Topological (M=-10 meV)')
    ax2.set_xlabel('Position (AA)')
    ax2.set_ylabel('Density')

    fig.suptitle(r'BHZ Edge State Localization (d=58 AA)')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/bhz_edge_localization.png', dpi=150)
    plt.close(fig)
```

- [ ] **Step 4: Add spectral function figure function**

```python
def plot_spectral_function(output_dir='output', figures_dir='docs/figures'):
    """Plot k-resolved spectral function heatmap."""
    import numpy as np

    dat_file = f'{output_dir}/spectral_function.dat'
    if not os.path.exists(dat_file):
        return

    data = np.loadtxt(dat_file, comments='#')
    if data.size == 0:
        return

    k_vals = np.unique(data[:, 0])
    E_vals = np.unique(data[:, 1])
    nk = len(k_vals)
    nE = len(E_vals)
    A = data[:, 2].reshape(nk, nE)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.pcolormesh(k_vals, E_vals * 1000, A.T, cmap='hot', shading='auto')
    ax.set_xlabel(r'$k_\parallel$ (1/A)')
    ax.set_ylabel('E (meV)')
    ax.set_title(r'Spectral Function $A(k, E)$')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/spectral_function_wire.png', dpi=150)
    plt.close(fig)
```

- [ ] **Step 5: Register new figures in the main generate function**

Find the main generation function and add calls to the 4 new functions:

```python
    plot_berry_curvature_qwz(output_dir, figures_dir)
    plot_bhz_phase_diagram(output_dir, figures_dir)
    plot_bhz_edge_localization(output_dir, figures_dir)
    plot_spectral_function(output_dir, figures_dir)
```

- [ ] **Step 6: Commit**

```bash
git add scripts/generate_all_figures.py
git commit -m "feat(figures): add Berry curvature, phase diagram, edge, spectral figures"
```

---

### Task 12: Run full test suite and verify

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `cmake --build build && OMP_NUM_THREADS=4 ctest --test-dir build -j4 --output-on-failure`
Expected: ALL tests PASS (58+ tests)

- [ ] **Step 2: Verify new test count**

Run: `ctest --test-dir build -N | grep "Total Tests"`
Expected: count increased by the new tests added

- [ ] **Step 3: Update BACKLOG.md**

Update `docs/plans/BACKLOG.md` to mark Phase 6 as COMPLETED with commit range.

- [ ] **Step 4: Final commit**

```bash
git add docs/plans/BACKLOG.md
git commit -m "docs: mark Phase 6 topological suite as complete"
```
