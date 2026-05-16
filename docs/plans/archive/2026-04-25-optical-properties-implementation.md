# Optical Properties Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the incorrect diagonal-P optical matrix elements with commutator-based velocity operators, add a standalone `opticalProperties` executable supporting all geometries (bulk/QW/wire) and all quantities (absorption/spontaneous emission/gain/ISBT) with optional spin resolution.

**Architecture:** New standalone program `opticalProperties` reads `input.cfg`, builds velocity matrices via the commutator `[r,H]` (shared with g-factor code), sweeps k_par, accumulates optical spectra, and writes output files. The existing `optical_spectra.f90` is rewritten to use commutator velocity instead of `pMatrixEleCalc`.

**Tech Stack:** Fortran 90, MKL SpBLAS (CSR sparse matrices), LAPACK/FEAST eigensolvers, pFUnit for testing.

**Dependency:** This plan assumes `build_velocity_matrices()` already exists in `hamiltonianConstructor.f90` (lines 1798-1834, implemented for the wire g-factor fix). It returns `vel_x` and `vel_y` for 2D grids. We extend it to also support 1D (QW) grids and add a `vel_z` variant.

---

## Phase 1: Infrastructure

### Task 1: Extend `build_velocity_matrices` for QW (1D confinement)

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:1798-1834`
- Test: `tests/unit/test_hamiltonian.pf`

**Context:** The existing `build_velocity_matrices(H_csr, grid, vel_x, vel_y)` only handles `grid%ndim == 2` (wire). For QW (`ndim == 1`), the grid is 1D along z with `grid%z(:)` coordinates. We need a new overload that returns `vel(3)` with all three velocity components.

**Step 1: Write the failing test**

Add to `tests/unit/test_hamiltonian.pf`:

```fortran
@assertThat(build_velocity_matrices_qw_1d)
subroutine test_build_velocity_1d()
  use hamiltonianConstructor
  use sparse_matrices
  use definitions
  ! Build a minimal 8x1 QW Hamiltonian (1 FD point, 8 bands)
  ! H is diagonal with band offsets. velocity should be zero (no neighbors).
  type(csr_matrix) :: H_csr, vel(3)
  type(spatial_grid) :: grid
  ! ... setup H_csr as 8x8 identity, grid as 1D with 1 point
  ! call build_velocity_matrices(H_csr, grid, vel)
  ! assertEqual(vel(1)%values, (0.0_dp, 0.0_dp))
end subroutine test_build_velocity_1d
```

**Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -L unit -R hamilton -V`
Expected: FAIL (no overload accepting 1D grid + vel(3) array)

**Step 3: Write minimal implementation**

In `hamiltonianConstructor.f90`, add an overloaded interface:

```fortran
interface build_velocity_matrices
  module procedure build_velocity_matrices_2d  ! existing (wire)
  module procedure build_velocity_matrices_1d  ! new (QW)
end interface
```

New subroutine `build_velocity_matrices_1d(H_csr, grid, vel)`:

```fortran
subroutine build_velocity_matrices_1d(H_csr, grid, vel)
  type(csr_matrix), intent(in)    :: H_csr
  type(spatial_grid), intent(in)  :: grid
  type(csr_matrix), intent(out)   :: vel(3)

  integer :: k, row, col, sp_row, sp_col, Ngrid
  real(kind=dp) :: z_diff

  Ngrid = grid%ny  ! QW: ny = fdStep = number of z-grid points

  call csr_clone_structure(H_csr, vel(1))
  call csr_clone_structure(H_csr, vel(2))
  call csr_clone_structure(H_csr, vel(3))

  do row = 1, H_csr%nrows
    do k = H_csr%rowptr(row), H_csr%rowptr(row + 1) - 1
      col = H_csr%colind(k)
      sp_row = mod(row - 1, Ngrid) + 1
      sp_col = mod(col - 1, Ngrid) + 1
      z_diff = grid%z(sp_row) - grid%z(sp_col)
      vel(3)%values(k) = cmplx(0.0_dp, -z_diff, kind=dp) * H_csr%values(k)
    end do
  end do
  ! vel(1) and vel(2) remain zero (no x/y grid in QW)
end subroutine
```

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build -L unit -R hamilton -V`
Expected: PASS

**Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 tests/unit/test_hamiltonian.pf
git commit -m "feat: add build_velocity_matrices overload for QW (1D) geometry"
```

---

### Task 2: Extend `optics_config` type in defs.f90

**Files:**
- Modify: `src/core/defs.f90:199-214`

**Context:** The existing `optics_config` needs new fields for spontaneous emission, spin resolution, and separate electron/hole carrier densities.

**Step 1: Add new fields to `optics_config`**

Edit `src/core/defs.f90` lines 199-214. Add after `isbt_enabled`:

```fortran
  ! Spontaneous emission
  logical          :: spontaneous_enabled = .false.
  ! Spin resolution
  logical          :: spin_resolved = .false.
```

No new test needed — the type change is structural and will be validated by compilation.

**Step 2: Build to verify**

Run: `cmake --build build`
Expected: BUILD SUCCESS (all existing code still compiles, new fields have defaults)

**Step 3: Commit**

```bash
git add src/core/defs.f90
git commit -m "feat: extend optics_config with spontaneous and spin_resolved fields"
```

---

### Task 3: Parse new optics fields in input_parser.f90

**Files:**
- Modify: `src/io/input_parser.f90:617-620`

**Context:** After the existing `isbt_enabled` line (line 620), add parsing for the two new fields.

**Step 1: Add parsing lines**

After line 620 (`read(data_unit, *, iostat=status) label, cfg%optics%isbt_enabled`), add:

```fortran
      ! Spontaneous emission
      read(data_unit, *, iostat=status) label, cfg%optics%spontaneous_enabled
      if (status /= 0) goto 101
      print *, trim(label), cfg%optics%spontaneous_enabled

      ! Spin resolution
      read(data_unit, *, iostat=status) label, cfg%optics%spin_resolved
      if (status /= 0) goto 101
      print *, trim(label), cfg%optics%spin_resolved
```

**Step 2: Update the regression configs**

Add the two new fields to test configs in `tests/regression/configs/` that contain `optics:` blocks. Add after `isbt_enabled` line:

```
spontaneous_enabled F
spin_resolved F
```

**Step 3: Build and run regression tests**

Run: `cmake --build build && ctest --test-dir build -L regression`
Expected: All regression tests PASS (new fields default to `.false.`)

**Step 4: Commit**

```bash
git add src/io/input_parser.f90 tests/regression/configs/
git commit -m "feat: parse spontaneous_enabled and spin_resolved in optics block"
```

---

### Task 4: Add `opticalProperties` CMake target

**Files:**
- Create: `src/apps/main_optics.f90`
- Modify: `src/CMakeLists.txt:71-85`

**Step 1: Create minimal main_optics.f90**

```fortran
program opticalProperties
  use definitions
  use input_parser
  implicit none
  type(simulation_config) :: cfg

  call parse_input(cfg)
  print '(a)', 'opticalProperties: input parsed successfully'
  print '(a,L)', '  optics enabled = ', cfg%optics%enabled

end program opticalProperties
```

**Step 2: Add CMake target**

In `src/CMakeLists.txt`, after line 79 (`target_link_libraries(gfactorCalculation ...)`), add:

```cmake
# Optical properties executable
add_executable(opticalProperties
    apps/main_optics.f90
)
target_link_libraries(opticalProperties PRIVATE 8bandkp_common OpenMP::OpenMP_Fortran)
```

Also update line 82 to include the new target:

```cmake
set_target_properties(8bandkp_common bandStructure gfactorCalculation opticalProperties
```

**Step 3: Build and run**

Run: `cmake --build build && ./build/src/opticalProperties`
Expected: Prints "opticalProperties: input parsed successfully"

**Step 4: Commit**

```bash
git add src/apps/main_optics.f90 src/CMakeLists.txt
git commit -m "feat: add opticalProperties executable skeleton with CMake target"
```

---

## Phase 2: Core Optics — QW Interband Absorption with Commutator

### Task 5: Rewrite `optics_accumulate` to use commutator velocity

**Files:**
- Modify: `src/physics/optical_spectra.f90:80-158`

**Context:** The current `optics_accumulate` calls `pMatrixEleCalc` which uses the incorrect diagonal-P approximation. Replace with CSR SpMV against commutator-based velocity matrices.

**Step 1: Change the subroutine signature**

Replace lines 80-82:

```fortran
subroutine optics_accumulate(optcfg, eigvals, eigvecs, k_weight, &
  & vel, numcb, numvb, fermi_level)
  type(optics_config), intent(in) :: optcfg
  real(kind=dp), intent(in) :: eigvals(:)
  complex(kind=dp), intent(in) :: eigvecs(:,:)
  real(kind=dp), intent(in) :: k_weight
  type(csr_matrix), intent(in) :: vel(3)       ! commutator velocity matrices
  integer, intent(in) :: numcb, numvb
  real(kind=dp), intent(in) :: fermi_level
```

**Step 2: Replace the velocity matrix element computation**

Replace the inner loop (lines 126-145) that calls `pMatrixEleCalc`:

```fortran
        ! Velocity matrix elements via commutator
        px = 0.0_dp; py = 0.0_dp; pz = 0.0_dp
        do dir = 1, 3
          call csr_spmv(vel(dir), eigvecs(:,i), Ytmp, ONE, ZERO)
          Pele = zdotc(dim, eigvecs(1,numvb+j), 1, Ytmp, 1)
          select case(dir)
          case(1); px = real(Pele * conjg(Pele), kind=dp)
          case(2); py = real(Pele * conjg(Pele), kind=dp)
          case(3); pz = real(Pele * conjg(Pele), kind=dp)
          end select
        end do
```

Where `Ytmp` is a `complex(kind=dp), allocatable` work array of size `dim`, and `dim = size(eigvecs, 1)`.

**Step 3: Remove the old imports**

The module no longer needs `use gfactorFunctions` (line 4). Replace with `use sparse_matrices, only: csr_matrix, csr_spmv`.

**Step 4: Build**

Run: `cmake --build build`
Expected: COMPILE FAIL — callers of `optics_accumulate` in `main.f90` need updated signatures (handled in Task 6).

**Step 5: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "refactor: rewrite optics_accumulate with commutator velocity matrices"
```

---

### Task 6: Update `main.f90` to pass velocity matrices to optics

**Files:**
- Modify: `src/apps/main.f90`

**Context:** `main.f90` currently calls `optics_accumulate` with profile/kpterms. After the signature change, it must pass `vel(3)` CSR matrices instead.

**Step 1: Add velocity matrix construction**

In the QW branch (confinement=1), after the Hamiltonian `HT` is constructed and before the k-sweep loop, build velocity matrices. Since QW uses dense matrices (`HT` is `complex(kind=dp), allocatable`), we need to convert to CSR first:

```fortran
  ! Build velocity matrices for optics
  type(csr_matrix) :: vel_opt(3)
  if (cfg%optics%enabled) then
    ! Convert dense HT to CSR for velocity computation
    call dense_to_csr(HT, HT_csr)
    call build_velocity_matrices(HT_csr, cfg%grid, vel_opt)
    ! Rebuild HT_csr at each k_par is not needed — velocity
    ! is k-independent (commutator with r, not k). Build once at k=0.
    call csr_free(HT_csr)
  end if
```

**Step 2: Update the optics_accumulate call site**

Replace the call in the k-loop:

```fortran
  call optics_accumulate(cfg%optics, eig(:,k), eigv(:,:,k), k_weight, &
    & vel_opt, cfg%numcb, cfg%numvb, fermi_level)
```

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -L regression`
Expected: Regression tests that use optics should still pass (or need slight numerical update due to improved velocity).

**Step 4: Commit**

```bash
git add src/apps/main.f90
git commit -m "refactor: pass commutator velocity matrices to optics in bandStructure"
```

---

### Task 7: Build `opticalProperties` main loop for QW

**Files:**
- Modify: `src/apps/main_optics.f90`

**Context:** Flesh out the standalone program. It reads config, initializes confinement, builds velocity matrices at k=0, sweeps k_par, accumulates spectra, finalizes.

**Step 1: Implement the full QW branch**

```fortran
program opticalProperties
  use definitions
  use parameters
  use hamiltonianConstructor
  use input_parser
  use optical_spectra
  use sparse_matrices
  use eigensolver
  use linalg, only: zheevx
  use strain_solver
  implicit none

  type(simulation_config) :: cfg
  real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
  type(csr_matrix) :: HT_csr, vel(3)
  complex(kind=dp), allocatable :: HT(:,:)
  complex(kind=dp), allocatable :: eigvecs(:,:), work(:)
  real(kind=dp), allocatable :: eigvals(:), rwork(:)
  real(kind=dp) :: k_par, dk, k_max, k_weight
  integer :: N, dim, ik, nk, info, M
  integer :: iounit, il, iuu, lwork, NB
  real(kind=dp) :: abstol, vl, vu
  integer, allocatable :: iwork(:), ifail(:)
  complex(kind=dp), allocatable :: eigv_all(:,:)
  real(kind=dp), allocatable :: eig_all(:)

  call parse_input(cfg)

  if (.not. cfg%optics%enabled) then
    print '(a)', 'Error: optics not enabled in input.cfg'
    stop 1
  end if

  ! --- QW branch ---
  if (cfg%confinement == 1) then
    print '(a)', '=== Optical Properties: QW ==='

    ! Initialize confinement (material parameters at each z-point)
    call confinementInitialization(cfg, profile, kpterms)

    ! Build Hamiltonian at k_par=0 for velocity matrix construction
    dim = 8 * cfg%fdStep
    N = dim
    allocate(HT(N, N))
    HT = (0.0_dp, 0.0_dp)
    call ZB8bandQW(HT, 0.0_dp, profile, kpterms)

    ! Convert to CSR and build velocity matrices
    call dense_to_csr(HT, HT_csr)
    call build_velocity_matrices(HT_csr, cfg%grid, vel)
    call csr_free(HT_csr)

    ! Initialize optics accumulation
    call optics_init(cfg%optics)

    ! k_par sweep (use same wavevector settings as bandStructure)
    nk = cfg%waveVectorStep
    k_max = cfg%waveVectorMax
    dk = k_max / max(nk - 1, 1)

    do ik = 1, nk
      k_par = (ik - 1) * dk

      ! Build H(k_par)
      HT = (0.0_dp, 0.0_dp)
      call ZB8bandQW(HT, k_par, profile, kpterms)

      ! Diagonalize (zheevx for all eigenvalues)
      ! ... standard LAPACK eigensolve boilerplate ...

      ! Compute k_weight (Simpson integration)
      k_weight = simpson_weight(ik, nk) * dk

      ! Accumulate optical spectra
      call optics_accumulate(cfg%optics, eigvals, eigvecs, k_weight, &
        & vel, cfg%numcb, cfg%numvb, 0.0_dp)  ! fermi_level from SC if available
    end do

    ! Finalize and write output
    call optics_finalize(cfg%optics)
    call optics_cleanup()

    ! Free velocity matrices
    call csr_free(vel(1)); call csr_free(vel(2)); call csr_free(vel(3))

  else if (cfg%confinement == 0) then
    print '(a)', 'Bulk optics: TODO (Phase 3)'
  else if (cfg%confinement == 2) then
    print '(a)', 'Wire optics: TODO (Phase 3)'
  end if

end program opticalProperties
```

Note: The `simpson_weight` function is in `src/core/utils.f90` as `simpson_wf`. The LAPACK boilerplate for zheevx should follow the pattern in `main.f90`.

**Step 2: Build**

Run: `cmake --build build`
Expected: BUILD SUCCESS

**Step 3: Test with a QW config**

Copy a QW regression config with optics enabled to `input.cfg` and run:

```bash
cp tests/regression/configs/qw_optics.cfg input.cfg
./build/src/opticalProperties
```

Expected: Output files in `output/` with absorption spectra.

**Step 4: Commit**

```bash
git add src/apps/main_optics.f90
git commit -m "feat: implement QW optical properties main loop with commutator velocity"
```

---

## Phase 3: All Geometries and Quantities

### Task 8: Add bulk optics branch

**Files:**
- Modify: `src/apps/main_optics.f90`
- Modify: `src/physics/optical_spectra.f90`

**Context:** Bulk has no spatial discretization — the Hamiltonian is 8x8. The commutator `[r,H]=0` for bulk (no spatial grid), so use the existing analytical `ZB8bandBulk(g='g')` for velocity matrix elements. For bulk, absorption is the joint density of states weighted by k-dependent matrix elements.

**Step 1: Add bulk velocity computation**

For bulk, `vel(d)` is the 8x8 dense matrix from `ZB8bandBulk(g='g', dir=d)`. Add a helper:

```fortran
subroutine build_bulk_velocity(params, vel_dense)
  type(paramStruct), intent(in) :: params(:)
  complex(kind=dp), intent(out) :: vel_dense(3, 8, 8)
  type(wavevector) :: wv_dummy
  integer :: d
  do d = 1, 3
    call set_wavevector_direction(wv_dummy, d)
    call ZB8bandBulk(vel_dense(d,:,:), wv_dummy, params, g='g')
  end do
end subroutine
```

**Step 2: Add bulk optics accumulation**

Bulk absorption is computed as a k-space integral over the 3D Brillouin zone. The matrix elements are k-dependent. Use the existing dense matrix-vector multiply for the 8x8 system.

**Step 3: Commit**

```bash
git add src/apps/main_optics.f90 src/physics/optical_spectra.f90
git commit -m "feat: add bulk optical properties with analytical velocity matrices"
```

---

### Task 9: Add wire optics branch

**Files:**
- Modify: `src/apps/main_optics.f90`

**Context:** Wire uses CSR Hamiltonian + FEAST eigensolver (same as g-factor). Velocity matrices from `build_velocity_matrices(HT_csr, grid, vel_x, vel_y)` plus `ZB8bandGeneralized(g='g3')` for vel_z.

**Step 1: Implement wire branch**

Follow the pattern from `main_gfactor.f90` lines 105-170:

```fortran
  ! confinement == 2
  call confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, &
    & profile_2d, kpterms_2d, cfg%FDorder)
  call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg, coo_cache)
  call build_velocity_matrices(HT_csr, cfg%grid, vel(1), vel(2))
  call ZB8bandGeneralized(vel(3), 1.0_dp, profile_2d, kpterms_2d, cfg, g='g3')
```

Then sweep kz, diagonalize with FEAST, accumulate.

**Step 2: Commit**

```bash
git add src/apps/main_optics.f90
git commit -m "feat: add wire optical properties with commutator velocity"
```

---

### Task 10: Add spontaneous emission accumulation

**Files:**
- Modify: `src/physics/optical_spectra.f90`

**Context:** Spontaneous emission uses occupation factor `f_c * (1 - f_v)` instead of `f_v - f_c`. Add new accumulation arrays and subroutine.

**Step 1: Add module-level arrays**

```fortran
  real(kind=dp), allocatable, save :: spont_te(:), spont_tm(:)
```

**Step 2: Add `optics_accumulate_spontaneous` subroutine**

Same structure as `optics_accumulate` but with:

```fortran
  occ_factor = f_c * (1.0_dp - f_v)
```

**Step 3: Write spontaneous output in `optics_finalize`**

```fortran
  if (optcfg%spontaneous_enabled .and. allocated(spont_te)) then
    ! Write output/spontaneous_TE.dat and output/spontaneous_TM.dat
  end if
```

**Step 4: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "feat: add spontaneous emission accumulation and output"
```

---

### Task 11: Rewrite gain computation to use commutator velocity

**Files:**
- Modify: `src/physics/optical_spectra.f90:660-840`

**Context:** The existing `compute_gain_qw` also uses `pMatrixEleCalc`. Replace with the same commutator-based velocity as absorption.

**Step 1: Change signature**

```fortran
subroutine compute_gain_qw(optcfg, eigvals, eigvecs, k_weight, &
  & vel, numcb, numvb, carrier_density)
```

**Step 2: Replace `pMatrixEleCalc` calls with `csr_spmv`**

Same pattern as Task 5.

**Step 3: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "refactor: rewrite gain computation with commutator velocity matrices"
```

---

### Task 12: Rewrite ISBT to use commutator velocity

**Files:**
- Modify: `src/physics/optical_spectra.f90:498-559`

**Context:** Replace `z_dipole` helper with commutator-based vel(3) SpMV. This works for QW (z-confinement) and wire (x,y-confinement).

**Step 1: Rewrite `compute_isbt_absorption`**

Replace:
```fortran
  z_ij = z_dipole(eigvecs, z_grid, dz, fdstep, state_i, state_j)
  z_ij_abs2 = real(z_ij * conjg(z_ij), kind=dp)
```

With:
```fortran
  call csr_spmv(vel(dir_isbt), eigvecs(:,state_j), Ytmp, ONE, ZERO)
  pele_ij = zdotc(dim, eigvecs(1,state_i), 1, Ytmp, 1)
  p_abs2 = real(pele_ij * conjg(pele_ij), kind=dp)
```

Where `dir_isbt = 3` for QW (z-dipole) or `dir_isbt = 1` for wire (transverse dipoles).

**Step 2: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "refactor: rewrite ISBT absorption with commutator velocity"
```

---

## Phase 4: Spin Resolution

### Task 13: Create `spin_projection.f90` module

**Files:**
- Create: `src/physics/spin_projection.f90`

**Context:** Implements the Clebsch-Gordan transformation from the 8-band k.p basis to the explicit spin-orbital basis. For each eigenstate, compute spin-up/down weights.

**Step 1: Create the module**

```fortran
module spin_projection
  use definitions
  implicit none
  private
  public :: spin_weights

contains

  ! Returns w_up = sum of |<orbital_up|psi>|^2 for a single eigenstate.
  ! w_dw = 1 - w_up.
  ! Basis ordering: HH_UP(1), LH_UP(2), LH_DW(3), HH_DW(4),
  !                 SO_UP(5), SO_DW(6), EL_UP(7), EL_DW(8)
  subroutine spin_weights(psi, Ngrid, w_up, w_dw)
    complex(kind=dp), intent(in) :: psi(:)  ! (8*Ngrid)
    integer, intent(in) :: Ngrid
    real(kind=dp), intent(out) :: w_up, w_dw

    integer :: n, b, idx
    complex(kind=dp) :: c(8)
    real(kind=dp) :: px_up, px_dw, py_up, py_dw, pz_up, pz_dw, ps_up, ps_dw

    w_up = 0.0_dp
    do n = 1, Ngrid
      do b = 1, 8
        idx = (b - 1) * Ngrid + n
        c(b) = psi(idx)
      end do

      ! X_up = |HH_UP>/sqrt(2) + |LH_DW>/sqrt(6) - i|SO_DW>/sqrt(3)
      px_up = abs(c(1)/sqrt(2.0_dp) + c(3)/sqrt(6.0_dp) &
        & - cmplx(0,1,dp)*c(6)/sqrt(3.0_dp))**2
      ! Y_up = i|HH_UP>/sqrt(2) - i|LH_DW>/sqrt(6) - |SO_DW>/sqrt(3)
      py_up = abs(cmplx(0,1,dp)*c(1)/sqrt(2.0_dp) &
        & - cmplx(0,1,dp)*c(3)/sqrt(6.0_dp) - c(6)/sqrt(3.0_dp))**2
      ! Z_up = -i*sqrt(2)*|LH_UP>/sqrt(3) + |SO_UP>/sqrt(3)
      pz_up = abs(-cmplx(0,1,dp)*sqrt(2.0_dp)*c(2)/sqrt(3.0_dp) &
        & + c(5)/sqrt(3.0_dp))**2
      ! S_up = |EL_UP>
      ps_up = abs(c(7))**2

      w_up = w_up + px_up + py_up + pz_up + ps_up
    end do
    w_dw = 1.0_dp - w_up
  end subroutine spin_weights

end module spin_projection
```

**Step 2: Add to CMake sources**

In `src/CMakeLists.txt`, add `physics/spin_projection.f90` to the `8bandkp_common` library sources.

**Step 3: Build**

Run: `cmake --build build`
Expected: BUILD SUCCESS

**Step 4: Commit**

```bash
git add src/physics/spin_projection.f90 src/CMakeLists.txt
git commit -m "feat: add spin_projection module for Clebsch-Gordan spin weights"
```

---

### Task 14: Add spin-resolved accumulation

**Files:**
- Modify: `src/physics/optical_spectra.f90`

**Context:** When `spin_resolved = .true.`, decompose the velocity matrix element `|<i|v|j>|^2` into spin-up and spin-down channels using weights from `spin_weights`.

**Step 1: Add spin-resolved accumulation arrays**

```fortran
  real(kind=dp), allocatable, save :: alpha_te_up(:), alpha_te_dw(:)
  real(kind=dp), allocatable, save :: alpha_tm_up(:), alpha_tm_dw(:)
```

**Step 2: In `optics_accumulate`, add spin decomposition**

```fortran
  if (optcfg%spin_resolved) then
    call spin_weights(eigvecs(:,numvb+j), Ngrid, w_up_j, w_dw_j)
    call spin_weights(eigvecs(:,i), Ngrid, w_up_i, w_dw_i)
    ! Spin-up channel
    alpha_te_up(ie) = alpha_te_up(ie) + occ_factor * (px + py) &
      & * w_up_i * w_up_j * lineshape * k_weight
    ! Spin-down channel
    alpha_te_dw(ie) = alpha_te_dw(ie) + occ_factor * (px + py) &
      & * w_dw_i * w_dw_j * lineshape * k_weight
  end if
```

**Step 3: Write spin-resolved output files**

In `optics_finalize`, when `spin_resolved`:
- `output/absorption_TE_up.dat`, `output/absorption_TE_dw.dat`
- `output/absorption_TM_up.dat`, `output/absorption_TM_dw.dat`

**Step 4: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "feat: add spin-resolved optical spectra accumulation"
```

---

## Phase 5: Documentation and Regression Tests

### Task 15: Update lecture 06-optical-properties.md

**Files:**
- Modify: `docs/lecture/06-optical-properties.md`

**Context:** Update the lecture to reflect the commutator-based approach, the standalone program, and the unified framework.

**Step 1: Rewrite the velocity operator section**

Replace the old description of diagonal-P matrix elements with:

- The Heisenberg equation of motion: `v_alpha = [r_alpha, H] / (i*hbar)`
- Discretized form on CSR matrices: `vel(i,j) = -i * (r_i - r_j) * H(i,j)`
- Why this is more correct than the diagonal-P approximation
- Connection to the wire g-factor fix

**Step 2: Add the unified accumulation framework**

Document the shared structure for all optical quantities:
- Same velocity matrix elements from commutator
- Same Voigt broadening
- Only occupation factor differs between absorption/spontaneous/gain/ISBT

**Step 3: Document the new program**

Add section on `opticalProperties` executable, input format, output files.

**Step 4: Commit**

```bash
git add docs/lecture/06-optical-properties.md
git commit -m "docs: update lecture 06 with commutator-based optical properties"
```

---

### Task 16: Add regression test config for optical properties

**Files:**
- Create: `tests/regression/configs/qw_optics_communicator.cfg`
- Create: `tests/regression/data/qw_optics_golden.dat` (after validation)

**Context:** Create a reference config that exercises the new optics path.

**Step 1: Create test config**

```
! GaAs/AlGaAs QW optical properties test
confinement 1
FDstep 100
FDorder 2
numLayers 3
materialN GaAs 0 50
materialN AlGaAs 50 200
materialN GaAs 200 250
...
optics:
T
linewidth_lorentzian 0.030
linewidth_gaussian 0.005
refractive_index 3.3
E_min 1.0
E_max 2.5
num_energy_points 200
temperature 300.0
carrier_density 0.0
gain_enabled F
gain_carrier_density 3.0e12
isbt_enabled F
spontaneous_enabled F
spin_resolved F
```

**Step 2: Generate golden output**

Run the program and save the output as reference data.

**Step 3: Add regression test script**

Create a test that runs `opticalProperties` with this config and compares output.

**Step 4: Commit**

```bash
git add tests/regression/configs/ tests/regression/data/ tests/regression/test_optics.sh
git commit -m "test: add regression test for commutator-based optical properties"
```

---

### Task 17: Update input reference documentation

**Files:**
- Modify: `docs/reference/input-reference.md` (if it exists, or the relevant docs file)

**Context:** Document the new `spontaneous_enabled` and `spin_resolved` fields in the optics block.

**Step 1: Add field documentation**

```
optics block fields (after isbt_enabled):
  spontaneous_enabled  logical  F       Enable spontaneous emission calculation
  spin_resolved        logical  F       Enable spin-up/down decomposition
```

**Step 2: Commit**

```bash
git add docs/reference/
git commit -m "docs: document new optics input fields"
```

---

## Validation Checklist

After all tasks are complete, verify:

1. **Build**: `cmake --build build` succeeds with no warnings
2. **Unit tests**: `ctest --test-dir build -L unit` — all pass
3. **Regression tests**: `ctest --test-dir build -L regression` — all pass
4. **QW absorption**: Compare GaAs/AlGaAs QW absorption with old code output — spectra should match at zone center but differ at finite k (commutator captures gradient terms)
5. **Bulk limit**: For large QW widths, absorption converges toward bulk
6. **Gain**: Negative alpha for inverted population
7. **Spin sum rule**: `alpha_up + alpha_dw ≈ alpha_total`
8. **Wire ISBT**: Nonzero transverse dipole transitions

## Summary of New Files

| File | Lines (est.) |
|---|---|
| `src/apps/main_optics.f90` | ~250 |
| `src/physics/spin_projection.f90` | ~80 |
| `tests/regression/configs/qw_optics_communicator.cfg` | ~40 |

## Summary of Modified Files

| File | Change |
|---|---|
| `src/physics/hamiltonianConstructor.f90` | Add `build_velocity_matrices_1d` overload |
| `src/core/defs.f90` | Add 2 fields to `optics_config` |
| `src/io/input_parser.f90` | Parse 2 new optics fields |
| `src/physics/optical_spectra.f90` | Rewrite accumulation with commutator velocity |
| `src/apps/main.f90` | Update optics call sites |
| `src/CMakeLists.txt` | Add `opticalProperties` target + `spin_projection.f90` |
| `docs/lecture/06-optical-properties.md` | Update to commutator formalism |
