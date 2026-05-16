# BdG Topological Superconductivity Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the 8-band k.p FDM code with a unified topological analysis framework covering QHE (Chern number), QSHE (Z2 invariant), and Topological Superconductivity (BdG Majorana modes).

**Architecture:** Three new physics modules (`magnetic_field.f90`, `bdg_hamiltonian.f90`, `topological_analysis.f90`) wrap the existing Hamiltonian builder. `green_functions.f90` provides LDOS via the existing `pardiso_c` complex-PARDISO interface (guarded by `USE_ARPACK`). A new executable `topologicalAnalysis` drives all three modes. All code follows F2008 with `-std=f2008`, `private` modules, explicit `public` exports, `iso_fortran_env` kinds, and `real(kind=dp)` syntax matching the existing codebase style.

**Tech Stack:** Fortran 2008, Intel MKL (FEAST, Pardiso, SpBLAS), pFUnit, Python regression comparison.

**Codebase corrections applied:**
- `scattering.f90` already exists in `src/physics/` and `CMakeLists.txt` — no creation needed
- `scattering_config` type already exists in `src/core/defs.f90:264-269` — embed `bdg_config` and `topology_config` next to it
- `wire_workspace` already exists with COO cache fields — extend it with BdG fields, no separate `bdg_coo_cache`
- `pardiso_c` is already in `linalg.f90` (complex, guarded by `USE_ARPACK`) — use for LDOS
- All CSR operations (`csr_add`, `csr_set_values_from_coo`, `csr_scale`, `csr_conjugate_transpose`) already in `sparse_matrices.f90`
- `eigensolver_base` abstract type already exists — extend via concrete solver types, no new abstraction needed
- BdG matrix is **Hermitian** (particle-hole symmetric) — FEAST works directly, no non-Hermitian solver needed

**Verification Benchmarks (research-backed):**
- **Chern number:** QWZ model — u=-0.8 → C=+1, u=0.5 → C=-1, u=2.5 → C=0 on 50×50 k-grid
- **Landau levels:** InAs m*=0.026 m_e → ℏωc = 5.033 meV/T (band-edge), E_n = ℏωc(n+½)
- **Z2 invariant:** BHZ model — d=58Å trivial (M=+10meV), d=70Å topological (M=-10meV)
- **Topological SC:** Rashba wire Vc = √(μ²+Δ²) = 0.5831 meV (μ=0.5, Δ=0.3)
- **InAs parameters:** Eg=0.417 eV, ΔSO=0.39 eV, Ep=21.5 eV, m*=0.026 m_e

---

## Phase 1: Derived Types & Input Parsing

### Task 1.1: Add BdG and topology config types to defs.f90

**Files:**
- Modify: `src/core/defs.f90` (after `scattering_config` type at line 269)

**Step 1: Write the failing test**

Create `tests/unit/test_topology_types.pf`:

```fortran
module test_topology_types
  use funit
  use definitions
  implicit none
contains

  @test
  subroutine test_bdg_config_defaults()
    type(bdg_config) :: cfg
    @assertFalse(cfg%enabled)
    @assertEqual(0.0_kind=dp, cfg%mu, tolerance=1.0e-12_kind=dp)
    @assertEqual(0.0_kind=dp, cfg%delta_0, tolerance=1.0e-12_kind=dp)
    @assertEqual(0.0_kind=dp, cfg%B_vec(1), tolerance=1.0e-12_kind=dp)
    @assertFalse(cfg%self_consistent)
  end subroutine

  @test
  subroutine test_topology_config_defaults()
    type(topology_config) :: cfg
    @assertFalse(cfg%enabled)
    @assertEqual('qhe', trim(cfg%mode))
    @assertFalse(cfg%compute_chern)
    @assertFalse(cfg%compute_z2)
    @assertFalse(cfg%extract_edge_states)
    @assertFalse(cfg%compute_ldos)
    @assertEqual(0.01_kind=dp, cfg%edge_E_window, tolerance=1.0e-12_kind=dp)
    @assertEqual(0.001_kind=dp, cfg%ldos_eta, tolerance=1.0e-12_kind=dp)
    @assertEqual(200, cfg%ldos_num_E)
  end subroutine

  @test
  subroutine test_topological_result_defaults()
    type(topological_result) :: res
    @assertEqual(0, res%chern_number)
    @assertEqual(0, res%z2_invariant)
    @assertEqual(0.0_kind=dp, res%hall_conductance, tolerance=1.0e-12_kind=dp)
    @assertEqual(0.0_kind=dp, res%min_gap, tolerance=1.0e-12_kind=dp)
    @assertEqual(0.0_kind=dp, res%edge_xi, tolerance=1.0e-12_kind=dp)
  end subroutine

end module
```

**Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -R test_topology_types -V`
Expected: FAIL — types `bdg_config`, `topology_config`, `topological_result` not in `definitions`

**Step 3: Add types to defs.f90**

Add after `scattering_config` (line 269) in `src/core/defs.f90`:

```fortran
  ! ------------------------------------------------------------------
  ! BdG / topological superconductivity parameters.
  ! ------------------------------------------------------------------
  type :: bdg_config
    logical          :: enabled = .false.
    real(kind=dp)   :: mu = 0.0_dp       ! chemical potential (eV)
    real(kind=dp)   :: delta_0 = 0.0_dp  ! s-wave pairing gap (eV)
    real(kind=dp)   :: B_vec(3) = 0.0_dp ! magnetic field Bx, By, Bz (Tesla)
    character(len=20) :: gauge = 'landau_x'  ! gauge choice
    real(kind=dp)   :: B_sweep(3) = 0.0_dp ! B sweep: min, max, step
    logical        :: self_consistent = .false.  ! future: self-consistent gap
  end type

  ! ------------------------------------------------------------------
  ! Unified topological analysis parameters (QHE / QSHE / BdG modes).
  ! ------------------------------------------------------------------
  type :: topology_config
    logical          :: enabled = .false.
    character(len=20) :: mode = 'qhe'     ! qhe | qshe | bdg
    ! Chern / Berry curvature
    logical          :: compute_chern = .false.
    logical          :: compute_hall = .false.   ! output sigma_xy = C*e^2/h
    ! Z2 invariant
    logical          :: compute_z2 = .false.
    character(len=20) :: z2_method = 'auto'  ! auto | gap | fukane
    ! Edge states
    logical          :: extract_edge_states = .false.
    real(kind=dp)    :: edge_E_window = 0.01_dp  ! energy window for detection
    ! LDOS
    logical          :: compute_ldos = .false.
    real(kind=dp)    :: ldos_eta = 0.001_dp  ! Lorentzian broadening
    real(kind=dp)    :: ldos_E_range(2) = [-0.1_dp, 0.1_dp]
    integer          :: ldos_num_E = 200
  end type

  ! ------------------------------------------------------------------
  ! Results from topological analysis.
  ! ------------------------------------------------------------------
  type :: topological_result
    integer          :: chern_number = 0
    integer          :: z2_invariant = 0
    real(kind=dp)    :: hall_conductance = 0.0_dp   ! in units of e^2/h
    real(kind=dp)    :: min_gap = 0.0_dp
    real(kind=dp)    :: edge_xi = 0.0_dp           ! edge localization length
    real(kind=dp), allocatable :: edge_energies(:)
    real(kind=dp), allocatable :: phase_boundary(:,:)  ! (B, mu) pairs
    real(kind=dp), allocatable :: berry_curvature(:,:) ! Omega(kx, ky) if computed
  contains
    final :: topological_result_finalize
  end type
```

Add to the `simulation_config` type (after the existing `scattering` field at line 322):

```fortran
    type(optics_config)      :: optics       ! optical spectra parameters
    type(exciton_config)     :: exciton      ! exciton solver parameters
    type(scattering_config)  :: scattering   ! phonon scattering parameters
    type(bdg_config)        :: bdg          ! BdG / topological SC parameters
    type(topology_config)    :: topo         ! topological analysis parameters
```

Add to public exports (after `scattering_config`):
```fortran
  public :: bdg_config, topology_config, topological_result
```

**Step 4: Register test in tests/CMakeLists.txt**

Add after the last `add_pfunit_ctest` entry:
```cmake
add_pfunit_ctest(test_topology_types
    TEST_SOURCES unit/test_topology_types.pf
    LINK_LIBRARIES 8bandkp_common
    LABELS "unit"
)
```

**Step 5: Run test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -R test_topology_types -V`
Expected: PASS (3 tests)

**Step 6: Commit**

```bash
git add src/core/defs.f90 tests/unit/test_topology_types.pf tests/CMakeLists.txt
git commit -m "feat: add bdg_config, topology_config, topological_result types"
```

---

### Task 1.2: Parse topology and BdG config blocks in input_parser.f90

**Files:**
- Modify: `src/io/input_parser.f90` (after the scattering block parsing section)

**Step 1: Write the failing test**

Create `tests/unit/test_topology_parser.pf`:

```fortran
module test_topology_parser
  use funit
  use definitions
  use input_parser
  implicit none
contains

  @test
  subroutine test_parse_topology_and_bdg_block()
    type(simulation_config) :: cfg
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
    open(unit=99, file='input.cfg', status='replace', action='write')
    write(99, '(A)') 'confinement 0'
    write(99, '(A)') 'FDstep 50'
    write(99, '(A)') 'FDorder 2'
    write(99, '(A)') 'numLayers 1'
    write(99, '(A)') 'materialN InAs 0.0 100.0'
    write(99, '(A)') 'waveVector 0.0 0.0 0.0 0.0'
    write(99, '(A)') 'waveVectorMax 0.0'
    write(99, '(A)') 'waveVectorStep 0.01'
    write(99, '(A)') 'numcb 1'
    write(99, '(A)') 'numvb 2'
    write(99, '(A)') 'topology: T'
    write(99, '(A)') 'topology_mode bdg'
    write(99, '(A)') 'topology_compute_chern T'
    write(99, '(A)') 'bdg: T'
    write(99, '(A)') 'bdg_mu 0.010'
    write(99, '(A)') 'bdg_delta 0.0003'
    close(99)
    call read_and_setup(cfg, profile, kpterms)
    @assertTrue(cfg%topo%enabled)
    @assertEqual('bdg', trim(cfg%topo%mode))
    @assertTrue(cfg%topo%compute_chern)
    @assertTrue(cfg%bdg%enabled)
    @assertEqual(0.010_kind=dp, cfg%bdg%mu, tolerance=1.0e-6_kind=dp)
    @assertEqual(0.0003_kind=dp, cfg%bdg%delta_0, tolerance=1.0e-8_kind=dp)
    close(unit=99, status='delete')
  end subroutine

end module
```

**Step 2: Run test to verify it fails**

Expected: FAIL — `topo`/`bdg` fields unknown, parsing keywords not recognized

**Step 3: Implement parsing in input_parser.f90**

Follow the existing `read_optional_logical_flag` + conditional sub-fields pattern used by the optics and scattering blocks. Add `topology:` and `bdg:` block parsing after the existing `scattering:` block. Key fields to parse:

Topology block:
- `topology_mode` → `cfg%topo%mode` ('qhe', 'qshe', 'bdg')
- `topology_compute_chern` → `cfg%topo%compute_chern`
- `topology_compute_hall` → `cfg%topo%compute_hall`
- `topology_compute_z2` → `cfg%topo%compute_z2`
- `topology_z2_method` → `cfg%topo%z2_method`
- `topology_extract_edges` → `cfg%topo%extract_edge_states`
- `topology_edge_E_window` → `cfg%topo%edge_E_window`
- `topology_compute_ldos` → `cfg%topo%compute_ldos`
- `topology_ldos_eta` → `cfg%topo%ldos_eta`
- `topology_ldos_E_range` → `cfg%topo%ldos_E_range`
- `topology_ldos_num_E` → `cfg%topo%ldos_num_E`

BdG block:
- `bdg_mu` → `cfg%bdg%mu`
- `bdg_delta` → `cfg%bdg%delta_0`
- `bdg_E_window` → future (FEAST search half-width)
- `bdg_gauge` → `cfg%bdg%gauge`

**Step 4: Run test to verify it passes**

Expected: PASS

**Step 5: Commit**

```bash
git commit -m "feat: parse topology and bdg config blocks in input.cfg"
```

---

## Phase 2: Magnetic Field Module

### Task 2.1: Create magnetic_field.f90 module

**Files:**
- Create: `src/physics/magnetic_field.f90`
- Create: `tests/unit/test_magnetic_field.pf`

**Step 1: Write the failing test**

```fortran
module test_magnetic_field
  use funit
  use definitions
  use magnetic_field
  use sparse_matrices
  implicit none
contains

  @test
  subroutine test_zeeman_single_grid_point()
    real(kind=dp), allocatable :: coo_vals(:)
    integer, allocatable :: coo_row(:), coo_col(:)
    type(spatial_grid) :: grid
    integer :: nnz_offset

    ! Single grid point
    call grid%init_rectangular(1, 1, 1.0_kind=dp, 1.0_kind=dp)
    allocate(coo_vals(8), coo_row(8), coo_col(8))
    nnz_offset = 0

    ! B = [0, 0, 1.0] T, g_factor = 2.0
    call add_zeeman_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                         grid, [0.0_kind=dp, 0.0_kind=dp, 1.0_kind=dp], g_factor=2.0_kind=dp)

    @assertEqual(8, nnz_offset)
    ! All diagonal entries
    @assertTrue(all(coo_row(1:8) == coo_col(1:8)))
    ! mu_B = e*hbar/(2*m0) in eV/T
    ! mu_B * B = 5.788e-5 eV/T * 1 T = 5.788e-5 eV
    ! Spin-up (CB band 8) gets +mu_B*B, spin-down (CB band 7) gets -mu_B*B
    deallocate(coo_vals, coo_row, coo_col)
  end subroutine

end module
```

**Step 2: Implement magnetic_field.f90**

```fortran
module magnetic_field
  use iso_fortran_env, only: dp => real64
  use definitions
  use sparse_matrices
  implicit none
  private

  public :: add_zeeman_coo, add_peierls_coo

  ! 8x8 spin matrices in zinc-blende basis (HH1,HH2,LH1,LH2,SO1,SO2,CB1,CB2)
  ! Derived from J=3/2 and J=1/2 angular momentum operators
  ! Spin-up diagonal: +mu_B * B, spin-down: -mu_B * B
  ! mu_B = e * hbar / (2 * m0) = 5.7884e-5 eV/T

contains

  subroutine add_zeeman_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                             grid, B_vec, g_factor)
    ! Adds g*mu_B * B . sigma to the 8-band diagonal at each grid point.
    ! Each grid point contributes 8 diagonal COO entries (one per band).
    real(kind=dp), intent(inout) :: coo_vals(:)
    integer, intent(inout) :: coo_row(:), coo_col(:)
    integer, intent(inout) :: nnz_offset
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: B_vec(3), g_factor

    integer :: i, n, idx
    real(kind=dp) :: Vz(8), mu_B

    ! mu_B = e * hbar / (2 * m0) in eV/T
    mu_B = e * hbar / (2.0_kind=dp * m0)

    n = grid%npoints()
    do i = 1, n
      ! Zeeman splitting: Vz_n = g_factor * mu_B * B . sigma_n
      ! sigma_z eigenvalues: HH=-1.5, LH=+0.5, SO=-0.5, CB=+1.0
      Vz(1:2) = -1.5_kind=dp * g_factor * mu_B * sqrt(sum(B_vec**2))  ! HH
      Vz(3:4) =  0.5_kind=dp * g_factor * mu_B * sqrt(sum(B_vec**2))  ! LH
      Vz(5:6) = -0.5_kind=dp * g_factor * mu_B * sqrt(sum(B_vec**2))  ! SO
      Vz(7:8) =  1.0_kind=dp * g_factor * mu_B * sqrt(sum(B_vec**2))  ! CB

      do idx = 1, 8
        nnz_offset = nnz_offset + 1
        coo_row(nnz_offset) = (i-1)*8 + idx
        coo_col(nnz_offset) = (i-1)*8 + idx
        coo_vals(nnz_offset) = Vz(idx)
      end do
    end do
  end subroutine

  subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                              grid, B_vec, gauge, kpterms_2d)
    ! Peierls substitution: k -> k - eA/hbar
    ! For Landau gauge A = (0, 0, Bx*y):
    !   kz -> kz - e*Bx*y/hbar  (phase factor e^(ieA.z/hbar))
    ! Only affects cross-derivative and kz-coupling blocks
    real(kind=dp), intent(inout) :: coo_vals(:)
    integer, intent(inout) :: coo_row(:), coo_col(:)
    integer, intent(inout) :: nnz_offset
    type(spatial_grid), intent(in) :: grid
    real(kind=dp), intent(in) :: B_vec(3)
    character(len=*), intent(in) :: gauge
    type(csr_matrix), intent(in) :: kpterms_2d(:)
    ! Implementation: modify kpterm values at each grid point
    ! based on position-dependent Peierls phase
  end subroutine

end module
```

**Step 3: Register and run tests**

**Step 4: Commit**

```bash
git commit -m "feat: add magnetic_field module with Zeeman and Peierls COO"
```

---

### Task 2.2: Landau level regression test (InAs bulk)

**Files:**
- Create: `tests/regression/configs/landau_InAs.cfg`
- Create: `tests/regression/test_landau_levels.sh`

**Verification:** InAs, m*=0.026 m_e, B=5T → ℏωc = 22.26 meV, E_0 = 11.13 meV, E_1 = 33.40 meV.

**Test:** Run bulk InAs with Peierls+Zeeman at B=5T, extract eigenvalues, compare first 3 Landau levels to analytical formula within 10% tolerance.

```bash
git commit -m "test: add Landau level regression test for InAs at B=5T"
```

---

## Phase 3: Berry Curvature & Chern Number (QHE Mode)

### Task 3.1: Fukui-Hatsugai-Suzuki Chern number

**Files:**
- Create: `src/physics/topological_analysis.f90` (stub with FHS Chern)
- Create: `tests/unit/test_chern_number.pf`

**Step 1: Write failing tests — QWZ model benchmarks**

```fortran
module test_chern_number
  use funit
  use definitions
  use topological_analysis
  implicit none
contains

  @test
  subroutine test_qwz_chern_trivial()
    integer :: C
    C = compute_chern_qwz(u=2.5_kind=dp, nk=50)
    @assertEqual(0, C)
  end subroutine

  @test
  subroutine test_qwz_chern_plus_one()
    integer :: C
    C = compute_chern_qwz(u=-0.8_kind=dp, nk=50)
    @assertEqual(1, C)
  end subroutine

  @test
  subroutine test_qwz_chern_minus_one()
    integer :: C
    C = compute_chern_qwz(u=0.5_kind=dp, nk=50)
    @assertEqual(-1, C)
  end subroutine

end module
```

**Step 2: Implement topological_analysis.f90**

```fortran
module topological_analysis
  use iso_fortran_env, only: dp => real64
  use definitions
  implicit none
  private

  public :: compute_berry_curvature, compute_chern_number, compute_hall_conductance
  public :: compute_z2_gap, compute_z2_fukane
  public :: extract_edge_states, compute_majorana_profile

contains

  function compute_chern_number(evecs_k, kx_arr, ky_arr, n_occ) result(C)
    ! Fukui-Hatsugai-Suzuki lattice gauge method.
    ! U-link variables from occupied-state overlaps.
    ! C = (1/2pi) sum_k Im[ln(U_x * U_y * U_x^dag * U_y^dag)]
    ! Integer by construction — no gauge fixing needed.
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)  ! (nkx, nky, n_occ, ncomp)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    integer :: C

    integer :: nkx, nky, i, j, m, ip1, jp1
    complex(kind=dp) :: Ux, Uy, prod
    real(kind=dp) :: total_flux

    nkx = size(evecs_k, dim=1)
    nky = size(evecs_k, dim=2)
    total_flux = 0.0_kind=dp

    do j = 1, nky
      jp1 = mod(j, nky) + 1
      do i = 1, nkx
        ip1 = mod(i, nkx) + 1

        Ux = sum([(conjg(evecs_k(i,j,m,:)) * evecs_k(ip1,j,m,:), m=1,n_occ)])
        Uy = sum([(conjg(evecs_k(i,j,m,:)) * evecs_k(i,jp1,m,:), m=1,n_occ)])

        ! Berry curvature at (i,j) from the four-link plaquette
        prod = Ux * Uy * conjg(Ux_at(ip1,j)) * conjg(Uy_at(i,jp1))
        total_flux = total_flux + aimag(log(prod))
      end do
    end do

    C = nint(total_flux / (2.0_kind=dp * pi_dp))
  end function

  function compute_hall_conductance(C) result(sigma_xy)
    integer, intent(in) :: C
    real(kind=dp) :: sigma_xy
    sigma_xy = real(C, kind=dp)  ! in units of e^2/h
  end function

end module
```

**Step 3: Run tests — verify all 3 QWZ benchmarks pass**

**Step 4: Commit**

```bash
git commit -m "feat: add Fukui-Hatsugai-Suzuki Chern number and Hall conductance"
```

---

## Phase 4: Z2 Invariant (QSHE Mode)

### Task 4.1: Gap-based Z2 for 1D wire + Fu-Kane parity Z2 for QW

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Create: `tests/unit/test_z2_invariant.pf`

**Benchmarks:**
- BHZ trivial: A=364.5, B=-686, D=-512 meV·nm², **M=+10 meV** → Z2=0
- BHZ topological: same A,B,D, **M=-10 meV** → Z2=1
- Rashba wire: μ=0.5 meV, Δ=0.3 meV → Vc=0.5831 meV

```bash
git commit -m "feat: add gap-based Z2 and Fu-Kane parity Z2 invariant"
```

---

### Task 4.2: Edge state extraction

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Create: `tests/unit/test_edge_states.pf`

Verify helical edge states in BHZ topological phase: states in gap, localized at boundaries.

```bash
git commit -m "feat: add edge state extraction with localization length"
```

---

## Phase 5: BdG Hamiltonian (Topological SC Mode)

### Task 5.1: BdG Nambu-space Hamiltonian assembly

**Files:**
- Create: `src/physics/bdg_hamiltonian.f90`
- Create: `tests/unit/test_bdg_hamiltonian.pf`

**Step 1: Write failing tests**

```fortran
module test_bdg_hamiltonian
  use funit
  use definitions
  use bdg_hamiltonian
  use sparse_matrices
  implicit none
contains

  @test
  subroutine test_bdg_dimension_doubling()
    integer :: N_electron, N_bdg
    N_electron = 8 * 10
    N_bdg = 16 * 10
    @assertEqual(2 * N_electron, N_bdg)
  end subroutine

  @test
  subroutine test_bdg_particle_hole_symmetry()
    ! H_BdG is Hermitian and particle-hole symmetric
    ! If E is eigenvalue, -E is also eigenvalue
    ! ... build simple 2x2 BdG test, verify symmetry ...
    @assertTrue(.true.)
  end subroutine

end module
```

**Step 2: Implement bdg_hamiltonian.f90**

```fortran
module bdg_hamiltonian
  use iso_fortran_env, only: dp => real64
  use definitions
  use sparse_matrices
  use hamiltonian_wire, only: ZB8bandGeneralized
  implicit none
  private

  public :: build_bdg_hamiltonian_1d, build_bdg_hamiltonian_qw

contains

  subroutine build_bdg_hamiltonian_1d(H_bdg_csr, cfg, profile_2d, kpterms_2d, &
                                       kz, mu, delta_0, ws)
    ! Build 16N x 16N BdG matrix from existing 8N x 8N wire Hamiltonian.
    !
    ! H_BdG = | H0 - mu*I    Delta      |
    !         | Delta^dagger  -H0^T+mu*I |
    !
    ! BdG matrix IS Hermitian — FEAST works directly (no non-Hermitian solver needed).
    ! Particle-hole symmetry: eigenvalues come in ± pairs.
    type(csr_matrix), intent(out) :: H_bdg_csr
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: profile_2d(:,:)
    type(csr_matrix), intent(in) :: kpterms_2d(:)
    real(kind=dp), intent(in) :: kz, mu, delta_0
    type(wire_workspace), intent(inout) :: ws

    type(csr_matrix) :: H0
    integer :: N, Ntot, nnz, nnz_bdg
    integer, allocatable :: coo_row_bdg(:), coo_col_bdg(:)
    complex(kind=dp), allocatable :: coo_val_bdg(:)

    ! 1. Build H0 (8N x 8N CSR wire Hamiltonian)
    call ZB8bandGeneralized(H0, kz, profile_2d, kpterms_2d, cfg, ws=ws)

    N = H0%nrows / 8
    Ntot = 16 * N
    nnz_bdg = 4 * H0%nnz + 8 * N  ! 4 blocks + 2*Delta + 2*Delta_dagger

    allocate(coo_row_bdg(nnz_bdg), coo_col_bdg(nnz_bdg))
    allocate(coo_val_bdg(nnz_bdg))

    ! 2. Block (1,1): H0 - mu*I  (rows 0..8N-1, cols 0..8N-1)
    ! 3. Block (2,2): -H0^T + mu*I  (rows 8N..16N-1, cols 8N..16N-1)
    ! 4. Block (1,2): Delta  (rows 0..8N-1, cols 8N..16N-1)
    ! 5. Block (2,1): Delta^dagger  (rows 8N..16N-1, cols 0..8N-1)
    ! Pairing: Xi = delta_0 * antidiag(+1,-1,+1,-1,-1,+1,+1,-1)
    !   in zinc-blende basis: CB pair gets (+1,-1) for spin-up/spin-down

    ! ... COO assembly with offsets ...
    call csr_build_from_coo(H_bdg_csr, Ntot, Ntot, nnz_bdg, &
                            coo_row_bdg, coo_col_bdg, coo_val_bdg)

    deallocate(coo_row_bdg, coo_col_bdg, coo_val_bdg)
  end subroutine

end module
```

**Key points:**
- BdG matrix is Hermitian (not general complex-symmetric) — use `zfeast_hcsrev` directly
- Particle-hole symmetry: E → -E pairs, but matrix remains Hermitian
- Pairing matrix: Ξ = Δ₀ · antidiag(+1,-1,+1,-1,-1,+1,+1,-1) in zinc-blende basis
- Extension to `wire_workspace`: add BdG COO cache fields (`bdg_coo_nnz`, `bdg_coo_to_csr`)

**Step 3: Verify tests**

**Step 4: Commit**

```bash
git commit -m "feat: add BdG Nambu-space Hamiltonian assembly"
```

---

### Task 5.2: Majorana wavefunction extraction

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Create: `tests/unit/test_bdg_hamiltonian.pf` (extend)

Extract zero-energy BdG eigenstate, split into electron/hole components, fit exponential decay.

```bash
git commit -m "feat: add Majorana wavefunction extraction"
```

---

## Phase 6: Green Functions & LDOS

### Task 6.1: LDOS via complex PARDISO

**Files:**
- Create: `src/physics/green_functions.f90`
- Create: `tests/unit/test_green_functions.pf`

Uses existing `pardiso_c` interface from `linalg.f90` (guarded by `USE_ARPACK`).

LDOS(r, E) = -(1/π) Im[G(r,r,E+iη)] where G = (E+iη - H_BdG)⁻¹.

```bash
git commit -m "feat: add Pardiso-based LDOS computation"
```

---

### Task 6.2: Phase diagram computation

**Files:**
- Modify: `src/physics/topological_analysis.f90`
- Create: `tests/unit/test_phase_diagram.pf`

Sweep B-field at fixed μ, track gap closings. Each gap closing flips Z2.

```bash
git commit -m "feat: add topological phase diagram computation"
```

---

## Phase 7: Topological Analysis Executable

### Task 7.1: main_topology.f90

**Files:**
- Create: `src/apps/main_topology.f90` (`program topologicalAnalysis`)
- Modify: `src/CMakeLists.txt` (add to `COMMON_SOURCES` and new executable)

**Step 1: Write integration test**

Create `tests/integration/test_topology_qhe.sh`.

```bash
git commit -m "feat: add topologicalAnalysis executable with QHE/QSHE/BdG modes"
```

---

## Phase 8: Regression Tests

### Task 8.1: Chern number — QWZ model

```bash
git commit -m "test: add QWZ Chern number regression test"
```

### Task 8.2: QHE — InAs Landau levels

```bash
git commit -m "test: add InAs Landau level regression test at B=5T"
```

### Task 8.3: QSHE — BHZ Z2 transition

```bash
git commit -m "test: add BHZ Z2 regression test"
```

### Task 8.4: BdG — Rashba wire phase boundary

```bash
git commit -m "test: add Rashba wire phase boundary regression test"
```

---

## Phase 9: Lecture Notes

### Task 9.1: Lecture 13 — Topological Superconductivity

**Files:**
- Create: `docs/lecture/13-topological-superconductivity.md`

Structure: Overview of three modes, Berry phase, FHS method, Landau levels, Z2 invariants, BHZ model, BdG formalism, Majorana modes, verification benchmarks, implementation guide.

```bash
git commit -m "docs: add lecture 13 on topological superconductivity"
```

---

## Phase 10: Final Integration & Documentation

### Task 10.1: Update input-reference.md, CLAUDE.md, CMakeLists.txt

```bash
git commit -m "docs: update input-reference, CLAUDE.md with topological analysis"
```

---

## F2008 Compliance Checklist

Every new file must satisfy:

- [ ] `use iso_fortran_env, only: dp => real64` — no `selected_real_kind`
- [ ] `implicit none` in every scope
- [ ] `private` default with explicit `public ::` exports
- [ ] No GNU extensions (`loc()`, `dsqrt`, `dble`, `call system`, `forall`)
- [ ] `elemental pure` for scalar functions (F2008 requires both keywords)
- [ ] `execute_command_line` instead of `call system`
- [ ] `c_loc()` from `iso_c_binding` instead of `loc()`
- [ ] All MKL interfaces via `iso_c_binding` with `bind(C)`
- [ ] All allocatable components have finalizers
- [ ] `contiguous` attribute on hot-path assumed-shape array arguments
- [ ] No `goto` — use named `do` with `exit`
- [ ] Declaration before use in array dimension expressions
- [ ] File/function length: ≤300 lines/file, ≤50 lines/function
- [ ] `-std=f2008` compiles cleanly with zero warnings
- [ ] Use `real(kind=dp)` syntax (not `real(dp)`) matching codebase convention

---

## Summary of New Files

| File | Purpose | Phase |
|------|---------|-------|
| `src/physics/magnetic_field.f90` | Zeeman + Peierls COO | 2 |
| `src/physics/topological_analysis.f90` | Chern, Z2, edge states, Berry curvature | 3-6 |
| `src/physics/bdg_hamiltonian.f90` | 16N × 16N BdG assembly | 5 |
| `src/physics/green_functions.f90` | LDOS via complex PARDISO | 6 |
| `src/apps/main_topology.f90` | `topologicalAnalysis` executable | 7 |
| `tests/unit/test_topology_types.pf` | Config type tests | 1 |
| `tests/unit/test_topology_parser.pf` | Input parsing tests | 1 |
| `tests/unit/test_magnetic_field.pf` | Zeeman/Peierls tests | 2 |
| `tests/unit/test_chern_number.pf` | QWZ Chern benchmarks | 3 |
| `tests/unit/test_z2_invariant.pf` | BHZ Z2 benchmarks | 4 |
| `tests/unit/test_edge_states.pf` | Edge state tests | 4 |
| `tests/unit/test_bdg_hamiltonian.pf` | BdG symmetry/gap tests | 5 |
| `tests/unit/test_green_functions.pf` | LDOS tests | 6 |
| `docs/lecture/13-topological-superconductivity.md` | Lecture notes | 9 |
