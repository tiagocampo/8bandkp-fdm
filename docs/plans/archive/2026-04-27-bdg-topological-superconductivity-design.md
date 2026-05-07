# Topological Band Theory & BdG Superconductivity Design

Date: 2026-04-27
Branch: feature/bdg-topological-sc (planned)
Author: Tiago de Campos

## Overview

Extend the 8-band k.p finite difference code with a **unified topological analysis framework** covering three classes of topological phenomena:

1. **Integer Quantum Hall Effect (QHE)** — Landau levels, Chern number, Hall plateaus (bare k.p + B-field)
2. **Quantum Spin Hall Effect / 2D Topological Insulators** — Fu-Kane Z2 invariant, helical edge states (bare k.p, no B-field)
3. **Topological Superconductivity / Majorana modes** — BdG formalism, Majorana wavefunctions (k.p + B-field + SC pairing)

All three share the same infrastructure: Berry curvature computation, edge state extraction, FEAST eigensolver, and spatial grid resolution. The framework adds:
- Full Nambu-space Hamiltonian (16N x 16N BdG matrix) for topological SC
- Magnetic field with Zeeman splitting and Peierls orbital effects (Landau gauge)
- Unified topological invariant computation (Chern number, Z2 by gap method and Fu-Kane parity)
- Edge state extraction (Majorana, chiral, helical)
- Berry curvature and Hall conductance
- Local density of states via sparse Green's function (Pardiso)

## Physics Background

### Three operational modes

The framework operates in three modes sharing the same Hamiltonian infrastructure:

| Mode | Hamiltonian | B-field | Topological invariant | Key physics |
|---|---|---|---|---|
| **QHE** | Bare k.p (8N) | Yes (Zeeman + Peierls) | Chern number | Landau levels, Hall plateaus σ_xy = Ce²/h |
| **QSHE / 2D TI** | Bare k.p (8N) | No | Z2 (Fu-Kane) | Band inversion, helical edge states |
| **Topological SC** | BdG doubled (16N) | Yes (Zeeman + Peierls) | Z2 (gap method) | Majorana zero modes |

All modes use the same Berry curvature machinery and edge state extraction.

### Quantum Hall Effect

A 2D electron gas in a perpendicular magnetic field develops quantized Landau levels. The Hall conductance is quantized as σ_xy = Ce²/h where C is the Chern number (TKNN invariant) of the occupied bands. The 8-band k.p with Peierls substitution captures the full Landau level structure including spin-orbit corrections (important for wide-gap semiconductors where SO coupling modifies the level degeneracies).

The Chern number is computed by integrating the Berry curvature over the 2D Brillouin zone using the Fukui-Hatsugai-Suzuki lattice gauge method.

**Chiral edge states** appear naturally in QHE mode: when C ≠ 0, the gap hosts |C| unidirectional edge channels. The edge state extraction identifies them from their spatial localization at the boundary.

**Connection to Weyl physics (chiral anomaly):** The chiral zeroth Landau level (n=0) propagates in only one direction along B — this is the 1D analogue of a Weyl fermion. Under parallel E·B, chiral charge is not conserved (chiral anomaly), producing negative magnetoresistance. The QHE mode captures this through:
- Chiral zeroth Landau level in the spectrum
- Longitudinal conductivity σ_zz(B) showing negative magnetoresistance (chiral anomaly signature)
- Edge states as the 2D analogue of Fermi arcs in 3D Weyl semimetals

### Quantum Spin Hall Effect / 2D Topological Insulators

In inversion-symmetric quantum wells with band inversion (e.g., HgTe/CdTe), the Z2 topological invariant distinguishes trivial (Z2 = 0) from topological (Z2 = 1) phases. The canonical example is the BHZ model (Bernevig-Hughes-Zhang 2006), derived from a k.p Hamiltonian equivalent to the one in this code.

For inversion-symmetric systems, the Z2 invariant is computed via the Fu-Kane parity formula:

```
(-1)^ν = Π_{i=1}^{4} Π_{m∈occ} ξ_m(Γ_i)
```

where ξ_m = ±1 are the parity eigenvalues of occupied bands at the four time-reversal invariant momenta (TRIM) Γ_i in the 2D BZ. This only requires eigenvalues at 4 k-points — no Berry curvature integration needed.

Key physics: the band inversion at Γ in HgTe/CdTe QWs flips the parity of the CB, changing Z2 from 0 to 1 above a critical well width, producing helical edge states.

### Topological Superconductivity

Three ingredients are needed:
1. **Spin-orbit coupling** — already present in the 8-band k.p Hamiltonian
2. **s-wave superconducting pairing** — proximity-induced from an SC contact (e.g., Al shell on InAs wire)
3. **Time-reversal symmetry breaking** — magnetic field (Zeeman + orbital)

When all three are present in a 1D wire, the system enters a topological phase (class D) hosting Majorana zero modes at the wire ends.

### BdG formalism

The Bogoliubov-de Gennes Hamiltonian in Nambu basis {c, c†}:

```
H_BdG = | H₀(k) - μI     Δ          |
        | Δ†            -H₀ᵀ(k) + μI |
```

Where:
- H₀(k) is the existing 8N x 8N k.p Hamiltonian
- μ is the chemical potential
- Δ is the pairing matrix (s-wave proximity gap)
- H₀ᵀ is the matrix transpose (hole sector)

The full matrix dimension is 16N (doubling from 8-band to include electron-hole Nambu space).

### s-wave pairing in the 8-band basis

The pairing matrix couples time-reversed partners via the operation |J, mⱼ⟩ ↔ (-1)^(J-mⱼ)|J, -mⱼ⟩.

In the zinc-blende basis ordering (HH↑, LH↑, LH↓, HH↓, SO↑, SO↓, CB↑, CB↓):

```
Ξ = Δ₀ · antidiag(+1, -1, +1, -1, -1, +1, +1, -1)
```

The full pairing matrix in the Nambu Hamiltonian is Δ = Ξ ⊗ I_N (diagonal in real space for uniform s-wave).

## Module Architecture

### New files

```
src/physics/
  bdg_hamiltonian.f90      ! Nambu-space Hamiltonian assembly (topological SC only)
  magnetic_field.f90        ! Zeeman + Peierls COO insertions (QHE + topological SC)
  topological_analysis.f90  ! Unified: Chern number, Z2 (gap + Fu-Kane), edge states, Hall conductance
  green_functions.f90       ! Sparse LDOS / spectral function via Pardiso
src/apps/
  main_topology.f90         ! New executable: unified topological analysis (all three modes)
```

### Module dependency graph

```
defs.f90
  <- magnetic_field.f90           (Zeeman 8x8 + Peierls correction — used by QHE + topological SC)
       <- bdg_hamiltonian.f90     (Nambu doubling — used only for topological SC mode)
  <- topological_analysis.f90     (unified: Chern, Z2 gap, Z2 Fu-Kane, edge states — used by all modes)
       <- green_functions.f90     (LDOS via Pardiso — used by all modes)
```

Key principle: **bdg_hamiltonian.f90 wraps the existing Hamiltonian builder, it does not replace it.** All current capabilities (strain, SC loop, optics) remain untouched.

### New derived types (defs.f90)

```fortran
type :: bdg_config
  logical  :: enabled = .false.
  real(dp) :: mu = 0.0_dp         ! chemical potential (eV)
  real(dp) :: delta_0 = 0.0_dp    ! s-wave pairing gap (eV)
  real(dp) :: B_vec(3) = 0.0_dp   ! magnetic field Bx, By, Bz (Tesla)
  character(len=20) :: gauge = 'landau_x'  ! gauge choice
  real(dp) :: B_sweep(3) = 0.0_dp ! B-field sweep: min, max, step
  logical  :: self_consistent = .false.  ! future: self-consistent gap
end type

type :: topology_config
  logical  :: enabled = .false.
  character(len=20) :: mode = 'qhe'  ! qhe | qshe | bdg
  ! Berry curvature / Chern number
  logical  :: compute_chern = .false.
  logical  :: compute_hall = .false.  ! output σ_xy = Ce²/h
  ! Z2 invariant
  logical  :: compute_z2 = .false.
  character(len=20) :: z2_method = 'auto'  ! auto | gap | fukane
  ! Edge states
  logical  :: extract_edge_states = .false.
  real(dp) :: edge_E_window = 0.01_dp  ! energy window for edge state detection
  ! LDOS
  logical  :: compute_ldos = .false.
  real(dp) :: ldos_eta = 0.001_dp  ! Lorentzian broadening
  real(dp) :: ldos_E_range(2) = [-0.1_dp, 0.1_dp]
  integer  :: ldos_num_E = 200
end type

type :: topological_result
  integer  :: chern_number = 0
  integer  :: z2_invariant = 0
  real(dp) :: hall_conductance = 0.0_dp  ! in units of e²/h
  real(dp) :: min_gap = 0.0_dp
  real(dp) :: edge_xi = 0.0_dp       ! edge state localization length
  real(dp), allocatable :: edge_energies(:)
  real(dp), allocatable :: phase_boundary(:,:)  ! (B, mu) pairs
  real(dp), allocatable :: berry_curvature(:,:) ! Ω(kx, ky) if computed
end type
```

## Component Designs

### 1. bdg_hamiltonian.f90

Central module that assembles the 16N x 16N BdG matrix.

**Public interfaces:**

```fortran
! Wire (2D confinement, sparse CSR)
subroutine build_bdg_hamiltonian_1d(H_bdg_csr, cfg, profile_2d, kpterms_2d, &
                                     mu, delta_0, ws)
  type(csr_matrix), intent(out) :: H_bdg_csr
  type(simulation_config), intent(in) :: cfg
  real(dp), intent(in) :: profile_2d(:,:)
  type(csr_matrix), intent(in) :: kpterms_2d(:)
  real(dp), intent(in) :: mu, delta_0
  type(wire_workspace), intent(inout) :: ws
end subroutine

! QW (1D confinement, dense)
subroutine build_bdg_hamiltonian_qw(H_bdg, cfg, profile, kpterms, &
                                     kpar, mu, delta_0)
  complex(dp), allocatable, intent(out) :: H_bdg(:,:)
  type(simulation_config), intent(in) :: cfg
  real(dp), intent(in) :: profile(:,:)
  real(dp), intent(in) :: kpterms(:,:,:)
  real(dp), intent(in) :: kpar, mu, delta_0
end subroutine
```

**Assembly algorithm (wire):**

1. Call `ZB8bandGeneralized` to build H₀ (8Ngrid x 8Ngrid CSR)
2. Add Zeeman COO corrections via `add_zeeman_coo`
3. Add Peierls COO corrections via `add_peierls_coo`
4. Shift diagonal by -μ: `H₀(i,i) → H₀(i,i) - μ`
5. Build the four BdG blocks:
   - Block (1,1): H₀ - μI → copy H₀ COO entries, offset row/col by 0
   - Block (2,2): -H₀ᵀ + μI → transpose COO entries, flip sign, offset by 8Ngrid
   - Block (1,2): Δ → diagonal COO entries at (i, i+8Ngrid) with pairing matrix values
   - Block (2,1): Δ† → conjugate of block (1,2)
6. Merge into single COO array, convert to CSR

**Workspace extension:** The `wire_workspace` type needs BdG-aware fields:
- `bdg_coo_cache` for the doubled matrix COO-to-CSR mapping
- Pre-allocated arrays for the four Nambu blocks

### 2. magnetic_field.f90

Adds magnetic field terms to the Hamiltonian via COO insertion.

**Public interfaces:**

```fortran
subroutine add_zeeman_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                           grid, B_vec, params)
  ! Adds g*μ_B·B·σ to COO arrays at each grid point
  ! Uses 8x8 spin matrices (SIGMA_X/Y/Z) projected onto the 8-band basis
  ! Position-dependent g-factor from material parameters
end subroutine

subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                            grid, B_vec, gauge, kpterms_2d)
  ! Adds k → k - eA/ℏ corrections to the kz-dependent k.p terms
  ! Landau gauge: A = B × r (specific form depends on gauge choice)
  ! Only affects cross-derivative and kz-coupling blocks
end subroutine
```

**Zeeman term details:**
- The 8x8 spin matrices in the zinc-blende basis already encode the J=3/2 and J=1/2 spin structure
- Material-dependent g*: InSb (g ≈ -50), InAs (g ≈ -15), GaAs (g ≈ -0.4), Al (g ≈ 2)
- For the wire, `params(grid%material_id(i))%gfactor` gives position-dependent g*

**Peierls term details:**
- For wire along z with B = Bx x̂, Landau gauge A = (0, 0, Bx·y)
- Correction: kz → kz - eBx·y/ℏ at each grid point
- Affects terms with kz in the Hamiltonian (Q, T, S blocks)
- Implemented as a position-dependent additive correction to the existing COO values
- Gauge choices: `landau_x`, `landau_y`, `symmetric` (config option, start with `landau_x`)

### 3. topological_analysis.f90

Unified module for topological classification across all three modes. Computes Chern numbers, Z2 invariants (by two methods), Berry curvature, edge state profiles, and Hall conductance.

**Public interfaces:**

```fortran
! --- Berry curvature and Chern number ---

subroutine compute_berry_curvature(evecs_k, kx_arr, ky_arr, n_occ, omega_k)
  ! Discretized Berry curvature Ω(kx,ky) via Fukui-Hatsugai-Suzuki method
  ! Works for any Hamiltonian (bare k.p or BdG)
  complex(dp), intent(in) :: evecs_k(:,:,:,:)  ! eigenvectors(kx,ky,band,grid)
  real(dp), intent(in) :: kx_arr(:), ky_arr(:)
  integer, intent(in) :: n_occ  ! number of occupied bands
  real(dp), allocatable, intent(out) :: omega_k(:,:)  ! Ω(kx,ky)
end subroutine

function compute_chern_number(evecs_k, kx_arr, ky_arr, n_occ) result(C)
  ! Chern number C = (1/2π) ∫ Ω dk
  ! Integer by construction (Fukui lattice gauge)
  complex(dp), intent(in) :: evecs_k(:,:,:,:)
  real(dp), intent(in) :: kx_arr(:), ky_arr(:)
  integer, intent(in) :: n_occ
  integer :: C
end function

function compute_hall_conductance(C) result(sigma_xy)
  ! σ_xy = C · e²/h  (in units of e²/h)
  integer, intent(in) :: C
  real(dp) :: sigma_xy
end function

function compute_longitudinal_conductance(evals_landau, T, mu) result(sigma_zz)
  ! Longitudinal conductivity σ_zz from Landau level spectrum
  ! Uses Boltzmann transport: σ_zz ∝ Σ_n v_z,n² · (-df/dE) at E_n
  ! Chiral anomaly signature: σ_zz increases with B (negative magnetoresistance)
  real(dp), intent(in) :: evals_landau(:,:)  ! Landau levels (B, n)
  real(dp), intent(in) :: T, mu
  real(dp) :: sigma_zz
end function

! --- Z2 invariants ---

function compute_z2_gap(evals_kz, kz_arr) result(z2)
  ! Wire / 1D: gap-based Z2 invariant (topological SC mode)
  ! Gap closing/reopening at E=0 signals topological phase transition
  real(dp), intent(in) :: evals_kz(:,:)   ! eigenvalues(kz, band)
  real(dp), intent(in) :: kz_arr(:)
  integer :: z2
end function

function compute_z2_fukane(evals_trim, evecs_trim, parity_op) result(z2)
  ! QW / 2D: Fu-Kane parity-based Z2 invariant (QSHE mode)
  ! (-1)^ν = Π_i Π_{m∈occ} ξ_m(Γ_i)  at 4 TRIM points
  ! Requires inversion-symmetric system
  real(dp), intent(in) :: evals_trim(:,:,:)     ! evals(TRIM_point, band, spin)
  complex(dp), intent(in) :: evecs_trim(:,:,:,:) ! evecs(TRIM, band, grid)
  integer, intent(in) :: parity_op(:)            ! 8-band parity eigenvalues
  integer :: z2
end function

! --- Edge state extraction ---

subroutine extract_edge_states(evals, evecs, grid, E_window, edge_profiles)
  ! General edge state extraction for any mode
  ! Identifies states within E_window of gap edges
  ! Computes spatial profile and localization length
  real(dp), intent(in) :: evals(:,:)
  complex(dp), intent(in) :: evecs(:,:,:)
  type(spatial_grid), intent(in) :: grid
  real(dp), intent(in) :: E_window
  type(edge_state), allocatable, intent(out) :: edge_profiles(:)
end subroutine

subroutine compute_majorana_profile(evec_zero, grid, psi_e, psi_h, xi)
  ! Specialized: extract Majorana wavefunction from BdG zero-energy eigenstate
  ! Splits into electron (psi_e) and hole (psi_h) components
  ! Fits exponential decay to get localization length xi
  complex(dp), intent(in) :: evec_zero(:)
  type(spatial_grid), intent(in) :: grid
  real(dp), allocatable, intent(out) :: psi_e(:), psi_h(:)
  real(dp), intent(out) :: xi
end subroutine
```

**Chern number (Fukui method) — used by QHE and Chern insulator modes:**
1. Define U-link variables on the k-mesh from occupied-state overlaps
2. Compute lattice Berry curvature F(k) = ln(U_x · U_y · U_x† · U_y†)
3. Chern number C = (1/2π) Σ_k F(k)
4. Integer by construction, no gauge-fixing needed
5. Hall conductance: σ_xy = C · e²/h

**Z2 Fu-Kane (parity method) — used by QSHE mode:**
1. Diagonalize H(k) at the 4 TRIM points: Γ=(0,0), M₁=(π/a,0), M₂=(0,π/a), M₃=(π/a,π/a)
2. Compute parity eigenvalues ξ_m of each occupied band at each TRIM
3. The 8-band parity operator in the zinc-blende basis: P = diag(+1, -1, -1, +1, -1, +1, +1, -1) for the valence/conduction band structure
4. (-1)^ν = Π_i Π_{m∈occ} ξ_m(Γ_i)
5. ν = 0 → trivial, ν = 1 → topological (Z2 TI)
6. Only requires 4 diagonalizations — extremely efficient

**Z2 gap method — used by topological SC mode (1D wire):**
1. For each kz, compute min|E| across all bands
2. The minimum gap = min_kz(min|E|)
3. As B (or μ) is swept, track gap closings
4. Each gap closing flips the Z2 invariant
5. More robust than Pfaffian computation for large sparse matrices

**Edge state extraction — unified across all modes:**
1. After computing the band structure, identify states within the gap
2. Compute spatial density |ψ(r)|² for each such state
3. Check if the density is localized at the boundary (edge state) vs. bulk
4. Fit exponential decay: |ψ(r)|² ~ exp(-2r/ξ) to get localization length ξ
5. For BdG mode: additionally split into electron/hole components

### 4. green_functions.f90

LDOS computation via Pardiso sparse LU factorization.

**Public interfaces:**

```fortran
subroutine compute_ldos(H_bdg_csr, E_grid, eta, grid, ldos_total, &
                         ldos_electron, ldos_hole)
  ! LDOS at each spatial point for each energy
  ! Uses Pardiso factorization of (E+iη - H_BdG)
  type(csr_matrix), intent(in) :: H_bdg_csr
  real(dp), intent(in) :: E_grid(:), eta
  type(spatial_grid), intent(in) :: grid
  real(dp), allocatable, intent(out) :: ldos_total(:,:)
  real(dp), allocatable, intent(out) :: ldos_electron(:,:)
  real(dp), allocatable, intent(out) :: ldos_hole(:,:)
end subroutine

subroutine compute_spectral_function(H_bdg_csr, E_arr, k_arr, eta, A_kE)
  ! k-resolved spectral function A(k,E) for ARPES-like visualization
  type(csr_matrix), intent(in) :: H_bdg_csr
  real(dp), intent(in) :: E_arr(:), k_arr(:), eta
  real(dp), allocatable, intent(out) :: A_kE(:,:)
end subroutine
```

**LDOS algorithm (Pardiso):**
1. For each energy E:
   a. Build shifted matrix: A = (E + iη)I - H_BdG
   b. Phase 1 (symbolic factorization) — done once
   c. Phase 2 (numerical factorization) — done per E if shift changes
   d. Phase 3 (solve): extract diagonal blocks of G = A⁻¹ by solving with unit sources
   e. LDOS(r) = -(1/π) Im[G(r,r)]
2. Split into electron/hole components from the Nambu structure (upper/lower 8N blocks)

**Optimization:** MKL Pardiso supports shift-invert mode. Factorize H_BdG once, then solve (E+iη - H)⁻¹ for multiple E with iterative refinement. This is much faster than full refactorization at each energy.

### 5. main_topology.f90

New executable `topologicalAnalysis` supporting all three modes.

**Mode selection** — determined by input config:
- `mode = qhe`: bare k.p + B-field → Landau levels + Chern number + Hall conductance
- `mode = qshe`: bare k.p, no B-field → Fu-Kane Z2 + helical edge states
- `mode = bdg`: BdG doubled + B-field → Majorana modes + gap-based Z2

**Computational workflow (unified):**

```
1. Parse input.cfg (extended with topology: and bdg: blocks)
2. Initialize grid, materials, kpterms (existing code)
3. Build Hamiltonian depending on mode:
   - QHE:  H₀(k) + Zeeman + Peierls  (8N)
   - QSHE: H₀(k)                      (8N)
   - BdG:  BdG[ H₀(k) + Zeeman + Peierls, μ, Δ ]  (16N)
4. k-space sweep:
   a. QHE:  2D mesh (kx, ky) for Chern number
   b. QSHE: TRIM points (4 k-points) for Fu-Kane + full (kx,ky) for edge states
   c. BdG:  1D kz sweep for gap-based Z2
5. Diagonalize with FEAST at each k-point
6. Post-process depending on mode:
   a. QHE:  Berry curvature → Chern number → Hall conductance σ_xy(B)
   b. QSHE: Fu-Kane Z2 invariant + edge state extraction
   c. BdG:  Z2 from gap closings + Majorana wavefunction extraction
7. Optional: LDOS computation at selected parameter points
8. Output files
```

**Eigensolver:** FEAST for both wire (CSR) and QW (dense). Search interval [-δE, +δE] around E=0 where δE is configurable (default: 0.1 eV). This avoids computing the many high-energy states (deep VB hole replicas in the BdG spectrum) that are physically irrelevant.

### 6. Input extension (input_parser.f90)

New `topology:` block for topological analysis (all modes):

```
topology: T                        ! enable topological analysis
topology_mode qhe                  ! qhe | qshe | bdg
topology_B 0.0 0.0 5.0             ! magnetic field Bx By Bz (Tesla)
topology_B_sweep 0.0 10.0 0.1     ! optional B sweep: Bmin Bmax Bstep
topology_gauge landau_x            ! gauge: landau_x | landau_y | symmetric
topology_compute_chern T           ! compute Chern number (QHE mode)
topology_compute_hall T            ! output Hall conductance (QHE mode)
topology_compute_z2 T              ! compute Z2 invariant (QSHE/BdG mode)
topology_z2_method auto            ! auto | fukane | gap
topology_extract_edges T           ! extract edge state profiles
topology_edge_E_window 0.01       ! energy window for edge detection (eV)
topology_compute_ldos F            ! enable LDOS computation
topology_ldos_eta 0.001            ! Lorentzian broadening (eV)
topology_ldos_E_range -0.1 0.1    ! LDOS energy range (eV)
topology_ldos_num_E 200            ! number of energy points for LDOS
```

BdG-specific parameters (only used when `topology_mode = bdg`):

```
bdg: T                             ! enable BdG mode (required if topology_mode=bdg)
bdg_mu 0.010                       ! chemical potential (eV)
bdg_delta 0.0003                   ! s-wave pairing gap (eV)
bdg_E_window 0.1                   ! FEAST search half-width around E=0 (eV)
```

## Eigensolver Strategy

| Geometry | Matrix format | FEAST interface | Notes |
|---|---|---|---|
| QW (1D) | Dense 16N x 16N | FEAST dense driver | Only states near E=0, avoids computing full spectrum |
| Wire (2D) | CSR 16Ngrid x 16Ngrid | FEAST CSR driver | Same as current wire workflow, doubled dimension |

FEAST contour interval set to [-bdg_E_window, +bdg_E_window]. For Majorana detection, a narrow window (±0.1 eV) captures the relevant states efficiently.

## Output Files

All output to `output/` directory, consistent with existing convention:

```
output/
  topo_spectrum.dat           ! E vs k for each B value (all modes)
  topo_phase_diagram.dat      ! B, mu, gap_size, invariant (Z2 or C)
  topo_edge_states.dat        ! spatial profile of edge states
  topo_berry_curvature.dat    ! Ω(kx, ky) if computed (QHE mode)
  topo_hall_conductance.dat   ! σ_xy vs B (QHE mode)
  topo_longitudinal_cond.dat  ! σ_zz vs B — chiral anomaly signature
  topo_ldos.dat               ! LDOS(E, r) — energy-resolved
  topo_spectral.dat           ! A(k, E) — k-resolved spectral function
```

Columnar text format for gnuplot visualization.

## Implementation Phases

### Phase 1: Magnetic field module (foundation for QHE + topological SC)
- `magnetic_field.f90`: Zeeman term + Peierls corrections (Landau gauge)
- `topology_config` type in defs.f90
- `topology:` block parsing in input_parser.f90
- Tests: unit tests for Zeeman COO insertion, Landau level verification in bulk

### Phase 2: Berry curvature and Chern number (QHE mode)
- `topological_analysis.f90`: Berry curvature (Fukui method), Chern number, Hall conductance
- `main_topology.f90`: QHE mode — Landau level spectrum + σ_xy(B)
- Tests: Chern number = 1 for lowest Landau level, known Hall plateaus

### Phase 3: Fu-Kane Z2 invariant (QSHE mode)
- `topological_analysis.f90`: Fu-Kane parity-based Z2, 8-band parity operator
- Edge state extraction for helical edge states
- `main_topology.f90`: QSHE mode — Z2 invariant + edge state profiles
- Tests: reproduce BHZ phase transition in HgTe/CdTe QW (Z2 flip at critical well width)

### Phase 4: BdG Hamiltonian (topological SC mode)
- `bdg_hamiltonian.f90`: Nambu doubling, μ shift, pairing matrix, workspace extension
- `bdg_config` type in defs.f90
- `main_topology.f90`: BdG mode — spectrum + gap-based Z2
- Tests: particle-hole symmetry verification, BCS gap in bulk limit

### Phase 5: Majorana modes and LDOS
- `topological_analysis.f90`: Majorana wavefunction extraction
- `green_functions.f90`: Pardiso-based LDOS and spectral function
- Phase diagram computation (Z2 vs B, μ)
- Tests: reproduce Lutchyn/Oreg phase boundary, Majorana localization length

### Phase 6 (future): Self-consistent gap
- Extend SC loop (DIIS) to self-consistent Δ equation
- Coupling constant V for the pairing interaction
- Proximity self-energy at the SC/semiconductor interface

## Verification Strategy

1. **Landau levels (QHE)**: With Peierls orbital effects, bulk InAs spectrum should show Landau level quantization En = ℏωc(n + 1/2) with ωc = eB/m*. Check spacing and degeneracy.

2. **Chern number = 1**: For the lowest occupied Landau level, the Chern number must be exactly 1. This is a strong numerical test of the Fukui method.

3. **Hall plateau structure**: As B is swept, σ_xy should show integer plateaus at multiples of e²/h with transitions at Landau level crossings.

4. **HgTe/CdTe band inversion (QSHE)**: The BHZ model predicts a topological phase transition at a critical well width d_c ≈ 6.3 nm. The 8-band code should reproduce this Z2 flip.

5. **Fu-Kane Z2**: For a trivial QW (e.g., GaAs/AlGaAs), Z2 = 0. For HgTe/CdTe above d_c, Z2 = 1. Test both cases.

6. **Helical edge states**: In the Z2 = 1 phase, edge states should cross the gap with counter-propagating spin-momentum locking. Verify spatial localization at the QW boundary.

7. **Bulk BdG limit**: For homogeneous material with uniform Δ, BdG spectrum should show the BCS gap 2Δ around the Fermi level. Analytically known.

8. **Particle-hole symmetry**: Verify H_BdG = -C·H_BdG·C⁻¹ where C is the particle-hole operator. Check numerically for assembled matrices.

9. **Topological SC phase boundary**: For a simple Rashba wire (2-band model), the topological transition occurs at Vz² = μ² + Δ² where Vz = gμ_BB/2. The 8-band code should reproduce this in the limit of weak spin-orbit coupling.

10. **Majorana localization**: In the topological phase, the zero-energy state should decay exponentially with length ξ = ℏv_F/Δ for a clean wire. Check the localization length against this prediction.

## Scope Boundaries

- **In scope**: s-wave pairing, Zeeman + Peierls, Landau levels, Chern number, Hall conductance, Fu-Kane Z2, helical edge states, BdG spectrum, gap-based Z2, Majorana modes, LDOS
- **Out of scope** (future work): p-wave pairing, Floquet topological SC, non-equilibrium (Keldysh) formalism, disorder, multi-terminal Josephson junctions, quantum anomalous Hall effect (requires magnetism), fractional QHE
- **Not modified**: Existing hamiltonianConstructor.f90, optical_spectra.f90, gfactor_functions.f90 — all current capabilities preserved
