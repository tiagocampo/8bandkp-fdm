# Quantum Wire (2D Confinement) Design Document

**Date**: 2026-04-02
**Author**: Tiago de Campos
**Status**: In Progress (Phases 0-2 complete)

## 1. Overview

Extend the 8-band zincblende k.p FDM code from 1D confinement (quantum wells) to 2D confinement (quantum wires). The wire axis is along x (free propagation direction k_x), with confinement in y and z. The Hamiltonian generalizes from 8N x 8N to 8 * Ny * Nz x 8 * Ny * Nz.

### Scope

- Electronic subband structure E_n(k_x) for quantum wires
- g-factor calculation for wire subbands
- Optical matrix elements (inter-subband transitions)
- Self-consistent Schrodinger-Poisson for doped wires
- Strain from lattice mismatch (continuum elasticity + Pikus-Bir)
- Arbitrary wire cross-sections (hexagonal, circular, rectangular, core-shell) via cut-cell immersed boundary on Cartesian FD grid
- Dimension-agnostic refactor: same infrastructure serves bulk, QW, and wire
- FEAST and ARPACK sparse eigensolvers

### Key References

1. **Stier & Bimberg**, "Modeling of strained quantum wires using eight-band k.p theory," PRB 55, 7726 (1997). DOI: 10.1103/PhysRevB.55.7726
2. **Bahder**, "Eight-band k.p model of strained zinc-blende crystals," PRB 41, 11992 (1990). DOI: 10.1103/PhysRevB.41.11992
3. **Foreman**, "Elimination of spurious solutions from eight-band k.p theory," PRB 56, R12748 (1997). DOI: 10.1103/PhysRevB.56.R12748
4. **Veprek, Steiger, Witzigmann**, "Reliable k.p band structure calculation for nanostructures using finite elements," J. Comput. Electron. 7, 521-529 (2008).
5. **Lew Yan Voon & Willatzen**, "The k.p Method: Electronic Properties of Semiconductors," Springer (2009). ISBN: 978-3-540-92871-3
6. **Winkler**, "Spin-Orbit Coupling Effects in Two-Dimensional Electron and Hole Systems," Springer (2003). DOI: 10.1007/b13586
7. **Trellakis et al.**, "The 3D nanometer device project nextnano," J. Comput. Electron. 5, 285-289 (2006).
8. **Birner**, "Modeling of semiconductor nanostructures and nextnano3," PhD thesis, TU Munich (2011).

## 2. Architecture

### 2.1 Module structure

```
src/
  core/
    defs.f90              -- extended: spatial_grid type, confinement=2, strain_config
    parameters.f90        -- extended: elastic constants, deformation potentials, lattice constants
  math/
    finitedifferences.f90 -- refactored: dim-agnostic FD operators returning sparse CSR
    sparse_matrices.f90   -- NEW: CSR builder, Kronecker products, cut-cell geometry
    geometry.f90          -- NEW: shape definitions, cut-cell volume/face fractions
  physics/
    hamiltonianConstructor.f90 -- refactored: ZB8bandGeneralized (dim-agnostic)
    strain_solver.f90     -- NEW: continuum elasticity (plane strain) + Pikus-Bir
    gfactor_functions.f90 -- extended: works with sparse Hamiltonian for wire subbands
    poisson.f90           -- extended: 2D Poisson (PARDISO or CG)
    charge_density.f90    -- extended: 2D charge density n(y,z), p(y,z)
    sc_loop.f90           -- extended: works with 2D Poisson + charge
  io/
    input_parser.f90      -- extended: wire geometry, strain params, 2D doping
    outputFunctions.f90   -- extended: 2D eigenfunction output
  apps/
    main.f90              -- unified: handles bulk/QW/wire via same path
    main_gfactor.f90      -- extended: g-factor for wire subbands
tests/
  unit/                   -- new: test_sparse_matrices, test_fd_2d, test_cutcell,
                           --       test_strain_solver, test_eigensolver, test_hamiltonian_2d
  integration/            -- new: wire rectangle, hexagon, strain, SC
  regression/             -- new: wire reference configs and data
```

### 2.2 Key new type: `spatial_grid`

Replaces scattered `z(:)`, `delta`, `dz` fields in `simulation_config` with a unified grid object:

```fortran
type spatial_grid
  integer :: ndim              ! 0=bulk, 1=QW, 2=wire
  integer :: ny, nz            ! grid points per confinement direction
  real(dp) :: dy, dz           ! grid spacings (Angstrom)
  real(dp), allocatable :: y(:), z(:)                    ! 1D coordinate arrays
  real(dp), allocatable :: coords(:,:)                   ! (2, ny*nz) flattened 2D grid
  integer, allocatable  :: material_id(:)                ! material at each point (ny*nz)
  ! Cut-cell immersed boundary fields:
  real(dp), allocatable :: cell_volume(:)                ! fractional volume inside domain [0,1]
  real(dp), allocatable :: face_fraction_y(:,:)          ! (ny*nz, 2) left/right y-face fractions
  real(dp), allocatable :: face_fraction_z(:,:)          ! (ny*nz, 2) bottom/top z-face fractions
  integer, allocatable  :: ghost_map(:,:)                ! (ny*nz, 4) nearest active neighbor indices
end type
```

For QW (confinement=1): `ny=1, nz=fdStep`, cut-cell fields all = 1.0. For bulk (confinement=0): `ny=nz=1`.

### 2.3 Module dependency graph (updated)

```
defs.f90
  <- parameters.f90           (extended: elastic constants, deformation potentials)
  <- sparse_matrices.f90      (NEW: CSR type, Kronecker products)
       <- geometry.f90        (NEW: shape definitions, cut-cell fractions)
  <- finitedifferences.f90    (refactored: dim-agnostic, returns sparse CSR)
       <- hamiltonianConstructor.f90  (refactored: ZB8bandGeneralized)
            <- gfactor_functions.f90  (extended: sparse g-factor)
       <- strain_solver.f90          (NEW: plane-strain elasticity + Pikus-Bir)
       <- poisson.f90                (extended: 2D Poisson)
       <- charge_density.f90         (extended: 2D charge density)
            <- sc_loop.f90           (extended: 2D SC loop)
  <- eigensolver.f90          (NEW: FEAST + ARPACK wrapper)
  <- outputFunctions.f90      (extended: 2D output)
  <- input_parser.f90         (extended: wire geometry, strain, solver config)
```

## 3. 2D Finite Difference Operators

### 3.1 Kronecker product decomposition

For a 2D grid with Ny x Nz points, operators are built as Kronecker products of 1D operators:

| Operator | Kronecker decomposition | Size |
|---|---|---|
| d^2/dy^2 | D2y (x) I_z | (Ny*Nz) x (Ny*Nz) |
| d^2/dz^2 | I_y (x) D2z | (Ny*Nz) x (Ny*Nz) |
| d/dy | D1y (x) I_z | (Ny*Nz) x (Ny*Nz) |
| d/dz | I_y (x) D1z | (Ny*Nz) x (Ny*Nz) |
| d^2/dydz | D1y (x) D1z | (Ny*Nz) x (Ny*Nz) |

Where D2y = buildFD2ndDerivMatrix(Ny, dy, order), I_z = identity(Nz).

Existing 1D FD stencil machinery from finitedifferences.f90 (orders 2-10) is reused. The new code only adds Kronecker product assembly and sparse storage.

### 3.2 Cross-derivative stencil (new for 2D)

The term D1y (x) D1z creates a 4-corner stencil at each grid point (2nd-order):

```
  (-1,-1)    (0,-1)    (+1,-1)
             |
  (-1, 0)  [center]  (+1, 0)
             |
  (-1,+1)    (0,+1)    (+1,+1)

Nonzeros from D1y (x) D1z: corners only -> 4 entries
```

This arises from the R coupling term: gamma_3 * k_y * k_z becomes gamma_3(y,z) * (d/dy)(d/dz).

### 3.3 Sparse CSR storage

New module `sparse_matrices.f90`:

```fortran
type :: csr_matrix
  integer :: nrows, ncols, nnz
  real(dp), allocatable :: values(:)      ! length nnz (complex stored as pairs or separate)
  integer, allocatable  :: colind(:)      ! length nnz
  integer, allocatable  :: rowptr(:)      ! length nrows+1
end type

interface kronecker_product
  module procedure kron_dense_dense    ! two dense -> sparse CSR
  module procedure kron_sparse_dense   ! sparse CSR (x) dense identity
  module procedure kron_dense_sparse   ! dense identity (x) sparse CSR
end interface
```

Key routines:
- `kron_dense_dense`: Kronecker product, output direct to COO then CSR. Skips zero entries.
- `apply_variable_coeff_2d`: Applies position-dependent coefficient profile to sparse FD operator, with cut-cell face fraction weighting.
- `sum_csr`: Sparse matrix addition for combining y and z Kronecker contributions.

### 3.4 Cut-cell boundary treatment

For points near the wire boundary, FD stencil weights are modified using face fractions from `spatial_grid`. Box-integration Laplacian:

```
Lap phi ~ (1/V) * [ (phi_{i+1} - phi_i)*fy_right/dy - (phi_i - phi_{i-1})*fy_left/dy
                  + (phi_{j+1} - phi_j)*fz_top/dz   - (phi_j - phi_{j-1})*fz_bot/dz ]
```

Where fy_right, fy_left, fz_top, fz_bot are face fractions in [0,1]. Interior points: all = 1.0 (standard stencil). Edge points: reduced fractions give natural BC enforcement.

This is the same box-integration discretization used by nextnano.

### 3.5 kpterms generalization

Currently `kpterms(N, N, 10)` holds ten NxN dense matrices. For 2D, ten CSR sparse matrices of size (Ny*Nz) x (Ny*Nz):

| Index | 1D (dense) | 2D (sparse CSR) |
|---|---|---|
| 1 | gamma_1 diagonal | gamma_1 diagonal (same structure) |
| 2 | gamma_2 diagonal | gamma_2 diagonal |
| 3 | gamma_3 diagonal | gamma_3 diagonal |
| 4 | P diagonal | P diagonal |
| 5 | A * d^2/dz^2 | A * (d^2/dy^2 + d^2/dz^2) via Kronecker sum |
| 6 | P * d/dz | P * d/dy + P * d/dz (both directions) |
| 7 | (gamma_1-2*gamma_2) * d^2/dz^2 | (gamma_1-2*gamma_2) * (d^2/dy^2 + d^2/dz^2) |
| 8 | (gamma_1+2*gamma_2) * d^2/dz^2 | (gamma_1+2*gamma_2) * (d^2/dy^2 + d^2/dz^2) |
| 9 | gamma_3 * d/dz | gamma_3 * d/dy + gamma_3 * d/dz |
| 10 | A diagonal | A diagonal |
| **11** (NEW) | N/A | gamma_3 * d^2/dydz (cross-derivative) |

For 1D QW (ny=1), Kronecker products collapse to existing 1D dense operators.

### 3.6 Memory and performance estimates

For a 100x100 grid, 2nd-order FD:

| Component | Estimate |
|---|---|
| Total DOF | 8 * 10,000 = 80,000 |
| NNZ per spatial operator | ~50,000 (5-point Laplacian) |
| Total Hamiltonian NNZ | ~4M (8 bands * cross-coupling) |
| CSR memory | ~100 MB |
| 1D FD matrix build | negligible (100x100 dense) |
| Kronecker assembly | ~0.1s |
| Eigensolve (FEAST, 20 eigenvalues) | ~5-30s |

## 4. Hamiltonian Assembly

### 4.1 8x8 block structure is dimension-independent

The zincblende 8x8 block topology is fixed. What changes with dimensionality is the internal structure of each block (scalar / dense / sparse):

```
        HH1  HH2  LH1  LH2  SO1  SO2   CB1  CB2
HH1  [  P+Q         R     -S           PZ    PM  ]
HH2  [  R*    P-Q         S           PZ    PM  ]
LH1  [ -S*    S*   P+Q        -R      PZ    PP  ]
LH2  [              -R*  P-Q           PZ    PP  ]
SO1  [  PZ*   PZ*   PZ*  PZ*  A+k^2         PZ'      ]
SO2  [  PM*   PM*   PP*  PP*             A+k^2  PZ'  ]
CB1  [                                Ec           ]
CB2  [                                      Ec      ]
```

### 4.2 Term definitions per confinement mode

| Term | Bulk (0) | QW (1) | Wire (2) |
|---|---|---|---|
| Q | scalar | dense NxN | sparse: -(gamma_1+gamma_2)*kx^2*I + kpterms_2d(7) |
| T | scalar | dense NxN | sparse: same pattern with (gamma_1+2*gamma_2) |
| R | scalar | scalar | sparse: gamma_2*kx^2*I + i*gamma_3*kx*(d/dy) + kpterms_2d(11) |
| S | scalar | dense (d/dz) | sparse: d/dy + d/dz contributions |
| PZ | scalar | dense (d/dz) | sparse: d/dy + d/dz contributions |
| PP/PM | scalar | scalar | scalar (kx only) |
| A | scalar | dense | sparse: Laplacian with A coefficient |

### 4.3 Assembly: COO insertion strategy

Two-phase build (same as existing g-factor code):

1. **Append phase**: For each nonzero 8x8 block (alpha, beta), compute Ngrid x Ngrid sparse block, append nonzeros as (row_global, col_global, value) triplets where row_global = (alpha-1)*Ngrid + block_row.

2. **Finalize phase**: Sort by (row, col), merge duplicates, build CSR rowptr.

The 8x8 topology has ~20 nonzero upper-triangular blocks. Each block is ~5*Ngrid nonzeros. Total COO entries: ~100*Ngrid.

### 4.4 Spurious solution mitigation

1. **Burt-Foreman operator ordering**: Asymmetric first-derivative discretization (forward/backward) per Foreman 1997.
2. **Cross-derivative symmetrization**: Split gamma_3*k_y*k_z as 0.5*(gamma_3*k_y*k_z + k_z*gamma_3*k_y) before discretization.
3. **Hermiticity verification**: After assembly, check max|H - H^dagger| < 1e-12.
4. **Energy filtering**: Discard eigenvalues outside physical range.

### 4.5 1D fast path

For small QW problems (N < dense_threshold), preserve existing ZHEEVX dense path:

```fortran
if (grid%ndim == 1 .and. N < dense_threshold) then
  call ZB8bandQW_dense(HT, kx, profile, kpterms, cfg)
  call zheevx(...)
else
  call ZB8bandGeneralized(HT_csr, kx, profile_2d, kpterms_2d, grid, cfg)
  call solve_sparse_evp(...)
end if
```

## 5. Strain Solver

### 5.1 Pipeline

```
Stage 1: Continuum elasticity        Stage 2: Pikus-Bir           Stage 3: Add to H
Solve for eps_ij(y,z) on the         Transform eps_ij into        Add delta_H_strain to
2D cross-section                     band edge shifts +           the 8x8 block Hamiltonian
                                      coupling terms               before eigenvalue solve
```

### 5.2 Plane-strain elasticity

For a wire along x with lattice mismatch, solve the plane-strain equations:

```
d(sigma_xx)/dy + d(sigma_xy)/dz = 0
d(sigma_xy)/dy + d(sigma_zz)/dz = 0
```

With anisotropic stress-strain (zincblende C11, C12, C44) and misfit strain epsilon_0. This is a system of 2 PDEs in 2 unknowns (displacement u_y, u_z) on the 2D cross-section.

Discretized with the same FD grid + cut-cell treatment as the Hamiltonian. Solved via MKL PARDISO (direct sparse). The stiffness matrix K is 2*Ngrid x 2*Ngrid sparse.

BCs: stress-free on outer boundary (natural in weak form), interface matching automatic from position-dependent C_ij.

### 5.3 New material parameters

Extended `paramStruct` in `parameters.f90`:

```fortran
! Elastic constants (Vurgaftman 2001 Table XII, GPa)
real(dp) :: C11, C12, C44
! Lattice constant (Angstrom)
real(dp) :: a0
! Deformation potentials (eV)
real(dp) :: ac    ! conduction band hydrostatic
real(dp) :: av    ! valence band hydrostatic
real(dp) :: b     ! shear (tetragonal)
real(dp) :: d     ! shear (rhombohedral)
real(dp) :: bp    ! CB shear coupling (Bahder a')
```

Sources: Vurgaftman 2001 (C_ij for 30+ materials), Bahder 1990 (deformation potentials), Winkler 2003.

### 5.4 Pikus-Bir strain Hamiltonian

Diagonal shifts (added to band offset profile):

```
delta_E_CB = ac * Tr(eps)
delta_E_HH = -P_eps + Q_eps
delta_E_LH = -P_eps - Q_eps
delta_E_SO = -P_eps

P_eps = -av * Tr(eps)
Q_eps = b/2 * (eps_zz - 0.5*(eps_yy + eps_xx))
```

Off-diagonal couplings (new sparse matrices added to Hamiltonian):

```
b * (eps_yy - eps_zz)   -> couples HH-LH (R block)
d * eps_yz               -> couples HH-LH (S block)
bp * (eps_yy + eps_zz)   -> small CB-SO coupling
```

### 5.5 1D QW strain: special case

For QW (confinement=1), strain reduces to simple biaxial algebra:

```
eps_xx = eps_yy = eps_0 = (a_well - a_barrier) / a_barrier
eps_zz = -2*C12/C11 * eps_0
```

No PDE solve needed. The strain solver returns this trivially when ndim=1. QW strain comes for free.

### 5.6 Piezoelectric effects (Phase 2, optional)

Zincblende: single piezoelectric constant e_14. Strain generates polarization P_y = 2*e_14*eps_yz. Creates internal field entering Poisson as source term. Can be added later without restructuring.

## 6. Eigenvalue Solvers

### 6.1 Unified interface

```fortran
type :: eigensolver_config
  character(len=10) :: method = 'FEAST'     ! 'FEAST' or 'ARPACK'
  integer  :: nev = 8
  real(dp) :: emin, emax                     ! search interval [FEAST] or shift [ARPACK]
  integer  :: max_iter = 100
  real(dp) :: tol = 1e-10_dp
  integer  :: feast_m0 = 0                   ! subspace estimate (0 = auto = 2*nev)
  integer  :: ncv = 0                        ! Lanczos vectors for ARPACK
  character(len=1) :: which = 'S'            ! 'S'=smallest, 'L'=largest
end type
```

### 6.2 FEAST (default)

MKL FEAST (`zfeast_hcsrev`) finds all eigenvalues in [Emin, Emax]. Best for wires where you want all subbands in a physical energy window.

Energy window heuristic:

```
emin_search = minval(profile(:,3)) - 0.1                    ! below lowest CB
emax_search = maxval(profile(:,3)) + 5 * hbar2_2m / dy^2   ! CB + kinetic cutoff
```

### 6.3 ARPACK (alternative)

MKL Extended Eigensolver (`mkl_sparse_z_ev`) wrapping ARPACK. Best when you know exactly how many states you want. Shift-and-invert mode for interior eigenvalues.

### 6.4 k_x sweep optimization

Sparsity pattern of H is independent of k_x. FEAST/ARPACK can reuse symbolic factorizations across k-points.

## 7. Self-Consistent 2D Loop

### 7.1 Structure (same as 1D, dimension-agnostic)

1. Solve Schrodinger: H(kx=0) psi = E psi
2. Compute charge density: n(y,z) = sum_n |psi_n|^2 * f(E_n, EF, T) (kx integral)
3. Solve Poisson: nabla^2 phi = -rho(y,z)/eps(y,z) (2D)
4. Update potential: V_new = mix(V_old, e*phi), DIIS acceleration
5. Converged? If not, go to 1.

### 7.2 2D Poisson solver

MKL PARDISO direct solve. Build nabla.(eps*nabla) with variable eps(y,z) and cut-cell BCs. Factorization reused across SC iterations (only RHS changes).

Optional FFT-based spectral Poisson for uniform-eps regular grids (FFTW, O(N log N)).

### 7.3 2D charge density

```
2D wire: n(y,z) = sum_n integral dk_x / (2pi) * |psi_n(y,z)|^2 * f(E_n(kx) - EF, T)
```

The k_|| integral (2D for QW) becomes kx integral (1D for wire) - simpler. Simpson's rule over computed k_x dispersion.

### 7.4 DIIS mixing

Dimension-agnostic - operates on flattened potential array V(:) of length Ngrid. No changes needed.

## 8. Input File Format

New fields are additive. Existing QW configs work unchanged.

### 8.1 Wire-specific fields

```ini
confinement = 2
wire_ny = 100
wire_nz = 100
wire_dy = 0.5
wire_dz = 0.5
wire_shape = hexagon               ! rectangle, circle, hexagon, polygon
wire_radius = 50.0
! OR polygon vertices:
wire_polygon = 6
wire_vertex = 0.0 50.0
wire_vertex = 43.3 25.0
...

numRegions = 2
region = core  InAs  0.0  30.0
region = shell GaSb 30.0  50.0
```

### 8.2 Strain fields

```ini
strain = .true.
strain_ref = substrate
strain_solver = pardiso
piezoelectric = .false.
```

### 8.3 Eigensolver fields

```ini
eigensolver = FEAST
eigenvalue_emin = -2.0
eigenvalue_emax = 2.0
eigenvalue_num = 20
```

### 8.4 Backward compatibility

| Field | bulk (0) | QW (1) | wire (2) |
|---|---|---|---|
| FDstep | ignored | N (z-grid) | ignored |
| wire_* | ignored | ignored | required |
| numLayers | 1 | any | ignored |
| numRegions | ignored | ignored | any |
| eigensolver | dense | dense or sparse | sparse only |

## 9. Testing Strategy

### Tier 1: Analytical validation (unit tests, pFUnit)

| Test | Validates | Reference |
|---|---|---|
| Free electron 2D box | 2D Laplacian eigenvalues | Exact: pi^2*hbar^2/2m*(ny^2/Ly^2 + nz^2/Lz^2) |
| 2D quantum harmonic oscillator | Kronecker product FD | Exact: hbar*omega*(2n+1) |
| Single-band wire | Subband energies vs. particle-in-a-box | Analytical |
| Cut-cell area fractions | Circle on grid: sum vs. pi*r^2 | Geometric |
| Cross-derivative stencil | d^2(xy)/dxdy = 1 | Exact |
| CSR Hermiticity | max|H - H^dagger| < 1e-12 | Numerical |
| Strain: uniform biaxial | QW algebraic vs. PDE solver | Analytical |
| Strain: hydrostatic pressure | Uniform eps -> correct delta_Ec, delta_Ev | Analytical |

### Tier 2: Regression against 1D results

| Test | Method |
|---|---|
| GaAs/AlGaAs QW, CB states | Set ny=1, compare against existing 1D results |
| GaSb/AlSb/InAs QW | Same |
| QW with strain | Compare strained vs. existing code |
| QW SC loop | DIIS convergence unchanged |

### Tier 3: Published benchmarks

| System | Reference | Expected |
|---|---|---|
| GaAs rectangular wire 10nm x 10nm | Effective mass | CB ground ~56 meV |
| InAs/GaSb core-shell r=10nm | Stier & Bimberg 1997 | Subband spacing |
| InAs wire g-factor | Winkler 2003 | g ~ -14.9 (bulk limit) |
| InAs/GaAs strained wire | k.p literature | Strain-modified band gap |

### Tier 4: Internal consistency

| Test | Checks |
|---|---|
| Convergence with grid | Halve dy,dz: eigenvalues converge at FD order rate |
| Convergence with FD order | FDorder 2->4->6: error decreases |
| Kramer's degeneracy | All eigenvalues doubly degenerate at kx=0 |
| Normalization | sum|psi_n(y,z)|^2 = 1 |
| Strain <-> no-strain | Lattice-matched: strain=0 gives same results |

### Test infrastructure

```
tests/
  unit/
    test_sparse_matrices.pf
    test_fd_2d.pf
    test_cutcell.pf
    test_strain_solver.pf
    test_eigensolver.pf
    test_hamiltonian_2d.pf
  integration/
    test_wire_rectangle.sh
    test_wire_hexagon.sh
    test_wire_strain.sh
    test_sc_wire.sh
  regression/
    configs/
      wire_gaas_rectangle.cfg
      wire_inas_gasb_coreshell.cfg
      wire_inas_hexagon.cfg
      wire_strained_inas_gaas.cfg
    data/
      wire_gaas_rectangle.ref
      wire_inas_gasb_coreshell.ref
      wire_inas_hexagon.ref
      wire_strained_inas_gaas.ref
    compare_output.py
```

## 10. Implementation Progress

### Completed: Phases 0-3 (23 commits on `feature/quantum-wire-2d`)

**Branch:** `feature/quantum-wire-2d`
**Tests:** 22/22 passing (12 unit + 10 regression)
**Files changed:** 18, ~3,700 lines added

#### Phase 0: Foundation Refactoring -- DONE

- `spatial_grid` type in `defs.f90` with ndim, nx/ny, dx/dy, material_id, cut-cell fields
- `csr_matrix` type + COO builder + merge sort in `sparse_matrices.f90`
- Kronecker product routines: `kron_dense_dense`, `kron_eye_dense`, `kron_dense_eye`
- CSR arithmetic: `csr_add`, `csr_scale`, `csr_apply_variable_coeff`
- `input_parser.f90` extended with wire geometry fields (wire_nx, wire_ny, wire_dx, wire_dy, wire_shape, wire_width, wire_height, numRegions)
- `confinementInitialization` accepts `spatial_grid` via overloaded interface

#### Phase 1: 2D FD + Sparse Assembly -- DONE

- 2D Laplacian via Kronecker sum: `Iy x D2x + D2y x Ix` (column-major flat index)
- 2D gradient via Kronecker sum: `Iy x D1x + D1y x Ix`
- Cross-derivative: `D1y x D1x`
- Cut-cell geometry engine in `geometry.f90` with Sutherland-Hodgman clipping for circle, hexagon, rectangle, polygon
- Cell volume and face fraction computation for immersed boundary
- `confinementInitialization_2d`: builds 11 CSR kpterms operators + profile_2d
- `csr_apply_variable_coeff` with cut-cell face fraction weighting
- `ZB8bandGeneralized`: sparse 8*(Nx*Ny) x 8*(Nx*Ny) Hamiltonian from COO insertion
- Band-offset profile on diagonal blocks

#### Phase 2: Eigensolver + Wire Band Structure -- DONE

- FEAST wrapper (`zfeast_hcsrev`) in `eigensolver.f90` with upper-triangle extraction
- Dense LAPACK fallback for small problems
- Energy window heuristic from profile band edges + kinetic cutoff
- `wire_coo_cache` for symbolic assembly reuse across k-points (O(NNZ) value update)
- kx sweep loop in `main.f90` (confinement=2 path)
- 2D eigenfunction output in `outputFunctions.f90`
- Wire regression test: GaAs 22nm x 22nm rectangle

#### Review fixes applied

All code review items resolved:
- C1: Flat index convention unified to `(iy-1)*nx + ix` (column-major)
- C2: Memory leak in ZB8bandGeneralized fixed (self-aliasing guard in csr_add)
- C3: COO overflow hard-stops replaced with warnings
- P1-A: O(N^2) insertion sort replaced with O(N log N) merge sort
- P1-B: Sparsity pattern reuse via `wire_coo_cache` + `csr_set_values_from_coo`
- P1-C: Self-aliasing guard in `csr_add` using `loc()` check
- P1-D: Non-square grid test (nx=4, ny=6) for ZB8bandGeneralized with Hermiticity verification
- P2 items: dead code removed, misleading names fixed, face fractions implemented, debug prints removed

#### Phase 3: Strain Solver -- DONE

- `strain_config` type in `defs.f90` (enabled, reference, solver, piezoelectric)
- 8 strain parameters added to `paramStruct`: C11, C12, C44, a0, ac, av, b_dp, d_dp
- All 33 materials in `parameters.f90` populated from Vurgaftman 2001 (Table XII/XIII) + Winkler 2003
- Alloys use Vegard interpolation between binary endpoints
- `strain_solver.f90` module (957 lines):
  - `compute_strain`: dispatcher for QW (biaxial algebra) and wire (plane-strain PDE via PARDISO)
  - `apply_pikus_bir`: diagonal band edge shifts from strain tensor (CB: ac*Tr(eps), HH/LH/SO via P_eps/Q_eps)
  - QW biaxial: eps_xx=eps_yy=eps_0, eps_zz=-2*C12/C11*eps_0
  - Wire plane-strain: 2*Ngrid DOF Navier-Cauchy PDE with cubic anisotropic elasticity, real symmetric indefinite PARDISO (mtype=-2)
  - Cut-cell face fraction weighting for stress-free BCs at wire boundary
  - Lattice-matched early exit (zero strain when all a0 identical)
- Input parsing: `strain`, `strain_ref`, `strain_solver`, `piezoelectric` fields in `input.cfg`
- Integration in `main.f90`: strain computed before eigenvalue solve for both wire and QW modes
- 6 unit tests: biaxial QW, hydrostatic, Pikus-Bir diagonal, disabled, lattice-matched, result_free
- av sign convention: positive in code (P_eps = -av*Tr(eps)), opposite to Vurgaftman Table XIII

#### Phase 3 review fixes

- Cross-derivative FD coefficients: corrected averaging factor (0.5 like Laplacian) and stencil denominator (4*dx*dy)
- PARDISO mtype changed from 2 (SPD) to -2 (symmetric indefinite) for positive semi-definite stiffness matrix with all-Neumann BCs
- PARDISO memory release added on error path (phase=-1)
- Removed O(N^2) duplicate search from COO entry (stencil produces no duplicates; merge sort handles merging)
- Fixed `return` -> pre-loop allocation guard in `apply_pikus_bir`

#### Current test coverage

| Test suite | Tests | Status |
|---|---|---|
| Unit: test_defs | 2 | Passing |
| Unit: test_finitedifferences | 7 | Passing |
| Unit: test_utils | 2 | Passing |
| Unit: test_parameters | 13 | Passing (original + 5 strain + 2 alloy/InSb verification) |
| Unit: test_hamiltonian | 1 | Passing |
| Unit: test_poisson | 2 | Passing |
| Unit: test_charge_density | 1 | Passing |
| Unit: test_sc_loop | 1 | Passing |
| Unit: test_geometry | ~12 | Passing |
| Unit: test_hamiltonian_2d | 11 | Passing (kron, Laplacian, cross-deriv, kpterms, Hermiticity, profile, CSR ops, full H) |
| Unit: test_eigensolver | ~8 | Passing (free electron, harmonic oscillator, CSR, FEAST) |
| Unit: test_strain_solver | 6 | Passing (biaxial QW, hydrostatic, Pikus-Bir, disabled, lattice-matched, free) |
| Regression: bulk_gaas | 1 | Passing |
| Regression: qw_alsbw_gasbw_inasw | 1 | Passing |
| Regression: gfactor_cb | 1 | Passing |
| Regression: bulk_inas | 1 | Passing |
| Regression: sc_gaas_alas_qw | 1 | Passing |
| Regression: qcse_gaas_algaas | 1 | Passing |
| Regression: qcse_gaas_algaas_ef | 1 | Passing |
| Regression: stark_shift | 1 | Passing |
| Regression: sc_mod_doped_gaas_algaas | 1 | Passing |
| Regression: wire_gaas_rectangle | 1 | Passing |

### Remaining Phases

| Phase | Description | Status |
|---|---|---|
| Phase 3 | Strain solver (plane-strain PDE + Pikus-Bir) | **DONE** |
| Phase 4 | Self-consistent 2D Schrodinger-Poisson | Not started |
| Phase 5 | G-factor + optical matrix elements | Not started |

## 11. Implementation Phasing

### Phase 0: Foundation Refactoring

Goal: Refactor core types and FD infrastructure to be dimension-agnostic. All existing tests pass.

- 0a. `spatial_grid` type in `defs.f90`
- 0b. `csr_matrix` type in new `sparse_matrices.f90`
- 0c. Kronecker product routines
- 0d. Refactor `input_parser.f90` to use `spatial_grid`
- 0e. Refactor `confinementInitialization` to accept `spatial_grid`
- 0f. Optional: sparse QW path via FEAST, benchmark vs. ZHEEVX

Tests: All 21 existing tests pass. New: test_sparse_matrices, test_spatial_grid.

### Phase 1: 2D FD + Sparse Assembly

Goal: Build 2D kpterms operators, assemble wire Hamiltonian. First eigenvalue from rectangular wire.

- 1a. 2D FD stencil construction (Kronecker products)
- 1b. Cut-cell geometry engine (`geometry.f90`)
- 1c. `apply_variable_coeff_2d` with cut-cell weighting
- 1d. 2D `confinementInitialization`
- 1e. 2D profile construction
- 1f. `ZB8bandGeneralized` assembly

Tests: test_fd_2d, test_cutcell, test_hamiltonian_2d.

### Phase 2: Eigensolver + Wire Band Structure

Goal: First wire band structure E_n(kx).

- 2a. FEAST wrapper (zfeast_hcsrev)
- 2b. ARPACK wrapper (mkl_sparse_z_ev)
- 2c. Energy window heuristic
- 2d. kx sweep loop in main.f90
- 2e. 2D eigenvector output + plotting
- 2f. Spurious solution filtering

Tests: test_eigensolver, regression: wire_gaas_rectangle, wire_inas_hexagon.

### Phase 3: Strain Solver -- DONE

Goal: Strained lattice-mismatched wires.

- 3a. Elastic constants + deformation potentials in `parameters.f90` -- DONE
- 3b. Plane-strain PDE solver via PARDISO -- DONE (mtype=-2, symmetric indefinite)
- 3c. Pikus-Bir module -- DONE
- 3d. Strain for QW (trivial biaxial) -- DONE
- 3e. Input fields for strain -- DONE (strain, strain_ref, strain_solver, piezoelectric)
- 3f. Strain visualization output -- DEFERRED (not needed for initial validation)

Tests: test_strain_solver (6 tests passing). Regression: wire_strained_inas_gaas DEFERRED to Phase 4.

### Phase 4: Self-Consistent 2D SP

Goal: Doped wires, gate-controlled carrier densities.

- 4a. 2D Poisson solver (PARDISO)
- 4b. 2D charge density
- 4c. 2D SC loop (generalize arrays)
- 4d. Per-region doping
- 4e. Optional: FFTW spectral Poisson

Tests: test_poisson_2d, test_sc_loop_2d, regression: test_sc_wire.

### Phase 5: G-Factor + Optics

Goal: Full optoelectronic characterization.

- 5a. Wire g-factor (Lowdin on sparse H)
- 5b. Inter-subband optical transition matrix elements
- 5c. Absorption/gain spectra
- 5d. `gfactorCalculation` wire path

Tests: regression: wire g-factor validation.

### Dependency graph

```
Phase 0 --- Phase 1 --- Phase 2 --- Phase 3 --- Phase 4 --- Phase 5
(types+FD)   (2D FD+     (eigenvalue  (strain)     (SC loop)    (g-factor,
             assembly)    solver)                                  optics)

                   QW backward compatibility maintained throughout
```

Each phase is a separate feature branch + PR. Tests from all previous phases must pass before merging.

## 12. FFTW Integration

Port patterns from the MagneticQD project (PostDoc codebase):

### 11.1 Spectral Poisson solver (Phase 4)

For uniform-eps wires on regular grids, solve nabla^2 phi = -rho/eps via 2D FFT:

1. FFT rho(y,z) -> rho_hat(ky, kz)
2. Divide by -(ky^2 + kz^2) in reciprocal space (avoid k=0 divergence)
3. Inverse FFT -> phi(y,z)

Cost: O(N log N). Uses same `fft_c2c_2d` pattern as MagneticQD `coulombFFTintegration.f90`.

Limitation: Only works for constant eps. For variable eps, fall back to PARDISO.

### 11.2 Coulomb integrals (future, if needed)

If many-body corrections or exciton calculations are added, the FFT-based Coulomb integral approach from MagneticQD can be ported directly.

## 13. Computational Limits

| Grid (Ny x Nz) | DOF (8*Ngrid) | NNZ | CSR Memory | Solve Time (est.) |
|---|---|---|---|---|
| 50 x 50 | 20,000 | ~1M | ~25 MB | ~5s |
| 100 x 100 | 80,000 | ~4M | ~100 MB | ~30s |
| 150 x 150 | 180,000 | ~9M | ~225 MB | ~2min |
| 200 x 200 | 320,000 | ~16M | ~400 MB | ~5min |

Single workstation (64 GB RAM): up to ~200x200 grid. Full k_x sweep with 100 points: ~1-10 hours.
