# Phase 6: Topological Suite Completion — Design Spec

Date: 2026-05-05
Branch: feature/bdg-topological-superconductivity
Approach: Extend existing modules (Approach A with selective splitting)

## Scope

All 6 BACKLOG items:
1. Fu-Kane Z2 invariant (full parity method on 8-band k.p QW)
2. QW BdG Hamiltonian (dense 16N x 16N variant)
3. Gap sweep / phase diagram (Z2 vs B, mu for wire and QW)
4. Berry curvature heatmap + BHZ phase transition figures
5. Unit tests (new test_z2_invariant.pf + extend 3 existing files)
6. Conductance (Kubo + Landauer-Büttiker) + spectral function A(k,E)

## Architecture Decision

Keep existing module boundaries. Fill stubs and add new functions to:
- `topological_analysis.f90` — Fu-Kane, Berry curvature, conductance
- `bdg_hamiltonian.f90` — QW BdG variant
- `green_functions.f90` — spectral function

Split `topological_analysis.f90` into sub-modules only if it exceeds ~1000 lines.

Figures via existing Python infrastructure (`generate_all_figures.py`).

---

## 1. Fu-Kane Z2 Invariant

### Algorithm

Compute Z2 from parity eigenvalues at 4 Time-Reversal Invariant Momenta (TRIM) in the 2D QW Brillouin zone:

1. Define TRIM: Gamma=(0,0), M1=(pi/a, 0), M2=(0, pi/a), M3=(pi/a, pi/a)
2. At each TRIM, diagonalize the 8N x 8N QW Hamiltonian (dense `zheevx`)
3. For occupied bands (below Fermi level), compute parity eigenvalue:
   - Band-space parity: `P_band = [+1,+1,+1,+1,+1,+1,-1,-1]` (inversion acts as +1 on valence, -1 on conduction)
   - Envelope parity: `int psi_n(z) * psi_n(-z) dz` (assumes symmetric QW — parity undefined for asymmetric structures)
   - Total: delta_n(k_i) = P_band(b) * envelope_parity
4. Product over occupied bands at each TRIM: delta_i = prod_n delta_n(k_i)
5. Z2 = prod_i (delta_i / delta_Gamma) — Z2=1 if negative (topological), Z2=0 if positive (trivial)

### Interface

```fortran
function compute_z2_fukane_qw(cfg, profile, kpterms, n_occ) result(z2)
  type(simulation_config), intent(in) :: cfg
  real(dp), contiguous, intent(in) :: profile(:,:)
  real(dp), contiguous, intent(in) :: kpterms(:,:,:)
  integer, intent(in) :: n_occ  ! number of occupied bands
  integer :: z2                 ! 0 = trivial, 1 = topological
end function
```

### Testing

- BHZ model with known Z2=0 (M > 0, trivial) and Z2=1 (M < 0, topological)
- Parity eigenvalue sign verification
- TRIM product consistency
- Gap closing edge case

---

## 2. QW BdG Hamiltonian

### Architecture

Dense 16N x 16N Nambu-space Hamiltonian for 1D quantum well at given k_parallel:

```
H_BdG = [[H(k_par) - mu*I,  Delta],
          [Delta^dagger,     -H*(k_par) + mu*I]]
```

- H(k_par) from existing `ZB8bandQW` (8N x 8N dense)
- Delta: antidiagonal pairing `(+1,-1,+1,-1,-1,+1,+1,-1)` replicated N times
- H*: complex conjugate (time-reversed partner)
- Zeeman only (no Peierls — QW is 1D in z)
- Dense storage (N=50-200 typical, 16N ~ 800-3200, well within LAPACK range)

### Interface

```fortran
subroutine build_bdg_hamiltonian_qw(H_bdg, cfg, profile, kpterms, k_par, &
                                     mu, delta_0, B_vec, g_factor)
  complex(dp), allocatable, intent(out) :: H_bdg(:,:)  ! 16N x 16N
  type(simulation_config), intent(in) :: cfg
  real(dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
  real(dp), intent(in) :: k_par, mu, delta_0
  real(dp), intent(in), optional :: B_vec(3), g_factor
end subroutine
```

### Testing

- Dimension check: size(H_bdg, 1) == 16 * N
- Hermiticity of H_BdG
- Particle-hole symmetry: eigenvalue pairs +/- E
- Recovery of wire BdG results for limiting case

---

## 3. Gap Sweep / Phase Diagram

### Architecture

Extend existing `compute_phase_diagram` with Z2-specific sweep:

```fortran
subroutine compute_z2_gap_sweep(cfg, B_min, B_max, nB, mu_min, mu_max, nMu, &
                                 gap_threshold, z2_map, gap_map, transitions)
  type(simulation_config), intent(in) :: cfg
  real(dp), intent(in) :: B_min, B_max, mu_min, mu_max, gap_threshold
  integer, intent(in) :: nB, nMu
  integer, allocatable, intent(out) :: z2_map(:,:)    ! nMu x nB
  real(dp), allocatable, intent(out) :: gap_map(:,:)  ! nMu x nB
  real(dp), allocatable, intent(out) :: transitions(:,:)  ! n_trans x 2 (B, mu)
end subroutine
```

- **Wire mode:** Z2 from M-sign (fast, analytical boundary at M=0)
- **QW mode:** Z2 from Fu-Kane parity at each grid point (expensive but correct)
- Phase boundary detection: track Z2 flips, linear interpolation for smooth boundary
- Output: `z2_phase_diagram.dat` (B, mu, z2, gap columns)

### Testing

- Known topological transition at M=0 for BHZ model
- Phase diagram array shape verification
- Gap closing at phase boundary

---

## 4. Berry Curvature + Figures

### Algorithm

Extract per-point Berry curvature from existing FHS (Fukui-Hatsugai-Suzuki) lattice gauge method:

```
Omega(kx, ky) = Im[ln(U_x * U_y(k+x) * U_x^-1(k+y) * U_y^-1)]
```

Where U_x = `<u_n(k)|u_n(k+x)> / |<u_n(k)|u_n(k+x)>|` is the normalized overlap.

Refactor `compute_chern_qwz` to also output the per-point curvature array (currently internal).

**Eigenvector source:** For QWZ, eigenvectors come from the QWZ model Hamiltonian on a (kx, ky) mesh. For QW, a k_parallel mesh sweep via `ZB8bandQW` provides eigenvectors at each mesh point. The calling code is responsible for providing the eigenvector array — the Berry curvature function is pure math on eigenvectors.

### Interface

```fortran
function compute_berry_curvature_lattice(evecs, kx_arr, ky_arr, n_occ) result(Omega)
  complex(dp), intent(in) :: evecs(:,:,:,:)  ! (basis, n_occ, nkx, nky)
  real(dp), intent(in) :: kx_arr(:), ky_arr(:)
  integer, intent(in) :: n_occ
  real(dp), allocatable :: Omega(:,:)  ! nkx x nky
end function
```

Works for both QWZ (2-band) and QW (8N-band).

### Figures (Python, via `generate_all_figures.py`)

| Figure | Content | Source |
|--------|---------|--------|
| `chern_berry_curvature_qwz.png` | Berry curvature heatmap for QWZ model at u=-0.8 (C=+1) | FHS curvature output |
| `bhz_z2_phase_transition.png` | Z2 vs (M, width) phase diagram | Gap sweep output |
| `bhz_edge_localization.png` | Edge state density vs position, trivial vs topological | Edge state extraction |

### Testing

- Berry curvature integrates to Chern number: `sum(Omega) * dk^2 / (2*pi) == C`
- Known curvature sign for QWZ model
- Curvature vanishes at TRIM points

---

## 5. Unit Tests

### New file

| File | Tests | Coverage |
|------|-------|----------|
| `test_z2_invariant.pf` | 4-5 | Fu-Kane on BHZ (Z2=0, Z2=1), parity eigenvalue sign, TRIM product, gap closing edge case |

### Extensions to existing files

| File | New tests | Coverage |
|------|-----------|----------|
| `test_edge_states.pf` (+3) | QW edge states via Fu-Kane, edge/bulk density ratio, decay length consistency |
| `test_green_functions.pf` (+3) | Spectral function normalization, k-resolved A(k,E) sum rule, LDOS-spectral consistency |
| `test_magnetic_field.pf` (+3) | Zeeman splitting magnitude, Peierls gauge invariance, Landau level count |

Total: ~12-14 new tests.

---

## 6. Conductance

### 6a. Kubo Linear Response

```fortran
function compute_conductance_kubo(berry_curvature, kx_arr, ky_arr, n_occ) result(sigma_xy)
  real(dp), contiguous, intent(in) :: berry_curvature(:,:)
  real(dp), contiguous, intent(in) :: kx_arr(:), ky_arr(:)
  integer, intent(in) :: n_occ
  real(dp) :: sigma_xy  ! in units of e^2/h
end function
```

- Hall conductance: `sigma_xy = (1/(2*pi)) * sum_k Omega(k) * dk^2` — this is the Chern number times e^2/h
- Longitudinal conductance from velocity-velocity correlation (requires velocity matrices from existing `build_velocity_matrices`)

### 6b. Landauer-Büttiker

```fortran
function compute_conductance_landauer(H_csr, grid, mu, eta) result(G)
  type(csr_matrix), intent(in) :: H_csr
  type(spatial_grid), intent(in) :: grid
  real(dp), intent(in) :: mu, eta
  real(dp) :: G  ! in units of 2e^2/h
end function
```

- Transmission: `T = Tr[Gamma_L * G^R * Gamma_R * G^A]`
- Lead self-energies from wide-band approximation
- Reuses PARDISO-based Green function infrastructure

**Priority:** Kubo first (reuses Berry curvature), Landauer second.

### Testing

- Kubo: quantized sigma_xy = C * e^2/h for known Chern numbers
- Kubo: sigma_xy = 0 for trivial insulator
- Landauer: G = 2e^2/h for single conducting channel
- Landauer: G = 0 for insulator (gap > E_F)

---

## 7. Spectral Function A(k,E)

### Algorithm

k-resolved spectral function — momentum-resolved LDOS:

**QW path (dense):**
1. At each k_parallel, diagonalize H(k_par) once (`zheevx`)
2. `A(k, E) = sum_n delta_eta(E - E_n(k))` where delta_eta is Lorentzian broadening
3. Fast: no matrix inversion needed, just eigenvalue evaluation

**Wire path (CSR):**
1. For each kz, call the existing wire Hamiltonian builder with that kz value (kz is a parameter, not a Bloch phase on the CSR)
2. For each E, solve `(E + i*eta - H(kz)) * x = e_j` via PARDISO
3. `A(kz, E) = -(1/pi) * Im[sum_j G_jj(kz, E)]`

### Interface

```fortran
subroutine compute_spectral_function(H_csr, grid, k_arr, E_arr, eta, A_kE)
  type(csr_matrix), intent(in) :: H_csr
  type(spatial_grid), intent(in) :: grid
  real(dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
  real(dp), intent(in) :: eta
  real(dp), allocatable, intent(out) :: A_kE(:,:)  ! size(nk, nE)
end subroutine

! Dense QW overload
subroutine compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, eta, A_kE)
  type(simulation_config), intent(in) :: cfg
  real(dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
  real(dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
  real(dp), intent(in) :: eta
  real(dp), allocatable, intent(out) :: A_kE(:,:)
end subroutine
```

Generic interface `compute_spectral_function` dispatches to CSR or dense path.

### Output

`spectral_function.dat` (k, E, A columns) for Python heatmap generation.

### Figure

`spectral_function_wire.png` — A(k,E) heatmap showing band dispersion and gap structure.

### Testing

- Sum rule: `int A(k, E) dE = 1` per k-point (normalized states)
- LDOS consistency: `int A(k, E) dk = LDOS(E)`
- Peak positions match eigenvalue dispersion
- Gap region has A(k,E) ~ 0

---

## Module Impact Summary

| File | Changes |
|------|---------|
| `topological_analysis.f90` | Fu-Kane implementation, Berry curvature extraction, Kubo conductance, gap sweep extension |
| `bdg_hamiltonian.f90` | New `build_bdg_hamiltonian_qw` dense variant |
| `green_functions.f90` | Spectral function computation (both dense and CSR) |
| `main_topology.f90` | New dispatch branches for Fu-Kane QW, QW BdG, gap sweep, conductance, spectral function |
| `input_parser.f90` | New config fields for gap sweep, conductance, spectral function parameters |
| `defs.f90` | New fields in `topology_config` and `topological_result` types |
| `generate_all_figures.py` | 4 new figure functions (Berry curvature heatmap, phase diagram, edge localization, spectral function) |
| `tests/unit/test_z2_invariant.pf` | New file, 4-5 tests |
| `tests/unit/test_edge_states.pf` | +3 tests |
| `tests/unit/test_green_functions.pf` | +3 tests |
| `tests/unit/test_magnetic_field.pf` | +3 tests |
