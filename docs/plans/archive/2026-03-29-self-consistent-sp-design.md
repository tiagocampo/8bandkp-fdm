# Self-Consistent Schrödinger-Poisson Design

**Date:** 2026-03-29
**Branch:** `feature/self-consistent-sp`
**Status:** Design approved, pending implementation

## 1. Overview

Extend the 8-band k.p FDM code with self-consistent Schrödinger-Poisson (SP) capability for both bulk and quantum well systems, with and without external electric fields. The implementation adds a Poisson solver, charge density calculation from 8-component k.p eigenstates, and an iterative self-consistency driver with DIIS acceleration.

## 2. Scope

| Mode | Self-consistency behavior |
|---|---|
| Bulk, no EF | Find Fermi level for charge neutrality with given doping |
| Bulk + EF | Apply uniform shift + find self-consistent Fermi level |
| QW, no EF | Full SP loop: solve k.p → charge density → Poisson → iterate |
| QW + EF | EF sets initial potential, SP loop screens it self-consistently |

## 3. References

### Methodology
- Tan, Snider, Chang, Hu, *J. Appl. Phys.* **68**, 4071 (1990) — canonical SP benchmark with Newton-Raphson Poisson
- Birner et al., *Acta Phys. Pol. A* **110**, 111 (2006) — nextnano 8-band k.p + Poisson implementation
- Trellakis, Galick, Pacelli, Ravaioli, *J. Appl. Phys.* **81**, 7880 (1997) — predictor-corrector iteration scheme
- Stern, *J. Comput. Phys.* **6**, 56 (1970) — iteration methods for self-consistent fields

### Mixing/convergence
- Pulay, *Chem. Phys. Lett.* **73**, 393 (1980) — original DIIS paper
- Pulay, *J. Comput. Chem.* **3**, 556 (1982) — improved SCF convergence acceleration
- Marks, *J. Chem. Theory Comput.* **17**, 5715 (2021), arXiv:2104.04384 — Anderson/Pulay/DIIS equivalence analysis

### k.p theory
- Winkler, *Spin-Orbit Coupling Effects in 2D Electron and Hole Systems*, Springer (2003)
- Chuang, Chang, *Semicond. Sci. Technol.* **12**, 252 (1997)
- Foreman, *Phys. Rev. B* **56**, R12748 (1997)

### Material parameters
- Vurgaftman, Meyer, Ram-Mohan, *J. Appl. Phys.* **89**, 5815 (2001) — III-V parameters including dielectric constants
- Adachi, *Properties of Group-IV, III-V and II-VI Semiconductors*, Wiley (2005)

### Benchmark validation
- Snider's 1D Poisson code (nd.edu/~gsnider) — reference GaAs/AlGaAs results
- Erenler, *Comput. Phys. Commun.* **200**, 247 (2016) — Aestimo 1D SP solver
- Pfeffer, *Phys. Rev. B* **59**, 3680 (1999) — InAs/AlSb band structure
- Sze, *Physics of Semiconductor Devices*, Wiley, 3rd ed. (2007) — bulk carrier statistics

## 4. Architecture

### New modules

```
src/physics/
  poisson.f90              -- Poisson solver (box-integration FD, Thomas algorithm)
  charge_density.f90       -- Compute n(z), p(z) from k.p eigenstates + k_∥ sampling
  sc_loop.f90              -- Self-consistent iteration driver (linear + DIIS mixing)
```

### Extended modules

```
src/core/defs.f90          -- New types: doping_spec, sc_config
src/core/parameters.f90    -- Add ε(0) dielectric constants to material database
src/io/input_parser.f90    -- Parse SC/doping input fields
src/apps/main.f90          -- Integrate SC loop for bandStructure
src/apps/main_gfactor.f90  -- Integrate SC loop for gfactorCalculation
```

### Data flow

```
Initial potential Φ₀(z)
    │
    ▼
┌─► Apply Φ(z) to profile → build H(k) → solve eigenproblem at each k_∥
│       │
│       ▼
│   Compute charge density: n(z) = Σ_s N_s(k_∥) |Ψ_s(z)|²   [charge_density.f90]
│       │
│       ▼
│   Solve Poisson: d/dz[ε(z) dΦ/dz] = -ρ(z)/ε₀               [poisson.f90]
│       │
│       ▼
│   Mix: Φ_new = (1-α)Φ_old + α·DIIS(Φ_Poisson)               [sc_loop.f90]
│       │
│       ▼
└─── Converged? |Φ_new - Φ_old|_∞ < tol ?
        │ Yes
        ▼
    Output final band structure / g-factor with self-consistent potential
```

**Key principle:** The existing `confinementInitialization` → `ZB8bandQW` pipeline remains unchanged. The SC loop modifies only the `profile` array before it enters the Hamiltonian, exactly like the current `externalFieldSetup_electricField` does.

## 5. New Derived Types

### `doping_spec` (in `defs.f90`)

```fortran
type doping_spec
    real(kind=dp) :: ND = 0.0_dp    ! donor concentration (cm^-3)
    real(kind=dp) :: NA = 0.0_dp    ! acceptor concentration (cm^-3)
end type doping_spec
```

### `sc_config` (in `defs.f90`)

```fortran
type sc_config
    integer       :: enabled = 0              ! 0=off, 1=on
    integer       :: max_iterations = 100     ! iteration cap
    real(kind=dp) :: tolerance = 1.0e-6_dp    ! convergence |ΔΦ|_∞ in eV
    real(kind=dp) :: mixing_alpha = 0.3_dp    ! linear mixing parameter
    integer       :: diis_history = 7         ! DIIS history length
    real(kind=dp) :: temperature = 300.0_dp   ! Kelvin
    integer       :: fermi_mode = 0           ! 0=charge_neutrality, 1=fixed
    real(kind=dp) :: fermi_level = 0.0_dp     ! fixed μ (eV), used if fermi_mode=1
    integer       :: num_kpar = 200           ! k_∥ sampling points
    real(kind=dp) :: kpar_max = 0.0_dp        ! max k_∥ (1/nm), 0=auto from Fermi
    character(len=3) :: bc_type = 'DD '       ! 'DD'=Dirichlet-Dirichlet, 'DN'=Dirichlet-Neumann
    real(kind=dp) :: bc_left = 0.0_dp         ! left BC potential (eV)
    real(kind=dp) :: bc_right = 0.0_dp        ! right BC potential (eV)
end type sc_config
```

### Extension to `simulation_config`

```fortran
type(doping_spec), allocatable :: doping(:)   ! per-layer doping, size numLayers
type(sc_config)                :: sc          ! self-consistency parameters
```

## 6. Input File Extensions

Backward-compatible. New fields appended after `EFParams`:

```
! Self-consistency
SC:            1                ! enable flag
SC_max_iter:   100
SC_tolerance:  1.0e-6
SC_mixing:     0.3
SC_diis:       7
SC_temperature: 300.0
SC_fermi_mode: 0               ! 0=charge_neutrality, 1=fixed
SC_fermi_level: 0.0            ! only used if fermi_mode=1
SC_num_kpar:   200
SC_kpar_max:   0.0             ! 0=auto
SC_bc:         DD              ! DD or DN
SC_bc_left:    0.0
SC_bc_right:   0.0
! Per-layer doping (repeated for each layer, after material lines)
doping:        0.0  0.0        ! ND, NA for layer 1
doping:        1e18  0.0       ! ND, NA for layer 2
```

Default `SC: 0` means no self-consistency — identical to current behavior.

## 7. Poisson Solver

**File:** `src/physics/poisson.f90`

### Box-integration discretization

Following Birner et al. (2006), integrate Poisson over each cell centered on grid point z_i:

$$\frac{\varepsilon_{i+1/2}(\Phi_{i+1}-\Phi_i) - \varepsilon_{i-1/2}(\Phi_i-\Phi_{i-1})}{h^2} = -\frac{\rho_i}{\varepsilon_0}$$

where $\varepsilon_{i\pm1/2} = (\varepsilon_i + \varepsilon_{i\pm1})/2$.

### Tridiagonal system

Resulting system AΦ = b solved by Thomas algorithm (O(N)):

- $b_i = \varepsilon_{i+1/2}/h^2$ (upper diagonal)
- $c_i = \varepsilon_{i-1/2}/h^2$ (lower diagonal)
- $a_i = -(b_i + c_i)$ (main diagonal)
- RHS: $-\rho_i/\varepsilon_0$

### Boundary conditions

- **Dirichlet:** replace row with identity (Φ₁ = V_left)
- **Neumann:** replace row with zero-flux (Φ₁ - Φ₂ = 0)

### Dielectric constants

Added to `parameters.f90` from Vurgaftman (2001) Table III:

| Material | ε_r |
|---|---|
| GaAs | 12.90 |
| AlAs | 10.06 |
| InAs | 15.15 |
| GaSb | 15.70 |
| InSb | 16.80 |
| AlSb | 12.04 |
| InP | 12.61 |
| GaP | 11.11 |
| AlP | 9.02 |

Alloy interpolation: ε(x) = (1-x)ε_A + xε_B.

### Main subroutine

```fortran
subroutine solve_poisson(phi, rho, epsilon, dz, N, bc_left, bc_right, bc_type)
```

## 8. Charge Density Calculation

**File:** `src/physics/charge_density.f90`

### Core formula (QW mode)

$$n(z) = \sum_s \int_0^{k_{max}} |\Psi_s(z, k_\parallel)|^2 \cdot f(E_s(k_\parallel) - \mu, T) \cdot \frac{k_\parallel}{2\pi} \, dk_\parallel$$

where $|\Psi_s(z)|^2 = \sum_{\nu=1}^{8} |\Psi_s^\nu(z)|^2$ (all 8 envelope components, per Winkler 2003).

### k_parallel sampling

- Solve k.p eigenproblem at k_∥ = 0, Δk, 2Δk, ..., k_max (isotropic in-plane)
- `kpar_max = 0` (auto): increase until f(E_s(k_max)) < 10⁻⁶
- `num_kpar` points (default 200), uniform spacing
- Integration via Simpson's rule (reusing `simpson` from `utils.f90`)

### Total charge

$$\rho(z) = q \left[N_D^+(z) - N_A^-(z) + p(z) - n(z)\right]$$

- n(z) = Σ CB subband contributions (electrons)
- p(z) = Σ VB subband contributions (holes = empty VB states)
- N_D⁺, N_A⁻: ionized dopant concentrations (full ionization initially)

### Spin degeneracy

Factor of 2 included in the occupation integral for each subband (both spin channels occupied identically in the 8-band basis).

### Bulk mode

No z-dependence. Charge density computed from 3D DOS via explicit k-sampling of the bulk E(k) dispersion. The SC loop simplifies to a Fermi level finder — no spatial Poisson solve needed.

### Main subroutines

```fortran
subroutine compute_charge_density_qw(n_electron, n_hole, eigenvectors, &
    & eigenvalues, kpar_grid, fermi_level, temperature, N, num_subbands, &
    & num_kpar, dz)

subroutine compute_charge_density_bulk(n_electron, n_hole, eigenvalues_k, &
    & k_grid, fermi_level, temperature, num_k, volume)
```

## 9. Self-Consistent Loop

**File:** `src/physics/sc_loop.f90`

### Algorithm: Linear + DIIS

Following Pulay (1980) and Marks (2021):

**Phase 1 — Linear mixing (warm-up, first m iterations):**
  Φ_{n+1} = (1-α)Φ_n + α·Φ_Poisson

**Phase 2 — DIIS (after m iterations):**
  1. Store m pairs {(Φ_i, r_i)} where r_i = Φ_Poisson - Φ_input
  2. Solve (m+1)×(m+1) least-squares system via LAPACK dgesv:
     min ||Σ c_j r_j||² subject to Σ c_j = 1
  3. Extrapolate: Φ_extrapolated = Σ c_j Φ_j
  4. Blend: Φ_{n+1} = (1-α)Φ_n + α·Φ_extrapolated

### Convergence criteria

| Criterion | Default threshold |
|---|---|
| |Φ_new - Φ_old|_∞ | 10⁻⁶ eV |
| |ρ_new - ρ_old|_∞ | 10⁻⁶ e·cm⁻³ |
| Iteration count | 100 (safety cap) |

### Fermi level iteration (charge neutrality mode)

Inner bisection loop within each SP iteration to find μ such that:

$$\int [n(z) - p(z)] dz = \int [N_D^+(z) - N_A^-(z)] dz$$

Bounds: μ ∈ [E_{CB,min} - 10k_BT, E_{VB,max} + 10k_BT], tolerance 10⁻⁸ eV. Converges in ~30 bisection steps.

### Integration with existing code

```fortran
! In main.f90 / main_gfactor.f90, after confinementInitialization:
if (cfg%sc%enabled == 1) then
    call self_consistent_loop(profile, cfg, kpterms, ...)
else
    ! existing single-shot solve
end if
```

### Logging

Each iteration prints: iteration number, |ΔΦ|_∞, |Δρ|_∞, Fermi level, total charge. DIIS coefficient norms logged for debugging.

## 10. Electric Field Corrections

### Current implementation audit

The existing `externalFieldSetup_electricField` applies a linear potential tilt:
```
profile(i,:) -= (Evalue * totalSize) * (z(i) + z(1)) / (2*z(1))
```

This is correct for a uniform field E in a domain [z(1), -z(1)] — it creates a linear ramp from -E·L at one boundary to 0 at the other. The formula will be verified and kept as-is.

### Bulk mode extension

Currently EF has no effect in bulk. To enable: add a uniform diagonal shift to the 8×8 bulk Hamiltonian proportional to field strength. For true self-consistency with spatial variation, the bulk case reduces to the Fermi level finder.

### EF + SC interaction

When both are enabled, EF sets the initial potential condition. The SC loop then screens this field self-consistently, modifying the profile at each iteration.

## 11. Bulk Self-Consistency

For bulk (confinement=0), the SP loop simplifies:

1. No z-grid, no spatial profile
2. Compute E(k) from 8×8 bulk Hamiltonian
3. Compute carrier density n, p from 3D DOS via k-sampling
4. Adjust Fermi level μ until charge neutrality: q(ND - NA + p - n) = 0
5. For EF mode: apply uniform shift to Hamiltonian, then find μ

No Poisson solve needed — the "potential" is just the Fermi level position.

## 12. Testing Strategy

### Tier 1: Unit tests (pFUnit)

**test_poisson.pf:**
- Uniform ε, uniform ρ → analytical Φ(z) = -ρz²/(2εε₀) + C₁z + C₂
- Two-layer ε → verify D = εE continuous across interface
- DD vs DN boundary conditions
- Thomas algorithm with known system

**test_charge_density.pf:**
- Gaussian wavefunction + known k-points → analytical integral comparison
- T=0 Fermi step function → verify cutoff
- T=300K Fermi tail → compare against explicit evaluation
- 8-component wavefunction → verify all components summed

**test_sc_loop.pf:**
- Linear mixing → verify geometric convergence rate (1-α)^n
- DIIS on 1D model → verify accelerated convergence
- Max-iterations safety cap
- Charge neutrality bisection → verify μ to tolerance

### Tier 2: Integration tests (shell scripts)

- `test_sc_qw_gaas_agaas.sh` — QW with SC, verify convergence log
- `test_sc_bulk_gaas.sh` — Bulk with doping, verify Fermi level output
- `test_sc_qw_with_ef.sh` — QW with EF + SC, verify screened potential

### Tier 3: Regression tests

**Benchmark 1: GaAs/AlGaAs modulation-doped QW**
- Reference: Tan, Snider, Chang, Hu, *J. Appl. Phys.* **68**, 4071 (1990)
- Cross-validated: nextnano (Birner et al. 2006), Aestimo (Erenler 2016), Snider's code
- Target: E₁ ≈ -3 meV, sheet density ≈ 6.6×10¹¹ cm⁻²
- Config: `tests/regression/configs/sc_qw_gaas_agaas_snider.cfg`

**Benchmark 2: GaAs/AlGaAs undoped QW with external field**
- Apply E = 50 kV/cm, run SC
- Cross-validate against Aestimo 1D (Erenler 2016)
- Config: `tests/regression/configs/sc_qw_ef_gaas_agaas.cfg`

**Benchmark 3: InAs/AlSb QW**
- Reference: Pfeffer, *Phys. Rev. B* **59**, 3680 (1999)
- Tests nonparabolic charge density (narrow-gap material)
- Config: `tests/regression/configs/sc_qw_inas_alsb.cfg`

**Benchmark 4: Bulk doped GaAs**
- ND = 10¹⁷ cm⁻³, T = 300K
- Reference: Sze, *Physics of Semiconductor Devices* (2007), Chapter 1
- Target: n ≈ ND, μ ≈ E_C - k_BT·ln(N_C/N_D)
- Config: `tests/regression/configs/sc_bulk_gaas_doped.cfg`

### Tier 4: Citation documentation

Every test file header cites validation references with expected tolerances. Example:
```fortran
! Validation: Tan et al., J. Appl. Phys. 68, 4071 (1990)
! Tolerance: eigenenergy ±0.5 meV, sheet density ±5%
```

## 13. Implementation Phases

### Phase 1: Material Parameters Extension (~50 lines)
- Add `eps0` to material parameter type
- Populate from Vurgaftman Table III
- Alloy interpolation
- Unit test for ε values

### Phase 2: Poisson Solver (~150 lines)
- Box-integration discretization
- Thomas algorithm
- DD and DN BCs
- Unit tests (4 cases)

### Phase 3: Charge Density Module (~200 lines)
- QW: 8-component density, Simpson integration, k_∥ sampling
- Bulk: 3D DOS from k-sampling
- Fermi-Dirac occupation
- Unit tests (4 cases)

### Phase 4: SC Loop Driver (~250 lines)
- Linear + DIIS mixing
- Fermi level bisection
- Convergence monitoring/logging
- Integration into main.f90 and main_gfactor.f90
- Extended defs.f90 types and input_parser.f90
- Unit tests (4 cases)

### Phase 5: Electric Field Fix + Bulk Support (~100 lines)
- Verify/fix externalFieldSetup_electricField
- Enable EF in bulk mode
- EF + SC interaction
- Integration tests

### Phase 6: Regression Tests and Validation (~300 lines)
- 4 benchmark configs + reference data
- ctest integration with `sc` label
- Citation documentation in test files
