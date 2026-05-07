# Optical Properties Redesign: Commutator-Based Velocity Operator

Date: 2026-04-25
Branch: feature/optical-properties-commutator

## Problem

The current optical properties implementation in `optical_spectra.f90` has several shortcomings:

1. **Incorrect matrix elements for confined systems**: `pMatrixEleCalc` uses the diagonal Kane P parameter (`kpterms(:,:,4)`) as a scalar approximation for `dH/dk_alpha`. This only captures the zone-center interband coupling and misses contributions from gradient operators in the QW and wire geometries — the same issue that breaks wire g-factors.

2. **Limited scope**: Only interband absorption (TE/TM) and ISBT for QWs. No spontaneous emission, no wire optics, no spin resolution.

3. **Embedded in bandStructure executable**: The optics accumulation is called from within the band-structure k-sweep loop, mixing two fundamentally different workflows.

4. **No bulk optics**: Bulk absorption/dielectric function not computed.

## Solution: Commutator-Based Velocity + Standalone Program

### Core principle

Compute the velocity operator from the Heisenberg equation of motion:

```
v_alpha = [r_alpha, H] / (i hbar)
```

Discretized for CSR matrices:

```
vel_alpha(i,j) = cmplx(0, -1) * (r_alpha(i) - r_alpha(j)) * H(i,j)
```

This element-wise scaling of the full Hamiltonian automatically captures all k.p coupling terms (PP, PM, Q, T, R, RC, A, S, SC) without separate perturbation modes. Same approach as the wire g-factor fix.

### Architecture

A new standalone `opticalProperties` executable handles all optical calculations. It builds velocity matrices via the commutator (shared with g-factor code) and accumulates spectra across a k_par sweep.

## Design Details

### 1. Shared Velocity Matrix Infrastructure

**File:** `src/physics/hamiltonianConstructor.f90`

**New subroutine:** `build_velocity_matrices(H_csr, grid, vel)`

- Input: full Hamiltonian CSR matrix, spatial grid with coordinates
- Output: array of 3 CSR matrices `vel(3)` (same sparsity pattern as H)
- For each nonzero entry `H(row, col)`:
  - Map row/col to spatial grid indices
  - `vel_x%val(k) = cmplx(0,-1,dp) * (x(grid_row) - x(grid_col)) * H%val(k)`
  - `vel_y%val(k) = cmplx(0,-1,dp) * (y(grid_row) - y(grid_col)) * H%val(k)`
  - `vel_z%val(k) = cmplx(0,-1,dp) * (z(grid_row) - z(grid_col)) * H%val(k)`

This subroutine is shared between g-factor and optics. The g-factor fix (branch `fix/wire-gfactor-commutator`) introduces it; optics consumes the same function.

**Geometry-specific grid coordinates:**

| Geometry | x-coords | y-coords | z-coords |
|---|---|---|---|
| Bulk | (no grid) | (no grid) | (no grid) |
| QW | 0 (all same) | 0 (all same) | z_grid(1:N) |
| Wire | x_grid(1:Ngrid) | y_grid(1:Ngrid) | scalar kz |

**Bulk special case:** The bulk Hamiltonian is 8x8 with no spatial discretization, so `[r,H]=0`. Use the existing analytical `ZB8bandBulk(g='g')` for bulk velocity. The commutator approach is only used for QW (confinement=1) and wire (confinement=2).

### 2. New Executable — `opticalProperties`

**Entry point:** `src/apps/main_optics.f90` (`program opticalProperties`)

**CMakeLists.txt:** Add new target alongside `bandStructure` and `gfactorCalculation`.

**Workflow:**

```
1. Parse input.cfg → simulation_config (with optics: block)
2. confinementInitialization → profile, kpterms
3. Build velocity matrices:
   bulk:  vel(d) = ZB8bandBulk(g='g', dir=d)           [analytical]
   QW:    build_velocity_matrices(H_csr, z_grid)        [commutator]
   wire:  build_velocity_matrices(H_csr, grid_2d)       [commutator]
4. k_par sweep loop:
   for each k_par:
     a. Construct H(k_par) — reuse existing Hamiltonian construction
     b. Diagonalize → eigvals, eigvecs
     c. optics_accumulate(eigvals, eigvecs, vel, k_weight, optcfg, ...)
5. optics_finalize() — apply prefactors, write output files
```

### 3. Optics Input Configuration

**Extended `optics_config` in `defs.f90`:**

```fortran
type :: optics_config
  real(dp) :: E_min, E_max
  integer  :: num_energy_points
  real(dp) :: linewidth_lorentzian, linewidth_gaussian
  real(dp) :: refractive_index
  real(dp) :: temperature

  ! Feature flags
  logical  :: absorption_enabled
  logical  :: spontaneous_enabled
  logical  :: gain_enabled
  logical  :: isbt_enabled
  logical  :: spin_resolved

  ! Gain parameters
  real(dp) :: gain_carrier_density  ! cm^-2

  ! Spin-resolved parameters
  real(dp) :: spin_polarization_c, spin_polarization_v
  real(dp) :: carrier_density_n, carrier_density_p  ! cm^-2 or cm^-3
end type
```

**Input file format** (`input.cfg` optics block):

```
optics:
  E_min E_max num_energy_points
  linewidth_lorentzian linewidth_gaussian
  refractive_index
  temperature
  absorption_enabled spontaneous_enabled gain_enabled isbt_enabled spin_resolved
  gain_carrier_density
```

### 4. Unified Accumulation Framework

All optical quantities share the same core structure:

```
S_alpha(E) = sum_{k_par} sum_{i,j} |<i|v_alpha|j>|^2 * Phi(E - E_ij) * W(k_par) * occ(i,j)
```

Where:
- `Phi(E - E_ij)` = Voigt lineshape (Lorentzian + Gaussian)
- `W(k_par)` = Simpson integration weight
- `occ(i,j)` = geometry-dependent occupation factor

**Occupation factors by quantity:**

| Quantity | States | occ(i,j) |
|---|---|---|
| Absorption | VB(i) → CB(j) | `f_v(E_i) - f_c(E_j)` |
| Spontaneous emission | CB(j) → VB(i) | `f_c(E_j) * (1 - f_v(E_i))` |
| Gain | VB(i) ↔ CB(j) | `f_v(E_i) - f_c(E_j)` (sign preserved) |
| ISBT | CB(i) → CB(j) | `f_c(E_i) - f_c(E_j)` |

For gain, the quasi-Fermi levels `mu_e` and `mu_h` are computed from the carrier density via bisection (reusing `charge_density` module's Fermi level finder).

**Velocity matrix element computation:**

```fortran
do dir = 1, 3
  call csr_spmv(vel(dir), state_j, Y)
  Pele = zdotc(dim, state_i, Y)
  p_abs2(dir) = real(Pele * conjg(Pele), dp)
end do
```

Polarization decomposition:
- TE: `p_abs2(1) + p_abs2(2)` (px + py)
- TM: `p_abs2(3)` (pz)

### 5. Spin-Resolved Optical Properties

**New file:** `src/physics/spin_projection.f90`

Implements the Clebsch-Gordan transformation from the 8-band k.p basis to the explicit spin-orbital basis (X_up, Y_up, Z_up, S_up, X_dw, Y_dw, Z_dw, S_dw).

The transformation is the same as the old code's basis manipulation in `sigmaElem`:

```
|X_up> = |HH_UP>/sqrt(2) + |LH_DW>/sqrt(6) - i|SO_DW>/sqrt(3)
|X_dw> = i|HH_DW>/sqrt(2) + i|LH_UP>/sqrt(6) + |SO_UP>/sqrt(3)
|Y_up> = i|HH_UP>/sqrt(2) - i|LH_DW>/sqrt(6) - |SO_DW>/sqrt(3)
|Y_dw> = |HH_DW>/sqrt(2) - |LH_UP>/sqrt(6) + i|SO_UP>/sqrt(3)
|Z_up> = -i*sqrt(2)*|LH_UP>/sqrt(3) + |SO_UP>/sqrt(3)
|Z_dw> = sqrt(2)*|LH_DW>/sqrt(3) + i|SO_DW>/sqrt(3)
|S_up> = |EL_UP>
|S_dw> = |EL_DW>
```

For each eigenstate, compute spin-up/down weights:

```
w_up = sum_{orbital} |<orbital_up|psi>|^2
w_dw = 1 - w_up
```

Spin-resolved spectra decompose `|<i|v_alpha|j>|^2` into spin channels. For interband transitions, cross-spin terms are small, so:

```
alpha_alpha_up ~ sum_{i,j} w_up(i) * w_up(j) * |<i|v_alpha|j>|^2 * occ
alpha_alpha_dw ~ sum_{i,j} w_dw(i) * w_dw(j) * |<i|v_alpha|j>|^2 * occ
```

Spin-resolved quasi-Fermi levels use separate spin-projected DOS, following the old code's approach.

### 6. ISBT for All Geometries

- **QW ISBT:** z-dipole `<CB_i|v_z|CB_j>` from commutator `[z, H]`. The FD stencil entries coupling neighboring z-points produce the velocity via `(z_i - z_j) * H(i,j)`. This correctly captures the envelope function transition dipole.

- **Wire ISBT:** Transverse dipoles from `[x, H]` and `[y, H]`. Captures confinement-induced transitions and material-inhomogeneity effects at interfaces.

- **Bulk ISBT:** No spatial variation, no subbands, no ISBT (correct).

### 7. Physical Prefactors (Finalization)

The accumulated raw spectra have units of `(eV * Angstrom)^2 * (1/Angstrom^dim)` where dim is the k-space dimensionality.

**Interband absorption coefficient:**

```
alpha(E) = C * (2/S) * sum |p_CV|^2 * occ * Phi(E)
C = 2*pi*e^2 / (n_r * c * eps0 * hbar^2 * E)
```

Where `|p_CV|^2 = (m0/hbar)^2 * |<CB|dH/dk|VB>|^2` and the final units are `cm^-1`.

**Spontaneous emission rate:**

```
R_sp(E) = C' * E * sum |p_CV|^2 * f_c * (1 - f_v) * Phi(E)
```

**Gain:** Same formula as absorption but with separate quasi-Fermi levels for electrons and holes. Negative alpha = gain.

### 8. Output Files

All files in `output/` directory:

```
output/absorption_TE.dat           ! E(eV)  alpha(cm^-1)
output/absorption_TM.dat           ! E(eV)  alpha(cm^-1)
output/absorption_ISBT.dat         ! E(eV)  alpha(cm^-1)
output/spontaneous_TE.dat          ! E(eV)  rate(arb. units)
output/spontaneous_TM.dat          ! E(eV)  rate(arb. units)
output/gain_TE.dat                 ! E(eV)  gain(cm^-1)
output/gain_TM.dat                 ! E(eV)  gain(cm^-1)
output/isbt_transitions.dat        ! i  j  E_ij(eV)  |z_ij|^2  f_ij
```

With spin resolution enabled, each gets `_up`/`_dw` variants:
```
output/absorption_TE_up.dat
output/absorption_TE_dw.dat
... etc.
```

## Files to Create/Modify

### New files
| File | Purpose |
|---|---|
| `src/apps/main_optics.f90` | Standalone program |
| `src/physics/spin_projection.f90` | Clebsch-Gordan transformation, spin weights |

### Modified files
| File | Changes |
|---|---|
| `src/physics/hamiltonianConstructor.f90` | Add `build_velocity_matrices()` (shared with g-factor fix) |
| `src/physics/optical_spectra.f90` | Rewrite: commutator-based, all geometries, all quantities |
| `src/core/defs.f90` | Extend `optics_config` type |
| `src/io/input_parser.f90` | Parse `optics:` block |
| `CMakeLists.txt` + `src/CMakeLists.txt` | Add `opticalProperties` target |
| `docs/lecture/06-optical-properties.md` | Update to reflect commutator approach |

### Files unchanged
| File | Reason |
|---|---|
| `src/physics/gfactor_functions.f90` | Continues to use its own velocity path; may share `build_velocity_matrices` later |
| `src/apps/main.f90` | Band structure program unchanged |
| `src/apps/main_gfactor.f90` | g-factor program unchanged |

## Implementation Phases

### Phase 1: Infrastructure (depends on g-factor commutator fix)
- `build_velocity_matrices()` in hamiltonianConstructor.f90
- Extend `optics_config` in defs.f90
- Parse optics: block in input_parser.f90
- CMake target for opticalProperties

### Phase 2: Core optics (QW interband absorption)
- `main_optics.f90` with QW k_par sweep
- Accumulation framework with commutator velocity
- Absorption (TE/TM) finalization and output
- Test against old code for GaAs/AlGaAs QW

### Phase 3: All geometries and quantities
- Bulk optics (analytical velocity)
- Wire optics (commutator velocity)
- Spontaneous emission accumulation
- Gain with quasi-Fermi levels
- ISBT for QW and wire

### Phase 4: Spin resolution
- `spin_projection.f90` module
- Spin-resolved accumulation
- Spin-resolved quasi-Fermi levels
- Output spin-up/down spectra

### Phase 5: Documentation
- Update `docs/lecture/06-optical-properties.md` with commutator formalism
- Add regression test configs for optical properties
- Update `docs/reference/input-reference.md` with optics block

## Validation

1. **Bulk GaAs**: Compare absorption edge with analytical Elliott formula
2. **QW GaAs/AlGaAs**: Compare absorption spectra with old code output
3. **QW ISBT**: Verify oscillator strength sum rule
4. **Wire**: Verify spectra converge to bulk for large diameters
5. **Gain**: Verify population inversion produces negative alpha
6. **Spin-resolved**: Verify spin-up + spin-down = total spectrum
