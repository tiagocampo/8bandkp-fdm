# Phase 5: Bulk Landau Levels (Confinement=3) — Design Spec

## Overview

New confinement mode (`confinement=3`) that computes orbital Landau levels for bulk semiconductors in arbitrary magnetic field directions. Uses x-direction finite-difference discretization with minimal coupling substitution to build 8Nx8N Hamiltonians, diagonalized with dense LAPACK.

## Physics

### Gauge Choice

Landau gauge: **A = (0, Bz·x, −By·x)**

- Captures orbital effects from B_y and B_z components
- B_x (parallel to discretization direction) contributes Zeeman spin splitting only
- Minimal substitution: Π_y(x) = k_y + eBz·x/ℏ, Π_z(x) = k_z − eBy·x/ℏ

### QW-to-Landau Mapping

| QW (confinement=1) | Landau (confinement=3) |
|---------------------|------------------------|
| FD direction: z | FD direction: x |
| In-plane: k_x, k_y (scalar) | Parameters: Π_y(x), Π_z(x) (position-dependent) |
| k_z (FD derivative) | k_x (FD derivative, unchanged by gauge) |

### Position-Dependent k.p Terms

The 8-band k.p terms that use k², k±, k_x²−k_y², k_x·k_y become position-dependent:

- Diagonal blocks (at x_i): k²(i) = Π_y(i)² + Π_z(i)²
- Off-diagonal FD-hopping (x_i ↔ x_{i±1}): evaluate at midpoint x_mid = (x_i + x_{i±1})/2
- FD stencil coefficients (kpterms 5-9): unchanged (A_x = 0, no Peierls phase on x-hopping)

### Zeeman Splitting

Full |B| = √(Bx² + By² + Bz²) used for spin splitting via `compute_zeeman_vz`. Already implemented pattern reused from `ZB8bandQW`.

### Boundary Conditions

Open boundaries (same as QW). Grid must span enough magnetic lengths l_B = √(ℏ/eB) for wavefunctions to decay before reaching edges. At B=5T: l_B ≈ 11.5 nm.

## Architecture

### New Files/Modifications

| File | Change |
|------|--------|
| `src/core/defs.f90` | Allow confinement=3; add `landau_width`, `landau_nx`, `landau_sweep` fields to `simulation_config`; decouple B-field from bdg%enabled |
| `src/physics/confinement_init.f90` | New `confinementInitialization_landau(material, params, x_grid, profile, kpterms, FDorder)` |
| `src/physics/hamiltonianConstructor.f90` | New `ZB8bandLandau(HT, wv, profile, kpterms, x_grid, cfg)` subroutine |
| `src/physics/magnetic_field.f90` | New `compute_gauge_shifts(x_grid, B_vec, ky, kz, Pi_y, Pi_z)` helper |
| `src/apps/main.f90` | Add `confinement==3` dispatch branch |
| `src/io/input_parser.f90` | Parse `landau_width`, `landau_nx`, `landau_sweep`; accept confinement=3; decouple B-field from bdg%enabled |
| `tests/unit/test_landau.pf` | Unit tests for gauge and assembly |
| `tests/regression/configs/landau_InAs.cfg` | InAs at B=5T regression config |
| `tests/regression/configs/landau_GaAs.cfg` | GaAs at B=5T regression config |
| `tests/regression/configs/landau_InAs_Bsweep.cfg` | B-sweep for fan diagram |

### Data Flow

```
input.cfg (confinement=3, b_field, material, landau_*)
  → input_parser → simulation_config
  → confinementInitialization_landau → kpterms(N,N,10), profile(N,3)
  → k-point loop (k_y, k_z) or B-sweep:
      → ZB8bandLandau:
          → compute_gauge_shifts → Π_y(i), Π_z(i), Π_y(mid), Π_z(mid)
          → build Q, T, S, R, A, PZ, PP, PM with per-point k-values
          → add profile + Zeeman
      → zheevx diagonalization
      → store eigenvalues
  → output files (eigenvalues vs sweep parameter)
```

### ZB8bandLandau Assembly

New subroutine following the same pattern as `ZB8bandQW` but with position-dependent k-values:

1. Compute gauge shifts: Π_y(i) = k_y + eBz·x_i/ℏ, Π_z(i) = k_z − eBy·x_i/ℏ
2. For each grid point pair (i, j):
   - If i == j: use Π_y(i), Π_z(i) for k², k±
   - If |i − j| ≤ FDorder: use midpoint values for off-diagonal FD terms
3. Build 8x8 block matrices: Q, T, S, SC, R, RC, PZ, PP, PM, A
4. Assemble 8Nx8N Hamiltonian (same block structure as QW)
5. Add band-edge profile (uniform for bulk)
6. Add Zeeman splitting (full |B|)

Key difference from `ZB8bandQW`: k², k±, k_x²−k_y², k_x·k_y are arrays indexed by grid point instead of scalars. The assembly loop structure changes from scalar multiplication to per-element computation.

### confinementInitialization_landau

Sets up kpterms for a single homogeneous material:
- Diagonal terms (1-4, 10): constant material parameters
- FD stencil terms (5-9): standard tridiagonal (order 2) or higher-order stencil

No material profile variation — all grid points use the same material. The profile array is constant: profile(:,1) = EV, profile(:,2) = EV − ΔSO, profile(:,3) = EC.

## Input Config

```
confinement: 3              ! Landau mode
material1: InAs             ! Single bulk material
landau_width: 2000.0        ! Domain width in Angstrom (default 2000 = 200 nm)
landau_nx: 100              ! Number of x-grid points (default 100)
b_field: 0.0 0.0 5.0        ! Bx By Bz in Tesla

! Sweep mode (replaces waveVector for Landau mode)
landau_sweep: ky             ! ky, kz, or B (default: ky)
waveVectorMax: 0.05
waveVectorStep: 0.001

! B-sweep (when landau_sweep: B)
b_sweep: 0.5 10.0 0.5       ! B_min B_max B_step (reuses bdg%B_sweep)

numcb: 4
numvb: 4
FDorder: 2
```

**Important:** When `confinement=3`, the parser must NOT set `cfg%bdg%enabled = .true.` upon parsing `b_field:`. The B_vec values are stored in `cfg%bdg%B_vec` (reusing existing storage) but BdG superconductivity remains disabled. A new boolean `cfg%landau_enabled` or checking `confinement==3` gates the Landau code path.

## Output Modes

Three sweep modes controlled by `landau_sweep` field (not `waveVector`, to avoid overloading the k-direction parser):

| Mode | landau_sweep | Sweep | Fixed | Output |
|------|-------------|-------|-------|--------|
| k_y degeneracy | `ky` | k_y | k_z=0, B fixed | E vs k_y (flat lines expected) |
| k_z dispersion | `kz` | k_z | k_y=0, B fixed | E vs k_z (Landau level dispersion) |
| B fan diagram | `B` | B | k_y=0, k_z=0 | E_n vs B (fan diagram) |

**Fan diagram output** (`output/landau_fan.dat`): columns `B E_0 E_1 ... E_n`. Effective mass extracted from slope: m* = ℏeB/(∂E/∂n).

## Testing

### Unit Tests (`test_landau.pf`)

1. `test_gauge_correction_zero_B` — Π_y = k_y, Π_z = k_z when B = 0
2. `test_gauge_correction_nonzero_B` — Verify Π_y(i) = k_y + eBz·x_i/ℏ for known B, x
3. `test_landau_hermiticity` — 8Nx8N Hamiltonian is Hermitian at nonzero B
4. `test_landau_recovers_bulk` — At B=0, Landau eigenvalues match bulk band structure

### Regression Tests

1. **InAs B=5T** (`landau_InAs.cfg`): CB Landau spacing ≈ ℏω_c = 22.26 meV, Zeeman ≈ 0.58 meV
2. **GaAs B=5T** (`landau_GaAs.cfg`): ℏω_c ≈ 8.67 meV
3. **InAs B-sweep** (`landau_InAs_Bsweep.cfg`): E_n ∝ B (linear fan diagram)

### Verification

Extend `scripts/verify_landau_levels.py` to read computed Landau eigenvalues and compare against analytical E_n = (n + ½)·ℏω_c + g·μ_B·B·s.

## Grid Sizing

Default: `landau_nx=100`, `landau_width=2000 Å` (200 nm). This covers ~17 magnetic lengths at B=5T (l_B ≈ 11.5 nm) with 2 nm grid spacing (Δx = 200/100 nm). User should increase `landau_nx` or `landau_width` for lower B fields (larger l_B) or higher accuracy.

Magnetic length formula: l_B = √(ℏ/(eB)) ≈ 256/√B(nm) where B is in Tesla.

## Scope Exclusions

- QW + Peierls (requires 2D discretization → use wire mode instead)
- Wire Peierls for non-BdG path (separate feature)
- Fu-Kane Z2 for QW (deferred to Phase 6)
- Higher-order FD stencils with gauge (future enhancement)
