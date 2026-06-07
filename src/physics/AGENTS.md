# Physics Engine

Owns the 8-band k·p Hamiltonian construction, all derived physical quantities (optics, g-factor, strain, topology, excitons, scattering), and self-consistent iteration. Does NOT own: material parameters (`src/core/parameters.f90`), FD stencils (`src/math/finitedifferences.f90`), eigensolvers (`src/math/eigensolver.f90`), I/O (`src/io/`), or input parsing (`src/io/input_parser.f90`).

## Entry Points

Called from `src/apps/` and `src/core/simulation_setup.f90`:
- `ZB8bandQW`, `ZB8bandBulk`, `ZB8bandLandau` (`hamiltonianConstructor.f90`) — dense Hamiltonian
- `ZB8bandGeneralized` (`hamiltonian_wire.f90`) — sparse wire Hamiltonian
- `build_bdg_hamiltonian_1d` (`bdg_hamiltonian.f90`) — Nambu-space Hamiltonian
- `self_consistent_loop`, `self_consistent_loop_wire` (`sc_loop.f90`) — SC iteration
- `gfactorCalculation`, `gfactorCalculation_wire` (`gfactor_functions.f90`) — g-factor
- `optics_accumulate`, `compute_isbt_absorption` (`optical_spectra.f90`) — optical spectra
- `compute_chern_qwz`, `compute_z2_fukane`, `compute_berry_curvature_lattice` (`topological_analysis.f90`) — topological invariants
- `compute_exciton_binding` (`exciton.f90`) — variational exciton
- `compute_phonon_scattering` (`scattering.f90`) — LO-phonon rates

## Module Inventory (17 files)

| File | Module | Lines | Role |
|------|--------|------:|------|
| `hamiltonian_blocks.f90` | `hamiltonian_blocks` | 260 | 52-entry k·p block table (KP_Q..KP_A). Single source of truth. |
| `confinement_init.f90` | `confinement_init` | 863 | kpterms FD operator matrices for QW/wire/Landau |
| `hamiltonianConstructor.f90` | `hamiltonianConstructor` | 758 | Dense bulk/QW/Landau Hamiltonian + velocity matrices |
| `hamiltonian_wire.f90` | `hamiltonian_wire` | 1357 | Sparse wire Hamiltonian (CSR/COO) + workspace cache |
| `bdg_hamiltonian.f90` | `bdg_hamiltonian` | 432 | BdG Nambu-space (16N×16N) with s-wave pairing |
| `magnetic_field.f90` | `magnetic_field` | 134 | Zeeman splitting + Peierls phase as COO insertions |
| `strain_solver.f90` | `strain_solver` | 1206 | Biaxial strain, plane-strain PDE, Bir-Pikus, Zeeman table |
| `gfactor_functions.f90` | `gfactorFunctions` | 1217 | Lowdin partitioning, spin matrices, optical matrix elements |
| `optical_spectra.f90` | `optical_spectra` | 1058 | Absorption (TE/TM), gain, spontaneous emission, ISBT |
| `charge_density.f90` | `charge_density` | 437 | n(z), p(z) from eigenstates; output in cm⁻³ |
| `poisson.f90` | `poisson` | 517 | 1D Thomas + 2D PARDISO Poisson; box-integration |
| `sc_loop.f90` | `sc_loop` | 1089 | SC Schrödinger-Poisson with DIIS; Fermi bisection |
| `topological_analysis.f90` | `topological_analysis` | 1433 | Chern, Z2, Berry curvature, Majorana, BHZ wire |
| `green_functions.f90` | `green_functions` | 501 | Spectral function, Landauer transmission, LDOS |
| `exciton.f90` | `exciton_solver` | 529 | Variational exciton (Bastard), Sommerfeld enhancement |
| `scattering.f90` | `scattering_solver` | 442 | LO-phonon Fröhlich intersubband scattering |
| `spin_projection.f90` | `spin_projection` | 71 | Spin-up/down weights, band character decomposition |

## Dependency DAG

```
Layer 0 (leaves):  hamiltonian_blocks, strain_solver, magnetic_field,
                   spin_projection, confinement_init, charge_density,
                   poisson, exciton, scattering
Layer 1:           hamiltonian_wire, optical_spectra, gfactorFunctions
Layer 2 (hubs):    hamiltonianConstructor, green_functions, sc_loop
Layer 3:           bdg_hamiltonian, topological_analysis
```

No cycles. `hamiltonianConstructor` depends on `hamiltonian_wire` (not vice versa).

## Contracts & Invariants

### Basis ordering (NEVER change)
```
Band 1: |3/2,+3/2⟩ HH↑    Band 2: |3/2,+1/2⟩ LH↑
Band 3: |3/2,-1/2⟩ LH↓    Band 4: |3/2,-3/2⟩ HH↓
Band 5: |1/2,+1/2⟩ SO↑    Band 6: |1/2,-1/2⟩ SO↓
Band 7: |S,+1/2⟩  CB↑     Band 8: |S,-1/2⟩  CB↓
```

### Band-major spatial layout
All modules: `idx = (band-1)*Ngrid + spatial_index`. Spatial recovered via `sp = mod(idx-1, Ngrid) + 1`.

### Single source of truth tables
- **k·p blocks**: `get_kp_block_table()` in `hamiltonian_blocks.f90` — 52 entries. Both dense and COO builders import this.
- **Bir-Pikus**: `compute_bp_scalar()` in `strain_solver.f90` — `elemental pure`. Never duplicate.
- **Strain table**: `get_strain_table()` in `strain_solver.f90` — band-pair topology.
- **Zeeman table**: `get_zeeman_table()` in `strain_solver.f90` — g-multipliers per band.

### kpterms sign convention
`kpterms(ii, jj, term_idx) = -result(ii, jj)` — stored with negative sign. All consumers must account for this.

### BdG Nambu structure
`H_BdG = [[H₀-μI, Δ],[Δ†, -H₀ᵀ+μI]]`. IS Hermitian. Pairing: `Δ = δ₀(iσ_y⊗I₄)`. Kramers: 1↔4, 2↔3, 5↔6, 7↔8. Layout: 1..8N = electron, 8N+1..16N = hole.

### Zeeman g-multipliers (0-based)
HH±3/2: ∓1.5, LH±1/2: ±0.5/∓0.5, SO±1/2: ∓0.5, CB±1/2: ±1.0/∓1.0

### Unit conventions
- Charge density output: cm⁻³ (converted from 1/nm³ via 1e21)
- Carrier density input to optics: cm⁻² (converted to Å⁻² internally)
- Poisson: box-integration (Birner et al. 2006), half-point dielectric via arithmetic average

## Patterns

### Adding a new physics observable
1. Add computation in a new or existing module in `src/physics/`
2. Call from the appropriate `main_*.f90` in `src/apps/`
3. Wire through `simulation_setup.f90` if needed at k-sweep level
4. Add `# COVERAGE:` annotations in corresponding test file
5. Add golden-reference test in `tests/regression/`

### Three Hamiltonian construction paths
- Dense QW/bulk: `hamiltonianConstructor` → `ZB8bandQW`/`ZB8bandBulk`/`ZB8bandLandau`
- Sparse wire: `hamiltonian_wire` → `ZB8bandGeneralized`
- BdG Nambu: `bdg_hamiltonian` → `build_bdg_hamiltonian_1d`/`_qw` (wraps wire path)

### Velocity matrices
Commutator-based: `v_α = -i[r_α, H]` element-wise. Generic interface `build_velocity_matrices` dispatches to `_2d` (wire) or `_1d` (QW).

## Anti-patterns

- Never hardcode k·p block indices — use `get_kp_block_table()` constants
- Never duplicate Bir-Pikus formulas — call `compute_bp_scalar()`
- Never bypass `wire_workspace` cache when building wire Hamiltonians across k-points
- Never change spin matrix phase convention — `SIGMA_X/Y/Z` in `gfactor_functions.f90` match `ZB8bandBulk`

## Pitfalls

- `kpterms` stored with negative sign — easy sign error if forgotten
- `gfactorFunctions` ≠ filename `gfactor_functions.f90` — module name has no underscore/different case
- `hamiltonianConstructor` ↔ `hamiltonianConstructor.f90` — module name matches filename but is inconsistent with other modules
- Wire strain PDE DOF ordering: comp=1 is u_y, comp=2 is u_z (not x,y,z)

## Related Context

- Material parameters: `src/core/parameters.f90`
- FD stencils: `src/math/finitedifferences.f90`
- Eigensolvers: `src/math/eigensolver.f90`
- Simulation orchestration: `src/core/simulation_setup.f90`
- Test infrastructure: `tests/integration/AGENTS.md`
- Lecture scripts: `scripts/AGENTS.md`
