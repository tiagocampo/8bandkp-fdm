# Phase 6 Completion Repair Design

Date: 2026-05-05
Branch: feature/bdg-topological-superconductivity
Scope: Full Phase 6 completion for the topological suite, including repairs for review findings and replacement of placeholder feature paths with tested implementations.

## Goal

Bring Phase 6 to a working, merge-ready state:

- Dense QW BdG is correct, tested, and reachable from `topologicalAnalysis`.
- Fu-Kane QW Z2 uses a real inversion-parity operator in band-major basis.
- Berry curvature and Kubo conductance are validated and executable from app modes.
- Spectral and LDOS paths handle dense and CSR cases correctly.
- Gap sweep produces real phase diagrams, not heuristic placeholders.
- Parser, docs, outputs, examples, and tests cover all exposed Phase 6 modes.

Correctness has priority over maintainability, readability, type-exactness, performance, and minimality.

## Non-Negotiable Constraints

- Do not change the global basis ordering: bands 1-4 valence, bands 5-6 split-off, bands 7-8 conduction.
- All QW and wire eigenvector consumers must use band-major indexing:

```fortran
irow = (iband - 1) * N + isite
```

- Do not modify material parameters in `parameters.f90` unless verified against published references.
- Do not silently return physically meaningful-looking zeros for unsupported or failed calculations.
- Every executable mode must either compute its advertised result or fail with a clear error.

## 1. Correctness Foundations

This checkpoint fixes cross-cutting defects before adding more Phase 6 physics.

### Required Changes

- Add a small internal helper for band-major row indexing where topological routines consume eigenvectors.
- Update `compute_majorana_profile` and `compute_z2_fukane_qw` to use band-major indexing.
- Restore `compute_zeeman_vz` import in `bdg_hamiltonian.f90`.
- Remove stale local symbols and imports:
  - unused `nnz_zeeman`
  - unused `zheevd`
  - unused `confinementInitialization` import in `main_topology.f90` if only `_2d` is used
  - unused dummy argument `n_occ` from `compute_conductance_kubo`, or document and use it
- Use `pi_dp` instead of repeated `acos(-1.0_dp)` in new or touched Phase 6 code.
- Rename `lorntz` to `lorentz`.
- Fix misleading wire pairing comments so they match actual bands 1-8.

### Validation Guards

- Berry curvature requires `nkx >= 2`, `nky >= 2`, `n_occ >= 1`, and enough occupied-vector columns.
- Spectral and LDOS broadening require `eta > 0`.
- Gap sweep dimensions must be positive and ranges must be finite.
- BdG Majorana mode runs require physically valid `delta_0`; default `delta_0 = 0` must not create zero-width eigensolver windows.
- Every LAPACK diagonalization checks `info`.
- Failed Majorana localization fits return a distinguishable failure state, preferably negative `xi`, and callers exclude failures from averages.

### Issues Closed

Review issues 1, 2, 4, 5, 6, 7, 8, 9, 14, 20, 21, 23, 24, C1, C2, C3, I1, I3, I4, and I8.

## 2. Fu-Kane QW Z2

The existing CB/VB-weight heuristic must be replaced with an inversion-parity operator.

### Algorithm

For a symmetric QW:

1. Determine the lattice constant `a` from the active material/profile data. Do not hardcode `6.0 A`.
2. Evaluate all four TRIM:
   - Gamma: `(0, 0)`
   - X: `(pi/a, 0)`
   - Y: `(0, pi/a)`
   - M: `(pi/a, pi/a)`
3. Build and diagonalize `ZB8bandQW` at each TRIM.
4. For each occupied Kramers pair, compute inversion expectation:

```text
<psi | P_band P_env | psi>
```

where:

- `P_band = +1` for bands 1-6 and `-1` for bands 7-8.
- `P_env` maps site `i` to mirrored site `N + 1 - i`.
- Eigenvector rows use band-major indexing.

5. Convert each pair expectation to a parity sign with an explicit tolerance.
6. Multiply signs over occupied Kramers pairs, not every occupied state.
7. Compute the Fu-Kane product using all four TRIM. The implementation must document the chosen convention and include Gamma correctly.

### Symmetry Requirement

Fu-Kane parity is only valid for inversion-symmetric QW profiles. Add a profile symmetry check. If the profile is asymmetric beyond tolerance, the routine must fail clearly or return a status indicating that Fu-Kane parity is not applicable.

### API Design

The current integer-only result is insufficient for failures. Introduce a status-aware internal routine, for example:

```fortran
subroutine compute_z2_fukane_qw_result(cfg, profile, kpterms, n_occ, z2, min_gap, status)
```

The existing `compute_z2_fukane_qw(...) result(z2)` may remain as a compatibility wrapper, but app dispatch should use the status-aware routine.

### App Wiring

`topology_mode: qshe` with `confinement=1` must call the corrected QW Fu-Kane path, set `topo_result%z2_invariant`, set `topo_result%min_gap`, and write both values.

### Tests

- Unit-test the parity operator on synthetic band-major vectors with `N > 1`.
- Test even and odd envelope parity.
- Test CB/VB band parity.
- Test controlled Fu-Kane products where Gamma changes the result.
- Test that asymmetric profiles do not silently report an invariant.
- Add a small end-to-end symmetric QW Fu-Kane test or a deterministic reduced Hamiltonian fixture plus a real regression config.

### Issues Closed

Review issues 2, 13, 17, Gamma/TRIM findings, CB/VB heuristic findings, Kramers double-counting findings, QSHE QW unreachable path, and Fu-Kane coverage gaps.

## 3. BdG QW And Majorana Paths

Dense QW BdG must become a tested analogue of the wire BdG path.

### BdG Assembly

Build electron and hole partners explicitly:

```text
H_e = H(+k_par)
H_h = H(-k_par)

H_BdG = [ H_e - mu I      Delta            ]
        [ Delta^dagger   -conjg(H_h)+mu I ]
```

Pairing blocks keep the 8-band antidiagonal pattern in band-major order.

### Zeeman Ownership

Use one Zeeman source of truth:

- Preferred: Hamiltonian constructors apply Zeeman from `cfg`.
- `build_bdg_hamiltonian_qw` must not add an independent second Zeeman diagonal after calling `ZB8bandQW`.
- If `B_vec` and `g_factor` remain in the public QW BdG signature, they must be converted into a local config before constructing `H_e` and `H_h`, or documented as deprecated and ignored when `cfg%bdg` already controls the field.

### App Wiring

`topology_mode: bdg` must support:

- `confinement=2`: existing wire path
- `confinement=1`: dense QW BdG path

Both paths must validate `delta_0`, `mu`, magnetic field, eigensolver window, and expected dimensions.

### Majorana Profile

`compute_majorana_profile` must accumulate electron and hole density in band-major order. Failed fits return a failure status or negative `xi`; averages include only successful fits and report the failed count.

### Tests

- QW BdG dimensions and Hermiticity.
- QW BdG particle-hole eigenvalue pairing at nonzero `k_par`.
- QW BdG Zeeman branch with a field, verifying no double splitting.
- Synthetic band-major Majorana profile with known spatial density.
- Small QW BdG regression config.

### Issues Closed

Review issues 1, 4, 7, 8, 23, C1, I1, I3, I4, and dense QW BdG Zeeman/particle-hole test gaps.

## 4. Berry Curvature, Kubo, And Conductance

The FHS Berry curvature implementation should remain, but conductance mode must compute its prerequisites instead of using default result fields.

### Berry Curvature

- Validate grid and occupied-state dimensions before allocation.
- Allocate `overlap(n_occ,n_occ)` once outside the `nkx * nky` loop and reuse it.
- Document units:
  - `Omega(i,j)` is continuum-normalized curvature.
  - `sum(Omega) * dkx * dky / (2*pi_dp)` yields the Chern number.
- If raw plaquette phase is needed for figures, expose it as `berry_flux`, not `Omega`.

### Kubo Conductance

`compute_conductance_kubo` integrates Berry curvature directly and returns `sigma_xy` in units of `e^2/h`.

Duplicate Chern-to-conductance helpers should be collapsed into one clear routine, for example:

```fortran
compute_hall_conductance_from_chern(C)
```

### App Conductance Mode

`topology_mode: conductance` must compute the required input:

- `method = kubo_chern`: compute QWZ Chern from `qwz_u`, then convert to conductance.
- `method = kubo_berry`: compute Berry curvature on a configured grid and integrate it.
- `method = landauer` or `longitudinal`: run the implemented Green-function/lead method.

No branch may fall back to default `chern_number = 0`.

### Landauer And Longitudinal Completion

Full Phase 6 requires both advertised conductance families:

1. `compute_longitudinal_conductance`
   - Implement Kubo-Greenwood longitudinal conductance using existing velocity matrix infrastructure.
   - Return `conductance_zz` in documented units.
   - Validate broadening and occupied/unoccupied state selection.

2. `compute_conductance_landauer`
   - Implement a minimal wire-compatible Landauer path with lead self-energies and transmission `T(E)`.
   - Start with a narrow, testable effective-chain helper, then connect it to the wire app path.
   - Do not silently route this method to Kubo.

Both methods must have analytic small-system tests. Method names, units, and limitations must be documented.

### Results

Resolve `hall_conductance` vs `conductance_xy` ambiguity:

- Prefer writing `conductance_xy` for Kubo/Landauer outputs.
- Keep `hall_conductance` only as a backward-compatible alias if needed.
- `topology_result.dat` must write the scalar that was actually computed.

### Tests

- Degenerate Berry grids reject.
- Berry curvature integrates to QWZ Chern.
- `compute_conductance_kubo` is directly tested.
- Conductance executable regression gives nonzero quantized output for a known Chern case.
- Landauer or longitudinal conductance has a small analytic test.

### Issues Closed

Review issues 3, 5, 10, 11, 12, 14, 21, 24, conductance default-zero, duplicate conductance semantics, and Kubo coverage gaps.

## 5. Spectral Function, LDOS, And Green Functions

Dense spectral broadening and CSR Green-function routines must be validated separately.

### QW Spectral Function

- Validate `eta > 0`.
- Check `zheev info` after every diagonalization.
- Keep full-band Lorentzian spectral sum unless a caller requires filtering.
- Document normalization: energy integral approaches number of states when the energy window covers the relevant spectrum.

### CSR LDOS

Fix shifted matrix assembly:

```fortran
a_val = -H%values
! add E + i*eta only to diagonal entries
```

If the CSR matrix lacks structural diagonal entries, either insert them before factorization or reject the matrix clearly. Do not add the energy shift to off-diagonal nonzeros.

### Spectral Modes

`topology_mode: spectral` must support:

- QW dense spectral function.
- Wire CSR/eigensolver or Green-function spectral function.
- Bulk 8x8 dense spectral function.

Because the backlog says `compute_spectral_function`, bulk, QW, and wire spectral paths are all in scope.

### Output

- Large heatmaps go to `spectral_function.dat`.
- Header includes dimensions, units, `eta`, and k/E grid ranges.
- `topology_result.dat` records spectral metadata only, not the full heatmap.

### Tests

- `eta = 0` rejection.
- QW spectral sum-rule over a sufficiently broad energy window.
- LDOS CSR small off-diagonal matrix test catches shift-on-every-nonzero regression.
- QW and wire executable/regression configs.

### Issues Closed

Review issues 6, 16 where parser affects spectral fields, 20, C2, LDOS shifted-matrix bug, missing wire spectral path, spectral output coverage, and spectral sum-rule gaps.

## 6. Gap Sweep And Phase Diagrams

The gap sweep must become evaluator-backed and physically meaningful.

### Shared Sweep Driver

Replace duplicated loops with a shared internal driver:

```text
for each mu in mu_grid
  for each B in B_grid
    evaluator(B, mu) -> z2, gap, status
```

Outputs:

- `z2_map(nMu,nB)`
- `gap_map(nMu,nB)`
- `transitions(:,2)` or richer transition records

### Evaluators

Implement these evaluators:

1. **BHZ analytic evaluator**
   - Fast model path for tests and figures.
   - Named explicitly as BHZ analytic, not generic wire physics.
   - Tests known transition at `M = 0`.

2. **Wire BdG evaluator**
   - Builds wire BdG Hamiltonian for each `(B, mu)`.
   - Solves near zero energy.
   - Computes gap as the lowest positive particle-hole gap or `min(abs(E))`, with convention documented.
   - Classifies finite-wire topology with a documented finite-system Majorana criterion: near-zero modes below `gap_threshold`, particle-hole pairing, and edge-localized density at both ends.
   - This finite-wire classifier must not be named a bulk invariant. If a periodic-wire Pfaffian/scattering invariant is added later, it can replace this evaluator behind the same status-aware interface.

3. **QW Fu-Kane evaluator**
   - Updates local config for `(B, mu)`.
   - Computes corrected Fu-Kane Z2.
   - Computes the direct occupied/unoccupied gap from the TRIM spectra.

### Transition Detection

- Scan neighboring cells in both B and mu directions.
- Record transitions when Z2 changes or gap crosses `gap_threshold`.
- Use `mu` and `gap_threshold`; they may not be ignored.
- `gap_map` may be zero only when computed physical gaps are zero.

### Output And Figures

- Write `z2_phase_diagram.dat` with columns `B`, `mu`, `z2`, `gap`.
- Write transition records with units.
- Figure scripts read these real outputs.

### Tests

- Shared sweep driver tested with a fake evaluator.
- BHZ analytic transition at `M=0`.
- Wire/QW tiny-grid smoke tests assert finite non-placeholder gaps.
- Regression for `topology_mode: sweep` checks file columns and non-placeholder gap values.

### Issues Closed

Review issues 18, 19, I5, placeholder `gap_map`, ignored `mu`, ignored `gap_threshold`, heuristic-blessing tests, and phase-diagram figure correctness.

## 7. Parser, App Dispatch, Output, Docs, And Tests

Phase 6 must be usable from `topologicalAnalysis`, not only from unit routines.

### Parser

- Preserve old `topology: T` compatibility.
- Replace fragile optional-field peek/backspace handling with name-aware optional parsing for Phase 6 topology fields.
- Validate `topology_config` and `bdg_config` after parsing.
- Accept only documented method names:
  - `kubo_chern`
  - `kubo_berry`
  - `landauer`
  - `longitudinal`
- Resolve control-flow ambiguity:
  - Either mode drives execution and `compute_*` flags gate subfeatures within that mode, or flags are removed from the user-facing docs.
  - The implementation must not parse flags that are then ignored.

### App Dispatch

Required executable behavior:

- `qhe`: Chern, Berry curvature, and Kubo conductance when requested.
- `qshe`: wire/BHZ and QW Fu-Kane paths reachable.
- `bdg`: wire and QW BdG paths reachable.
- `spectral`: QW, wire, and bulk paths where implemented.
- `conductance`: computes prerequisites, never uses default result fields.
- `sweep`: evaluator-backed real sweep.

### Output

`topology_result.dat` writes all scalar fields actually computed:

- `chern_number`
- `z2_invariant`
- `conductance_xy`
- `conductance_zz`
- `min_gap`
- Majorana count
- Majorana fit failures
- edge localization summaries

Large arrays use dedicated files:

- `berry_curvature.dat`
- `spectral_function.dat`
- `z2_phase_diagram.dat`
- `gap_map.dat` if separate

Headers must include units and grid dimensions.

### Documentation And Examples

Update `docs/reference/input-reference.md` with:

- new topology modes
- conductance methods
- spectral grids
- gap sweep grids
- QW Fu-Kane symmetry requirement
- output files and units

Add example configs under `tests/regression/configs/` for:

- QW Fu-Kane
- QW BdG
- spectral QW
- spectral wire
- conductance
- sweep

### Test Requirements

No test may assert known placeholder behavior. Add coverage for:

- `compute_z2_fukane_qw`
- Fu-Kane parity operator and product
- QW BdG with Zeeman
- QW BdG particle-hole symmetry at nonzero `k_par`
- `compute_conductance_kubo`
- new executable modes
- old and new topology parser block formats
- spectral sum rule
- CSR LDOS shifted matrix
- gap sweep with nonzero gaps

### Maintainability

If `topological_analysis.f90` grows further during implementation, split focused code into smaller modules:

- `topology_chern`
- `topology_z2`
- `topology_phase_diagram`

Keep a stable facade if existing callers expect `use topological_analysis`.

Review dead exports such as `add_zeeman_coo`. If retained, mark as legacy or fix its indexing before any future use.

### Issues Closed

Review issues 3, 10, 11, 12, 15, 16, 19, 22, parser/order fragility, agent-native gaps, dropped result fields, documentation gaps, and executable-mode coverage gaps.

## Implementation Checkpoints

1. Foundation fixes and tests.
2. Fu-Kane QW parity rewrite and QSHE app dispatch.
3. Dense QW BdG, Zeeman ownership, and QW BdG app dispatch.
4. Berry curvature validation and Kubo/conductance completion.
5. Spectral/LDOS fixes and wire/bulk spectral support.
6. Real gap sweep evaluators and phase diagram outputs.
7. Parser, docs, examples, regression tests, cleanup.

Each checkpoint must build and run focused tests before proceeding. Full `ctest` must pass before Phase 6 is considered complete.

## Success Criteria

- Fresh build succeeds from a clean CMake/Ninja tree.
- Unit tests cover every listed P0/P1 defect and every new Phase 6 API.
- Regression tests exercise every new executable mode.
- No Phase 6 output is a placeholder unless explicitly labeled as an analytic model output.
- `docs/reference/input-reference.md` matches parser behavior.
- `topology_result.dat` and dedicated data files contain the computed results with units.
- All known review findings listed in this spec are fixed, removed by design, or documented as intentionally unsupported with executable errors.
