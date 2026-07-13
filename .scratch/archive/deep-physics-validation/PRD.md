**Status**: COMPLETE (2026-07-05)

# PRD: Deep Physics Validation (Phase 20)

## Problem Statement

The solver's four executables produce eigenvalues, g-factors, optical spectra, and strained band structures — but the validation stops at eigenvalue correctness. The post-eigenvalue physics pipeline (commutator-based velocity matrices, Lowdin g-factor partitioning, optical accumulation) is tested only indirectly through regression golden files and qualitative lecture-script checks. No test validates that the velocity matrices satisfy the Thomas-Reiche-Kuhn sum rule, that QW dense Hamiltonians are Hermitian, or that wire strain band-edge shifts match analytical Bir-Pikus predictions. The commutator-based velocity operator is the most subtle algebraic operation in the codebase, underpinning g-factor, optical absorption, gain, and ISBT — yet it has no algebraic identity check.

## Solution

Implement four validation capabilities that close the deepest gaps in the testing hierarchy:

1. **QW and wire Hermiticity with block-level diagnostics** — assert H = H† for all geometry+physics combinations currently lacking direct checks, with diagnostics that report which 8×8 block pair violates Hermiticity
2. **TRK sum rule via resolvent identity** — validate commutator-based velocity matrices against the double commutator ½⟨0|[x,[H,x]]|0⟩, computed as a linear system solve (no full eigenbasis needed)
3. **Verification ladder rungs 7-8** — extend the existing 6-rung ladder to cover g-factor (Roth analytical formula) and optical observables (Kane Ep, absorption edge)
4. **Wire strain quantitative validation** — replace the 60% tolerance profile test with 5-10% band-edge shift and HH-LH splitting checks against analytical Bir-Pikus in the wide-core limit

## User Stories

### Hermiticity (Idea 8)

1. As a physicist running QW bandStructure, I want the 8N×8N dense Hamiltonian to be verified Hermitian to machine precision, so that I trust the eigenvalues from LAPACK zheevd are real
2. As a physicist running QW with strain enabled, I want the strained QW Hamiltonian verified Hermitian, so that the Bir-Pikus insertion doesn't break H = H†
3. As a physicist running QW with Zeeman splitting, I want the magnetic QW Hamiltonian verified Hermitian, so that the Zeeman block insertion is algebraically correct
4. As a physicist running wire with strain enabled, I want the strained CSR Hamiltonian verified Hermitian, so that the COO strain insertion preserves CSR Hermiticity
5. As a physicist running wire with Zeeman splitting, I want the magnetic CSR Hamiltonian verified Hermitian, so that the Zeeman COO insertion is algebraically correct
6. As a physicist running wire g3 mode (dH/dkz for velocity matrices), I want the g3 Hamiltonian verified Hermitian, so that the velocity operator derivation is on a sound footing
7. As a developer debugging a Hermiticity failure, I want the test to report which 8×8 block pair (i,j) vs (j,i) violates the condition and by how much, so that I can pinpoint the sign or conjugation error in the k.p block table
8. As a developer adding a new k.p term, I want the block-level Hermiticity check to catch missing conjugate pairs, so that I don't ship a non-Hermitian Hamiltonian

### TRK Sum Rule (Idea 7)

9. As a physicist validating the bulk k.p model, I want the velocity matrix elements (dH/dk) to satisfy the TRK identity against d²H/dk², so that I know the k.p coupling is algebraically self-consistent
10. As a physicist validating the QW commutator velocity matrices, I want the z-velocity to satisfy the TRK identity against the double commutator, so that I know the FD-discretized commutator is correct
11. As a physicist validating the wire commutator velocity matrices, I want the transverse velocities (vx, vy) to satisfy the TRK identity, so that I know the CSR element-wise commutator is correct
12. As a developer modifying the velocity matrix construction, I want an automated test that catches any corruption of the commutator algebra, so that I don't silently break g-factor, optical absorption, and ISBT simultaneously
13. As a physicist reviewing the code, I want the TRK completeness ratio (sum/hbar²/2m₀) reported as a diagnostic, so that I can see how much oscillator strength the 8-band model concentrates vs. the free-electron value

### Verification Ladder Extension (Idea 5)

14. As a physicist validating g-factorCalculation, I want bulk CB g-factor compared against the Roth analytical formula for multiple materials, so that I know the Lowdin partitioning code path is correct
15. As a physicist validating QW g-factorCalculation, I want QW g-factor compared against Roth with confinement correction, so that I know the QW-specific g-factor code path is correct
16. As a physicist running opticalProperties on a QW, I want the absorption onset energy verified against Eg + ECB1 + EVB1 (self-consistent with ladder rung 3), so that I know the optical pipeline produces physically meaningful spectra
17. As a physicist validating bulk optical properties, I want the Kane matrix element Ep extracted from dH/dk and compared against Vurgaftman tabulated values, so that I know the interband coupling strength is correct
18. As a developer adding new optical features, I want the ladder rungs to catch regressions in the optical accumulation pipeline, so that I don't break existing spectra

### Wire Strain Quantitative (Idea 6)

19. As a physicist validating wire strain, I want the interior band-edge shift compared against analytical Bir-Pikus at 5-10% tolerance (not 60%), so that the strain solver produces quantitatively correct physics
20. As a physicist validating wire strain, I want the HH-LH splitting under strain compared against the analytical Bir-Pikus prediction, so that the valence band deformation is correct
21. As a developer modifying the Navier-Cauchy solver, I want a quantitative regression test that catches degradation in strain accuracy, so that PDE solver changes don't silently worsen physics output

### Developer Stories

22. As a developer, I want all Phase 20 tests registered under appropriate ctest labels (unit, verification, strain-validation), so that they integrate with the existing CI structure
23. As a developer, I want the minor carry-over items (Richardson absorption_edge gap, PR14 return None fix) resolved, so that no loose ends remain from previous phases

## Implementation Decisions

### Hermiticity: all 6 gaps + block-level diagnostics

Six geometry+physics combinations lack direct Hermiticity checks:

| ID | Geometry | Physics | Builder | Test location |
|----|----------|---------|---------|---------------|
| H1 | QW dense | Base H | ZB8bandQW | test_hamiltonian.pf |
| H2 | QW dense | + Strain | ZB8bandQW + apply_strain_table_dense | test_hamiltonian.pf |
| H3 | QW dense | + Zeeman | ZB8bandQW + B_vec | test_hamiltonian.pf |
| H4 | Wire CSR | + Strain | ZB8bandGeneralized + insert_strain_coo | test_hamiltonian_2d.pf |
| H5 | Wire CSR | + Zeeman | ZB8bandGeneralized + B_vec | test_hamiltonian_2d.pf |
| H6 | Wire CSR | g3 mode | ZB8bandGeneralized(g='g3') | test_hamiltonian_2d.pf |

Block-level diagnostics for QW: extract 8×8 block submatrices from the 8N×8N dense Hamiltonian (block(i,j) = rows [(i-1)*N+1 : i*N], cols [(j-1)*N+1 : j*N]), check each block pair against its conjugate transpose, report the maximum violation with block indices.

Tolerances: 1e-12 for dense, 1e-10 for CSR — consistent with existing Hermiticity tests in the codebase.

### TRK sum rule: resolvent approach with double commutator reference

The TRK sum for state |0⟩ is:

```
S_α = Σ_{n≠0} |⟨n|v_α|0⟩|² / (E_n - E_0) = ½ ⟨0|[x_α,[H,x_α]]|0⟩
```

**LHS (resolvent):**
1. b = v_α · ψ_0 (CSR SpMV)
2. Solve (H - E₀I + εP₀) y = b (PARDISO for CSR, LAPACK for dense)
3. c = v_α · y (CSR SpMV)
4. S = Re[ψ₀ᴴ · c]

Where P₀ = |0⟩⟨0| is the projector that regularizes the singular system (E₀ is an eigenvalue). The system is consistent because ⟨0|v_α|0⟩ = 0 for any eigenvector of H.

**RHS (double commutator):**
- Bulk (k-space): d²H/dk² evaluated analytically from the k.p structure
- QW/wire (real space): -(x_i - x_j)² H_ij summed over nonzero entries weighted by ψ_0

**Comparison:** Assert |LHS - RHS| / max(|RHS|, ε) < 1e-10. This is a zero-parameter, material-independent algebraic identity.

**Diagnostic:** Report S/hbar2O2m0 as the "TRK completeness ratio."

**Geometries:**
- Bulk: 8×8 dense at k=0, all 3 velocity directions (dH/dkx, dH/dky, dH/dkz), initial states: CB (band 7), HH (band 1)
- QW: 8N×8N dense, z-velocity only (commutator), initial state: CB1 (ground conduction subband)
- Wire: 8N×8N CSR, transverse velocities (vx, vy) only (commutator), initial state: CB1

**Singularity handling:** Shift (H - E₀I) by adding ε|0⟩⟨0| to the diagonal, where ε is a small number (1e-6 times the spectral gap). This removes the null eigenvalue without affecting the solution in the orthogonal subspace.

### Verification ladder rungs 7-8

**Rung 7 (g-factor):** Python script `verify_8band_rung7_gfactor.py` runs `gfactorCalculation` for:
- Bulk GaAs, InAs, InSb, GaSb: CB g-factor vs Roth formula (roth_gfactor from star_helpers.py), 2% tolerance
- QW InAs/GaAs: g-factor vs Roth with confinement correction, 5% tolerance

Materials parameters (Eg, Ep, DeltaSO) from Vurgaftman 2001 via star_helpers constants or paramDatabase lookup.

**Rung 8 (optical):** Python script `verify_8band_rung8_optical.py` runs `opticalProperties` for:
- Bulk GaAs: Extract Kane Ep from velocity matrix elements at k=0 (computed numerically from dH/dk), compare against Vurgaftman Ep = 28.8 eV, 1% tolerance
- QW GaAs/AlGaAs: Absorption onset energy vs Eg + ECB1 + EVB1 (self-consistent with rung 3 subband energies), 10 meV tolerance (Lorentzian broadening)

### Wire strain quantitative

**Configuration:** New TOML config `wire_inas_gaas_strain_wide.toml` with 100Å InAs core, 200Å GaAs domain, 40×40 grid, 5Å spacing. Strain reference = GaAs.

**Script:** New Python script `verify_strain_wire_quantitative.py` that:
- Runs bandStructure with strain and without strain
- Extracts interior band-edge shifts (CB bottom, VB top) from strained vs unstrained
- Compares against `bir_pikus_biaxial_001()` from star_helpers.py
- Checks HH-LH splitting against analytical Bir-Pikus prediction
- Tolerance: 5-10% (accounts for finite-size correction to the wide-core limit)

### Minor carry-over

- Richardson convergence helper: investigate why `absorption_edge` is absent from stored JSON results despite code being complete
- Replace `return None` in `test_bulk_zeeman.py:44` with `raise RuntimeError` for consistency

### Modules modified

| Module | Change | Issue |
|--------|--------|-------|
| `tests/unit/test_hamiltonian.pf` | Add QW Hermiticity tests (H1-H3) + block diagnostics | #1 |
| `tests/unit/test_hamiltonian_2d.pf` | Add wire Hermiticity tests (H4-H6) | #2 |
| `tests/support/csr_test_helpers.f90` | No changes needed (infrastructure exists) | — |
| `tests/unit/test_trk_sum_rule.pf` | New file: TRK sum rule tests for bulk, QW, wire | #3, #4 |
| `tests/support/trk_helpers.f90` | New file: double commutator computation, resolvent TRK solver | #3 |
| `tests/regression/configs/wire_inas_gaas_strain_wide.toml` | New wide-core wire strain config | #5 |
| `tests/integration/verify_strain_wire_quantitative.py` | New wire strain quantitative script | #5 |
| `tests/integration/verify_8band_rung7_gfactor.py` | New g-factor ladder rung | #6 |
| `tests/regression/configs/` (2-4 new configs) | g-factor test configs for rung 7 | #6 |
| `tests/integration/verify_8band_rung8_optical.py` | New optical ladder rung | #7 |
| `tests/regression/configs/` (2 new configs) | Optical test configs for rung 8 | #7 |
| `tests/integration/convergence_helpers.py` | Fix absorption_edge JSON gap | #8 |
| `validation/bulk/test_bulk_zeeman.py` | Replace return None with RuntimeError | #8 |

### Implementation slices (all AFK)

| Issue | Slice | Blocked by |
|-------|-------|------------|
| #1 | QW Hermiticity + block diagnostics | None |
| #2 | Wire Hermiticity with perturbations | None |
| #3 | TRK sum rule helper + bulk test | None |
| #4 | TRK sum rule QW + wire tests | #3 |
| #5 | Wire strain quantitative validation | None |
| #6 | Verification ladder rung 7 (g-factor) | None |
| #7 | Verification ladder rung 8 (optical) | None |
| #8 | Minor carry-over fixes | None |

## Testing Decisions

### Test philosophy

Hermiticity and TRK tests verify algebraic identities (H = H†, TRK = double commutator) that hold to machine precision regardless of material parameters or grid spacing. These are zero-parameter oracle tests — no external reference needed, no tolerance tuning required.

Ladder rung and wire strain tests verify physics against analytical references (Roth, Bir-Pikus, Vurgaftman Ep). These use the existing tolerance tiers from star_helpers.py (TOL_ANALYTICAL = 1%, TOL_NUMERICAL = 5%).

### Modules tested

- **Hermiticity:** `ZB8bandQW` (dense QW builder), `ZB8bandGeneralized` (CSR wire builder) — tested via direct H = H† assertion
- **TRK:** `build_velocity_matrices_1d`, `build_velocity_matrices_2d`, PARDISO linear solve — tested via resolvent vs double commutator identity
- **g-factor:** `gfactorCalculation` (bulk/QW Lowdin), `gfactorCalculation_wire` — tested via Roth formula comparison
- **Optical:** `opticalProperties` executable pipeline — tested via absorption edge and Kane Ep
- **Wire strain:** Navier-Cauchy solver + Bir-Pikus insertion — tested via wide-core analytical limit

### Prior art

- **Hermiticity:** `test_bulk_hermitian` in test_hamiltonian.pf (manual loop, 1e-12), `csr_hermitian_error` in csr_test_helpers.f90 (CSR, 1e-10)
- **TRK:** No prior art — new test type. Pattern: `csr_spmv` + `zdotc` from optical_spectra.f90, PARDISO from bdg_hamiltonian.f90
- **Ladder rungs:** `verify_8band_rung1_bulk_k0.py` through `verify_8band_rung6_strain_qw.py`
- **Wire strain:** `verify_strain_wire_profile.py` (60% tolerance), `verify_strain_rung5_bulk.py` (1% Bir-Pikus, uses `bir_pikus_biaxial_001`)
- **Standard-star g-factor:** `verify_star_inas_bulk.py` (Roth formula, 1%), `verify_star_insb_bulk.py` (Roth, 1%)

## Out of Scope

- **Comparing TRK sum against ℏ²/2m₀** — the 8-band model does not satisfy the free-electron TRK due to remote band contributions folded into Kane parameters. Only the double commutator identity (internal consistency) is tested.
- **Wire TRK for z-direction** — wire vel(3) is dH/dkz (k-derivative, not commutator), which has a different algebraic structure. Deferred.
- **Optical accumulation unit tests** — `optics_accumulate` and `compute_gain_qw` have no unit tests; adding them is a separate effort.
- **Spontaneous emission tests** — zero coverage but out of scope for this phase.
- **Spin-resolved absorption tests** — zero coverage but out of scope.
- **Phase 18 (architectural cleanup)** — deferred.
- **Phase 19 (extended cross-code validation)** — blocked on kdotpy API.

## Further Notes

- This work implements Phase 20 from `docs/plans/BACKLOG.md`, sourced from `docs/ideation/2026-05-08-core-kp-validation-ideation.md` (Ideas 5, 6, 7, 8 partial).
- The TRK resolvent approach avoids the FEAST incomplete-spectrum problem: it needs only ONE eigenvector (not the full basis), computed via any eigensolver.
- The double commutator reference is the same algebraic identity regardless of geometry — only the computation method differs (k-space derivative for bulk, real-space entry loop for QW/wire).
- Block-level Hermiticity diagnostics reuse the 8×8 band structure (bands 1-4 = valence, 5-6 = split-off, 7-8 = conduction) from `hamiltonian_blocks.f90`.
- Wire strain wide-core config (100Å InAs core) is intentionally small enough for fast regression testing (~10s for PARDISO strain solve + FEAST eigensolve) while large enough that the interior strain is within 5-10% of the biaxial limit.
- All Python integration tests use shared infrastructure from `tests/integration/star_helpers.py` (`run_exe`, `roth_gfactor`, `bir_pikus_biaxial_001`, physical constants).
