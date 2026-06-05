---
date: 2026-05-17
topic: core-kp-validation
focus: Validating core k.p physics against published literature, independent codes, and internal consistency checks
mode: repo-grounded
---

# Ideation: Core k.p Physics Validation

## Grounding Context

**Codebase Context (updated June 2026):** Fortran 2018 k.p FD code with 4 executables. 111 tests: 35 unit (pFUnit), 44 regression, 16 verification (8 ladder rungs + 2 strain rungs + 2 wire strain), 7 standard-star benchmarks (S1-S7), 6 convergence, 2 strain-validation, 1 coverage. 88 regression configs. 15 lecture scripts (L00-L14). Verification ladder (8 rungs), standard-stars (S1-S7: GaAs/InAs/InSb bulk, GaAs/AlGaAs QW, InAs/GaSb QW, InAs/GaAs strained QW, InAs wire), CSR structure + Krylov snapshot tests, coverage matrix with validation_universe.yml.

**Remaining Gaps:** Richardson extrapolation only on FD stencil coefficients, not physics observables. Wire strain: profile test but no physics observable validation. Optics wire validation thin, no exciton binding energy regression. SC loop tests are "does it converge" checks. Standard-star suite missing g-factor, optical spectra, topological invariants. No cross-code comparison capability.

**Past Learnings:** CB effective mass from 8-band model ≠ Vurgaftman (29% deviation is expected model physics, not a bug). Lecture-test review found 13 systemic failures (vacuous assertions, material mismatch, FEAST window misses, hardcoded PASS). Richardson convergence: max-rate strategy, order-adaptive tolerance. FD coefficient bugs at order-10 were silently wrong. CSR test infrastructure had tautological assertions. Automated review has 80%+ false-positive rate on physics code. Last-layer-wins overlapping material layer pitfall in SC configs.

**External Context:** kdotpy (Wuerzburg, SciPost Phys. Codebases 47, 2025) is closest open-source 8-band k.p comparator — Python, plane-wave discretization, same Vurgaftman parameters. No shared inter-code benchmark suite exists (CECAM 2024 workshop identified gap). Gershoni 1993 is canonical optical benchmark. WannierTools has BdG with Fu-Kane analytical benchmarks. Richardson extrapolation theoretically proven for FD eigenvalue problems (E(h) = E_exact + c_2 h^2 + c_4 h^4 + ...).

## Previously Ranked Ideas (May 2026 Status Update)

### 1. 8-Band Verification Ladder — IMPLEMENTED
**Status:** 4 rungs (bulk k=0, bulk dispersion, QW subbands, wire eigenvalues) + 2 strain rungs. All passing in CI.

### 2. CSR Structure Testing — IMPLEMENTED
**Status:** csr_spmv, csr_helpers, csr_structural unit tests + Krylov snapshot tests for SC loop, optics wire, gfactor wire.

### 3. Standard-Star Benchmark Systems — IMPLEMENTED
**Status:** S1-S7 (GaAs bulk, InAs bulk, InSb bulk, GaAs/AlGaAs QW, InAs/GaSb QW, InAs/GaAs strained QW, InAs wire).

### 4. Executable Lecture-Test Pairs — IMPLEMENTED
**Status:** 14 lecture scripts (L00-L14) serving dual pedagogy + validation roles.

### 5. Validation Coverage Matrix — IMPLEMENTED
**Status:** validation_universe.yml (49 required + 12 aspirational cells), 61 COVERAGE annotations, ctest -L coverage.

### 6. Richardson Extrapolation Convergence Fixture — IMPLEMENTED
**Status:** convergence_helpers.py (416 lines) + 6 ctest convergence tests covering CB1, subband spacing, effective mass, g-factor, absorption edge, Fermi level, subband shift, charge integral, exciton binding energy.

### 7. Strain Validation Across Geometries — IMPLEMENTED
**Status:** Bulk + QW rungs + wire profile test + wire strain quantitative (R10a-c qualitative + R11 quantitative Bir-Pikus gap shift at 15%).

## Ranked Ideas (May 2026 — Continued Ideation)

### 1. kdotpy Cross-Code Comparison — IMPLEMENTED
**Status:** 12/12 pipeline tests passing. Bulk k=0, dispersion, QW subbands/dispersion/convergence, strain bandedge (cross-code). Zeeman, Landau, wire, g-factor, strain QW, SC QW (analytical/single-code). Merged via `33859cc` + `0f5adc7`. Extended items (wire cross-code, Landau fan, Berry curvature) blocked on kdotpy API.

### 2. Richardson Extrapolation on Physics Observables — IMPLEMENTED
**Status:** convergence_helpers.py (416 lines) + 6 ctest convergence tests. Observables: CB1, subband spacing, effective mass, g-factor, absorption edge, Fermi level, subband shift, charge integral, exciton binding energy. Merged via `0f5adc7`.

### 3. Exciton Binding Energy Validation — IMPLEMENTED
**Status:** Covered by convergence_exciton ctest (Richardson extrapolation on exciton binding energy). Merged via `0f5adc7`.

### 4. SC Loop Physics-Output Regression — PARTIALLY IMPLEMENTED
**Status:** convergence_sc ctest validates Fermi level convergence. SC convergence test checks quantitative output. Full physics-output regression (charge profile integration, Bastard subband shift comparison) not yet done.

### 5. Verification Ladder Extension to Optical and g-Factor — IMPLEMENTED
**Status:** Rungs 7-8 delivered. verify_8band_rung7_gfactor.py (R7.1 bulk Roth 4 materials at 2%, R7.2 QW g-factor consistency) + verify_8band_rung8_optical.py (R8.1 Kane self-consistency, R8.2 QW absorption edge 10meV, R8.3 TE/TM ordering). Commits `b8f16a2` through `3596737`.

### 6. Wire Strain Physics Observable Validation — PARTIALLY IMPLEMENTED
**Status:** verify_strain_wire_quantitative.py delivers 4 checks: R10a (gap increases under compressive strain, qualitative), R10b (HH above LH, qualitative), R10c (measurable VB shift), R11 (quantitative Bir-Pikus gap shift at 15% tolerance). 25x25 grid with narrowed FEAST window. Full quantitative Bir-Pikus comparison for individual CB/VB shifts requires finer grids than CI-feasible.

### 7. TRK Sum Rule / Optical Invariant Gate — IMPLEMENTED
**Status:** trk_helpers.f90 (6 routines: double commutator bulk, velocity matrix, band curvature, TRK sum explicit, TRK sum CSR, double commutator CSR) + test_trk_sum_rule.pf (7 tests: GaAs CB/HH, InAs CB/HH, cubic symmetry, QW CSR, wire CSR). Commits `b8f16a2` through `3596737`.

### 8. Hermiticity Smoke Test for All Geometries — IMPLEMENTED
**Status:** 20 hermiticity tests across 6 test files: test_hamiltonian.pf (5: bulk GaAs/InAs, QW, QW+strain, QW+Zeeman), test_hamiltonian_2d.pf (6: 2D kpterms, heterostructure, strain, Zeeman, g3), test_bdg_hamiltonian.pf (4: QW, Peierls, zero field, no-Peierls), test_landau.pf (1), test_edge_states.pf (1: BHZ), test_hamiltonian_blocks.pf (1). Covers bulk, QW, wire, BdG, Landau, BHZ, blocks. Commits `b8f16a2`.
**Rationale:** The verification ladder catches wrong eigenvalues at specific k-points and materials. A dedicated Hermiticity test catches the bug at the matrix level with clear diagnostics, protecting every geometry with a single cheap unambiguous test.
**Downsides:** CSR Hermiticity test requires CSR-to-dense conversion (infrastructure exists in csr_test_helpers). BdG Nambu-space assembly is complex (16N x 16N).
**Confidence:** 90%
**Complexity:** Low
**Status:** Unexplored

## Rejection Summary

### Original Rejections (May 2026)
| # | Idea | Reason Rejected |
|---|------|-----------------|
| 1 | Per-material g-factor gap dashboard | Subsumed by Standard-Star Benchmarks (implemented) |
| 2 | Negative-validation / model limitation tests | Subsumed by Standard-Stars and 8-Band Ladder |
| 3 | Spurious mode sentinel (ellipticity) | Better as brainstorm topic for new feature |
| 4 | Parameter provenance / Vurgaftman audit | Lower priority than physics validation |
| 5 | Sengupta matrix diff (block-level comparison) | Subsumed by 8-Band Ladder and CSR Structure Testing |
| 6 | Topological test oracle (integrality checks) | Better as brainstorm topic for topological track |
| 7 | Cross-executable consistency tests | Largely assured if 8-Band Ladder passes |
| 8 | Validation-as-publication (citable tests) | Subsumed by Executable Lecture-Test Pairs |
| 9 | Forced-degradation testing | High implementation burden for uncertain value |
| 10 | Fast physics gate (30-second smoke tests) | Off-topic: DX/performance, not physics validation |
| 11 | One-command test scaffolding | Off-topic: DX tooling |
| 12 | Dimensional consistency checks | Off-topic: general code quality |
| 13 | Cross-code validation pipeline (Nextnano) | Too expensive; requires Nextnano installation |
| 14 | k.p phantoms (independent solver) | Too expensive; analytical benchmarks cover similar ground |
| 15 | Temporal regression (git history) | Too expensive; high infrastructure cost |
| 16 | Published-plot-to-tolerance pipeline | Too expensive; PDF curve extraction is separate project |
| 17 | Inverse parameter recovery | Subsumed by Standard-Stars effective mass validation |
| 18 | Proof-load benchmarks (extreme FD orders) | Subsumed by Richardson Extrapolation Fixture |
| 19 | Delta reports with drift attribution | Nice-to-have infrastructure, not essential |
| 20 | Physics-labeled regression assertions | DX improvement, not physics validation |
| 21 | Sweep-to-config auto-benchmark factory | Too ambitious / vague |
| 22 | Physics-module coverage sentinel | Subsumed by Validation Coverage Matrix |
| 23 | Analytical shadow suite | Subsumed by 8-Band Ladder rungs 1-2 |

### Continued Ideation Rejections (May 2026)
| # | Idea | Reason Rejected |
|---|------|-----------------|
| 24 | Topological Phase Diagram Sweep | Lower leverage than Richardson/TRK/kdotpy; single-point invariant tests adequate for now |
| 25 | Mutation Testing / Vacuous Assertion Detection | Meta-level test-quality idea, better as brainstorm topic |
| 26 | Kill Golden Files / Analytical Invariant Tests | Too radical; some observables lack clean analytical invariants |
| 27 | COVERAGE Annotation Lint | DX tooling, not physics validation |
| 28 | Aspirational-to-Required Promotion | Process improvement, subsumed by broader ideas |
| 29 | Geometry-Limit Consistency | Subsumed by Hermiticity + Richardson + kdotpy combination |
| 30 | Optical Gauge Consistency | Narrow, subsumed by TRK sum rule |
| 31 | Systematic Error Budget | Documentation exercise, not core validation |
| 32 | Regression Config Physics Assertion Layer | Too expensive; standard-stars already cover key systems |
| 33 | Eigenstate Orthonormality | Subsumed by Hermiticity (LAPACK returns orthonormal vectors from Hermitian H) |
| 34 | Hamiltonian Test Vectors | Subsumed by kdotpy comparison |
| 35 | Blind Analysis for Benchmarks | Process methodology, difficult to enforce |
| 36 | Negative-Space Testing | Conflicts with CLAUDE.md boundary on parameter modification |
| 37 | 30-Second Smoke Test Tier | DX improvement, not physics validation |
| 38 | Spurious-Mode Sentinel | Tests model pathology, better as feature brainstorm |
| 39 | CECAM Benchmark Suite | Publication-oriented, subsumed by kdotpy comparison as first entry |
| 40 | Parameter-Set Invariance (W-variant) | Narrow, subsumed by kdotpy |
| 41 | FD Order Convergence as Regression | Duplicate of Richardson on Observables |
| 42 | Topological Perturbation Robustness | Narrow to topology module |
