---
date: 2026-05-17
topic: core-kp-validation
focus: Validating core k.p physics against published literature, independent codes, and internal consistency checks
mode: repo-grounded
---

# Ideation: Core k.p Physics Validation

## Grounding Context

**Codebase Context (updated May 2026):** Fortran 2018 k.p FD code with 4 executables. 91 tests: 31 unit (pFUnit), 44 regression, 13 verification (4 ladder rungs + 2 strain rungs), 7 standard-star benchmarks (S1-S7), 1 coverage matrix. 87 regression configs. 14 lecture scripts (L00-L14). Verification ladder (4 rungs), standard-stars (S1-S7: GaAs/InAs/InSb bulk, GaAs/AlGaAs QW, InAs/GaSb QW, InAs/GaAs strained QW, InAs wire), CSR structure + Krylov snapshot tests, coverage matrix with validation_universe.yml (49 required + 12 aspirational cells, 61 COVERAGE annotations).

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

### 6. Richardson Extrapolation Convergence Fixture — PARTIALLY IMPLEMENTED
**Status:** Unit test exists for FD stencil coefficients only. NOT extended to physics observables (Eg, subband spacing, g-factor, absorption edge).

### 7. Strain Validation Across Geometries — PARTIALLY IMPLEMENTED
**Status:** Bulk + QW rungs done. Wire has profile test but no physics observable validation (strain_shift, subband spacing under strain).

## Ranked Ideas (May 2026 — Continued Ideation)

### 1. kdotpy Cross-Code Comparison
**Description:** Build an automated comparison pipeline against kdotpy (Wuerzburg, Python-based, open-source 8-band k.p solver using plane-wave expansion). Translate parameters.f90 material data to kdotpy format, run identical bulk and QW configurations in both codes, compare eigenvalues, effective masses, and g-factors. Start with the 7 standard-star materials where expected values are already documented. JSON/CSV interchange format. The comparison does not need to cover wire, topology, or optics — just the shared physics (bulk + QW eigenvalues and masses).
**Warrant:** `external:` kdotpy (Beugeling et al., arXiv:2407.12651, SciPost Phys. Codebases 47, 2025) is an independent Python-based 8-band k.p implementation using different discretization. No shared inter-code benchmark exists for k.p codes — CECAM 2024 workshop identified this gap. `direct:` CLAUDE.md states "no cross-code comparison capability."
**Rationale:** The deepest validation gap is not "we test too few things" but "we test against references that cannot distinguish 8-band-correct from 8-band-wrong." When the 8-band CB mass deviates 29% from Vurgaftman, the test tolerance accommodates that deviation — but also accommodates a wrong implementation. An independent 8-band code provides ground truth that 2-band analytical formulas cannot. "Our code agrees with kdotpy to X ppm on these benchmarks" is publishable in a way that "our code matches the 2-band Kane formula to within 30%" is not.
**Downsides:** Requires kdotpy installation, parameter mapping, and config translation. One-time setup effort for bulk and QW.
**Confidence:** 80%
**Complexity:** High
**Status:** Unexplored

### 2. Richardson Extrapolation on Physics Observables
**Description:** Build a convergence fixture that runs the same config at FD orders 2, 4, 6, 8 and multiple grid spacings, extracts a named physics observable (subband energy, g-factor, absorption edge), fits E(h) = E_exact + c_2 h^2 + c_4 h^4 + ..., and asserts the observed convergence rate matches the theoretical order. For standard-stars, report a Grid Convergence Index (GCI) — quantified discretization uncertainty alongside every benchmark number. The extrapolated continuum value becomes a self-generated reference that any future config can be compared against.
**Warrant:** `direct:` CLAUDE.md explicitly states "Richardson extrapolation only on FD stencil coefficients, NOT on physics observables." The order-10 D1 coefficient bug produced silently wrong results caught only by operator-level convergence testing. `external:` Roache 1998 GCI methodology is the community standard for FD code verification (ASME V&V 20, AIAA guidelines). Theoretical basis: FD eigenvalue approximations satisfy E(h) = E_exact + c_2 h^2 + c_4 h^4 + ... (J. Comput. Appl. Math. 2009).
**Rationale:** The codebase supports FD orders 2-10 but no test uses more than one order at a time on physics observables. Richardson turns this latent capability into a validation tool: every observable gets a convergence certificate. If the convergence rate is wrong, the FD stencils or matrix assembly have a bug — regardless of whether the numbers look plausible. The GCI turns standard-star snapshots into certified measurements with uncertainty bounds, making them scientifically publishable.
**Downsides:** Multiplies runtime per config by 3-5x. Some observables may not converge cleanly at high orders. Some configs already take minutes at one resolution.
**Confidence:** 95%
**Complexity:** Medium
**Status:** Unexplored

### 3. Exciton Binding Energy Validation
**Description:** Add validation for the exciton module (exciton.f90), which is wired into production executables (bandStructure, opticalProperties) but has zero test coverage — no unit test, no regression test, no verification script, no standard-star benchmark. Validate against the 2D hydrogen model (E_b = Rydberg/4 for infinite barriers, with finite-barrier corrections) and published QW exciton binding energies (e.g., Miller et al., Phys. Rev. B 1985 for GaAs/AlGaAs).
**Warrant:** `direct:` Searching for exciton tests yields nothing. The module is called in main_optics.f90:548-575 but no test validates the output. A figure in the lecture docs hardcodes E_b = 10 meV while the code outputs 4.82 meV — this discrepancy is undetected by any automated test. The validation_universe.yml has no exciton observable at all.
**Rationale:** Shipped code with zero validation is the highest-risk category. The 2D hydrogen model provides an analytical reference. Even a minimal test (exciton binding energy within 50% of R/4 for a GaAs/AlGaAs QW) would be a vast improvement over the current state of no test at all.
**Downsides:** Exciton physics is sensitive to barrier height and well width. Finite-barrier corrections reduce binding energy from the infinite-barrier limit. Non-parabolicity corrections from the 8-band model may cause deviations from the 2D hydrogen model.
**Confidence:** 90%
**Complexity:** Low
**Status:** Unexplored

### 4. SC Loop Physics-Output Regression
**Description:** After self-consistent convergence, assert quantitative physics outputs: Fermi level matches charge-neutrality expectation, subband energy shift relative to flat-band matches analytical estimate (Bastard 1981 for QW subband shift under electric field), charge density profile integrates to the doping level. Currently tests only check "did the loop converge." The SC benchmark verifier checks gap within a 120x range (0.05-6.0 eV) that nearly any wrong answer falls within.
**Warrant:** `direct:` CLAUDE.md explicitly states "SC loop convergence tests are 'does it converge' checks, no physics-output regression." A bug in charge_density.f90 (wrong Fermi integration, incorrect k_parallel sampling) or poisson.f90 (wrong boundary condition application) would still produce a convergent SC loop — converging to an unphysical steady state.
**Rationale:** SC is the most complex multi-module interaction (Hamiltonian + Poisson + charge density + DIIS + k-parallel sampling). "Converges" is necessary but not sufficient — the converged state must be physically correct. A single GaAs/AlGaAs reference with frozen Fermi level, subband energies, and charge profile would catch regressions in any of the five coupled modules.
**Downsides:** SC results are grid- and k-point-dependent. Reference values must be established at a specific configuration. Different doping configurations (uniform vs modulation-doped) may need different physics invariants.
**Confidence:** 85%
**Complexity:** Medium
**Status:** Unexplored

### 5. Verification Ladder Extension to Optical and g-Factor
**Description:** Extend the existing 4-rung verification ladder (band structure observables) to optical and g-factor observables. Rung 1: bulk Kane interband matrix element P vs Vurgaftman/Winkler tabulated values. Rung 2: QW intersubband transition energies vs effective-mass approximation (known to be approximate, tests the deviation direction). Rung 3: bulk g-factor vs Roth formula (already analytically verified in lecture scripts). Rung 4: QW g-factor vs published experimental data for InGaAs/GaAs. The ladder structure provides diagnostic precision: if rung 1 passes but rung 3 fails, the bug is in QW-specific g-factor calculation.
**Warrant:** `direct:` Standard-star suite is missing g-factor and optical spectra benchmarks. The optical module and g-factor module are tested by lecture scripts but lack the systematic rung structure that makes the band-structure ladder so diagnostically valuable. `external:` Roth formula is exact for 2-band g-factor; Winkler 2003 Table 6.2 lists expected accuracy for 8-band g-factors.
**Rationale:** Optical and g-factor calculations depend on the full k.p machinery — Hamiltonian construction, eigensolution, velocity matrices, and post-processing. They are integration tests for the entire pipeline. The verification ladder structure provides diagnostic precision and fills cells in the coverage matrix automatically.
**Downsides:** 8-band g-factor has known 20-30% shortfall vs experiment (model limitation). Optical spectra benchmarks require linewidth assumptions.
**Confidence:** 85%
**Complexity:** Medium
**Status:** Unexplored

### 6. Wire Strain Physics Observable Validation
**Description:** Replace the 60% tolerance profile-shape test with quantitative band-edge shift checks. In the wide-core limit (wire diameter >> lattice constant), the strain field should converge to bulk biaxial strain values — analytically computable from Vurgaftman deformation potentials (a_v, a_c, b, d) via Chuang strained-band formulas. Add tests: (a) strained wire band gap shift vs bulk biaxial prediction in wide-core limit; (b) HH-LH splitting under strain matches analytical Bir-Pikus prediction; (c) subband spacing under strain is physically consistent (compression increases confinement, tension decreases it).
**Warrant:** `direct:` CLAUDE.md states "Wire strain: profile test exists but no physics observable validation." The verify_strain_wire_profile.py tolerance of 60% on strain magnitude is enormous — even a qualitatively wrong PDE solution could pass. `external:` Vurgaftman 2001 Tables XIV-XV provide deformation potentials; Chuang Chapter 4 derives strained band edges analytically.
**Rationale:** Wire strain is the gateway to wire g-factor, wire optics, and wire topological analysis. If the strain field is wrong, every downstream observable built on top of it inherits the error. The wide-core analytic limit provides a rigorous anchor that the current 60% tolerance cannot deliver.
**Downsides:** Wide-core limit requires large wire configs (expensive). Strain field non-uniformity at finite diameter means the limit is approximate.
**Confidence:** 80%
**Complexity:** Medium
**Status:** Unexplored

### 7. TRK Sum Rule / Optical Invariant Gate
**Description:** Implement the Thomas-Reiche-Kuhn f-sum rule (sum of oscillator strengths = 1 for a complete eigenbasis) as a parameter-free, material-independent check on the commutator-based velocity matrices. Extend to Kramers degeneracy (E(+k) = E(-k) without B-field) and time-reversal symmetry. These are algebraic identities that hold to machine precision for a correct Hamiltonian and eigensolver.
**Warrant:** `direct:` Velocity matrices v_alpha = -i[r_alpha, H] are tested only indirectly through optical absorption peaks. The commutator construction on CSR is the most subtle algebraic operation in the codebase. `reasoned:` TRK is an algebraic identity of the complete eigenbasis, independent of material parameters, geometry, or grid spacing. No external reference needed.
**Rationale:** Currently every validation compares against published numbers. Sum rules and symmetry checks are zero-parameter oracles that catch entire classes of bugs (wrong velocity matrix, missing band contributions, eigenstate normalization errors) that literature comparisons miss. A code change that silently corrupts the velocity operator corrupts g-factor, optical absorption, and ISBT simultaneously — the TRK check catches this at the source.
**Downsides:** Finite-basis truncation means the sum rule holds approximately, not exactly. Tolerance setting requires physics judgment (how many bands must be included for the sum to converge).
**Confidence:** 90%
**Complexity:** Medium
**Status:** Unexplored

### 8. Hermiticity Smoke Test for All Geometries
**Description:** Assert H = H† to machine precision for bulk, QW (dense), wire (CSR), and BdG (16N x 16N Nambu) Hamiltonians at representative k-points and multiple materials. Eigenvalue reality check: all eigenvalues purely real. Material-independent, zero-tolerance pass/fail. Clear diagnostics: "block (3,5) differs from conjugate of (5,3) by 2.3e-4" rather than cryptic eigenvalue mismatch.
**Warrant:** `direct:` Wire CSR assembly has no Hermiticity test. Current bulk Hermiticity test only checks GaAs at one arbitrary k-point. `reasoned:` Hermiticity is an algebraic identity checkable to machine precision. Any non-Hermitian term (wrong sign, missing conjugate, incorrect k.p coupling) is caught immediately.
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
