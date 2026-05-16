---
date: 2026-05-08
topic: core-kp-validation
focus: Validating core k.p physics against published literature through tests, regression checks, and lecture documentation
mode: repo-grounded
---

# Ideation: Core k.p Physics Validation

## Grounding Context

**Codebase Context:** Fortran 2018 k.p FD code with 4 executables (bandStructure, gfactorCalculation, opticalProperties, topologicalAnalysis). 20+ regression test configs in `tests/regression/configs/` with golden-output comparison via `compare_output.py`. `docs/reference/benchmarks.md` documents 6 validation regimes but lacks systematic publishable benchmarks. 15 unit + 24 regression tests pass; pFUnit-based unit testing. Known gaps: sparse matrix assembly tested only via eigenvalues (not CSR structure); lecture docs don't cover all physics modules. Known 8-band model limitation: 20-30% g-factor shortfall vs experiment (model, not bug).

**Past Learnings:** G-factor Roth formula gives exact analytical match (GaAs -0.315, InAsW -14.858). SC-SP had hidden unit-mismatch bug (Angstroms vs nm) only visible under strongest feedback mode. Topological modules had 5 index bugs — need structural CSR tests. gfortran -O3 can corrupt derived-type allocatable components via temporaries.

**External Context:** Nextnano (Birner thesis) validates via Fock-Darwin, HgTe/CdTe crossover (~6.5 nm), InAs/GaSb vs TB, SiGe QCL vs experiment. Sengupta et al. (arXiv:1409.4376) provides explicit 8-band k.p FD coefficient matrices. Canonical benchmarks: bulk effective masses (GaAs 0.067 m0), band gaps (InAs 0.354 eV), spin-orbit splitting (GaAs 0.34 eV), BHZ edge states. Spurious solution detection via ellipticity criterion. Richardson extrapolation for convergence testing. No open-source competitor covers bulk+QW+wire+optics+topology.

## Ranked Ideas

### 1. 8-Band Verification Ladder
**Description:** Build a hierarchy of 8-band benchmarks from simplest to most complex, each with published/analytical reference: (1) 8-band bulk at k=0 — eigenvalues match band gap + spin-orbit splitting from `parameters.f90` by construction; (2) 8-band bulk dispersion — effective masses (m*_e, m*_hh, m*_lh, m*_so) near Gamma compared against Vurgaftman 2001 tables for all 25+ materials; (3) 8-band QW — subband energies vs Bastard analytical results and nextnano cross-validation; (4) 8-band wire — 2D confinement eigenvalues vs published wire results. Each rung validates the full 8-band Hamiltonian at increasing complexity.
**Warrant:** `direct:` benchmarks.md documents 6 validation regimes but lacks systematic progression from bulk k=0 through dispersion through confinement. `external:` Vurgaftman 2001 provides tabulated effective masses and band gaps for all III-V materials; Bastard 1981 gives analytical QW subband spacings.
**Rationale:** Validates the 8-band model — the core of the codebase — at every level. When a regression fails, the hierarchy tells you exactly where: bulk at k=0 (parameters), bulk dispersion (Hamiltonian construction), QW (confinement + FD), or wire (2D + CSR assembly).
**Downsides:** Requires parabolic fitting for effective masses from computed dispersion. Wire benchmark literature is thin.
**Confidence:** 90%
**Complexity:** Medium
**Status:** Unexplored

### 2. CSR Structure Testing
**Description:** Add unit tests that validate sparse matrix assembly independently of eigenvalue comparison. Two complementary approaches: (a) golden CSR structure files (row pointers, column indices, nnz count, spot-check element values) for canonical configs, and (b) Krylov-action snapshots (SpMV chain at fixed seed vector) for exponential sensitivity to structural differences.
**Warrant:** `direct:` "Known gap: sparse matrix assembly tested only via eigenvalues, not CSR structure" — explicitly documented gap. The topological/magnetic bug report confirms CSR indexing bugs passed eigenvalue tests.
**Rationale:** Eigenvalue identity is necessary but insufficient — two different matrices can share the same spectrum. Structural tests break the circularity and protect the wire and topology executables.
**Downsides:** Golden CSR files must be regenerated when Hamiltonian construction intentionally changes.
**Confidence:** 85%
**Complexity:** Low-Medium
**Status:** Unexplored

### 3. Standard-Star Benchmark Systems
**Description:** Define 3-4 canonical material systems (GaAs bulk, InAs/GaSb QW, HgTe/CdTe topological, InSb extreme spin-orbit) and validate each against multiple independent published observables — not just band structure but also effective mass, g-factor, and optical absorption edge. A GaAs "standard system passes" means band gap matches Vurgaftman, effective mass matches Kane prediction, AND g-factor matches Roth formula — all with the same parameter set.
**Warrant:** `direct:` benchmarks.md lacks systematic published benchmarks. `external:` Nextnano validates InAs/GaSb against tight-binding AND HgTe/CdTe crossover against BHZ — multiple properties from same material (Birner thesis, TUM 2011).
**Rationale:** Turns "does the code match itself?" into "does it reproduce Vurgaftman Table X for GaAs?" Multi-property validation per material makes the g-factor shortfall quantitatively traceable.
**Downsides:** Requires literature lookup for multiple observables per material. Some properties may lack clean published values.
**Confidence:** 90%
**Complexity:** Medium
**Status:** Unexplored

### 4. Executable Lecture-Test Pairs
**Description:** Refactor each lecture doc into a triple: (1) derivation of expected physics with equations, (2) committed regression config that exercises exactly that physics, (3) script that runs the config and plots output alongside the analytical derivation. Running `make lecture-roth` reproduces the derivation with code output.
**Warrant:** `direct:` Lecture docs exist in docs/lecture/ but don't cover all physics modules — the BdG/topological module has a design spec but no lecture. `reasoned:` Connecting derivation to code doubles as validation: if output diverges from derivation, either the code or the derivation is wrong.
**Rationale:** Every new lecture creates a new validation point at marginal cost. The topological module becomes accessible to anyone who didn't write it. No open-source k.p competitor offers executable lectures.
**Downsides:** Significant documentation effort. Plot scripts need maintenance as output formats change.
**Confidence:** 80%
**Complexity:** Medium
**Status:** Unexplored

### 5. Validation Coverage Matrix
**Description:** Build a heat map with axes: observable type (eigenvalues, g-factor, optical peak, Chern number, Z2, LDOS) x geometry (bulk, QW, wire) x material family (III-As, III-Sb, II-VI, nitride). Each cell is green (published/analytical benchmark), yellow (regression test, no published reference), or red (no test). Drives test-creation priorities.
**Warrant:** `direct:` 24 regression tests exist but their cross-product coverage is undocumented. `reasoned:` Coverage matrix is the standard first step in systematic testing — makes gaps visible before they become bugs.
**Rationale:** The highest-leverage meta-investment. Without it, test creation is ad hoc. With it, you identify that (e.g.) optical-properties-in-wire-mode has no published benchmark and prioritize filling it.
**Downsides:** Requires upfront audit of all 20+ configs. Must be maintained as tests are added.
**Confidence:** 85%
**Complexity:** Low
**Status:** Unexplored

### 6. Richardson Extrapolation Convergence Fixture
**Description:** Build a reusable test fixture (pFUnit or Python wrapper) that runs any regression config at multiple FD orders (2, 4, 6, 8) and grid spacings, then asserts the observed convergence rate matches the theoretical order. Outputs a convergence certificate per observable. Once built, applies to every existing and future config with zero additional test code.
**Warrant:** `external:` Richardson extrapolation is the standard method for verifying FD codes (Roache 1998). `direct:` FD stencils support orders 2-10 and CLAUDE.md flags stencil coefficients as "require approval" — critical infrastructure lacking isolated convergence testing.
**Rationale:** The FD stencils are the mathematical foundation. If a coefficient is wrong, every result is silently corrupted. Wrong coefficients produce wrong convergence rates, detectable even when answers look close. The fixture compounds because new modules inherit it for free.
**Downsides:** Running configs at multiple FD orders is computationally expensive. Some configs may not converge cleanly at all orders.
**Confidence:** 85%
**Complexity:** Medium
**Status:** Unexplored

### 7. Strain Validation Across Geometries
**Description:** Dedicated strain validation suite covering all three geometries against published results: (a) bulk — hydrostatic + biaxial strain shifts (Eg shift, HH-LH splitting) vs Vurgaftman deformation potentials (a_v, a_c, b, d) and Chuang strained-band formulas; (b) QW — strained InGaAs/GaAs subband energies and HH-LH crossover vs Bastard/Winkler analytical results and nextnano comparisons; (c) wire — strain in nanowire cross-section vs published FEM or analytical strain profiles, verifying Bir-Pikus Hamiltonian modification produces correct band-edge shifts in 2D confinement.
**Warrant:** `direct:` CLAUDE.md flags Bir-Pikus sign convention as "NEVER change" — the `P_eps = -av * Tr(eps)` sign flip is a known source of subtle errors. The gfortran -O3 segfault in QW strain path confirms the strain code path is fragile and undertested. `external:` Vurgaftman 2001 Tables XIV-XV provide deformation potentials; Chuang Chapter 4 derives strained band edges analytically; Winkler 2003 Section 4.4 covers Bir-Pikus in 8-band context.
**Rationale:** Strain is one of the most error-prone parts of the k.p implementation (sign conventions, unit consistency, geometry-dependent tensor components). It modifies the Hamiltonian directly — if wrong, everything downstream is silently corrupted. Validating across all three geometries ensures Bir-Pikus terms are correct regardless of confinement dimensionality.
**Downsides:** Wire strain comparisons are scarce in literature — may need synthetic reference data from independent FEM solver.
**Confidence:** 85%
**Complexity:** Medium
**Status:** Unexplored

## Rejection Summary

| # | Idea | Reason Rejected |
|---|------|-----------------|
| 1 | Per-material g-factor gap dashboard | Subsumed by Standard-Star Benchmarks (#3) |
| 2 | Negative-validation / model limitation tests | Subsumed by Standard-Stars and 8-Band Ladder |
| 3 | Spurious mode sentinel (ellipticity) | Important but more of a code feature than validation strategy; better as brainstorm topic |
| 4 | Parameter provenance / Vurgaftman audit | Lower priority than physics validation; useful follow-up |
| 5 | Sengupta matrix diff (block-level comparison) | Subsumed by 8-Band Ladder and CSR Structure Testing |
| 6 | Topological test oracle (integrality checks) | Strong but narrow; better as brainstorm topic for topological track |
| 7 | Cross-executable consistency tests | Largely assured if 8-Band Ladder passes (shared Hamiltonian construction) |
| 8 | Validation-as-publication (citable tests) | Subsumed by Executable Lecture-Test Pairs |
| 9 | Forced-degradation testing | Clever but high implementation burden for uncertain value |
| 10 | Fast physics gate (30-second smoke tests) | Off-topic: DX/performance, not physics validation |
| 11 | One-command test scaffolding | Off-topic: DX tooling |
| 12 | Dimensional consistency checks | Off-topic: general code quality |
| 13 | Cross-code validation pipeline (Nextnano) | Too expensive; requires Nextnano installation and config translation |
| 14 | k.p phantoms (independent solver) | Too expensive; analytical benchmarks cover similar ground |
| 15 | Temporal regression (git history) | Too expensive; high infrastructure cost |
| 16 | Published-plot-to-tolerance pipeline | Too expensive; PDF curve extraction is a separate project |
| 17 | Inverse parameter recovery | Subsumed by Standard-Stars effective mass validation |
| 18 | Proof-load benchmarks (extreme FD orders) | Subsumed by Richardson Extrapolation Fixture |
| 19 | Delta reports with drift attribution | Nice-to-have infrastructure, not essential |
| 20 | Physics-labeled regression assertions | DX improvement for debugging, not physics validation |
| 21 | Sweep-to-config auto-benchmark factory | Too ambitious / vague |
| 22 | Physics-module coverage sentinel | Subsumed by Validation Coverage Matrix |
| 23 | Analytical shadow suite (Roth + free-particle + infinite-barrier) | Subsumed by 8-Band Ladder rung 1 and rung 2 |
