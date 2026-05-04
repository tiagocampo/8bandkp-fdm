# Plans Review — Implementation Status

Review date: 2026-05-03

| # | Group | Files | Status | Notes |
|---|-------|-------|--------|-------|
| 1 | 2026-03-28-codebase-improvement | design | COMPLETE | Archived |
| 2 | 2026-03-29-cmake-testing | design | COMPLETE | Archived |
| 3 | 2026-03-29-gfactor-improvements | design | COMPLETE | Archived |
| 4 | 2026-03-29-gfactor-validation | standalone | INCOMPLETE | Missing regression tests for 5/6 configs |
| 5 | 2026-03-29-pr-review-fixes | design | COMPLETE | Archived (Makefile target + z(1) guard added) |
| 6 | 2026-03-29-self-consistent-sp | design | COMPLETE | Archived (bulk EF shift, delta-doping, gfactor SC wired; bulk SC uses QW path with confDir=z) |
| 7 | 2026-03-30-sc-validation-report | standalone | COMPLETE | Archived |
| 8 | 2026-04-02-quantum-wire | design | COMPLETE | Archived (piezoelectric excluded: ZB [001]=zero by symmetry; FFTW Poisson not needed, PARDISO used) |
| 9 | 2026-04-03-linalg-backend | design | COMPLETE | Archived |
| 10 | 2026-04-04-documentation-overhaul | design | COMPLETE | Archived (wire_optical_spectrum figure added; 7 renamed figures have better names) |
| 11 | 2026-04-04-pr-review-fixes | design | COMPLETE | Archived |
| 12 | 2026-04-06-bulk-chapter-overhaul | design | COMPLETE | Archived |
| 13 | 2026-04-06-documentation-deep-rewrite | standalone | COMPLETE | Archived |
| 14 | 2026-04-06-quickstart-fixes | design | COMPLETE | Archived |
| 15 | 2026-04-07-strain-unification | design + plan | COMPLETE | Archived |
| 16 | 2026-04-12-qw-documentation-overhaul | design + 2 plans | COMPLETE | Archived (anticrossing prose updated with 34 meV gap from realistic InAsW(8nm)/GaSbW(8nm) config) |
| 17 | 2026-04-17-fd-variable-coeff-fix | design + plan | COMPLETE | Archived (staggered-grid deviation justified) |
| 18 | 2026-04-18-ch01-03-review-fix | design + plan | COMPLETE | Archived |
| 19 | 2026-04-18-ch04-06-review-fix | plan | COMPLETE | Archived |
| 20 | 2026-04-18-ch04-strain-figures | design + plan | COMPLETE | Archived |
| 21 | 2026-04-19-ch03-ch04-figures | plan | COMPLETE | Archived |
| 22 | 2026-04-19-qw-tutorials | design | COMPLETE | All optics figure functions + 5 new Group #50 figures. ISBT cross-validated |
| 23 | 2026-04-20-remaining-work | design + plan | INCOMPLETE | Task 8: core-shell benchmark config missing |
| 24 | 2026-04-20-strain-sign-and-buffer-fix | standalone | COMPLETE | Archived |
| 25 | 2026-04-21-benchmark-matrix | standalone | COMPLETE | Archived (reference doc) |
| 26 | 2026-04-21-docs-physics-revamp | design + plan | INCOMPLETE | 10/12 tasks missing, design archived |
| 27 | 2026-04-21-failure-ledger | standalone | COMPLETE | Archived (reference doc) |
| 28 | 2026-04-21-figure-provenance | standalone | COMPLETE | Archived (reference doc) |
| 29 | 2026-04-21-nextnano-bulk-qw-comparison | standalone | COMPLETE | Archived (reference doc) |
| 30 | 2026-04-21-nextnano-qcse-optics-comparison | standalone | COMPLETE | Archived (reference doc) |
| 31 | 2026-04-21-nextnano-wire-comparison | standalone | COMPLETE | Archived (reference doc, 5 follow-ups pending) |
| 32 | 2026-04-21-validation-summary | standalone | COMPLETE | Archived (reference doc) |
| 33 | 2026-04-21-wire-validation-summary | standalone | COMPLETE | Archived (reference doc) |
| 34 | 2026-04-25-optical-properties | redesign + implementation | COMPLETE | Archived |
| 35 | 2026-04-25-wire-gfactor-commutator | design | COMPLETE | Archived |
| 36 | 2026-04-26-hamiltonian-performance-refactor | design + plan | COMPLETE | Archived |
| 37 | 2026-04-26-modern-fortran-migration | design + plan | MOSTLY COMPLETE | validate_simulation_config missing |
| 38 | 2026-04-27-bdg-topological-sc | design + impl plan | MOSTLY COMPLETE | QW BdG, 4 unit test files, 2 analysis funcs missing |
| 39 | 2026-04-27-pr11-review-fixes | design + plan | COMPLETE | Archived |
| 40 | 2026-04-27-wire-fastpath-followup | doc + plan | COMPLETE | Archived |
| 41 | 2026-04-28-modern-fortran-phase2 | plan | COMPLETE | Archived |
| 42 | 2026-04-28-pr11-cache-safety-fixes | plan | COMPLETE | Archived |
| 43 | 2026-04-29-modern-fortran-remaining | plan | COMPLETE | Archived |
| 44 | 2026-04-29-phase4-exploratory | standalone | COMPLETE | Archived |
| 45 | 2026-04-29-pr12-fixes | plan | COMPLETE | Archived (zdotc decision documented in CLAUDE.md) |
| 46 | 2026-05-01-bhz-remaining-work | plan | COMPLETE | Archived |
| 47 | 2026-05-02-magnetic-field-landau | design + impl plan | INCOMPLETE | Peierls not in bulk, ExternalField float fix skipped |
| 48 | 2026-05-02-physics-figures-implementation | plan | INCOMPLETE | Same blockers as group 47 |
| 49 | 2026-05-02-topological-suite-verification | spec + plan | INCOMPLETE | Fu-Kane, gap-sweep, Berry curvature figures missing |
| 50 | 2026-05-03-physics-figures-extended | plan | COMPLETE | All 5 phases done: bulk E(k), QW subbands, wavefunctions, wire geometry, Zeeman fan |

---

## Detailed Findings

(Filled in as each group is reviewed)

### 4. 2026-03-29-gfactor-validation

**Status: INCOMPLETE**

All code-level implementation is done (g-factor functions, sigma matrices, self-term exclusion, near-zero guards). Missing:

- **5/6 test configs lack reference data and automated regression tests:**
  - `gfactor_bulk_gaas_cb.cfg` — no reference data, no CMake test
  - `gfactor_bulk_gaasw_cb.cfg` — no reference data, no CMake test
  - `gfactor_bulk_inasw_cb.cfg` — no reference data, no CMake test
  - `gfactor_bulk_gaas_vb.cfg` — no reference data, no CMake test
  - `gfactor_qw_vb.cfg` — no reference data, no CMake test
- Only `gfactor_qw_cb.cfg` has full regression coverage
- Document values are stale (QW CB g-factor values differ from current reference data)

### 6. 2026-03-29-self-consistent-sp

**Status: COMPLETE**

All items resolved in PR #13 (Phase 2 physics wiring):
- **Bulk EF shift**: Uniform diagonal shift added to `ZB8bandBulk` via `cfg%Evalue` and `cfg%sc_potential_shift` (TDD test: `test_bulk_ef_shift`)
- **Delta-doping**: Extended `doping_spec` with Gaussian profile (`delta<N>: NS FWHM POS` syntax), implemented in `build_doping_charge` (TDD test: `test_delta_doping_gaussian`)
- **gfactor SC**: `self_consistent_loop` (QW) and `self_consistent_loop_wire` (wire) wired into `main_gfactor.f90`
- **Bulk SC**: Uses QW path (`confDir='z'`, single material + doping) — no separate bulk SC routine needed. Warning printed if user tries `confDir='n'` with SC enabled.
- **b_field parsing**: Bulk Landau level magnetic field support with peek/backspace pattern

### 22. 2026-04-19-qw-tutorials

**Status: COMPLETE**

Fortran code fully done. All 15 optics figure functions implemented in generate_all_figures.py. All PNG files generated and verified in docs/figures/. Verified on 2026-05-04: absorption edge positions, ISBT peak energies, gain blue-shift behavior all physically correct.

ISBT cross-validation: z-dipole vs commutator velocity consistency confirmed (see test_isbt_dipole_velocity_consistency in test_optical_qw.pf, rel_diff = 0.007%).

### 26. 2026-04-21-docs-physics-revamp

**Status: INCOMPLETE**

Design doc archived (defines scope only). Plan has 12 tasks; Tasks 1-2 (control docs, provenance) DONE. Missing Tasks 3-12:

- **Tasks 3-5**: Figure fixes (broken refs, missing optics/gain/ISBT figures, PNG generation)
- **Tasks 6-7**: Solver repair (ISBT dipole sign, gain quasi-Fermi integration)
- **Tasks 8-9**: Wire hardening (FEAST convergence, branch tracking edge cases)
- **Tasks 10-11**: Chapter rebuilds (Ch06 optics, Ch08 wire with real data)
- **Task 12**: Final verification sweep

### 23. 2026-04-20-remaining-work

**Status: INCOMPLETE**

13/14 tasks done. Missing:

- **Task 8 (Stier & Bimberg core-shell wire benchmark)**: Integration test `test_wire_core_shell.sh` exists but regression config `wire_inas_gaas_core_shell.cfg` was never created — test cannot function without it

### 37. 2026-04-26-modern-fortran-migration

**Status: MOSTLY COMPLETE**

Phases 1-2 done (forall→do concurrent, goto→named blocks, dsqrt→sqrt, external→interfaces, elemental pure, private default, do concurrent in hot paths, finalizers, type-bound methods). Missing:

- **`validate_simulation_config` subroutine** (Plan Task 10): Post-parse validation routine not implemented in input_parser.f90
- **`config_validation_result` helper type**: Not present in defs.f90
- **Known `contiguous` gaps** documented in CLAUDE.md: ZB8bandQW profile/kpterms, utils.f90 dns(:,:), spin_projection.f90 psi(:)

### 38. 2026-04-27-bdg-topological-sc

**Status: MOSTLY COMPLETE**

Core physics fully implemented: Chern number (Fukui-Hatsugai-Suzuki), Z2 (gap + Fu-Kane parity), BdG 1D wire assembly, LDOS via complex PARDISO, Majorana extraction, phase diagrams, Berry curvature, topologicalAnalysis executable, lecture 13. Design doc archived. Missing:

- **`build_bdg_hamiltonian_qw`**: Only wire CSR path exists; no dense QW BdG variant — limits BdG analysis to 2D wires
- **`compute_longitudinal_conductance`** (sigma_zz chiral anomaly): Not implemented
- **`compute_spectral_function`** (k-resolved A(k,E)): Not implemented
- **4 unit test files missing**: test_magnetic_field.pf, test_z2_invariant.pf, test_edge_states.pf, test_green_functions.pf — integration/regression tests partially cover this

### 47. 2026-05-02-magnetic-field-landau

**Status: INCOMPLETE**

Phases 1-2 done (wire Peierls verified, QW Zeeman added). Phase 4 done (Majorana phase diagram). Missing:

- **Phase 3: Peierls not integrated into ZB8bandBulk**: `add_peierls_coo` exists but is NOT called from bulk Hamiltonian. Lecture 13 confirms "PENDING: Peierls substitution not yet integrated into bulk Hamiltonian"
- **ExternalField float parsing**: `input_parser.f90:492` still reads integer; `b_field` config workaround used instead
- **Phase 5: Fu-Kane Z2 for QW**: Stub returns 0; intentionally deferred but not implemented

### 48. 2026-05-02-physics-figures-implementation

**Status: INCOMPLETE**

Tasks 1.1-1.2, 3.1-3.3, 4.1-4.3 done. Missing (same blockers as group 47):

- **Task 2.1: Peierls in ZB8bandBulk**: Not done — blocks Landau level bulk verification
- **Task 2.2: Fix ExternalField parsing**: Not done — workaround in place
- **Task 4.4: Final verification**: Landau levels still PENDING

### 49. 2026-05-02-topological-suite-verification

**Status: INCOMPLETE**

Phases 1, 3-5 largely done (QWZ Chern verification, Peierls implementation, Landau config/figures, BdG parsing, Majorana figure, lecture 13 updated). Missing:

- **Phase 2: Fu-Kane Z2 for QW**: `compute_z2_fukane` stub returns 0 — not implemented
- **Gap-sweep function**: `compute_z2_gap_sweep` does not exist; wire uses M-sign heuristic
- **Berry curvature heatmap figures**: `chern_berry_curvature_*.png` not generated
- **BHZ phase transition figures**: `bhz_z2_phase_transition.png` and `bhz_edge_localization.png` not generated

### 50. 2026-05-03-physics-figures-extended

**Status: COMPLETE**

All 5 phases delivered in Phase 4 (commit 79688ab):
- Phase 1: Bulk E(k) for GaAs, InAs, InSb (3-panel figure)
- Phase 2: QW wavefunctions with band edge profile
- Phase 3: Wire material map + radial potential cut
- Phase 4: Analytical Zeeman fan diagram (spin splitting vs B)
- Phase 5: Figures integrated into `generate_all_figures.py` (not separate script)
