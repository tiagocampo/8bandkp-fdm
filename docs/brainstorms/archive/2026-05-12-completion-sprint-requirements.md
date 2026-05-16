---
date: 2026-05-12
topic: completion-sprint
focus: Complete remaining validation and infrastructure gaps in 3 phases on the current feature branch
mode: repo-grounded
---

# Requirements: Completion Sprint (Phases 14-16)

## Context

Phases 1-13 delivered: full 4-executable feature set, topological suite, Landau mode, CSR structure testing, 8-band verification ladder (4 rungs), standard-star benchmarks (S1-S7), 14 executable lecture-test pairs, validation coverage matrix, and strain validation across geometries.

Nine items remain from the original review backlog and ideation-driven work. This doc scopes them into 3 phases ordered by effort (quick wins first).

## Phasing Strategy

Stay on `feature/bdg-topological-superconductivity`. Quick wins first, then validation tightening, then integration + Rashba physics.

---

## Phase 14: Quick Wins Sprint

Low-effort items that can be completed in a single focused session each. No physics investigation needed — mechanical changes, data generation, and doc commits.

### P14.1: Contiguous attribute gaps (Backlog #37)

Add `contiguous` attribute to the two remaining hot-path array arguments:
- `dns(:,:)` in `src/core/utils.f90:19`
- `psi(:)` in `src/physics/spin_projection.f90:14`

### P14.2: Coverage matrix orphan annotations

Fix 3 observable name mismatches between test annotations and `validation_universe.yml`:
- `test_gfactor_no_optics.sh`: `observable=gfactor` should be `observable=g*_cb`
- `test_qw_bandstructure.sh`: `observable=E_sub` should be `observable=subband_spacing`
- `test_qw_character_and_output_dir.sh`: `observable=state_character` should match an existing universe cell (likely `CB_ground_state`)

### P14.3: Strain validation docs commit

Commit the two untracked planning documents:
- `docs/brainstorms/2026-05-11-strain-validation-requirements.md`
- `docs/plans/2026-05-11-strain-validation-plan.md`

Tolerances stay relaxed (justified by physics: Bastard formula invalid for narrow wells, wire free-surface relaxation).

### P14.4: gfactor regression golden data (Backlog #4)

Generate golden reference data for 3 gfactor configs that currently lack automated regression tests:
- `gfactor_bulk_gaas_vb.cfg` — VB g-factor for GaAs bulk
- `gfactor_bulk_gaasw_cb.cfg` — CB g-factor for GaAsW bulk
- `gfactor_qw_vb.cfg` — VB g-factor for GaAs/AlGaAs QW

Each needs: (1) run executable to generate output, (2) commit reference data in `tests/regression/data/<test_name>/`, (3) register in `tests/CMakeLists.txt` with `regression` label.

---

## Phase 15: Validation Tightening

Low-to-medium effort items completing the validation infrastructure.

### P15.1: CSR Krylov snapshot completion (Phase 8 gap)

Add Krylov snapshot tests for 3 remaining code paths:
- **SC loop** — snapshot the initial Hamiltonian before any SC iteration
- **Optics wire mode** — CSR assembly path used by `opticalProperties` in wire mode
- **G-factor wire mode** — CSR assembly path used by `gfactorCalculation` in wire mode

Each needs: (1) add canonical small-grid config to `generate_krylov_reference.f90`, (2) generate reference data, (3) add test case to `test_krylov_snapshots.pf`.

### P15.2: Standard-star assertion tightening (Phase 10 gap)

- **S4 absorption onset** (`verify_star_gaas_algaas_qw.py`): Run `opticalProperties` to establish `ONSET_REF`, replace range check with `compare_value`
- **S7 wire g-factor** (`verify_star_inas_wire.py`): Establish `G_WIRE_REF` from a reference run, replace 50% range with `compare_value` at 10% tolerance
- **Wrapper centralization** (S3/S4/S5/S6): Replace local `run_bandstructure`/`run_gfactor`/`run_optical` wrappers with `star_helpers.run_exe`
- **benchmarks.md Roth formula** (KD2): Fix incorrect formula at line 206

### P15.3: Re-scope docs physics revamp (Backlog #26)

Audit the 12 original tasks from `docs/plans/2026-04-21-docs-physics-revamp-plan.md` against what lecture-test pairs already delivered. Produce a short list of genuinely remaining gaps (likely ISBT dipole sign, gain quasi-Fermi integration). Either close the item or update the plan with the reduced scope.

---

## Phase 16: Integration + Rashba Physics

Medium-effort items requiring physics investigation.

### P16.1: Integration tests (Backlog #8)

- **Wire hexagon integration test**: End-to-end test running a hexagonal wire config through the solver, verifying geometry + eigenvalues. Config exists (`tests/regression/configs/`); needs shell wrapper and golden data.
- **SC wire integration test**: Self-consistent Poisson-Schrodinger loop on a wire geometry. Needs a dedicated config and golden data.

### P16.2: Rashba BdG physics calibration (Backlog — Rashba)

Fix the 4 known issues in the Rashba BdG sweep:

1. **`mu` must match a subband energy.** Run `bandStructure` on the wire config first to find subband positions. Set `mu` to the first CB subband energy.
2. **FEAST window calibration.** Size the window to the actual subband spacing around the chosen `mu`, not a fixed default.
3. **Grid spacing.** Use `wire_dx ~ wire_width / (wire_nx - 1)` for realistic confinement energies.
4. **False positive regression test.** `regression_topology_rashba_phase` should verify that FEAST actually finds eigenvalues when BdG is enabled, not just that gap < threshold.

Deliverable: a working `sweep_rashba_bdg.py` (or equivalent within `lecture_13_topological.py`) that demonstrates the Majorana phase transition.

---

## Success Criteria

- All 9 backlog items resolved or explicitly re-scoped
- Full test suite passes (`ctest --test-dir build -j4`)
- `ctest --test-dir build -L coverage` passes with zero orphan annotations
- REVIEW.md updated: all items marked COMPLETE or closed with documented rationale
- BACKLOG.md "Remaining Backlog" section empty or converted to a "Known Limitations" section

## Dependencies and Assumptions

- No new Fortran source changes needed for Phase 14 (only test infrastructure and docs)
- Phase 16.2 requires understanding wire subband physics — may need 2-3 iterations to calibrate
- The `sweep_rashba_bdg.py` script may have been removed during lecture refactoring; if so, the Rashba work happens within `lecture_13_topological.py`
- Items #5 (integration tests) and #7 (Rashba) could swap order if Rashba proves harder than expected
- Branch is `feature/bdg-topological-superconductivity` — these phases complete its story before merge to main
- **Performance improvements are in-scope**: when touching Fortran source for any reason, apply opportunistic optimizations (loop restructuring, cache-friendly access patterns, unnecessary allocations). Do not open separate performance work items, but fix what you see while in the code.
