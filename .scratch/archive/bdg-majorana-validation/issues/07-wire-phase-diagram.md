# Issue 07 — Wire minigap colormap + slim projected Pfaffian witness (U10)

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U10)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-07-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md`
**Phase**: PR-C, after PR-B (Issue 05 done)

## What to build

The wire-rung topological witness in three pieces:

1. **1D minigap curve (AE3 carry-over from PR40).** `min_gap(B)` at fixed `μ`, the open→close→reopen curve. Already produced by PR40; this slice just asserts the regression test `regression_wire_bdg_topological` is green at `μ=0.6601 eV` + transverse B (B_crit≈2.8 T) and that the removed Gershgorin auto-window fallback is no longer used.

2. **2D minigap colormap (A-2D extension).** At each `(B, μ)` grid point, compute `min_gap` (via Issue 00's evaluator) and shade the "topological region" where `min_gap < 5·δ₀` (the near-zero regime). The colormap replaces the heuristic `compute_z2_gap` output that previously labeled the all-zero PASS. Honest empirical colormap, not a heuristic invariant — the invariant is "did the gap close" (mechanical), the shaded region is its visual expression.

3. **Slim projected Pfaffian witness.** At one `(B, μ) = (B_crit, μ=0.6601 eV)` point, evaluate the projected Pfaffian sign with two subspace choices and assert agreement:
   - **S1 (empirical):** diagonalize the wire BdG matrix once at a reference `B` (e.g., `B = 2·B_crit`); pick the two lowest single-particle states at `kz=0`; project onto them. 4-dim Nambu.
   - **S2 (analytical):** project onto bands 7-8 (the conduction-band edge defined by the k·p block table SSOT); no diagonalization needed. 4-dim Nambu.
   - Witness: Pfaffian signs from S1 and S2 agree within sign tolerance. Disagreement is documented as a strong-SOC regime flag — a finding, not a bug.

The full wire Pfaffian sweep (U13's periodic/Bloch path) is deferred; this slice's Pfaffian evaluates at one point, not a grid.

The phase-diagram code lives in `src/apps/main_topology.f90` (the per-point evaluators inside `compute_wire_bdg_gap_sweep` and `compute_qw_fukane_gap_sweep`); the wire's slim Pfaffian wrapper calls into `src/physics/topological_analysis.f90`. The `detect_z2_transitions` / `is_z2_transition` helpers stay — the `z2_map`/`gap_map`/`transitions` contract is shape-preserving.

## Acceptance criteria

- [ ] `regression_wire_bdg_topological` asserts the open→close→reopen curve (AE3) at `μ=0.6601 eV` + transverse B.
- [ ] `regression_wire_bdg_topological` asserts the removed Gershgorin auto-window fallback is no longer used on the BdG path (the μ-in-gap run warns and returns the `min_gap=-1` sentinel).
- [ ] `validation_rejects_bad_topology`: `validate_semantic` rejects explicitly-set Gershgorin-scale (>1 eV) BdG solver windows.
- [ ] The 2D minigap colormap at fixed μ, B swept is non-flat (varies across `(B, μ)`) — guards the flat-diagram regression.
- [ ] `wire_bdg` sweep path uses the Pfaffian (not the heuristic) at every grid point for the gap_map / z2_map columns.
- [ ] Slim Pfaffian witness at `(B_crit, μ=0.6601 eV)`: S1 and S2 signs agree within tolerance.
- [ ] S1 disagreement with S2 (if any) is documented as a strong-SOC regime flag in the issue report — not a defect.
- [ ] Wire spectrum matches the dense-QW control (Issue 05) at comparable parameters (the bisection-control property).
- [ ] The previous "1.22 T (auto Gershgorin fallback)" PASS line is removed from `docs/lecture/figures/rashba_majorana_phase_diagram.txt`.
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

### PR40 work
- `regression_wire_bdg_topological` already exists (PR40 introduced it).
- Wire open→close→reopen curve at B_crit ≈ 2.8 T, μ=0.6601 eV + transverse B.

### `src/apps/main_topology.f90::compute_wire_bdg_gap_sweep` + `compute_qw_fukane_gap_sweep`
- Per-point evaluators inside `main_topology`; build-and-solve per grid point.
- Use `eval_bdg_point` from Issue 00's evaluator.

### `src/physics/topological_analysis.f90::detect_z2_transitions` / `is_z2_transition`
- Shape-preserving contract: `z2_map`, `gap_map`, `transitions`.

### `tests/integration/validation_universe.yml`
- Existing `minigap` cell (line 478) covers wire/InAs/required/Lutchyn_Oreg2010.
- Need to confirm coverage extends to the 2D colormap.

### `src/core/defs.f90::validate_semantic`
- Line 855. Topology-mode validation at lines 915–928. BdG mode validation at 940–961.
- Need to add: reject explicitly-set Gershgorin-scale (>1 eV) BdG solver windows.

### `docs/lecture/figures/rashba_majorana_phase_diagram.txt`
- All zeros data (the legacy artifact that the regression test must fail to reproduce post-fix).

## Constraints from CLAUDE.md + ADRs

- **ADR 0001 (fat derived type)**: no new types.
- **ADR 0003 (build-and-solve in app)**: build-and-solve in `main_topology`; slim Pfaffian in `topological_analysis.f90`.
- **ADR 0007 (hole-block canonical)**: the slim Pfaffian evaluates on the correct BdG matrix.
- **CLAUDE.md Boundaries**: NO approval gate (no Hamilton-construction code changes; the canonical hole block is already in place from Issue 03).
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.

## File ownership (exhaustive)

### New files (you create)
- `tests/regression/configs/regression_wire_bdg_topological_2d.toml` — new config for 2D colormap (B sweep, μ fixed at 0.6601 eV; OR B fixed at 2·B_crit, μ swept).
- `tests/regression/data/regression_wire_bdg_topological_2d_*.dat` — golden reference (regenerate from successful run).
- `tests/regression/test_wire_bdg_topological_2d.sh` — shell driver for the 2D colormap.
- `tests/integration/test_validation_rejects_bad_topology.sh` — assert `validate_semantic` rejects Gershgorin-scale BdG solver windows.
- `tests/integration/test_slim_wire_pfaffian_witness.py` — slim Pfaffian witness at (B_crit, μ=0.6601 eV).

### Modified files
- `src/apps/main_topology.f90` — `wire_bdg` sweep path uses Pfaffian (not heuristic) at every grid point for `gap_map`/`z2_map` columns; add the slim Pfaffian witness call at one `(B, μ)` point (S1 + S2 strategies).
- `src/physics/topological_analysis.f90` — add slim projected Pfaffian wrapper (S1 empirical + S2 analytical bands 7-8).
- `src/core/defs.f90` — `validate_semantic` rejects explicitly-set Gershgorin-scale (>1 eV) BdG solver windows.
- `docs/lecture/figures/rashba_majorana_phase_diagram.txt` — regenerate from new simulation (or remove if superseded).
- `tests/integration/validation_universe.yml` — confirm 2D colormap covered (existing `minigap` cell may suffice).

### NOT modified (out of scope)
- `src/physics/bdg_hamiltonian.f90` (Issue 03 owns).
- `src/physics/bdg_observables.f90` (Issue 00 owns).
- `src/math/pfaffian.f90` (Issue 01 owns).
- `scripts/lecture_13_topological.py` (Issue 08 owns).
- `docs/lecture/13-topological-superconductivity.md` (Issue 08 owns; Issue 07 only touches the phase_diagram.txt figure).

## Tests required

### `tests/regression/test_wire_bdg_topological_2d.sh`
- Run `wire_bdg` sweep in 2D (B swept, μ fixed at 0.6601 eV; OR B fixed, μ swept).
- Assert the colormap is non-flat (varies across the grid).
- Compare against golden reference.

### `tests/integration/test_validation_rejects_bad_topology.sh`
- Construct TOML with BdG solver window = 5.0 eV (Gershgorin-scale).
- Run `topologicalAnalysis`.
- Assert non-zero exit code (rejection).

### `tests/integration/test_slim_wire_pfaffian_witness.py`
- Build wire BdG at (B_crit, μ=0.6601 eV).
- S1: diagonalize BdG at reference B=2·B_crit, project onto 2 lowest single-particle states at kz=0. Compute Pf of 4-dim projected ω·H.
- S2: project onto bands 7-8 (conduction-band edge). Compute Pf of 4-dim projected ω·H.
- Assert sign(s1) == sign(s2).
- If sign disagrees: document as strong-SOC regime flag (don't fail).

## TDD discipline (mandatory)

1. Verify `regression_wire_bdg_topological` is currently GREEN (PR40 work).
2. Write the slim Pfaffian witness test FIRST (RED).
3. Confirm test fails (function doesn't exist).
4. Implement slim Pfaffian in `topological_analysis.f90`.
5. Confirm test passes.
6. Implement the 2D colormap test.
7. Implement `validate_semantic` rejection test.
8. Wire `wire_bdg` sweep path to use Pfaffian (not heuristic).
9. Regenerate the colormap golden reference from a successful run.
10. Confirm all new tests pass.

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
```

## Out of scope (must NOT do)

- Do NOT modify `src/physics/bdg_hamiltonian.f90`.
- Do NOT modify `src/physics/bdg_observables.f90`.
- Do NOT modify `src/math/pfaffian.f90`.
- Do NOT modify `src/apps/main_topology.f90`'s basic build-and-solve (only the sweep-path Pfaffian routing).
- Do NOT add new TOML fields.
- Do NOT add `sweep_model` enum value.
- The full wire Pfaffian sweep (U13) is deferred — do NOT implement periodic/Bloch BdG.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-07-report.md`

The report must include:
- Status (DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT)
- Files modified
- The detected B_crit value (numerically, should match Issue 05's 2.0 T ± tolerance — note: wire rung may differ from QW rung)
- The S1 vs S2 sign agreement result
- The 2D colormap non-flat result (numerically)
- The validate_semantic rejection test result
- Test results (full unit + regression suite pass count; pre-existing failures expected)
- Commit SHA
- Self-review findings
- Concerns

## Risk notes

- The wire B_crit may differ from the QW B_crit (different geometries, different effective couplings). Document both.
- The S1 vs S2 sign disagreement is allowed as a finding. Don't fail the test; document it.
- The 2D colormap golden reference must be regenerated from a SUCCESSFUL run (otherwise the regression test fails on a stale reference). Be careful about the regeneration timing.
- The `validate_semantic` rejection test requires a TOML that explicitly sets a Gershgorin-scale BdG solver window. Verify the TOML parser accepts such an explicit value (it should, by design — only the validate_semantic should reject it).

---

## Outcome (as executed)

- **2D minigap colormap**: non-flat across (B, μ) grid (test green). Generated plot replaces the all-zero `rashba_majorana_phase_diagram.txt` legacy artifact.
- **`validate_semantic` rejection**: Gershgorin-scale (>1 eV) BdG solver windows rejected (test green). Defense-in-depth guard per Issue 07 AC.
- **Slim projected Pfaffian witness**: S1 (empirical) and S2 (analytical bands 7-8) Pfaffian signs agree at (B_crit, μ=0.6601 eV).
- **`wire_bdg` sweep path**: uses Pfaffian (not heuristic) at every grid point for gap_map/z2_map columns.
- **Slim Pfaffian wrapper**: ~140 lines in `topological_analysis.f90` (slightly over 100-line target; acceptable).
- **Regression test**: `regression_wire_bdg_topological` green; open→close→reopen curve at μ=0.6601 eV + transverse B verified; Gershgorin auto-window fallback confirmed not used.
- **Plots**: 2D colormap regenerated and copied to `docs/lecture/figures/`.
- **Plot generation**: 2D colormap plot generated by `verify_wire_bdg_topological.py` and copied to `docs/lecture/figures/`.
- **Deferred**:
  - S1 (full diagonalization) sweep: needs U13 periodic/Bloch BdG.
  - Slim Pfaffian wrapper ~140 lines: slightly over 100-line target; acceptable.