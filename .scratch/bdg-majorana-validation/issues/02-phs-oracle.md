# Issue 00 — Pure per-point BdG evaluator extracted from app glue (U2)

> **File-numbering note**: This file is `02-phs-oracle.md` in `.scratch/bdg-majorana-validation/issues/`. The content is **Issue 00** (pure evaluator) — sourced from `.superpowers/sdd/issue-02-brief.md`. The brief filenames were off by 2 from the issue numbers; this corrected numbering aligns the file's content with its issue number in the dependency graph.

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U2)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-02-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md` (full Understand-phase analysis)
**Phase**: PR-A, foundational seam (Issue 00 must land before Issue 05 / 07 / 08)

## What to build

A pure-function seam for the per-point BdG gap/evaluator that today lives inline in `main_topology.f90`. The new evaluator takes BdG eigenpairs (not a matrix-build step) plus a small parameter record, and returns `(minigap, near_zero_count, invariant_flag)`. Build-and-solve stays in the caller (the app or test harness owns it). The three `sweep_model` per-point call sites in `main_topology` (`compute_qw_fukane_gap_sweep`, `compute_wire_bdg_gap_sweep`, the `bhz_analytic` dispatch) all call this one physics routine. The triplicated `2·min|E|` gap convention and the hardcoded `0.001·δ₀` near-zero literals are lifted into one parametrized, unit-testable call.

The pure-function seam is the foundational enablement for every later slice (the polarization routine, the Pfaffian wrapper, the LDOS observable all consume a uniform gap/evaluator contract).

## Acceptance criteria

- [ ] A pure `eval_bdg_point(eigenvalues, params) → (minigap, near_zero_count, invariant_flag)` exists in a physics module. Recommended location: NEW `src/physics/bdg_observables.f90` (keeps `bdg_hamiltonian.f90` untouched by Issue 00, so Issue 03 owns `bdg_hamiltonian.f90` exclusively in its PR). Alternative: colocate inside `bdg_hamiltonian.f90` — but this conflicts with Issue 03's exclusive ownership. **Use the new module.**
- [ ] `main_topology`'s three per-point call sites delegate to this evaluator; no inline `2·min|E|` or `0.001·δ₀` literal remains in `main_topology.f90`.
- [ ] The near-zero count uses a single named threshold on the params record (default `0.001·δ₀`, overridable per call).
- [ ] The evaluator is callable with a pure array argument and a small params record — no `simulation_config`, no filesystem.
- [ ] The heuristic invariant flag matches the existing `compute_z2_gap` output on the same input (behavior preserved during extraction).
- [ ] Unit test `tests/unit/test_bdg_evaluator.pf` (new) is green: synthetic ±E-symmetric spectrum returns `minigap = 2·min|E|` exactly; near-zero count matches `0.001·δ₀` behavior on the same spectrum.
- [ ] No regression in existing topology regression tests.
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

### Inline gap/threshold literals to extract

In `src/apps/main_topology.f90` (1320 lines):
- `2.0_dp * minval(abs(eigvals_bdg))` at lines **534, 703, 1313** (three sites — Issue 00's seam replaces all three).
- `0.001_dp * cfg%bdg%delta_0` near-zero threshold at lines **550, 570, 623** (three sites — also lift into params record).

### Existing `compute_z2_gap` (topological_analysis.f90:260-287) behavior to preserve

Pure 1D heuristic — counts eigenvalues inside `[-gap_threshold, +gap_threshold]`, returns 1 if `n_in_gap ≥ 2`, else 0. The extracted invariant_flag must match this on identical input.

### `bdg_zero_energy_gap` (topological_analysis.f90:891-900)

`pure function` returns `minval(abs(eigenvalues))` — related but simpler (just the smallest |E|). Issue 00's `eval_bdg_point` is the umbrella that subsumes this for BdG.

### `run_bdg_qw` minigap (line 703)

`result%min_gap = 2.0_dp * minval(abs(eigvals_bdg))` — applied to the full 16N spectrum (dominated by eV-scale valence bands). Issue 00's evaluator routes this through the seam; Issue 05's verifier exercises the fix.

## Constraints from CLAUDE.md + ADRs

- **ADR 0001 (fat derived type)**: pure functions only; no polymorphic types.
- **ADR 0003 (build-and-solve in app)**: per-point evaluator stays in physics; build-and-solve stays in `main_topology`.
- **CLAUDE.md Engineering Principles (DRY/SSOT)**: evaluator does not duplicate `compute_z2_gap` heuristic — it CONSUMES the existing function or re-implements the same logic at the new seam. Either is acceptable; the AC "behavior preserved during extraction" means same output on same input.
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line (no `&` continuation — known build failure mode).
- **CLAUDE.md Boundaries**: NO approval gate for Issue 00 (does not touch `bdg_hamiltonian.f90` constructor code; only adds a new file + modifies `main_topology.f90`'s per-point literal substitution).

## File ownership (exhaustive)

### New files (you create)
- `src/physics/bdg_observables.f90` — the evaluator module. Pure functions only.
- `tests/unit/test_bdg_evaluator.pf` — pFUnit test for the evaluator.

### Modified files
- `src/apps/main_topology.f90` — replace the three `2·min|E|` literals and three `0.001·δ₀` literals with calls to the new evaluator.
- `src/physics/AGENTS.md` — Module Inventory row for new `bdg_observables.f90`.
- `src/CMakeLists.txt` — add `bdg_observables.f90` to `COMMON_SOURCES`.
- `tests/CMakeLists.txt` — `add_pfunit_ctest(test_bdg_evaluator ...)` registration.

### NOT modified (out of scope for this issue)
- `src/physics/bdg_hamiltonian.f90` (Issue 03 owns this; preserve untouched).
- `src/core/defs.f90` (Issue 06/07 may add enum value; not Issue 00).
- Any docs (no AGENTS.md doc changes beyond module inventory line).

## Suggested API (illustrative — adjust if the implementer prefers a different shape)

```fortran
module bdg_observables
  implicit none
  private
  public :: bdg_eval_params_t, eval_bdg_point

  type :: bdg_eval_params_t
    real(dp) :: delta_0          ! SC gap magnitude (eV) — used for near-zero scale
    real(dp) :: near_zero_frac   ! default 0.001_dp; |E| < near_zero_frac*delta_0 counts as near-zero
    real(dp) :: invariant_tol    ! default 1.0e-10_dp; tolerance for invariant_flag comparison
  end type

  type :: bdg_eval_result_t
    real(dp) :: minigap          ! 2*minval(|E|)
    integer  :: near_zero_count  ! count of |E| < near_zero_frac*delta_0
    integer  :: invariant_flag   ! 1 if topological (near_zero_count >= 2), 0 otherwise
  end type

  pure function eval_bdg_point(eigenvalues, params) result(r)
    real(dp), intent(in) :: eigenvalues(:)
    type(bdg_eval_params_t), intent(in) :: params
    type(bdg_eval_result_t) :: r
    ! ... single implementation of the three per-point steps ...
  end function
end module
```

## Tests required

### `tests/unit/test_bdg_evaluator.pf`
1. ±E-symmetric spectrum → `minigap == 2*minval(|E|)` to roundoff.
2. Near-zero count matches `0.001·δ₀` behavior: 3 eigenvalues inside the band → `near_zero_count == 3`; 0 inside → 0.
3. Invariant flag: `near_zero_count == 2` → `1`; otherwise `0`.
4. Overrideable threshold: `params%near_zero_frac = 0.01` counts more eigenvalues.
5. Empty spectrum → `minigap = 0`, `near_zero_count = 0`, `invariant_flag = 0` (defensive; not a hot path).
6. Pure-function assertion: `eval_bdg_point(e, p)%minigap == eval_bdg_point(e, p)%minigap` trivially; the function must not retain state.

## TDD discipline (mandatory)

1. Write `tests/unit/test_bdg_evaluator.pf` FIRST (RED).
2. Confirm test fails because module doesn't exist (or signature mismatch).
3. Implement `src/physics/bdg_observables.f90` (GREEN).
4. Refactor inline literals in `main_topology.f90` (KEEP GREEN).
5. Run the full unit suite (`ctest --test-dir build -L unit --output-on-failure`) to confirm no regressions.
6. Run the full ctest (`ctest --test-dir build -L unit -L regression --output-on-failure`) for a final pass.

## Build commands (use these exact forms)

```bash
# Verify build (no rebuild needed if files unchanged)
cmake --build build

# Run pFUnit unit tests (with OMP cap to avoid timeouts)
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure

# Run regression tests (slower; cap OMP)
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure
```

## Out of scope (must NOT do)

- Do NOT modify `src/physics/bdg_hamiltonian.f90`. Issue 03 owns this exclusively.
- Do NOT modify `src/core/defs.f90`. Issue 06/07 owns.
- Do NOT add new TOML fields. ADR 0002 forbids it for this slice.
- Do NOT touch `compute_z2_gap` itself (only route call sites through the new seam; the function is reused as-is or replaced by an equivalent inline check inside `eval_bdg_point`).
- Do NOT change the heuristic invariant_flag logic (preserve behavior per AC).

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-02-report.md`

The report must include:
- Status (DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT)
- Files created / modified
- TDD evidence (RED → GREEN) with exact commands and output
- Test results (full unit suite pass count)
- Self-review findings
- Commit SHAs (run `git log --oneline -5` to capture)
- Any concerns

---

## Outcome (as executed)

- **Module**: `src/physics/bdg_observables.f90` (96 lines; preferred new-module location per AC to keep `bdg_hamiltonian.f90` untouched for Issue 03's exclusive ownership).
- **Three call sites lifted**: lines 534, 703, 1313 of `main_topology.f90`; lines 550, 570, 623 (`0.001·δ₀` near-zero literals) lifted to params record.
- **Heuristic invariant flag**: preserved — matches `compute_z2_gap` on identical input.
- **Unit test**: `tests/unit/test_bdg_evaluator.pf` (new) green.
- **Full unit suite**: 35/35 GREEN pre-Issue-00 → 36/36 GREEN post-Issue-00 (new test added).
- **No regression** in existing topology regression tests.
- **TDD evidence**: RED (module doesn't exist) → GREEN (after implementation + main_topology refactor).