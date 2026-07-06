**Status**: COMPLETE (2026-07-05)

# PR #27 Review Fixes — Batch 2 PRD

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all findings from the comprehensive 8-agent workflow review of PR #27 that remain after Batch 1 (commits `18e2239..5a5fa7f`).

**Architecture:** Three commit-sized modules — thread-safety fix (hamiltonian_blocks, strain_solver, main), validation fixes (defs, input_parser), parser cleanup (input_parser, hamiltonian_wire). Each module is independently committable and testable.

**Tech Stack:** Fortran 2018, pFUnit 4.x, CMake/Ninja, TOML (toml-f library), bash ctest

---

## Problem Statement

A comprehensive 8-agent workflow review (Codex P1 investigation, correctness core/physics, architecture, tests, performance, reliability, synthesis — 468 tool calls) of PR #27 identified 3 critical and 6 important issues that remain after Batch 1 fixes. The most severe is a thread-unsafe SAVE cache in `get_kp_block_table()` that races under OpenMP when multiple threads call `ZB8bandQW` simultaneously. Additionally, the `fermi_mode` enum silently falls back to `charge_neutrality` on typos (violating ADR 0002's "no silent corrections" principle), the `band_idx` range check skips bulk gfactor, 20 parser error exits use deprecated `stop 1`, polygon/b_field parse errors are silently ignored, and wire hot-path code uses heap allocation where stack suffices.

The review also confirmed that both Codex P1 findings are false positives: (1) `b_field` copy at lines 699-700 executes before the early return at line 703, and (2) `strainSubstrate` survives through `paramDatabase` because that routine never touches the field.

## Solution

Fix all confirmed findings in 3 logical commits: thread-safety first (it's the only one that can produce silent wrong results), then validation gaps, then parser cleanup. Add rejection tests for the validation fixes. Respond to Codex P1 comments explaining false positives.

## User Stories

1. As a **physicist running QW band structure with OpenMP**, I want the k.p block table to be thread-safe, so that multi-threaded runs produce correct results without data races.
2. As a **physicist running bulk g-factor**, I want `band_idx=2` with `num_cb=2` to be rejected at validation time, so that I get a clear error instead of an out-of-bounds crash.
3. As a **physicist configuring self-consistent calculations**, I want a typo in `fermi_mode` to produce an error, not silently fall back to `charge_neutrality`, so that I know my config is wrong.
4. As a **developer debugging parser errors**, I want all parser error exits to use `error stop` with descriptive messages, so that I get consistent diagnostics via stderr.
5. As a **user setting polygon vertices**, I want a type mismatch in vertex coordinates to produce an error, not silently default to 0.0, so that my wire geometry is correct.
6. As a **user setting b_field components**, I want a type mismatch in `Bx/By/Bz` to produce an error, not silently default to 0.0, so that my magnetic field is what I specified.
7. As a **physicist running wire calculations**, I want the k.p table lookup to use stack allocation, so that the kz-sweep hot path avoids unnecessary heap alloc/dealloc cycles.
8. As a **code reviewer**, I want the Codex P1 comments on the PR to be addressed with evidence, so that future reviewers know they were investigated.

## Implementation Decisions

### Commit 1: Thread-safety fix

**C1 — Pre-initialize block table caches before OpenMP parallel regions:**

The SAVE caches in `hamiltonian_blocks.f90` (`kp_table_cached`/`kp_table_cache`) and `strain_solver.f90` (`strain_table_cached`/`strain_table_cache`, `zeeman_table_cached`/`zeeman_table_cache`) use a lazy-init pattern that is not thread-safe. When `ZB8bandQW` is called inside `!$omp parallel` regions in `main.f90` (lines 616 and 675), multiple threads race on these SAVE variables.

Decision: Pre-initialize all three caches in serial code before the OpenMP parallel regions. Add initialization calls in `main.f90` before the QW parallel region (~line 614) and before the Landau parallel region (~line 673). The initialization is: calling `get_kp_block_table()`, `get_strain_table()`, and `get_zeeman_table()` once each, which populates the SAVE caches deterministically before any thread forks.

The SAVE cache pattern is kept (not eliminated) because:
- It works correctly after pre-initialization
- The cache is read-only after init (no mutation)
- Eliminating the cache would change more files and is a larger refactor

### Commit 2: Validation fixes + rejection tests

**C2 — Extend band_idx range check to all confinement modes:**

Remove the confinement guard (`confinement == 'wire' .or. confinement == 'qw'`) from S1 check in `validate_semantic`. The check `band_idx + 1 <= num_cb` is universally valid for gfactor (bulk, QW, wire all access `cb_state(:, bandIdx+1)`). The check is already gated by `app_name == 'gfactor'`, so it won't affect non-gfactor executables.

**C3 — Reject unknown fermi_mode values:**

Change the `case default` in `parse_sc` from `cfg%sc%fermi_mode = 0` to `error stop` with a message listing valid values (`charge_neutrality`, `fixed`). This matches the pattern used for other enums in the parser (`wave_vector mode`, `topology mode`, `sweep_model`, `conductance_method`).

Decision: Fix in parser (not `validate()`) because the parser is where string→integer conversion happens. After the parser, `fermi_mode` is an integer — `validate()` cannot know if the original string was valid.

**Rejection tests:**

Add two shell test cases to `test_validate_rejects_bad_configs.sh`:
- Bulk gfactor with `band_idx=2, num_cb=2` → expect error stop with bandIdx message
- SC config with `fermi_mode = 'typo'` → expect error stop with valid values message

### Commit 3: Parser cleanup

**I1 — Replace 20x `stop 1` with `error stop`:**

Replace all `print *, 'Error: ...'; stop 1` patterns with single `error stop 'Error: ...'` lines. The 20 occurrences are at lines: 40, 106, 155, 161, 175, 196, 202, 217, 264, 288, 306, 325, 339, 362, 379, 393, 914, 931, 948, 962 of `input_parser.f90`.

Decision: Use single `error stop` (not `print` + `error stop`) because `error stop` already sends the message to stderr. The `print *` before `stop 1` is a legacy pattern from the old `.cfg` parser.

**I2 — Add `check_optional_stat` for polygon vertices and b_field components:**

In `parse_wire` (~lines 297-300), polygon vertex `x`/`y` coordinates extract `stat` but never check it. Add `check_optional_stat(stat, 'x', 'wire.geometry.polygon')` (and same for `y`) after each `get_value`.

In `parse_b_field` (~lines 450-455), `Bx/By/Bz` components currently have `if (stat /= 0) ... = 0.0_dp` — silent fallback. Replace with `check_optional_stat(stat, 'Bx', 'b_field')` (and same for `By`, `Bz`). Key absence is fine (default 0.0 via `get_value`), but type mismatch must error.

**I4 — Change `allocatable table(:)` to `table(52)` in hot path:**

In `hamiltonian_wire.f90` lines 1019 and 1087, change `type(kp_entry), allocatable :: table(:)` to `type(kp_entry) :: table(52)`. The `get_kp_block_table()` function always returns exactly 52 entries — using an allocatable local triggers heap alloc/dealloc per call in the kz-sweep hot path.

### Codex P1 responses

Post inline comments on the two Codex P1 findings:
1. **b_field copy** (input_parser.f90:703): Explain that lines 699-700 perform the copy before the early return at line 703. The flow is correct — not a bug.
2. **strainSubstrate** (input_parser.f90:871): Explain that `paramDatabase` never touches `strainSubstrate`, so the value survives. Fragile but correct. Noted as suggestion S1 for future hardening.

## Testing Decisions

**What makes a good test:** Tests verify external behavior (correct error messages on invalid configs, correct exit codes) not implementation details. Thread-safety cannot be directly tested in CI — the fix is verified by code review (pre-init in serial) and the existing OpenMP regression tests (which would show different results if the cache were corrupted).

**Modules tested:**
- C1 (thread-safety): Verified by existing QW/Landau regression tests running with OpenMP. Pre-init ensures no race.
- C2 (band_idx): New rejection test — bulk gfactor TOML with band_idx=2, num_cb=2, expect non-zero exit + "bandIdx" in stderr.
- C3 (fermi_mode): New rejection test — TOML with fermi_mode='typo', expect non-zero exit + "fermi_mode" in stderr.
- I1/I2 (parser cleanup): Verified by running existing regression tests (no behavioral change for valid configs).
- I4 (stack alloc): Verified by existing wire regression tests (no behavioral change).

**Prior art:**
- `tests/integration/test_validate_rejects_bad_configs.sh` — existing rejection test pattern for `error stop` branches
- `tests/unit/test_defs.F90` — pFUnit tests for validation passing paths

## Out of Scope

- **I3 — Deep copy of 12 CSR matrices per kz-point:** Performance optimization, not a bug. Deferred to backlog. Wire functionality is correct.
- **I5/I6 — Missing rejection tests for ~20 validate() branches and all validate_semantic() branches:** Substantial work (~2h) better done in a dedicated PR.
- **S1-S12 suggestions:** Fragile coupling (strainSubstrate), god module (simulation_setup), block table spot-checks, convergence tolerance, thread_workspace was_freed guard — all deferred to backlog.
- **Codacy 11 high issues:** Static analysis warnings, not blocking. Review separately.

## Further Notes

- All fixes target the `feat/richardson-observables-convergence` branch via PR #27.
- Batch 1 fixes (commits `18e2239..5a5fa7f`) already addressed: B_vec copy positioning, strain_substrate scoping, doping extraction, fd_step threshold, V8 tolerance, I7-I12 semantic checks, stale confinement references, input-reference docs, CLAUDE.md updates, double-free guards, stop 1→error stop in simulation_setup/executables.
- After these 3 commits, re-run `ctest --test-dir build -j4` to verify all 107 tests pass.
- The grill session confirmed all decisions with the project owner before writing this PRD.
