**Status**: closed (decision locked 2026-07-14; code executed + verified 2026-07-25) — ticket 01 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed
Claimed-by: Claude (session 2026-07-14)
Resolved-by: Claude (session 2026-07-14) via wayfinder
Executed-by: Claude (session 2026-07-25) via wayfinder work-through — code landed in working tree on `feat/bdg-u2-actual-ship` (uncommitted). Build clean; ctest -L unit 52/52 PASS; legacy `test_bdg_eval_params_factory_uses_pfaffian_floor` deleted (asserted the dropped field), its coverage reasserted on the correct type by `test_bdg_pfaffian_params_default_uses_ssot`.
Blocked by: —
Related: design.md §"Solution" items 1-3; §"Data Flow"; §"Error Handling" (the one new error stop)

# 01 — SSOT wire-up: `bdg_pfaffian_params_t` from default to call site

## Question

How do we turn the *declared-but-not-consumed* `bdg_default_pfaffian_floor` SSOT at `src/physics/bdg_observables.f90:45` into an actually-consumed SSOT that flows from declaration through the seam sibling down to the 3 hard-coded literal sites at `topological_analysis.f90:1680, :1698, :1781`?

## Background

This ticket addresses **finding H1** (the cross-axis worst issue, flagged independently by both Standards and Spec review sub-agents) plus findings S1, S2, S3, S4 from the PR #42 review. The current state: the SSOT is declared as a module parameter and advertised in 5 places (its own docstring, `tests/unit/test_bdg_evaluator.pf`, lecture doc, parent plan status_note, AGENTS.md), but the seam sibling `eval_bdg_pfaffian_witness_csr` accepts it via `params%pfaffian_floor` and silently discards it. The 3 hard-coded `1.0e-12_dp` literals in `topological_analysis.f90` are unchanged.

## Sub-decisions LOCKED (resolution 2026-07-14)

**Sub-decision 1 — Use type `bdg_pfaffian_params_t` with single field `pfaffian_floor: real(dp)` (option a, design.md locks).**
- Type lives in `src/physics/bdg_observables.f90`, adjacent to `bdg_eval_params_t`.
- Default initial value: `pfaffian_floor = bdg_default_pfaffian_floor` so the structure constructor `bdg_pfaffian_params_t()` (no args) returns the SSOT — which is what the pre-written happy-path test at `tests/unit/test_bdg_evaluator.pf:184-191` (`test_bdg_pfaffian_params_default_uses_ssot`) expects.
- Projection-band indices 7-8 stay as a comment in `wire_pfaffian_witness_sweep`; the canonical SSOT remains in `hamiltonian_blocks.f90`.
- Rationale: per design.md and CLAUDE.md YAGNI; no second field is concretely needed.

**Sub-decision 2 — Drop `pfaffian_floor` from `bdg_eval_params_t`.**
- Field has been dead since the seam sibling was added in `52357b4` (params accepted, silently discarded — receipt is the "NOT consumed" comment at `bdg_observables.f90:143-145`).
- The 4 lambda-style one-liners at `main_topology.f90:530, :552, :826, :1360` only read `delta_0` and `near_zero_frac`; the field is unused.
- `bdg_eval_params_with_delta` factory drops its `p%pfaffian_floor = bdg_default_pfaffian_floor` line.
- Rationale: SRP (S3 review finding); cross-domain fields pollute `bdg_eval_params_t`.

**Sub-decision 3 — Do NOT mark `eval_bdg_pfaffian_witness_csr` or `eval_bdg_kitaev_majorana` as `pure` in this ticket. Document the dependency.**
- The blocker is NOT internal allocations (those are F2018-pure-compatible). The blocker is the call chain: `eval_bdg_pfaffian_witness_csr` → `wire_pfaffian_witness_sweep` → `complex_pfaffian` (in `src/math/pfaffian.f90`, not pure). Same for `eval_bdg_kitaev_majorana` → `kitaev_majorana_number` (also not pure).
- A pure wrapper requires `pure` on `complex_pfaffian` and its helpers (`zskpfa_laplace`, `zskpfa_reduction`, `real_pfaffian`, etc.) — that's `src/math/pfaffian.f90`, a different module. Ticket 01 owns `bdg_observables.f90` and `wire_pfaffian_witness_sweep` plumbing only.
- Decision: leave both seam siblings non-pure. Document in wrapper docstring with forward-pointer to Codacy triage (ticket 04) for the future-swap path.
- `eval_bdg_point` stays `pure` (already pure, no callee chain).
- Rationale: YAGNI; scope-creep into `pfaffian.f90` is a separate PR.

**Sub-decision 4 — Add factory `bdg_pfaffian_params_with_floor` with one error stop.**
- Factory signature: `bdg_pfaffian_params_with_floor(pfaffian_floor) result(p)` with optional `pfaffian_floor` arg.
- Validation: `if (pfaffian_floor <= 0.0_dp) error stop 'bdg_pfaffian_params_t: pfaffian_floor must be > 0'`.
- Default (no arg): uses `bdg_default_pfaffian_floor` SSOT.
- Pattern matches defensive validation style of `validate_semantic` BdG checks in `defs.f90:949-985`.
- Rationale: per design.md §"Error Handling". The factory is the validation site; the structure constructor (`bdg_pfaffian_params_t()` with default initial value) does not validate but always returns the safe SSOT.

## Wiring contract (for whoever implements)

- **Type:** `bdg_pfaffian_params_t` in `src/physics/bdg_observables.f90`, single field `pfaffian_floor: real(dp) = bdg_default_pfaffian_floor`.
- **Factory:** `bdg_pfaffian_params_with_floor(floor)` — optional arg, validates, returns typed record. Public export.
- **Public list update:** add `bdg_pfaffian_params_t` and `bdg_pfaffian_params_with_floor` to the `public ::` exports.
- **Seam sibling signature change:** `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params)` where `params :: bdg_pfaffian_params_t`. Docstring beef-up: cite User Story 1 (seam contract), the L3-import exception (forward-pointer to ticket 04), and the U13 future-swap path (Issue 05).
- **`wire_pfaffian_witness_sweep` signature change:** add `pfaffian_floor: real(dp), intent(in)` arg. Use it at `topological_analysis.f90:1781`. The other 2 literals (`:1680`, `:1698`) sit inside `wire_pfaffian_witness`, which ticket 02 deletes — they will not be parameterized in ticket 01's scope.
- **`main_topology.f90:1371` call site:** pass `bdg_pfaffian_params_t(pfaffian_floor = bdg_default_pfaffian_floor)` to `eval_bdg_pfaffian_witness_csr`. Import list gains `bdg_pfaffian_params_t` and `bdg_default_pfaffian_floor`.
- **`bdg_eval_params_t` cleanup:** drop `pfaffian_floor` field; drop the `p%pfaffian_floor = bdg_default_pfaffian_floor` line in `bdg_eval_params_with_delta`.
- **`src/physics/AGENTS.md`:** `bdg_observables.f90` row + DAG block — deferred to ticket 07.

## Out of scope for this ticket (locked)

- The 2 hard-coded literals at `topological_analysis.f90:1680, :1698` — both inside the dense `wire_pfaffian_witness` that ticket 02 deletes. Wire-up contract only parameterizes the surviving 1 site (`:1781`).
- Test amplification (SSOT error-stop tests, s2_sign tightening, meta-test) — ticket 03.
- Codacy triage — ticket 04.
- `pure` attribute on `src/math/pfaffian.f90` helpers — separate module, separate PR.

## Files to touch

| File | Change |
|---|---|
| `src/physics/bdg_observables.f90` | Add `bdg_pfaffian_params_t` + `bdg_pfaffian_params_with_floor` factory; drop `pfaffian_floor` from `bdg_eval_params_t` and from `bdg_eval_params_with_delta`; update `eval_bdg_pfaffian_witness_csr` signature to `(H_bdg_csr, Nbdg, bdg_pfaffian_params_t) → s2_sign`; do **NOT** mark seam siblings pure (sub-decision 3); tighten docstring (User Story 1 seam contract citation, L3-import exception forward-pointer to ticket 04, U13 forward reference) |
| `src/physics/topological_analysis.f90` | `wire_pfaffian_witness_sweep` accepts `pfaffian_floor: real(dp), intent(in)` arg; uses it at `:1781`. The 2 literals at `:1680` and `:1698` sit inside `wire_pfaffian_witness` (dense) — ticket 02 retires that subroutine, so they are NOT parameterized here |
| `src/apps/main_topology.f90` | `:1371` call site passes `bdg_pfaffian_params_t(pfaffian_floor = bdg_default_pfaffian_floor)` to seam sibling; import list gains `bdg_pfaffian_params_t` and `bdg_default_pfaffian_floor` |
| `src/physics/AGENTS.md` | (deferred to ticket 07) `bdg_observables.f90` inventory row + DAG block updated for the new type |

## Constraints

- TDD per CLAUDE.md. Test for `bdg_pfaffian_params_t(pfaffian_floor <= 0) → error stop` lands in ticket 03, but the type and its factory MUST exist before ticket 03 writes the test. The happy-path placeholder test `test_bdg_pfaffian_params_default_uses_ssot` already lives at `tests/unit/test_bdg_evaluator.pf:184-191` (added by prior session) — it MUST pass after ticket 01 lands.
- ADR 0001: dispatch by enum-tagged plain record or procedure choice. `bdg_pfaffian_params_t` is a plain record with one field; valid.
- All modules use `private` default with explicit `public ::` exports. The new type's default `pfaffian_floor = bdg_default_pfaffian_floor` keeps the factory one-liner.
- Sub-decision 3: do NOT mark `eval_bdg_pfaffian_witness_csr` / `eval_bdg_kitaev_majorana` `pure` — `pfaffian.f90` callees are not pure and ticket 01 doesn't own that module.

## Verification

- `cmake --build build` clean (no new warnings).
- `ctest -L unit` runs the existing 52 tests + new "happy path" placeholder test (now compiles and passes) = 53 tests, all PASS.
- `tests/unit/test_bdg_evaluator.pf` existing 2 SSOT-pin tests still GREEN (added in commit `c25b337` per design.md §"Components").
- `tests/unit/test_bdg_evaluator.pf` happy-path placeholder `test_bdg_pfaffian_params_default_uses_ssot` (line 184-191) GREEN — this is the receipt that the type + default initial value landed.
- Manual: `grep -nE '1\.0e-12_dp' src/physics/topological_analysis.f90` shows **2 surviving literals** at `:1680` and `:1698` (inside the dense `wire_pfaffian_witness` — ticket 02's responsibility). The `:1781` site in `wire_pfaffian_witness_sweep` reads `pfaffian_floor`.
- Manual: `grep -nE 'bdg_eval_params_t' src/apps/main_topology.f90 | grep pfaffian_floor` returns nothing (the field is gone from the type).

## Interaction

- **Ticket 02** deletes `wire_pfaffian_witness` and the 2 surviving literals. If ticket 02 lands first, only 1 literal at `:1781` remains; ticket 01's wire-up contract still applies to `wire_pfaffian_witness_sweep` only. If ticket 01 lands first (recommended): the `:1781` site consumes `pfaffian_floor`, the `:1680`/`:1698` literals stay as orphan-literals-pending-deletion until ticket 02 retires the dense subroutine.
- **Ticket 03** writes tests that pin the new type's error stop. The tests will fail to compile until ticket 01's type + factory exist. **Strict ordering required.**
- **Ticket 05** re-runs `ctest -L unit` and counts tests. Adds ≥1 (the happy-path placeholder, now compilable); likely adds more from ticket 03.
- **Ticket 04 (Codacy triage):** if `pure` placement on the seam siblings is a flagged `ErrorProne` (likely candidate), the wiring contract here already documents the trade-off (sub-decision 3); ticket 04 records the flag status but no further code change is required from ticket 01's resolution.