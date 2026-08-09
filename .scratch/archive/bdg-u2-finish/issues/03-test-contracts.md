**Status**: closed — ticket 03 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed
Claimed-by: Claude (session 2026-07-13)
Resolved-by: Claude (session 2026-07-13) via wayfinder
Blocked by: 01 ✓ (closed 2026-07-14)
Related: design.md §"Solution" items 5, 6, 7; §"Testing"; §"Spec Amendment"

# 03 — Test contracts: SSOT error-stops + s2_sign tightening + meta-test

## Question

After ticket 01 introduces `bdg_pfaffian_params_t` and ticket 02 retires the dense `wire_pfaffian_witness`, what test amplification is needed to (a) pin the new type's error-stop contract, (b) tighten the loosened `s2_sign /= 0` assertion back to a specific sign, (c) pin the User Story 5 spec contract as a self-documenting meta-test, and (d) dedupe the "Pf_real = -0.0225" comment block?

## Background

**S7 finding (loosened test):** the strict assertion `@assertTrue(s1 == s2 .and. s1 /= 0)` was loosened to `@assertTrue(s2 /= 0)` because the S1+S2 strict sign-agreement is structurally unattainable on the 16×16 synthetic BdG fixtures (the antisymmetric part of `h_proj @ omega_local` vanishes by Nambu particle-hole symmetry — `Asym(1,3) = -i*(E_min + (-E_min))/2 = 0`). Per the previous map's verification subagent: a focused reconstruction on paper + 200 random BdG seeds + 2 curated fixtures confirmed the loosening is correct, but the loosened assertion is too permissive — it lets a `s2_sign = +1` (trivial) fixture pass. Tighten to `@assertEqual(-1, s2_sign)` on the canonical fixture, which is the documented slim-witness "topological" sign.

**S6 finding (duplicated comment):** the "Pf_real = -0.0225" analysis block in `test_wire_pfaffian_witness.pf` duplicates the same explanation in `test_bdg_pfaffian_witness_csr.pf`. Move it to a shared fixture doc (top of `test_bdg_pfaffian_witness_csr.pf`).

**H2 finding (EOF newline):** `tests/unit/test_bdg_evaluator.pf` lacks a trailing newline. Trivial.

**Spec amendment pinning:** the loosened test was correct, but the rationale needs to live in a self-documenting test name so future agents know they violated User Story 5 if they tighten further.

## Sub-decisions LOCKED (resolution 2026-07-13)

**Sub-decision 1 — Error-stop tests on the new factory (LOCKED).** `test_bdg_evaluator.pf` gains 2 tests:
- `test_bdg_pfaffian_params_rejects_zero_floor` — `bdg_pfaffian_params_with_floor(0.0_dp)` → `error stop`. Pattern matches existing `validate()` rejection tests.
- `test_bdg_pfaffian_params_rejects_negative_floor` — `bdg_pfaffian_params_with_floor(-1.0_dp)` → `error stop`.
Both must be single-line `@assertEqual` per `feedback_pfunit_macro_single_line`. **Rationale:** the factory is the validation site (per ticket 01 sub-decision 4); tests pin the error-stop contract so future refactors cannot silently change it. The structure constructor (`bdg_pfaffian_params_t()` no-arg, default initial value) intentionally does NOT validate — it always returns the safe SSOT — so error-stops target the factory exclusively.

**Sub-decision 2 — Re-route the 2 SSOT-pin tests through the new type (LOCKED).** Existing 2 SSOT-pin tests (added in commit `c25b337` per ticket 02 §3 of the previous map) test `bdg_default_pfaffian_floor` via `bdg_eval_params_t%pfaffian_floor`. After ticket 01 drops the field from `bdg_eval_params_t` and introduces `bdg_pfaffian_params_t`, the tests must:
- Construct via factory `bdg_pfaffian_params_t()` (default = SSOT, via default initial value `pfaffian_floor = bdg_default_pfaffian_floor`).
- Pass the new type to `eval_bdg_pfaffian_witness_csr` (signature change per ticket 01 sub-decision 1).
- Verify the SSOT flows through (e.g., set a custom `pfaffian_floor = 1.0e-6_dp` and verify a near-zero-Pf fixture evaluates to `s2_sign = 0` instead of `s2_sign = -1`). **Rationale:** the previous tests pinned the value of the constant (`bdg_default_pfaffian_floor == 1.0e-12_dp`) but not its flow into the seam sibling. The re-routed tests pin both, which is what the H1 finding was actually about (declared-but-not-consumed SSOT).

**Sub-decision 3 — Tighten `_nondiagonal_nonzero_sign` to `@assertEqual(-1, s2_sign)` (LOCKED).** In `test_bdg_pfaffian_witness_csr.pf`'s `_nondiagonal_nonzero_sign` test, change `@assertTrue(s2_sign /= 0)` → `@assertEqual(-1, s2_sign)`. The canonical fixture's `Pf_real < 0` (per "Pf_real = -0.0225" comment block) maps to `s2_sign = -1` by the construction in `wire_pfaffian_witness_sweep` (`if real(pf_val) > 0: sign = 1, else: sign = -1`). **Rationale:** the loosened `s2 /= 0` was correct as a relaxation of the unattainable strict `s1 == s2` sign agreement (per the previous map's verification subagent: 200 random BdG seeds + 2 curated fixtures confirmed strict is structurally unattainable on 16×16 synthetic BdG due to Nambu particle-hole antisymmetric vanishing), but it is too permissive — it lets a `s2 = +1` (trivial-phase) fixture pass as GREEN. Tightening to `-1` pins the documented slim-witness "topological" sign on the canonical fixture while preserving the loosened-test rationale for non-canonical fixtures (covered by the meta-test in sub-decision 5).

**Sub-decision 4 — Add `test_seam_sibling_consumes_pfaffian_floor` (LOCKED).** New 5th test:
- Construct `bdg_pfaffian_params_t(pfaffian_floor = 1.0e-6_dp)`.
- Call `eval_bdg_pfaffian_witness_csr(canonical_fixture, ...)` with the custom floor.
- Verify the seam sibling accepts without error (proves the SSOT arg is consumed, not discarded — the H1 receipt test).
- A near-zero-Pf fixture with the looser custom floor should produce `s2_sign = 0` instead of `s2_sign = -1`, proving the floor flows into the sign extraction.
**Rationale:** this is the test that would have caught H1 (declared-but-not-consumed SSOT) if it had existed before PR #42 review. Per sub-decision 2 it covers both "the floor arg is consumed" AND "the value matters" — the two halves of the SSOT contract.

**Sub-decision 5 — Add `test_pfaffian_witness_spec_user_story_5_contract` (LOCKED).** Meta-test that pins the User Story 5 invariant as a self-documenting test name. Asserts:
- `s2_sign ∈ {-1, 0, +1}` (range contract).
- `s2_sign /= 0` on non-diagonal synthetic fixture (non-zeroness contract — preserves the loosened-test rationale).
- `s2_sign == -1` on canonical fixture (specific-sign contract — preserves the sub-decision 3 tightening).
The test name includes `spec_user_story_5_contract` so future audits can `grep` for it. **Rationale:** the meta-test name tells future agents which User Story they violated if they tighten further (e.g., reintroducing `@assertTrue(s1 == s2)` would surface this test name in the failure log with a pointer back to User Story 5). Doc propagation of the meta-test name → BACKLOG row 82 + L11 summary + parent plan status_note is ticket 07's responsibility; the test name itself is a derivative artifact of ticket 06's spec-amendment text.

**Sub-decision 6 — Move "Pf_real = -0.0225" analysis block to top-of-file in `test_bdg_pfaffian_witness_csr.pf` (LOCKED).** The S6 finding: the analysis block in `test_wire_pfaffian_witness.pf` duplicates the same explanation in `test_bdg_pfaffian_witness_csr.pf`. After ticket 02 deletes `test_wire_pfaffian_witness.pf` entirely (option A per ticket 02 sub-decision 3), the dedup reduces to: ensure the canonical home is `test_bdg_pfaffian_witness_csr.pf`'s top-of-file comment block, and that any future test that needs the analysis can grep for the constant `Pf_real = -0.0225` in the canonical home. **Rationale:** the analysis is the canonical explanation of why `s2_sign` tightens to `-1` on this fixture; keeping it in one place (with cross-references in any new test that depends on it) prevents drift.

**Sub-decision 7 — EOF newline on `tests/unit/test_bdg_evaluator.pf` (LOCKED).** Trivial H2 finding. The current file lacks a trailing newline (visible in `git diff` as `\ No newline at end of file`). **Rationale:** POSIX text-file convention; some tools (pFUnit preprocessor, POSIX `cat`) misbehave on missing trailing newlines. Single-character fix.

## Out of scope for this ticket

- The 2 Codacy high `ErrorProne` flags — ticket 04.
- Doc propagation of the meta-test name — ticket 07 (BACKLOG row 82 + L11 summary + parent plan status_note cite the spec-amendment text; the test name is a derivative artifact).
- Unit-count reconciliation — ticket 05 runs `ctest -L unit` AFTER this ticket's tests land and records the new total.

## Files to touch

| File | Change |
|---|---|
| `tests/unit/test_bdg_evaluator.pf` | Re-route 2 SSOT tests through `bdg_pfaffian_params_t` factory; add 2 error-stop tests; EOF newline (H2) |
| `tests/unit/test_bdg_pfaffian_witness_csr.pf` | Tighten `_nondiagonal_nonzero_sign` assertion to `@assertEqual(-1, s2_sign)` (S7); add `test_seam_sibling_consumes_pfaffian_floor` (5th test); add `test_pfaffian_witness_spec_user_story_5_contract` (meta-test); top-of-file "Pf_real = -0.0225" comment block |
| `tests/unit/test_wire_pfaffian_witness.pf` | Dedup "Pf_real = -0.0225" comment to cross-reference; verify 0 references to deleted dense subroutine (after ticket 02) |
| `tests/CMakeLists.txt` | (No change) — the new tests live in existing files; no new `add_pfunit_ctest` blocks needed |

## Constraints

- pFUnit `@assertEqual`/`@assertTrue` must be single-line (no `&` continuation) per `feedback_pfunit_macro_single_line`.
- Per ADR 0001: test types are plain records; no test-fixture class hierarchies.
- The meta-test name must include `spec_user_story_5_contract` so it can be grep'd from any future audit.

## Verification

- All test changes compile (no pFUnit macro errors).
- `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 -L unit --output-on-failure` shows new tests PASS and existing **51** stay PASS (baseline is 51 after ticket 02 closed one ctest target — ctest counts `add_pfunit_ctest` targets, not `@test` subroutines; the new tests here live in existing files so no new target, total stays 51).
- `wc -l tests/unit/test_bdg_evaluator.pf` shows +1 line (EOF newline only — SD-1 dropped, see Execution).
- `wc -l tests/unit/test_bdg_pfaffian_witness_csr.pf` shows the SD-3/4/5 additions + the SD-6 header reconciliation.
- `test_wire_pfaffian_witness.pf` — N/A: deleted by ticket 02 (option A). SD-6's dedup reduces to canonical-home already present in `test_bdg_pfaffian_witness_csr.pf`.

## Execution (2026-07-25, wayfinder work-through)

Claimed and executed on `feat/bdg-u2-actual-ship` (uncommitted, same working tree as tickets 01/02).

**SD-1 — resolved as INFEASIBLE (not executed).** pFUnit 4.16 (installed at `$HOME/.local/pfunit/PFUNIT-4.16`) cannot catch `error stop`; this is a documented codebase invariant (`tests/unit/test_parameters.pf:223` — "error stop cannot be caught by pFUnit 4.x"). The repo's rejection tests route to **integration shell scripts** via exit-code checks (`tests/integration/test_topology_validate_rejects.sh`), which target TOML-config-driven executables — but `bdg_pfaffian_params_with_floor` is not config-reachable today (the floor is SSOT-defaulted with no TOML kn supported; adding one would be YAGNI / scope creep). A literal `@test` calling the factory with `pfaffian_floor <= 0` would **`error stop` the whole `ctest` process on every run — turning GREEN into RED**, a defect not a test. Decision (HITL 2026-07-25): drop SD-1 from this ticket; the error-stop contract remains enforced by the factory source at `src/physics/bdg_observables.f90:153-154` (`if (pfaffian_floor <= 0.0_dp) error stop '...'`). The two tests SD-1 named are not authored. Last line of `test_bdg_evaluator.pf` gains the SD-7 trailing newline only.

**SD-2 — already executed by ticket 01** (verified, not re-touched). `test_bdg_evaluator.pf:174-182` (`test_bdg_pfaffian_params_default_uses_ssot`) routes through `bdg_pfaffian_params_t()` factory; `test_bdg_pfaffian_witness_csr.pf:67-88` (`..._delegation_contract`) pins byte-identical delegation to `wire_pfaffian_witness_sweep` with the floor threaded through. Together they pin arg-consumption + delegation; SD-4 adds the value-matters half.

**SD-3 — executed.** `test_bdg_pfaffian_witness_csr.pf` `_nondiagonal_nonzero_sign`: `@assertTrue(s2_sign /= 0)` → `@assertEqual(-1, s2_sign)` (single-line). Canonical fixture `build_nondiagonal_bdg_csr` (Pf_real = -0.0225 < 0 ⇒ -1 by the `real(pf_val) > 0 ⇒ +1 else ⇒ -1` branch in `wire_pfaffian_witness_sweep`).

**SD-4 — executed** as `test_seam_sibling_consumes_pfaffian_floor`. Both halves of the SSOT contract on one canonical fixture: default floor (1.0e-12_dp) ⇒ `@assertEqual(-1, s2_sign_default)`, loose floor (`bdg_pfaffian_params_with_floor(1.0_dp)`) ⇒ `@assertEqual(0, s2_sign_loose)`. The sign flip −1→0 as floor crosses |Pf|≈0.0225 proves the value drives the branch. (Chose a *large* floor exceeding |Pf| rather than engineering a fragile tiny-Pf fixture: same semantic — `best_pf > pfaffian_floor` false ⇒ default `s2_sign = 0` — robust to arithmetic drift.)

**SD-5 — executed** as `test_pfaffian_witness_spec_user_story_5_contract` (name grep-able). Pins range (s2∈{-1,0,+1}) on both `build_zero_bdg_csr` (⇒0) and `build_nondiagonal_bdg_csr` (⇒-1) fixtures, plus `@assertEqual(-1, s2_nd)` for the specific-sign contract. Single-line `@assertEqual`/`@assertTrue` throughout.

**SD-6 — executed by reconciliation.** The "Pf_real = -0.0225" analysis block is already canonical-home in `test_bdg_pfaffian_witness_csr.pf` (at the `_nondiagonal..._sign` test L100 and the fixture L220; ticket 02 deleted the duplicate `test_wire_pfaffian_witness.pf`). Reconciled the now-stale module-header comment that claimed "`params` is therefore unused inside the seam sibling at this revision" — SD-4's test explicitly contradicts that post-ticket-01. Header now states the floor is threaded through and cross-references SD-4.

**SD-7 — executed.** `test_bdg_evaluator.pf` trailing newline appended (was `\ No newline at end of file` in the diff).

## Interaction

- **Ticket 01** is blocked-by: this ticket's tests reference `bdg_pfaffian_params_t` which doesn't exist until ticket 01. Strict ordering.
- **Ticket 05** re-counts `ctest -L unit` after this ticket's tests land. The new total is the source-of-truth for doc propagation.
- **Ticket 06** writes the spec-amendment text. The meta-test name `test_pfaffian_witness_spec_user_story_5_contract` cross-references the spec's User Story 5 numbering.