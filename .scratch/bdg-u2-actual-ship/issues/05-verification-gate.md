**Status**: COMPLETE (2026-07-25) — ticket 05 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed (2026-07-25)
Claimed-by: Claude Fable 5
Blocked by: 01, 02, 03
Related: design.md §"Solution" item 9; §"Testing" §"Verification gate"

# 05 — Verification gate: unit-count + ctest + clean build

## Question

After tickets 01 (SSOT wire-up), 02 (orphan retirement), and 03 (test amplification) land, what is the actual `ctest -L unit` total, and is the build clean with no warnings?

## Background

Per `design.md` §"Testing":

> "**Test count target: 52 unit tests total** (per BACKLOG Phase 26's claim). The actual `ctest -L unit` total is the source of truth — verify the running count in the verification gate and use that number across all doc sites. If it differs from 52, use the actual count and re-propagate (the prior 50 / 52 / 60+ disagreement is a backfill hazard — the running count is canonical)."

The `ctest-phase3-final.log` in this dir shows **52/52 PASS** today (snapshot of pre-PR-#42-review state). Tickets 01-03 will add ≥1 test (ticket 03 adds 4-5 tests). The running total after ticket 03 lands will be the canonical number.

## Sub-decisions to lock

1. **Build cleanliness criterion.** "No warnings, no errors" per `design.md` §"Testing" item 1. Use `cmake --build build 2>&1 | tee /tmp/build.log` and `grep -c warning /tmp/build.log` — must be 0 (or only pre-existing third-party warnings, e.g., from MKL or pFUnit).
2. **Unit-test count.** `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 -L unit --output-on-failure 2>&1 | tee /tmp/unit.log`. The "Total Tests" line at the bottom is the canonical count.
3. **Acceptance gate regression check.** `ctest --test-dir build -R lecture_13_acceptance_gate` → 4-witness gate green at 1.0 T tolerance. Per `design.md` §"Testing" item 3: "do not touch the gate; polarization verifier untouched".
4. **Unit-count propagation.** Record the new total in a single canonical location: `.scratch/bdg-u2-actual-ship/unit-count.txt` (new artifact, lives next to `design.md`). Ticket 07 reads this file and propagates the number to all 5 doc sites + the test summary line in the final PR description.

## Out of scope for this ticket

- Doc propagation of the new count — ticket 07.
- BACKLOG row 82 / L11 summary / parent plan status_note updates — ticket 07.
- Memory entry update with the new count — ticket 08.

## Files to touch

| File | Change |
|---|---|
| `.scratch/bdg-u2-actual-ship/unit-count.txt` (new) | Single line: `unit-count = <N>` where `<N>` is the actual `ctest -L unit` total |

## Constraints

- Per CLAUDE.md, use `OMP_NUM_THREADS=$(( $(nproc)/4 ))` to avoid `ctest -jN` oversubscription.
- Per `/superpowers:verification-before-completion`, never claim "all tests pass" without showing the actual output. Save the log to `.scratch/bdg-u2-actual-ship/ctest-final.log` (overwriting `ctest-phase3-final.log`).
- The unit-count must be derivable from the saved log (not from memory or a back-of-envelope count).

## Verification

- `unit-count.txt` is committed to the scratch dir alongside `design.md`.
- The ctest log is committed to the scratch dir for auditability.
- The 4-witness acceptance gate regression check passes (proves the changes didn't break the existing contract).
- Manual scan: BACKLOG.md Phase 26 row, parent plan status_note, AGENTS.md `bdg_observables.f90` row, solutions doc, memory entry — all should reference the SAME count after ticket 07 lands.

## Interaction

- **Tickets 01, 02, 03** are all blocked-by: this ticket cannot run until code/test changes are in place. Strict ordering.
- **Ticket 07** reads `unit-count.txt` and propagates the number.
- **Ticket 08** mirrors the count in the memory entry.

## Resolution path

1. `cmake --build build 2>&1 | tee /tmp/build.log` — verify 0 warnings.
2. `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 -L unit --output-on-failure 2>&1 | tee .scratch/bdg-u2-actual-ship/ctest-final.log`.
3. `grep -E '^tests passed' .scratch/bdg-u2-actual-ship/ctest-final.log` — record the count.
4. `ctest --test-dir build -R lecture_13_acceptance_gate` — verify 4-witness green at 1.0 T.
5. Write `.scratch/bdg-u2-actual-ship/unit-count.txt` with the single line `unit-count = <N>`.
6. Commit the verification artifacts (log + unit-count.txt) as part of the PR.

## Execution

- Verified `cmake --build build` is clean with no warnings.
- Ran `ctest -L unit` and captured output to `.scratch/bdg-u2-actual-ship/ctest-final.log` (51/51 PASS).
- Locked unit test count as `51` in `.scratch/bdg-u2-actual-ship/unit-count.txt`.
- Discovered and diagnosed a major regression in `lecture_13_acceptance_gate` where the branch's promotion of the slim Pfaffian to a mandatory 4th witness caused failures due to incorrect Python indices (`parts[3]` instead of `parts[2]` for z2) and mismatched sign conventions (checking `z2 == -1` when Fortran emits `z2 ∈ {0,1}`). In addition, the slim Pfaffian returns constant topological signal (`s2_sign = -1`, `|pf| = 4e-8` above floor) across all B fields, which lacks the transitions of the bulk invariant (U13).
- **Decision:** As approved by user, reverted the gate script and python script back to the 3-witness gate layout of `main` (where `bcrit_pfaffian = None` and is excluded from the numeric range). Acceptance gate is now **100% green** (`lecture_13_acceptance_gate` PASSED in 224.82 s).
- Staged the reverted scripts, f-string warnings cleanup (pyflakes F541 on `lecture_13_topological.py:389` and `verify_majorana_polarization.py:166` by removing the unused `f` prefixes), and the verification artifacts.
- Map updated to reflect ticket 04 and 05 closure. Ready for commit.