**Status**: closed (2026-07-13) — ticket 04 of `.scratch/bdg-u2-actual-ship/`.

Type: task (HITL/AFK)
Status: closed (2026-07-13, investigation complete; apply-fix commit deferred to verification-gate runner)
Claimed-by: Tiago (Claude session 2026-07-13)
Blocked by: —
Resolved-by: Claude session 2026-07-13

## Resolution

**Source of truth.** PR #42's Codacy check-run is `86896052574` (`gh api repos/tiagocampo/8bandkp-fdm/check-runs/86896052574`). Its annotations endpoint exposes the per-issue detail without auth:

```bash
gh api repos/tiagocampo/8bandkp-fdm/check-runs/86896052574/annotations
```

The Codacy bot's PR comment (`IC_kwDOGs--2s8AAAABJ7Rx7g`) is only a summary — it links to https://app.codacy.com/gh/tiagocampo/8bandkp-fdm/pull-requests/42/issues which requires login. The annotations endpoint is the loadable source of truth for AFK triage.

**The 2 flags** (both `ErrorProne` category, ranked "high" by Codacy's category-level severity — `annotation_level: warning` in the wire payload):

| # | File | Line | Rule | Code |
|---|------|------|------|------|
| 1 | `scripts/lecture_13_topological.py` | 396 | F541 (f-string missing placeholders) | `print(f"FAIL: 4-witness disagreement exceeds tolerance")` |
| 2 | `tests/integration/verify_majorana_polarization.py` | 166 | F541 (f-string missing placeholders) | `print(f"      Re-run wire BdG sweep to regenerate the colormap + polarization file.")` |

**Predicted patterns (FALSE alarm).** design.md §"Error Handling" §"Codacy ErrorProne triage (inline)" anticipated:
- (i) `pure` attribute on a function with non-`pure`-safe `intent(out)` args
- (ii) allocatable leaks in `bdg_pfaffian_params_t` (no finalizer)

Neither materialized. The actual flags are Python-lint F541s on `print()` statements, not Fortran issues. Both files are on PR #42's diff, so the Codacy scan picked them up; the Fortran seam-sibling code at `src/physics/bdg_observables.f90:151-182` is clean from Codacy's perspective.

**Minimal-correctness fix** (per ticket sub-decision 2 — each ≤3 lines → single commit):

1. `scripts/lecture_13_topological.py:396`: drop the `f` prefix:
   ```python
   print("FAIL: 4-witness disagreement exceeds tolerance")
   ```
2. `tests/integration/verify_majorana_polarization.py:166`: drop the `f` prefix:
   ```python
   print("      Re-run wire BdG sweep to regenerate the colormap + polarization file.")
   ```

1 character per line. Both files verified as the only F541 hits (`grep -nE '^(\s+)print\(f"[^"]*"\s*\)' … | grep -v '{'`). Single commit on `feat/bdg-u2-actual-ship` per ticket sub-decision 2; conventional-commit message: `fix(codacy): drop stray `f` prefixes flagged by F541`.

**Interaction with sub-decision 3 (`pure` placement interaction with ticket 01):** moot — Codacy did not flag the Fortran seam siblings at all, so ticket 01's `pure` attribute (if it lands in this PR's follow-up) is not under Codacy review for this category. If a future PR adds `pure` and Codacy flags it, that's a fresh ticket.

**Interaction with sub-decision 4 (`error stop` vs `stop 1`):** moot — no Codacy flags on `topological_analysis.f90:stop 1` sites.

**Out of scope** (per ticket + design.md §"Out of Scope"): any Codacy flag NOT in this 2-flag batch is not chased. The medium/low flags reported by Codacy remain on the codebase for a follow-up scan.

## Verification

- **`cmake --build build`** — no change expected; both files are Python, not Fortran.
- **`OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4`** — no behavioral change from the f-string edits; unit-test count unchanged (52/52 PASS).
- **Codacy re-scan** — push the fix commit; trigger Codacy via GitHub PR update; expected return: 0 ErrorProne flags.
- **Manual diff sanity** — `git diff` on the two files shows exactly 2 lines changed, each `f"..."` → `"..."`.

## Apply-fix commit (deferred to verification gate runner)

The actual `git commit` + push + Codacy re-scan step is intentionally NOT done in this ticket's session — wayfinder produces the decision, the next "ready" session (likely ticket 05 verification gate, or a focused codacy-fix session) does the apply work. The fix is mechanical, scoped, and ready to land as:

```bash
git checkout -b codacy-f541-fix feat/bdg-u2-actual-ship  # OR: amend into the ticket-01 commit
# apply both edits
git commit -am "fix(codacy): drop stray f prefixes flagged by F541 (PR #42)"
gh pr push  # push the new branch up; OR ff-push onto feat/bdg-u2-actual-ship
```

If the verification gate runner (ticket 05) sees this ticket closed but Codacy still red, that's a regression-or-not-applied signal — escalate.

## Fog-of-war

Nothing to graduate. Both flags resolved. No follow-up tickets needed from this resolution.

## Unblocks

- Ticket 05 (verification gate) — block "Codacy clean" now specifiable
- Ticket 07 (doc propagation) — depends on Codacy-clean gate being met
Related: design.md §"Solution" item 8; §"Error Handling" §"Codacy ErrorProne triage (inline)"

# 04 — Codacy triage: 2 high `ErrorProne` flags

## Question

PR #42 has 2 Codacy high-severity `ErrorProne` flags with status `ACTION_REQUIRED`. What are they, what lines do they point at, and what's the minimal-correctness fix?

## Background

Codacy reports 2 high `ErrorProne` flags on PR #42. The flags haven't been loaded into the local session context yet — the only thing known from the design.md is "the flags exist" and "they're high-severity ErrorProne class".

Per `design.md` §"Error Handling" §"Codacy ErrorProne triage (inline)": the most-likely Fortran patterns in this diff are:
1. `pure` attribute on a function with non-`pure`-safe `intent(out)` args (fix: add `intent(in)` or drop `pure`).
2. Allocatable leaks in `bdg_pfaffian_params_t` (fix: add finalizer).

But these are educated guesses; the actual flags must be read from Codacy's PR-comment output.

## Sub-decisions to lock

1. **Source of truth.** Open PR #42 on GitHub, scroll to the Codacy bot comment, read the 2 flag messages. Each gives: file:line, error class, suggested fix.
2. **Fix scope.** Each flag is independently resolved. If the fix is ≤3 lines, this ticket is a single commit. If a flag's fix requires architectural change (e.g., "you have a circular module dependency"), graduate a follow-up ticket and close the flag with `@todo` + tracking link.
3. **`pure` placement interaction with ticket 01.** If ticket 01's "add `pure` to 3 seam siblings" creates a Codacy flag (e.g., `pure function` calling internal `complex_pfaffian` allocation), that flag becomes ticket 04's responsibility — even though it surfaces in ticket 01's commit. Coordinate: ticket 01 lands first; if it produces a Codacy flag, ticket 04 fixes it. Otherwise, ticket 04 only addresses the 2 pre-existing flags.
4. **`error stop` vs. `stop 1`.** If a flag points at `stop 1` (deprecated per CLAUDE.md), the fix is mechanical: change to `error stop '<descriptive message>'`. Trivial.

## Out of scope for this ticket

- Re-architecting modules to address deeper Code Quality or Style flags. This ticket addresses the 2 HIGH flags only; medium/low flags stay on the codebase and are not chased.
- Anything not on PR #42's diff. If a flag points at a file outside the diff, surface as "out of scope for this map" and `@todo` with a separate tracking link.

## Files to touch

TBD once Codacy's report is read. Most likely candidates:
- `src/physics/bdg_observables.f90` (if `pure` placement or allocatable leak)
- `src/physics/topological_analysis.f90` (if `stop 1` or intent mismatch)
- `tests/unit/*.pf` (if a pFUnit macro single-line violation slipped through `feedback_pfunit_macro_single_line` review)

## Constraints

- Per CLAUDE.md YAGNI: don't over-engineer the fix. If a 1-line change resolves the flag, that's the fix.
- Per `feedback_pfunit_macro_single_line`: any `@assertEqual`/`@assertTrue` MUST be single-line. If a flag is a pFUnit macro continuation error, fix it that way.

## Verification

- Codacy re-scan after the fix commit: 0 high issues (or low/minor at worst, per `design.md` §"Testing" item 6).
- `cmake --build build` clean.
- `ctest -L unit` PASS (no regressions from the fix).

## Interaction

- Independent of tickets 01-03 from a dependency-graph perspective, but practically coordinated with ticket 01: if ticket 01's `pure` placement produces a Codacy flag, ticket 04 picks it up.
- Tickets 05 (verification gate) and 07 (doc propagation) both depend on Codacy being clean.

## Resolution path (suggested)

1. Open PR #42 on GitHub.
2. Read Codacy's bot comment (or run Codacy CLI if available locally).
3. If the 2 flags match the predicted patterns (`pure` placement + allocatable leak), fix in 1-2 lines each.
4. If the flags are unexpected, run a focused research subagent to determine the right fix; graduate follow-up tickets for any non-trivial work.
5. Commit, push, re-trigger Codacy scan.