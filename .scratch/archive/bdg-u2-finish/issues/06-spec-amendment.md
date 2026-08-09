**Status**: closed (2026-07-25) — ticket 06 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed
Claimed-by: —
Blocked by: 03
Executed: 2026-07-25 (wayfinder work-through)

## Execution (2026-07-25)

- **SD-1 line range:** confirmed `spec.md:95-99` held the *old* "strict assertion to be GREEN" wording (the wording that predates the discovered invariant). Replaced by content (line numbers drifted by ±0 since the spec is in the archived, gitignored dir — located by content match on "5. As a researcher running cross-builder").
- **SD-2 verbatim preservation:** the new paragraph is byte-identical to `design.md` §"Spec Amendment" line 185 (the blockquote body), appended with the SD-3 citation.
- **SD-3 cross-reference:** appended `(test: \`test_pfaffian_witness_spec_user_story_5_contract\` in \`tests/unit/test_bdg_pfaffian_witness_csr.pf\`)`. Meta-test existence confirmed via `grep` — present at `test_bdg_pfaffian_witness_csr.pf:147` (landed by ticket 03).
- **SD-4 sibling stories:** stories 1–4 and 6–8 untouched; only the User Story 5 paragraph changed.

**Files touched:** `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (User Story 5 paragraph only).

**Note:** tracking is split. `.scratch/bdg-u2-actual-ship/` (the *live* dir — map, tickets, ctest-final.log) **is tracked in git** (added by commit 8e9128d). `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (the *archived* spec) is **untracked** (`??` in `git status`) — so this ticket's single source edit is not captured by `git diff` by design, and verification is by content read + design.md text comparison (per ticket §"Verification" bullet 2: "copy the new User Story 5 text, grep it against the diff" → here, grep the spec text against design.md §"Spec Amendment" blockquote body; the only diff is the leading `> ` marker, content byte-identical ✓). The map and ticket-06 execution-record edits *are* tracked and appear in `git status`.
Related: design.md §"Solution" item 12; §"Spec Amendment (User Story 5 verbatim text)"

# 06 — Spec amendment: User Story 5 verbatim text

## Question

What is the verbatim text that replaces `.scratch/archive/bdg-evaluator-pfaffian/spec.md` lines 95-99 (User Story 5) to encode the discovered invariant?

## Background

Per `design.md` §"Spec Amendment (User Story 5 verbatim text)":

> "Replace `.scratch/archive/bdg-evaluator-pfaffian/spec.md:95-99` with:
>
> > 5. As a researcher running cross-builder BdG unit tests, I want the Pfaffian sibling's witness assertions to be GREEN today (independent of the deferred Bloch-periodic BdG construction U13), where the witness is the S2-projected Pfaffian sign `s2 ∈ {-1, 0, +1}` (range contract pinned by `s2 >= -1 .and. s2 <= 1` and non-zeroness pinned by `s2 /= 0` on a non-diagonal synthetic fixture). Full S1+S2 strict sign-agreement is the Majorana-basis Pfaffian problem (Issue 05), deferred to a follow-up scope — for U2 the S1 path returns `(0, -1)` on the canonical fixtures, which is a structurally-valid gap-closure signal in the S1 eigenspace, not a defect."

The text is locked by the design.md; this ticket is **mechanical execution** (copy verbatim, paste into spec.md, commit).

## Sub-decisions to lock

1. **Line range confirmation.** Read `.scratch/archive/bdg-evaluator-pfaffian/spec.md` lines 95-99. Verify the existing User Story 5 text matches what design.md says to replace. If the line numbers have drifted (e.g., from prior commits), locate by content (`5\. As a researcher running cross-builder`).
2. **Verbatim preservation.** The text in `design.md` §"Spec Amendment" must NOT be paraphrased or "improved" during paste. Future agents grep for the exact wording in the 5 doc-propagation sites (ticket 07). Drift = backfill hazard (per `codebase-doc-drift-event-3.md`).
3. **Cross-reference additions.** The new User Story 5 should reference the meta-test name `test_pfaffian_witness_spec_user_story_5_contract` (added in ticket 03) as the executable form of the spec contract. Add a parenthetical at the end of User Story 5: "(test: `test_pfaffian_witness_spec_user_story_5_contract` in `tests/unit/test_bdg_pfaffian_witness_csr.pf`)".
4. **Sibling story updates.** User Stories 1, 2, 3, 4 in the same spec are unchanged. Verify by reading the surrounding lines before/after the replacement.

## Out of scope for this ticket

- Doc propagation to BACKLOG.md, parent plan, AGENTS.md, solutions doc, memory — ticket 07 (except for the verbatim wording record, which propagates from this ticket's diff).
- Spec changes outside User Story 5 — separate scope.

## Files to touch

| File | Change |
|---|---|
| `.scratch/archive/bdg-evaluator-pfaffian/spec.md` | Replace lines 95-99 (User Story 5) with the verbatim text from `design.md` §"Spec Amendment (User Story 5 verbatim text)" |

## Constraints

- Verbatim preservation per design.md.
- Per CLAUDE.md Known Issues, the spec lives in the *archived* dir (`archive/bdg-evaluator-pfaffian/`), not the live dir — because the live state of U2 is documented in this map (`.scratch/bdg-u2-actual-ship/`).

## Verification

- `git diff .scratch/archive/bdg-evaluator-pfaffian/spec.md` shows the User Story 5 paragraph changed; no other paragraph changed.
- Manual scan: copy the new User Story 5 text, grep it against the diff to verify byte-identical.
- The new User Story 5 references the meta-test name from ticket 03.

## Interaction

- **Ticket 03** is blocked-by: this ticket's spec amendment cites the meta-test name `test_pfaffian_witness_spec_user_story_5_contract` which doesn't exist until ticket 03 lands. Strict ordering.
- **Ticket 07** propagates the new User Story 5 wording to 5 doc sites (BACKLOG row 82, BACKLOG L11 summary, parent plan status_note, AGENTS.md, solutions doc).
- **Ticket 08** mirrors the wording in the memory entry.