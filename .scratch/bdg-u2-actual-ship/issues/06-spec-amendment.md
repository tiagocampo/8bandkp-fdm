**Status**: pending — ticket 06 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: pending
Claimed-by: —
Blocked by: 03
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