# BdG evaluator seam SSOT — archived spec-of-record — Wayfinder handoff

> **Read this first** (and then go to the live dir). This is the
> *archived* handoff for `.scratch/archive/bdg-evaluator-pfaffian/` — the
> 7-ticket original map + the spec-of-record that landed via PR #42
> (2026-07-13). The spec-of-record shipped; the **follow-up ticket work
> (codacy flags + PR-#42 review findings + spec amendment + doc propagation)
> lives on the live dir `.scratch/bdg-u2-actual-ship/`** — Read its
> `HANDOFF.md` next for current state and the open frontier.

## Status footer

**SHIPPED via PR #42 (2026-07-13) — spec-of-record only.** 7 tickets resolved
(original map: evaluator-API shape, dense-witness retirement, sign-agreement
test scope, BHZ-only heuristic rename, slim-Pfaffian band-major projection,
wire-polarization emitter/verifier, test-coverage audit). The archive moved
from `.scratch/bdg-evaluator-pfaffian/` to `.scratch/archive/bdg-evaluator-pfaffian/`
after the original map closed. **Follow-up ticket work (issues 01–08 in the
live dir) is still in progress on `feat/bdg-u2-actual-ship`** — PR #42 is
**open**, not merged.

## What shipped here

- A 3-face seam on `bdg_observables.f90`: `eval_bdg_point`,
  `eval_bdg_pfaffian_witness_csr` (S2-projected Pfaffian sign),
  `eval_bdg_kitaev_majorana` (wraps `kitaev_majorana_number`).
- The spec-of-record with 8 user stories (User Story 5 carries the
  S2-only invariant contract — verbatim, amendment from the live dir's
  `design.md` §"Spec Amendment").
- 7 original locked decisions — see `map.md` + `issues/01-…07-*.md`.

## What shipped in the *follow-up* (live dir, not here)

- SSOT `bdg_pfaffian_params_t` type + factory (ticket 01).
- Orphan dense `wire_pfaffian_witness` + helpers retired (ticket 02).
- Test contracts tightened (`s2 ∈ {-1, 0, +1}`, meta-test, factory) (ticket 03).
- Codacy F541 flags fixed (ticket 04).
- 3-witness gate verification GREEN; 51/51 unit PASS (tickets 05).
- User Story 5 amended into this archived `spec.md:95` (ticket 06).
- 5 doc propagation sites + 2 HANDOFF files (this one + the live one) (ticket 07).

## Out of scope for this archive (live-dir or onward-scope)

- **U10** phase-diagram — separate handoff at `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`.
- **Issue 05** (Majorana-basis strict S1==S2) — deferred to U13.
- **U13** (periodic/Bloch BdG construction) — CLAUDE.md Known Issues; carries the 4th slim-Pfaffian gate row.

## Standing preferences (the live-dir handoff carries the live reminders)

See `.scratch/bdg-u2-actual-ship/HANDOFF.md` §"Standing preferences" for the
current convention reminders (pFUnit single-line macros, no Co-Authored-By
trailer, error stop not stop 1, ctest-target-vs-`@test`-count gotcha,
doc-drift verbatim propagation). The archived versions are superseded by
the live ones for active work.

## Files to Read first

1. **`.scratch/bdg-u2-actual-ship/HANDOFF.md`** — the live follow-up state (read this next, not this file, for current work).
2. `.scratch/archive/bdg-evaluator-pfaffian/spec.md` — the spec-of-record (User Story 5 at line 95).
3. `.scratch/archive/bdg-evaluator-pfaffian/map.md` — the original 7-ticket map.
4. `.scratch/archive/bdg-evaluator-pfaffian/issues/01-…07-*.md` — the original locked decisions.
5. `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` — durable pattern record.

## Source

- Live follow-up: `.scratch/bdg-u2-actual-ship/` (map + 8 tickets + design.md).
- Spec-of-record: `.scratch/archive/bdg-evaluator-pfaffian/spec.md`.
- Solution doc: `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`.
- Next-scope sibling: `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`.
- Memory entry (live): `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`.
