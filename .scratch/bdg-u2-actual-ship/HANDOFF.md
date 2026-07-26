# U2 actual-ship — follow-up PR #42 review findings + Codacy — Wayfinder handoff

> **Read this first.** Future `/wayfinder` session entry point for the U2
> finish follow-up on `feat/bdg-u2-actual-ship`. This is the *live* scratch;
> the **spec-of-record** lives archived at
> `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (the User Story 5 verbatim
> text + the 7 original tickets). This handoff names the destination, fixes
> the state, marks the loose, and points at the standing preferences so the
> next chart starts grounded. The archived handoff
> (`.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md`) points back here.

## Status footer

**IN PROGRESS (2026-07-25).** U2 finish design follow-up on
`feat/bdg-u2-actual-ship` — closing PR #42 review findings (9 Standards + 3
Spec, two-axis code review) + the 2 Codacy high `ErrorProne` flags (pyflakes
F541, both fixed in `8e9128d`). Tests green: `ctest -L unit` 51/51 PASS,
`lecture_13_acceptance_gate` GREEN (3-witness gate). PR #42 is **open**;
the branch has not merged.

## Destination

Land the follow-up PR #42 close-out: every PR-#42 review finding closed +
Codacy flags fixed + `User Story 5` amended in the archived spec + the 5
doc-propagation sites carry the verbatim wording + the 2 HANDOFF files
exist + the memory entry refreshed. Single-FSM extension of PR #42 (no
squash, multi-commit TDD), then umbrella-merge `feat/bdg-u2-actual-ship`.

A "shipped" U2 finish looks like:

- Tickets 01–07 executed (SSOT wire-up, orphan retirement, test contracts,
  Codacy triage, verification gate, spec amendment, doc propagation) and
  ticket 08 (memory) closed.
- `bdg_observables.f90` seam: `eval_bdg_pfaffian_witness_csr` retyped to
  `bdg_pfaffian_params_t`; `bdg_pfaffian_params_with_floor` factory carries
  the `error stop`; orphan dense `wire_pfaffian_witness` + helpers removed.
- `ctest -L unit` 51/51 PASS (the 52→51 drop is the ctest-target count after
  retiring `test_wire_pfaffian_witness.pf`); `cmake --build build` clean.
- 3-witness lecture-13 acceptance gate GREEN (`WITNESS_LABEL="3-witness"`);
  the 4-witness slim-Pfaffian row regression reverted in `8e9128d`; that row
  is reserved for U13 (Bloch-Pfaffian).
- 5 doc sites (BACKLOG preamble + Phase-26 row + parent-plan U2 block +
  AGENTS.md + solutions doc) carry the User Story 5 verbatim wording and
  the 51 unit count; no 52, no "4-witness live", no "SHIPPED via PR #42".

## State of the tree (what already shipped on `feat/bdg-u2-actual-ship`)

These are the locked decisions this map **must not re-litigate**. If a
ticket re-derives one, close it and point to the source.

| Locked decision | Ticket / file | Commits / evidence |
|---|---|---|
| SSOT `type :: bdg_pfaffian_params_t` + factory `bdg_pfaffian_params_with_floor` (`error stop` on `pfaffian_floor <= 0`); `bdg_eval_params_t%pfaffian_floor` field dropped | `issues/01-ssot-wire-up.md` | `d908b56` |
| Orphan dense `wire_pfaffian_witness` + `s1_project`/`s2_construct` retired; CSR sweep `wire_pfaffian_witness_sweep` retained | `issues/02-orphan-retirement.md` | `d908b56` |
| `test_bdg_pfaffian_witness_csr.pf` carries 6 `@test`: delegation + zero + range + nondiagonal (`@assertEqual(-1,s2_sign)`) + floor-consumption + User Story 5 meta-test; SD-1 (factory error-stop `@test`) dropped — pFUnit 4.16 can't catch `error stop` | `issues/03-test-contracts.md` | `d908b56` |
| 2 Codacy F541 flags = drop `f`-prefix on `lecture_13_topological.py:396` + `verify_majorana_polarization.py`; both applied | `issues/04-codacy-triage.md` | `8e9128d` |
| Verification gate: 51/51 unit PASS, 3-witness acceptance-gate GREEN, `unit-count.txt = 51` | `issues/05-verification-gate.md` | `8e9128d` (gate regression fix + close) |
| User Story 5 verbatim (from `design.md` §"Spec Amendment") amended into `archive/bdg-evaluator-pfaffian/spec.md:95` + meta-test citation | `issues/06-spec-amendment.md` | this session (untracked spec) |

Frontier (2026-07-25): **ticket 07** (doc propagation — in progress) then
**ticket 08** (memory entry refresh, blocked by 07).

## Loose (chart as tickets only if a gap survives)

- **PR #42 merge.** Not a code ticket — review sign-off + umbrella-merge
  of `feat/bdg-u2-actual-ship` once tickets 07–08 land. Track outside this
  map (PR timeline) unless a fresh review finding reopens.
- **U13 slim-Pfaffian → Bloch-Pfaffian.** Reserved scope (CLAUDE.md Known
  Issues; map "Out of scope"). The 3-witness gate stays until U13 lands;
  the slim sweep is the interim diagnostic only.

## Standing preferences (every session on this map should consult)

- **Skills.** `/superpowers:test-driven-development` (the seam-sibling
  tests are TDD-doubled); `/superpowers:verification-before-completion`
  before closing any ticket (51/51 unit + gate GREEN); `/mattpocock-skills:wayfinder`
  to claim/work the map.
- **Convention reminders.** pFUnit `@assertEqual`/`@assertTrue` single-line
  (no `&` continuation — `feedback_pfunit_macro_single_line`); `error stop
  '<descriptive message>'` not `stop 1`; no Co-Authored-By trailer on commits
  (`feedback_no_coauthor_trailer`); ctest counts *targets* not `@test`
  subroutines — deleting a `.pf` shifts the total by the target, not the
  `@test` count (`feedback_ctest_counts_targets`).
- **Doc-drift discipline.** `codebase-doc-drift-prevention` + `codebase-doc-drift-event-3`:
  every behavior-touching change updates the spec + plan-status footer +
  BACKLOG-row + AGENTS-row + solutions-doc + memory, *verbatim* across all
  5 sites; drift is the rework root cause (PR #27, U8 reviews, PR #41 audit).
- **Engineering principles.** Repo root `CLAUDE.md` → "Engineering
  Principles": SSOT (`bdg_pfaffian_params_t`), YAGNI (no speculative config
  knob), KISS (slim Pfaffian is an interim diagnostic, not a mandatory
  witness — reverts are cheap).

## Out of scope (rule-of-scope decisions for this map)

- **U9** (BdG spectral + LDOS), **U10** (bulk minigap + (B, μ) phase
  diagram), **U11** (lecture 13 full revamp), **U13** (periodic/Bloch BdG
  construction) — carried forward (`design.md` §"Out of Scope"; U10 has its
  own handoff at `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`).
- **Issue 05** (Majorana-basis Pfaffian projection yielding strict S1==S2)
  — separate scope, deferred; User Story 5 encodes the carve-out.
- **`topological_analysis.f90` module split** (~1640 lines) — BACKLOG
  Phase 18 architecture-deepening item.
- **`pure` attribute on `src/math/pfaffian.f90` helpers** — ticket 01
  SD-3 defers (different module).
- **Retroactive editing of historical PR rows** (e.g., BACKLOG Phase-25 row
  describing PR #41's "1 expected unit fail" — true at PR #41 merge) —
  out of scope; doc-drift prevention fixes *current* state, not history.

## Files to Read first

When starting a Wayfinder chart session on this map, Read in order:

1. **This file** (`.scratch/bdg-u2-actual-ship/HANDOFF.md`) — destination + state.
2. `.scratch/bdg-u2-actual-ship/map.md` — the 8 tickets + frontier.
3. `.scratch/archive/bdg-evaluator-pfaffian/spec.md` — the spec-of-record (User Story 5 verbatim).
4. `.scratch/bdg-u2-actual-ship/design.md` — the 14-step as-built retrospective.
5. `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md` — the recover-and-execute plan.
6. `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` — durable pattern record.
7. `src/physics/AGENTS.md` §"Dependency DAG" + the `bdg_observables.f90` row.

## Source

- Parent plan: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U2.
- Original U2 plan: `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md`.
- Status footer: PR #42 (open) on `feat/bdg-u2-actual-ship` HEAD `8e9128d`;
  `ctest -L unit` 51/51 PASS; 3-witness acceptance-gate GREEN.
- Archived companion: `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md`
  (SHIPPED-via-PR-#42 spec-of-record, points back here for live follow-up state).
