**Status**: pending — ticket 07 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: pending
Claimed-by: —
Blocked by: 01, 03, 05, 06
Related: design.md §"Solution" items 10, 11, 13; §"Components" rows 9-13; §"Cross-References"

# 07 — Doc propagation: HANDOFF files + 5 doc sites

## Question

After tickets 01 (SSOT wire-up), 03 (test contracts), 05 (verification gate), and 06 (spec amendment) land, what doc surfaces need to be touched to (a) preserve the doc-drift-prevention invariant (verbatim wording in 5 propagation sites), and (b) create HANDOFF.md files that disambiguate the live `.scratch/bdg-u2-actual-ship/` dir from the archived `.scratch/archive/bdg-evaluator-pfaffian/` dir?

## Background

Per `design.md` §"Solution" items 10, 11, 13 and the previous map's `codebase-doc-drift-event-3.md` lesson: doc-drift is the rework root cause. Five propagation sites must carry the new User Story 5 wording verbatim:

1. `docs/plans/BACKLOG.md` row 82 (Z2-gap heuristic row)
2. `docs/plans/BACKLOG.md` L11 summary
3. `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` U2 close-out block (lines 54+)
4. `src/physics/AGENTS.md` `bdg_observables.f90` inventory row (line 44)
5. `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` §"Strict S1+S2 sign agreement" section

Plus two new HANDOFF.md files:
- `.scratch/bdg-u2-actual-ship/HANDOFF.md` (live-state, status footer "IN PROGRESS")
- `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` (archived-state, status footer "SHIPPED via PR #42"; points future agents to the live dir)

Plus the unit-count propagation: ticket 05 wrote `unit-count.txt` with the canonical count; this ticket propagates the count to the same 5 doc sites.

## Sub-decisions to lock

1. **Verbatim wording source.** Read the new User Story 5 text from `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (per ticket 06's diff). Propagate byte-identical text to all 5 sites.
2. **Unit-count propagation source.** Read `.scratch/bdg-u2-actual-ship/unit-count.txt` (per ticket 05). Replace any existing count reference (currently 50 / 52 / 60+ disagreement per design.md §"Problem Statement") with the canonical count.
3. **HANDOFF.md live state.** `.scratch/bdg-u2-actual-ship/HANDOFF.md` documents:
   - Status footer: "IN PROGRESS (2026-07-13). U2 finish design work — closing PR #42 review findings + Codacy flags."
   - Destination: PR #42 review findings closed.
   - 8 tickets enumerated (link to map.md).
   - Frontier: ticket 01 (lowest-numbered open).
4. **HANDOFF.md archived state.** `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` documents:
   - Status footer: "SHIPPED via PR #42 (2026-07-13). 7 tickets resolved; archive moved from `.scratch/bdg-evaluator-pfaffian/` to `.scratch/archive/bdg-evaluator-pfaffian/` after ticket 07 close-out."
   - Next destination (U10 phase diagram): see `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`.
   - This is the same pattern the existing `bdg-u10-pfaffian-phase-diagram/HANDOFF.md` follows (read first by future sessions).
5. **AGENTS.md DAG update.** The DAG block at `src/physics/AGENTS.md:46-58` (per previous map's ticket 02 §"Premise correction") needs to show `bdg_observables.f90:bdg_pfaffian_params_t` and the new `pure` placement on the seam siblings. The inventory row (line 44) cites the new type.

## Out of scope for this ticket

- The spec.md User Story 5 amendment — ticket 06 (this ticket propagates the wording, doesn't author it).
- Memory entry — ticket 08.
- BACKLOG.md row content beyond the verbatim User Story 5 wording + unit-count propagation (other BACKLOG rows may have tangential updates, but those are out of scope here).

## Files to touch

| File | Change |
|---|---|
| `docs/plans/BACKLOG.md` | Row 82 (Z2-gap heuristic) + L11 summary: replace any User Story 5 wording with the verbatim new text; replace any unit-count reference with `unit-count.txt` value |
| `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` | U2 close-out block (lines 54+): replace User Story 5 wording with verbatim new text; update unit count |
| `src/physics/AGENTS.md` | `bdg_observables.f90` inventory row (line 44): cite `bdg_pfaffian_params_t`; update DAG block at lines 46-58 |
| `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` | §"Strict S1+S2 sign agreement" section: reference Issue 05 / U13 + the new User Story 5 wording; update unit count |
| `.scratch/bdg-u2-actual-ship/HANDOFF.md` (new) | Live-state doc per the format in `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md` |
| `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` (new) | Archived-state doc per the same format |

## Constraints

- Verbatim preservation per `codebase-doc-drift-prevention.md` + `codebase-doc-drift-event-3.md`.
- The 5 doc sites must cite the SAME wording — drift surfaces as diff and is caught by the next doc-drift audit.
- The HANDOFF.md files must follow the format established by `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md` (Destination, State of the tree, Loose, Out of scope, Standing preferences, Files to Read first, Source).

## Verification

- `grep -l '<new User Story 5 wording fragment>' docs/ src/physics/ .scratch/` returns 5 doc sites + 1 spec.md + 1 HANDOFF.md (live) + 1 HANDOFF.md (archived) = 8 hits.
- `cat .scratch/bdg-u2-actual-ship/unit-count.txt` and `grep -E 'unit[- ]count|Total Test' docs/plans/BACKLOG.md docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` show consistent numbers.
- Manual scan: each of the 5 doc sites has the new User Story 5 wording verbatim, no paraphrase.

## Interaction

- **Tickets 01, 03, 05, 06** all blocked-by: this ticket cites the new SSOT (ticket 01), the meta-test name (ticket 03), the unit count (ticket 05), and the spec wording (ticket 06). All four must be in place first.
- **Ticket 08** mirrors the wording in the memory entry.

## Resolution path

1. Read the new User Story 5 from `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (per ticket 06).
2. Read `.scratch/bdg-u2-actual-ship/unit-count.txt` (per ticket 05).
3. Edit each of the 5 doc sites: replace existing User Story 5 wording (if any) with verbatim new text; replace unit-count references with the canonical count.
4. Create `.scratch/bdg-u2-actual-ship/HANDOFF.md` per the live-state format.
5. Create `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` per the archived-state format.
6. Run the verification greps to confirm consistency.
7. Commit all 7 file changes in one doc-propagation commit.