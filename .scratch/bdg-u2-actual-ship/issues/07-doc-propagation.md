**Status**: closed (2026-07-25) — ticket 07 of `.scratch/bdg-u2-actual-ship/`.

Type: task (AFK)
Status: closed
Claimed-by: —
Blocked by: 01, 03, 05, 06
Executed: 2026-07-25 (wayfinder work-through)

## Execution (2026-07-25)

Re-charted-then-executed after the premise-correction block above. Files touched:

| File | Change |
|---|---|
| `docs/plans/BACKLOG.md` | Preamble (L5-11): 52→51, "U2 closed"→"PR #42 open follow-up". Phase-26 detailed block (L670-682): SSOT type framing, orphan retirement, dense-test retirement, 3-witness gate revert, 51 count, HEAD `8e9128d`. Appended User Story 5 verbatim paragraph (L684). Phase-26 table row (L758): IN PROGRESS, 3-witness, 51, `bdg_pfaffian_params_t`. Trailer (L760): 51 + 6 seam-sibling tests + "PR #42 open". |
| `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` | U2 close-out block (L54-85): corrected framing (open PR #42, SSOT type, orphan retirement, 3-witness revert, 51 count) + embedded User Story 5 verbatim paragraph with meta-test citation. |
| `src/physics/AGENTS.md` | `bdg_observables.f90` row (L44): named `type :: bdg_pfaffian_params_t` + factory + User Story 5 contract fragment (range + non-zeroness + Issue-05 carve-out) byte-identical to spec. (DAG block already documented the L3-symbol trade-off; no DAG change needed beyond the row.) |
| `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` | "Problem"/"Solution"/"Why"/"Latent bug"/"Strict S1+S2"/"When to use"/"Source" sections reconciled: S2-only (not S1⊗S2), orphan retirement, `bdg_pfaffian_params_t` SSOT, 3-witness (4th reserved for U13), dense-variant line ranges marked retired, User Story 5 verbatim + meta-test + 51 count. |
| `.scratch/bdg-u2-actual-ship/HANDOFF.md` (new) | Live-state handoff: IN PROGRESS footer, PR #42 open, frontier, locked decisions table, standing preferences, files-to-read. |
| `.scratch/archive/bdg-evaluator-pfaffian/HANDOFF.md` (new) | Archived-state handoff: SHIPPED-via-PR-#42 footer (spec-of-record only; points to live dir for follow-up). |

**Correction principle applied (do NOT re-litigate):** HANDOFF files are *pointer* docs and do NOT embed the User Story 5 verbatim paragraph (DRY — a 7th verbatim copy would itself be a drift hazard; the U10 HANDOFF template they follow doesn't embed prose either). The "5 propagation sites" in SD-1 = the 5 doc files; the 2 HANDOFFs are SD-3/SD-4 status/footer docs that *point* at the spec + solutions doc. So the verbatim lives in exactly 6 places (spec + 5 doc sites) and is greppable via the range-contract fragment `s2 >= -1 .and. s2 <= 1`.

**Out-of-scope honored:** Historical Phase-25 row (PR #41 "1 expected unit fail") left as accurate history — doc-drift prevention fixes *current* state, not retroactive history. U10 HANDOFF locked-decisions table noted as having its own stale references (4-witness, dense `wire_pfaffian_witness`) but updating U10's handoff is U10 scope (this map's ticket-07 §"Out of scope"), not here.

## Premise correction (2026-07-25, before execution)

The frozen sub-decisions below were authored 2026-07-13 against then-current state. Since then three factual deltas landed that the frozen SDs would re-introduce as drift if executed verbatim. Locked corrections (do NOT re-litigate):

- **Unit count = 51, NOT 52.** `unit-count.txt` = 51; `ctest-final.log` footer = "100% tests passed, 0 tests failed out of 51". The frozen SD-2 ("replace 50/52/60+ disagreement with the canonical count") — the canonical count is **51**. BACKLOG rows 7/681/755/757 and parent-plan line 63 currently say "52/52" and must be overwritten to 51. The 52→51 drop was the ctest-target-vs-`@test`-subroutines count after deleting `test_wire_pfaffian_witness.pf` (map ticket-02 execution note; `feedback_ctest_counts_targets` memory).
- **Acceptance gate = 3-witness, NOT 4-witness live.** Commit 8e9128d (ticket-05 regression fix) reverted the 4th slim-Pfaffian row back to a 3-witness gate — the slim Pfaffian site-by-site sweep is an interim diagnostic reserved for the U13 Bloch-Pfaffian (map "Out of scope"; `project_bdg_allzero_mu_misparam` / handoff SD). The frozen SD-3/SD-4 framing ("4-witness design") must NOT be propagated to BACKLOG row 755 (currently says "4-witness acceptance gate (slim Pfaffian row live, colormap-extracted, ticket 05)") — correct it to "3-witness acceptance gate (slim-Pfaffian row reserved for U13)". The lecture-script comment strings that say "4-witness design; slim Pfaffian row reserved" are intentional design-aspiration markers and stay.
- **Status = PR #42 OPEN + follow-up branch, NOT "SHIPPED via PR #42".** PR #42 is still open (per `project_bdg_u2_ship_pr42`); `feat/bdg-u2-actual-ship` follow-up branch has not merged. The frozen SD-4 HANDOFF footer "SHIPPED via PR #42 (2026-07-13)" is itself drift as of today. The archived HANDOFF footer reads "spec-of-record shipped via PR #42 (2026-07-13); follow-up ticket work on `feat/bdg-u2-actual-ship` (PR #42 still open)".
- **AGENTS.md SD-5 is outstanding, not already-done.** `src/physics/AGENTS.md:39,44` cite the `bdg_default_*` constants and the seam-sibling delegations + the cross-level `wire_pfaffian_witness_sweep` trade-off, but do NOT name the `bdg_pfaffian_params_t` type. SD-5's "name `bdg_pfaffian_params_t` in inventory row + DAG" is real outstanding work. (Type exists at `bdg_observables.f90:36,62`.)

These corrections follow the locked-decision convention of this map (corrections recorded IN the ticket, not silently executed) and the `codebase-doc-drift-prevention` discipline.

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