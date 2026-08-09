# P6 — Pick-D scope gaps (T1a/T1b closure markers + BACKLOG row 26 contradiction) (IMPORTANT)

## Question (task, AFK)

The Pick-D commit `1a7767e` (status-stamps for BACKLOG.md / REVIEW.md / parent plan §U10 footer) covered T2, T3, T4 BLOCKED, T4r REFUTED, T5, T6 — but missed two items:

1. **T1a/T1b closure markers absent from parent plan §U10 footer** — `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md:467-483` lists T2 through T6 but does not include T1a (commit `1ee9ba2`, closed 2026-07-26) or T1b (commit `0a62143`, closed 2026-07-26). T1a/T1b were the foundational enablers — without them, T2/T3/T5/T6 wouldn't have a seam sibling or per-B proxy producer. Per `codebase-doc-drift-prevention`, every behavior PR must update plan status footers.

2. **BACKLOG.md row 26 internal contradiction** — `docs/plans/BACKLOG.md:759` has:
   - Row 26 cell: "PR #42 MERGED 2026-07-26 (squash `257b3c4`) → U2 closed"
   - Row 26 description text: "Branch `feat/bdg-u2-actual-ship` follow-up to PR #42 (pending merge), HEAD `4986d5a`"
   - Post-line 761: "U2 close-out on `feat/bdg-u2-actual-ship` (PR #42 open)"
   
   The cell is correct (verified via `git log main --oneline`); the description and post-line are stale.

## Acceptance

- [ ] Add a T1a block to `2026-06-14-001-feat-bdg-majorana-validation-plan.md` parent plan §U10 footer (before T2):
  ```
  **T1a (2026-07-26) closed (commit `1ee9ba2`):** optional `best_pf_abs` out-arg threaded through `wire_pfaffian_witness_sweep` and `eval_bdg_pfaffian_witness_csr`. 51/51 unit green.
  **T1b (2026-07-26) closed (commit `0a62143`):** per-B `min-|Pf|` proxy producer wired through `eval_wire_bdg_gap` → `compute_wire_bdg_gap_sweep` → `run_gap_sweep` → `write_wire_slim_pfaffian_witness`; new `regression_wire_slim_pfaffian_witness` ctest entry. 51/51 unit + 6 BdG/Pfaffian tests green.
  ```
- [ ] Update `BACKLOG.md:759` row 26 description: remove "(pending merge)" or replace with "(merged 2026-07-26 squash `257b3c4`)".
- [ ] Update `BACKLOG.md:761` post-line: drop "PR #42 open" → replace with "Phase 26 (U2 close-out) complete; PR #42 merged 2026-07-26".

## Cross-references

- T1a ticket: `tickets/T1a_pfaffian_magnitude_seam.md`
- T1b ticket: `tickets/T1b_per_b_min_pf_proxy.md`
- Memory: `project_bdg_u10_t5_execution.md` (Pick-D split context)
- Memory: `project_docs_overhaul.md` (codebase-doc-drift-prevention rule)

## Blocks / blocked-by

- **Blocks**: nothing.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
