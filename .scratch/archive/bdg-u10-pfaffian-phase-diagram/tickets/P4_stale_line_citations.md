# P4 — Stale line-number citations in test comments (IMPORTANT)

## Question (task, AFK)

Three pre-T1b line-number citations are stale post-T1b's seam plumbing (which added ~37 lines to `main_topology.f90`):

1. **`tests/integration/test_wire_bdg_topological_phase.sh:87`** — references `main_topology.f90:1420-1423` for the pre-T3 remap. Actual pre-T3 lines were `1381-1383` (verified via `git show main:src/apps/main_topology.f90`). Post-T1b, the line numbers shifted; the current `1420-1423` is the T1b seam plumbing region with no remap.

2. **`tests/integration/test_wire_slim_pfaffian_witness.sh:33`** — references `lecture_13_topological.py:226` for the regex. Actual regex at `lecture_13_topological.py:240`. Off by 14 lines (T6 added the `_PFAFFIAN_DEGENERACY_TOL` constant + the file-present-but-all-at-floor branch above the existing regex).

3. **`tests/regression/configs/wire_inas_gaas_bdg_topological_phase.toml:8`** — header comment says "widens the B grid to nB=21 for finer B_crit transition detection". Actual config has `gap_sweep_nB = 6` (line 65). The parent plan §U10 status footer correctly describes this as "B-grid (nB=6)".

## Acceptance

- [ ] Update `test_wire_bdg_topological_phase.sh:87` to either reference the pre-T3 function name (`eval_wire_bdg_gap`) or the verified pre-T3 line range (`main_topology.f90:1381-1383`).
- [ ] Update `test_wire_slim_pfaffian_witness.sh:33` to reference `lecture_13_topological.py:240`.
- [ ] Update `wire_inas_gaas_bdg_topological_phase.toml:8` to read "widens the B grid to nB=6 (vs canonical nB=5)" or similar.

## Cross-references

- T1b ticket: `tickets/T1b_per_b_min_pf_proxy.md` (the line-shift event)
- T6 ticket: `tickets/T6_lutchyn_oreg_proxy.md` (the line-shift event for lecture script)
- T2 ticket: `tickets/T2_mu_window.md` (the fixture authorship)

## Blocks / blocked-by

- **Blocks**: nothing.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
