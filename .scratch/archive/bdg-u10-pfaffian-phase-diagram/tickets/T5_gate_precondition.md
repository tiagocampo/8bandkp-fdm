# T5 — Gate precondition update z2==-1 native (acceptance surfaces)

## Question (AFK, task-ticket)

The acceptance-gate precondition `gate_row_colormap_present` (`verify_majorana_polarization.py:47-76`) reads `z2 == -1` from `output/z2_phase_diagram.dat`. Post-T3, the column is Pfaffian-native (`-1/+1/0`). The Python code already reads `int(float(parts[2]))` and checks `z2 == -1` — no change needed to the read logic; it should naturally match now.

But three downstream artifacts assume the OLD `z2 ∈ {0, 1}` convention:
- `lecture_13_topological.py:178-196` reads `z2` column looking for the 1D bcrit alignment; its `bcrit_2d` extraction doesn't filter on z2 value, so safe.
- `lecture_13_topological.py:222-256` `bcrit_pfaffian` extraction reads `output/wire_slim_pfaffian_witness.dat` (per T1). Need to update the "deferred to U13" status string to "approximation (open-chain projected; full Bloch-Pfaffian deferred U13)" per Q4 acceptance.
- `tests/integration/test_wire_bdg_topological_2d.sh` and `test_topology_sweep.sh` — check `awk` for row count and non-negative gap; no z2 filter, safe.

Open: is there a `tests/integration/verify_*.py` script that *currently* filters `z2 == 1` (old topological)? If yes, update filter to `z2 == -1` (new native). Search during execution.

## Acceptance

- [ ] `verify_majorana_polarization.py:74` `gate_row_colormap_present` reads `z2 == -1` (no change — verified).
- [ ] `lecture_13_topological.py:246` status string updated for approximation note.
- [ ] Any other `verify_*.py` or `lecture_*.py` that filtered old convention updated.
- [ ] Lecture 13 status table regenerated from machine-readable output (deferred to U11 — track via U11 dependency, not in this ticket).

## Blocks / blocked-by

- **Blocks**: U11 (lecture revamp will use the new schema).
- **Blocked by**: T3 (native schema) + T1 (proxy file). **Both closed** — T5 is now unblocked and claimable.
- **Unaffected by T4 closure** (2026-07-27): T5's gate logic reads `z2 == -1` natively (post-T3 column-write semantics) and does not depend on the colormap being non-flat. T5's work is independent of T4.

## Type
task (AFK).

## Branch
`feat/bdg-u10-pfaffian-phase-diagram`.

## Claimed by
Tiago de Campos (2026-07-27 — fresh session, post-T4 BLOCKED handoff).

## Resolution (2026-07-27, full cleanup per scope resolution)

T5 closed with **two commits** on `feat/bdg-u10-pfaffian-phase-diagram`:

1. `1a7767e` — `docs(bdg): U10 T4 BLOCKED + T5 active status footers (Pick-D bookkeeping)`: BACKLOG.md Phase 26 row + parent plan §U10 footer + REVIEW.md row 79.
2. `d17f75b` — `feat(bdg): U10 T5 — gate precondition accepts z2==-1 native (full cleanup)`: lecture-13 status strings + acceptance-gate comment blocks + slim-Pfaffian-projection doc-drift.

### Acceptance items (all closed)

- [x] `verify_majorana_polarization.py:74` `gate_row_colormap_present` reads `z2 == -1` (no change — verified T3 standing decision; gate works natively post-T3).
- [x] `lecture_13_topological.py` lines 238-246 WARN + line 255 BCRIT print + line 327 reconciliation table cell: status strings updated to Q4 approximation phrase ("approximation (open-chain projected; full Bloch-Pfaffian deferred U13)"). Outdated "no Fortran producer emits" / "reserved for U13" / "ce-doc-review P0 finding" comments collapsed.
- [x] Searched `tests/integration/verify_*.py` and `lecture_*.py` for stale `z2 == 1` (old topological) filter: only `test_edge_states.pf:80`, `test_z2_invariant.pf:74`, `test_phase_diagram.pf:178` match — these pin BHZ to `{0,1}` intentionally per T3 standing decision (would cause 7+ unit-test churn to migrate; out of scope for T5).
- [x] Lecture 13 status table regen deferred to U11 (per ticket body).

### Adjacent scope (per resolution C-full cleanup)

- `tests/integration/test_lecture_13_acceptance_gate.sh` comment blocks (49-56, 79-86, 95-111) + status label in deferred branch: rewritten to describe live auto-detect behavior (line 81 regex test unchanged). Logic preserved.
- `tests/integration/test_slim_pfaffian_witness_projection.py` lines 15, 43-44, 107-109, 167, 178-179: stale "RESERVED FOR U13" / "producer not yet in place" / "current code path does not emit the file" rewritten to reflect T1b emission state. Test logic + deferred-else fallback branch preserved (per scope resolution δ).

### Verification

- 51/51 unit green (33.9 s)
- 4/4 wire_bdg regression green (578 s) — T1b + T2 + T3 baseline preserved
- 2/2 slim_pfaffian regression green (317 s) — T1b emission intact
- `python3 -c "import lecture_13_topological"` → OK (syntax clean)

### Out of scope (locked by this resolution)

- BHZ unit-test churn (`test_edge_states.pf:80`, `test_z2_invariant.pf:74`, `test_phase_diagram.pf:178` pinning BHZ to `{0,1}`) — out-of-scope per T3 standing decision.
- L13 reconciliation-table 4-witness label sweep beyond line 327 — T11 dependency, not T5.
- `bcrit_pfaffian` cross-check logic vs. `bcrit_curve` / `bcrit_2d` (±0.5 T tolerance) — **T6** ticket.
- `test_slim_pfaffian_witness_projection.py` parse_slim_pf_witness logic change — test already in TDD-green state (T1b baseline); test logic preserved per δ.
