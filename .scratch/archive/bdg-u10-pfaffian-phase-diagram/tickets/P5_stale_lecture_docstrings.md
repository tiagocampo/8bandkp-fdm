# P5 — Stale docstrings in `lecture_13_topological.py` (T5 scope gap) (IMPORTANT)

## Question (task, AFK)

T5 commit `d17f75b` claimed to be a "full cleanup" but missed several docstring updates:

1. **`scripts/lecture_13_topological.py:25-29`** (top-of-file docstring) still says "3 active B_crit witnesses (wire_curve, wire_2d, qw_dense)" + "slim Pfaffian row reserved for U13 and excluded from the numeric range until that lands". Contradicted by the actual auto-detect logic: gate runs as 3-witness when slim Pfaffian row is missing/truncated, extends to 4-witness when emitted as a numeric value (T1b producer at `output/wire_slim_pfaffian_witness.dat`).

2. **`scripts/lecture_13_topological.py:33-35`** (same docstring) — same stale phrasing.

3. **`scripts/lecture_13_topological.py:336-341`** (`render_reconciliation_table` docstring) still says "(deferred to U13) since no Fortran producer emits" / "reserved for U13". T1b landed the producer; U10 T6 closed the lecture script's consumption.

## Acceptance

- [ ] Update `lecture_13_topological.py:25-29` to describe the actual auto-detect logic: "3 active B_crit witnesses by default (wire_curve, wire_2d, qw_dense); extends to 4-witness when `output/wire_slim_pfaffian_witness.dat` emits a numeric argmin (T1b proxy stand-in; full Bloch-Pfaffian deferred U13)".
- [ ] Update `:33-35` similarly.
- [ ] Update `render_reconciliation_table` docstring at `:336-341` to remove "reserved for U13" wording.
- [ ] Apply the Q4 approximation phrase canonical from T5 to all three sites.

## Cross-references

- T5 ticket: `tickets/T5_gate_precondition.md` (the closure; "full cleanup" scope)
- Memory: `project_bdg_u10_t5_execution.md` (Q4 approximation phrase canonical)
- T6 ticket: `tickets/T6_lutchyn_oreg_proxy.md` (consumption side)

## Blocks / blocked-by

- **Blocks**: nothing.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
