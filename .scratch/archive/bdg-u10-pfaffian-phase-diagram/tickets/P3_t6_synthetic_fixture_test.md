# P3 — T6 saturation branch not covered by automated test (IMPORTANT, defense-in-depth)

## Question (task, AFK)

T6 closed 2026-07-27 with a memory claim: "Three failure-mode branches verified end-to-end via synthetic file injection." Inspection of `scripts/lecture_13_topological.py:230-258` shows the branches fire on real data, but no automated test pins them. The acceptance gate (`test_lecture_13_acceptance_gate.sh`) exercises the canonical path (file saturated, all rows `4.0E-08`), but it doesn't assert the specific WARN message or `bcrit_pfaffian = None` return — only that the gate passes overall.

If a future maintainer silently breaks the degenerate branch (e.g., changes the threshold from `1e-6` to `0.5`), the canonical run would no longer be "saturated" relative to the threshold, and the gate would emit a numeric `bcrit_pfaffian` — a regression the gate wouldn't catch because the range check still passes.

## Acceptance

- [ ] Lift `section_wire_rung(exe, workdir=None)` into a parameterizable function (currently called only from `main()` at `lecture_13_topological.py:415` with no module-level pytest harness).
- [ ] Add three synthetic-fixture pytest cases (or shell-test cases that drive the python via `python3 -c`):
  - **Case A (file absent):** `WORKDIR` has no `wire_slim_pfaffian_witness.dat`. Assert `bcrit_pfaffian is None` and the "approximation (open-chain projected; full Bloch-Pfaffian deferred U13)" label is emitted.
  - **Case B (file saturated):** `WORKDIR/wire_slim_pfaffian_witness.dat` has all rows at the same magnitude (e.g., `B=0 |Pf|=4.0E-08`, `B=1 |Pf|=4.0E-08`, ...). Assert `bcrit_pfaffian is None` and the degenerate WARN is emitted.
  - **Case C (file varying):** `WORKDIR/wire_slim_pfaffian_witness.dat` has a non-degenerate profile (e.g., monotonically increasing `|Pf|` from `1e-10` to `1e-3`). Assert `bcrit_pfaffian` is a numeric value (the argmin B).
  - **Case D (file all-zero):** added per P2 — all rows `B=0 |Pf|=0.0`. Assert `bcrit_pfaffian is None`.
- [ ] Either as a new pytest file `tests/integration/test_lecture_13_synthetic_fixtures.py` or as shell + python-heredoc pattern in a new `tests/integration/test_lecture_13_degenerate_branches.sh`.

## Why this matters

T6 was a manual verification, not a regression net. The branch is shipping with the saturated regime as its canonical behavior — without an automated pin, future changes to the threshold or the predicate structure can silently flip the gate's behavior.

## Cross-references

- P1: `tickets/P1_slim_pf_parser_format.md` (related: parser + writer contract)
- P2: `tickets/P2_pmax_zero_short_circuit.md` (must land first; P3 tests the all-zero case D)
- T6 ticket: `tickets/T6_lutchyn_oreg_proxy.md` (the original closure)

## Blocks / blocked-by

- **Blocks**: nothing.
- **Blocked by**: P2 (all-zero branch fix must land first).

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
