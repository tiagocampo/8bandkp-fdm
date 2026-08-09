# P2 — `lecture_13_topological.py` all-zero witness short-circuit (CRITICAL)

## Question (task, AFK)

`scripts/lecture_13_topological.py:245` has the degeneracy check:

```python
if pmax > 0 and (pmax - pmin) / pmax < _PFAFFIAN_DEGENERACY_TOL:
    bcrit_pfaffian = None
    print(f"WARN: ...")
else:
    bcrit_pfaffian = min(pf_mags, key=lambda B: pf_mags[B])
```

When `pf_mags` contains only zero entries (all-closure regime, every `s2_sign = 0` so `best_pf_abs = 0`), `pmax == 0` short-circuits to False, falling into the `else` branch. `min(pf_mags, key=lambda B: pf_mags[B])` returns the **first B key** (Python dict preserves insertion order) — not a phase boundary, but a meaningless first-occurrence pick.

The T6 docs claim "Three failure modes verified end-to-end (file absent / file saturated / file varying)". The **all-zero case** (closure-only regime) is NOT in that list. The current gate would treat the closure regime as a valid B_crit witness and add it to the 4-witness list.

## Acceptance

- [x] Add an explicit `elif pmax == 0: bcrit_pfaffian = None` branch with the same WARN message as the saturated case (collapse to "bcrit_pfaffian marked 'approximation (open-chain projected; full Bloch-Pfaffian deferred U13)'").
- [x] **OR** tighten the predicate to drop the `pmax > 0` guard and let `(pmax - pmin) / pmax` evaluate to `0 / 0 = NaN` → handled by the degenerate branch. Verify NaN propagation is correct.
- [x] Preferred fix: **add the explicit `if pmax == 0:`** branch for clarity (avoids NaN handling). Implemented as `if pmax == 0:` (independent guard) rather than `elif`, so each degenerate case stands alone with its own WARN; both branches reach the same fallback (`bcrit_pfaffian = None`).
- [x] Add the all-zero case to the T6 failure-mode verification list (was: file absent / file saturated / file varying → becomes: file absent / file empty-match / file saturated / file all-zero / file varying). Note: ticket says "file saturated / file all-zero / file varying" but the actual list is now four: file absent / file empty-match (T6) / file saturated (T6) / file all-zero (P2) / file varying. Updated T6 ticket verification accordingly.
- [x] Update the in-script comment block (lines 230-249) to acknowledge the all-zero case. Block now enumerates all four degenerate regimes with the P2 all-closure description.

## Test

`tests/integration/test_pfaffian_degeneracy_detection.py` — Python regression
test pinning all five failure modes (absent / empty-match / saturated / all-zero /
varying). Imports `lecture_13_topological`, monkey-patches `REPO` to a tmpdir
with a synthetic `wire_slim_pfaffian_witness.dat`, stubs `_run_verifier` so
`section_wire_rung()` reaches the Pfaffian block, and asserts `bcrit_pfaffian
is None` for the four degenerate cases plus numeric-minimum for the varying case.
Pre-fix: all-zero test FAIL (`bpf == 0.0` = first B key). Post-fix: all 5 PASS.

## Why this matters

The gate's 4-witness range check (`test_lecture_13_acceptance_gate.sh:86`) would treat a meaningless first-B value (e.g., `0.000` T) as a valid Pfaffian B_crit. Range check could pass trivially if other witnesses are similar, masking real failures. This is the kind of bug that ships a false-positive "phase diagram OK" result to a user running lecture 13.

## Cross-references

- T6 ticket: `tickets/T6_lutchyn_oreg_proxy.md` (saturated-degenerate detection; P2 is the missing all-zero sibling)
- Lecture script context: `scripts/lecture_13_topological.py:230-258`
- Gate auto-detection: `tests/integration/test_lecture_13_acceptance_gate.sh:79-87`

## Blocks / blocked-by

- **Blocks**: P3 (T6 synthetic-fixture test) — P2 must land first to give P3 a fourth failure mode to verify.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Closed 2026-08-08 (single session, this PR-43 follow-up).

## Resolution (2026-08-08)

Implemented as a separate `if pmax == 0:` guard (not `elif` on the saturated
branch) so each failure mode has an independent, explicit WARN. The two
degenerate cases (all-zero vs saturated) differ in their diagnostic context
(`max=0.000e+00` vs `rel_var=N.Ne-NN`), which is useful for distinguishing
closure-regime (s2_sign=0) from FEST-floor (s2_sign=±1 but numerically
degenerate). Test pinned via `test_pfaffian_degeneracy_detection.py`.
