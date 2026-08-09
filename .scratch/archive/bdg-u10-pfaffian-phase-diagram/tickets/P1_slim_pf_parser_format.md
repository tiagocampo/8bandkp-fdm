# P1 — `test_slim_pfaffian_witness_projection.py` parse function format mismatch (CRITICAL)

## Question (task, AFK)

`tests/integration/test_slim_pfaffian_witness_projection.py:82-95` defines `parse_slim_pf_witness` which looks for `line.startswith("slim_pf_sign")` rows and returns int signs. The T1b producer (`write_wire_slim_pfaffian_witness` in `src/io/outputFunctions.f90`, lines ~833-874) emits `B=<val> |Pf|=<val>` rows with magnitudes, not signs.

The test passes when the file is absent (deferred branch at line 109). Now that `regression_wire_slim_pfaffian_witness` and `regression_wire_bdg_topological_2d` both emit the file, **the test fails whenever invoked** with `FAIL: slim Pfaffian witness returned 0 across all 0 witness rows`. Confirmed by direct invocation 2026-08-08 against the file left behind from a prior regression run.

**Silent failure mode:** the test is not registered as a ctest entry — it's only run by hand. So `ctest -l regression` doesn't surface it; only direct invocation (or a future ctest registration) would. The branch currently builds clean + all unit/regression tests pass, hiding the bug.

## Acceptance

- [ ] `parse_slim_pf_witness` updated to match the actual emitted format: `re.finditer(r"B=([\d.eE+-]+)\s+\|Pf\|=([\d.eE+-]+)", text)` (the same regex the lecture script uses at `lecture_13_topological.py:240`).
- [ ] Return type changed: list of `(B, magnitude)` tuples or list of dicts `{B: float, Pf: float}` instead of `list of int signs`.
- [ ] Update the parse contract (docstring at `test_slim_pfaffian_witness_projection.py:82-90`) to reflect the new semantics: per-B min-|Pf| magnitude proxy, not per-(B, μ) sign.
- [ ] Downstream assertions at lines 116, 152, 158 updated to match the new return shape.
- [ ] Test passes when run against the canonical fixture's emitted `output/wire_slim_pfaffian_witness.dat`.
- [ ] Test still passes when file is absent (deferred branch preserved).
- [ ] Register the test as a ctest entry under the `regression` label so future regressions surface immediately.

## Why this matters

If unfixed, the test silently fails whenever a developer invokes it directly post-regression-run, and the failure message is misleading ("0 across all 0 witness rows" implies the Pfaffian is broken when in fact it's a parser contract mismatch). Worse: if a future developer registers this test as a ctest entry without fixing the parser, every PR will be red.

## Cross-references

- Memory: `project_bdg_u10_post_t6_docdrift.md` (doc-drift cleanup context)
- T1b ticket: `tickets/T1b_per_b_min_pf_proxy.md` (producer landed in commit `0a62143`)
- T1a ticket: `tickets/T1a_pfaffian_magnitude_seam.md` (seam extension)
- Lecture script regex precedent: `scripts/lecture_13_topological.py:240`

## Blocks / blocked-by

- **Blocks**: P3 (T6 synthetic-fixture test) — both P1 and P3 touch the same parse path.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
