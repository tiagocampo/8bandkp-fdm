# P8 — Defense-in-depth polish batch (5 SUGGESTIONS from PR review)

## Question (task, AFK)

PR review surfaced 5 SUGGESTION-level polish items not blocking merge. Bundle into a single ticket so they're tracked but don't block PR #43 review.

### Sub-items

1. **`# z2 semantics:` header emission not pinned by any test** — `src/io/outputFunctions.f90:789` emits `# z2 semantics: {-1,+1,0} = Pfaffian native (topological, trivial, closure)` (wire_bdg path) and the BHZ heuristic header at line 794 (other paths). No test asserts these header lines exist. A regression that silently dropped the annotation would not break any test.
   - **Acceptance:** extend `test_wire_bdg_topological_phase.sh` to assert the `# z2 semantics:` line is present and well-formed in the emitted `output/z2_phase_diagram.dat`.

2. **`test_lecture_13_acceptance_gate.sh:108` stale comment** — still reads "append the slim-Pfaffian value only when U13 has landed a live numeric row (PFAFFIAN_DEFERRED == 0)". The T5 update added the "T1b producer" context to the comment block above but didn't refresh this line.
   - **Acceptance:** update line 108 to "append the slim-Pfaffian value when the row is numeric (T1b proxy stand-in; full Bloch-Pfaffian at PHS-invariant k deferred to U13)".

3. **`pf_mag intent(out)` not explicitly initialized** in `eval_wire_bdg_gap` at `src/apps/main_topology.f90:1319`. Fortran does NOT initialize `intent(out)` locals. Currently the only early returns use `error stop` (program termination), but fragile to future maintainer edits.
   - **Acceptance:** add explicit `pf_mag = 0.0_dp` near the top of `eval_wire_bdg_gap`.

4. **Unused `HERE` variable** in 2 new shell tests: `tests/integration/test_wire_bdg_topological_phase.sh:12` and `tests/integration/test_wire_slim_pfaffian_witness.sh:14`. Dead code carried over from `test_wire_bdg_topological.sh` pattern (where it IS used).
   - **Acceptance:** drop the unused `HERE="$(cd "$(dirname "$0")" && pwd)"` lines, or use them where the pattern intends.

5. **Trailing newlines missing** in `tests/integration/test_wire_bdg_topological_phase.sh` and `tests/regression/configs/wire_inas_gaas_bdg_topological_phase.toml`. POSIX tooling warns on missing trailing newlines; canonical files in the same dirs all end with newlines.
   - **Acceptance:** add trailing newlines to both files.

### Bonus item: bash regex tightening

6. **`test_lecture_13_acceptance_gate.sh:86` regex `[0-9.+-]+`** accepts malformed numerics like `++3`, `1.+-`, `.`. Currently the lecture script uses `f"{bcrit_pfaffian:.3f}"` which produces well-formed numerics, but tightening to `^[0-9]+(\.[0-9]+)?$` is defense-in-depth.
   - **Acceptance:** tighten the bash regex.

### Bonus item: YAGNI on optional arg

7. **`bdg_observables.f90:211,219` unused `best_pf_abs_local`** in the `else` branch — the optional `best_pf_abs` argument to the seam could be made required (no `optional`) since the only production caller (`eval_wire_bdg_gap`) always wants it. Compiler optimizes it out anyway, but the optional is over-engineered.
   - **Acceptance:** make `best_pf_abs` required in `eval_bdg_pfaffian_witness_csr`. Update the unit test contract.

## Cross-references

- PR review pass 2026-08-08 (full branch diff, +676/-93)

## Blocks / blocked-by

- **Blocks**: nothing.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
