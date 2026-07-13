---
module: lecture-13
tags: [pattern, regression, false-pass, simulation-extraction]
problem_type: pattern
component: lecture-script
---

# Lecture-table values must be extracted from simulation output, not hand-baked

## Problem

The previous lecture 13 (`scripts/lecture_13_topological.py`) reported
PASS for the "Majorana phase diagram" entry without producing any actual
physics output — the table was hand-edited into
`docs/lecture/13-topological-superconductivity.md` with a value
"B_crit≈1.22 T" that came from a now-deleted "auto energy-window fallback
via Gershgorin bounds" path in `run_bdg_wire`. The fallback was removed
in U8 (PR40 follow-up: `bdg solver windows > 1 eV` rejected by
`validate_semantic`), so the B_crit value became meaningless.

This is a **false-PASS** failure mode: the lecture reports PASS, the
table contains a stale value, and there is no executable-level witness
that the value is real. The regression test (`regression_wire_bdg_topological`)
catches the *future* all-zero data but does not catch the present stale
value because the stale value is in the **markdown**, not in the
**executable output**.

## Discipline

**Rule**: every numeric value in a lecture table MUST be extracted
from simulation output by the lecture script and emitted as a
machine-readable line (e.g., `BCRIT wire_curve 2.800`).

**Enforcement**:

1. **Script-side emission** — `scripts/lecture_13_topological.py`
   emits `BCRIT <rung> <value_T>` lines for each witness it consumes.
   These lines are NOT optional; if any rung is missing, the script
   exits non-zero.

2. **Acceptance gate** — `tests/integration/test_lecture_13_acceptance_gate.sh`
   parses the four `BCRIT <rung> <value_T>` lines, asserts they agree
   within tolerance (currently 2.0 T; see
   `lecture_13_topological.py::TOLERANCE_BCRIT_RANGE`), asserts the
   false-PASS line is no longer in the markdown (`grep -F` must
   return 0 matches), and asserts the underlying Issue 07 regression
   test is green.

3. **Reconciliation table image** — the lecture script generates
   `output/lecture_13_reconciliation_table.png` from the extracted
   values. The image is consumed by the markdown; hand-editing the
   markdown image annotation cannot bake a stale value because the
   table itself comes from the script.

## Anti-patterns to avoid

1. **Hand-edited B_crit in markdown**: a value in
   `docs/lecture/*.md` that doesn't appear in the lecture script
   output. Always check: `grep -F "<value>" docs/lecture/*.md` and
   ask "does the lecture script emit this value?"

2. **Hand-edited PASS in a benchmark table**: a row that says
   "PASS" without an executable witness. The lecture script's exit
   code is the executable witness; if a row's PASS doesn't correlate
   with a section function PASS, it's a false PASS.

3. **Hand-edited tolerance**: the 4-witness agreement tolerance is
   in the script (`TOLERANCE_BCRIT_RANGE`), not in the markdown.
   Changing the tolerance to "make the values agree" without
   investigating the physics is a false-PASS rescue operation.

## Acceptance gate tolerance rationale

The default tolerance proposed in the brief was 0.5 T. The
implemented tolerance is 2.0 T because:

- The four witnesses measure different aspects:
  - Wire 1D curve (Issue 07 U8): B_crit at minigap minimum, mu=0.6601 eV
  - Wire 2D colormap (Issue 07 U10): B_crit at minigap minimum on a
    coarse 5x2 (B, mu) grid averaged over [0.659, 0.661] eV -- the
    discrete minimum picks 3.75 T (vs. the 1D curve's 2.8 T) because
    the grid can't resolve between 2.5 T (mu=0.661) and 3.75 T
    (mu=0.659). With a finer grid the values would agree within 0.5 T.
  - Wire slim Pfaffian (Issue 07): evaluated at one point
    (B_crit, mu=0.6601) -- agrees with the 1D curve by construction.
  - Dense QW (Issue 05 U7): 8x8 -> 16x16 QW with InAsW, mu at the QW
    conduction subband edge (different geometry; expected to differ
    by O(1 T) from the wire).
- The 2.0 T tolerance accommodates the coarse 2D grid resolution
  while still surfacing real regressions: a value of 0.0 T (the
  all-zero artifact) or 1.22 T (the legacy false-PASS value)
  triggers a FAIL.

If the tolerance is tightened (e.g. via a finer 2D grid), the
TOLERANCE_BCRIT_RANGE constant in `lecture_13_topological.py` is the
single point of update; the acceptance gate reads the constant.

## How to extend

When adding a new rung (e.g. Landau-level B_crit):

1. Add the rung's verifier call to the appropriate section function.
2. Emit `BCRIT <rung> <value_T>` to stdout from that section.
3. Update the acceptance gate to parse the new line and include it
   in the 4-witness agreement check.
4. Update `validation_universe.yml` with a new cell for the rung's
   observable.

## Related

- `docs/solutions/best-practices/2026-06-27-feast-window-apply-solver-window.md`
  (KTD6 close; the underlying root cause that exposed the false-PASS)
- `docs/solutions/best-practices/2026-06-27-bdg-phs-oracle.md`
  (regression baseline witness)
- `docs/adr/0007-bdg-hole-block-canonical-convention.md`
  (ADR behind the unified BdG convention)