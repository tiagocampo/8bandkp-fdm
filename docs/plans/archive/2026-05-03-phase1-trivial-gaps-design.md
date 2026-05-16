# Phase 1: Close the Trivial Gaps — Design

Date: 2026-05-03
Consolidates: REVIEW.md groups #45, #5, #16, #10, #8

## Goal

Close 5 INCOMPLETE/MOSTLY COMPLETE groups with minimal effort, update tracking
docs, and archive their plan files. No new physics — just closing gaps and
documenting deliberate decisions.

## Tasks

### Task 1: Close #45 — zdotc subroutine form (effort: 0)

The PR12 plan recommended converting `zdotc` to subroutine ABI. This was
deliberately NOT done: on this platform (gfortran + MKL LP64), `zdotc` returns
`complex(dp)` by value. CLAUDE.md already documents this decision. No code
changes needed.

Files to archive: `2026-04-29-pr12-fixes-plan.md`

### Task 2: Close #5 — Makefile target + z(1) guard (effort: low)

Two minor items from the 2026-03-29 PR review:

**(a) Add `test-fd-sign` Makefile target**
- `tests/fd_sign_compare.f90` exists but has no build/run integration
- Add target to root Makefile (or verify CMake integration)
- Verify the test builds and runs

**(b) Add z(1) guard in `externalFieldSetup_electricField`**
- `src/physics/hamiltonianConstructor.f90`: subroutine has no z(1)==0 check
- Caller-side guard exists in `input_parser.f90` (defense at call site)
- Add defense-in-depth: `if (size(z) > 0 .and. z(1) == 0.0_dp)` check with
  `print *, "ERROR"` + `stop 1`

Files to archive: `2026-03-29-pr-review-fixes-design.md`

### Task 3: Close #16 — Anticrossing prose quantitative (effort: low)

`docs/lecture/02-quantum-well.md` Section B.5 (lines 796-811) describes the
InAsW/GaSbW anticrossing in vague terms: "computed from the eigenvalue sweep"
and "typically 10--20 meV". The gapfix plan called for embedding actual
computed values.

Steps:
1. Run `bandStructure` with `qw_inas_gasb_broken_gap_kpar.cfg`
2. From eigenvalue output, find the k_parallel where e1 and lh1 are closest
3. Compute the hybridization gap in meV
4. Replace lines 799 and 808-809 with explicit values:
   - "At $k_\parallel = X.XX$ Å$^{-1}$"
   - "The computed hybridization gap is **XX.X meV**"

Files to archive: `2026-04-12-qw-documentation-overhaul-design.md`,
`2026-04-12-qw-docs-phase1-plan.md`, `2026-04-12-qw-phase1-gapfix-plan.md`

### Task 4: Close #10 — Figures renamed + wire_optical.png (effort: low)

**(a) Renamed figures — close as-is**
Seven design-specified figures have better names in the implementation. No
action needed; the current names are more descriptive.

**(b) Generate `wire_optical_spectrum.png`**
- Config `wire_gaas_optical_window.cfg` exists
- Run `opticalProperties` to produce absorption data
- Add `fig_wire_optical_spectrum` function to `generate_all_figures.py`
- Plot absorption vs energy for the wire geometry
- Save to `docs/figures/wire_optical_spectrum.png`

Files to archive: `2026-04-04-documentation-overhaul-design.md`

### Task 5: Close #8 — Remove piezoelectric from backlog (effort: 0)

Piezoelectric is excluded by design:
- ZB [001] QW: zero piezoelectric by symmetry
- Wires: negligible compared to band offsets
- The parsed field (`cfg%strain%piezoelectric`) stays in code (doesn't hurt)
- Just document the decision in tracking docs

Files to archive: `2026-04-02-quantum-wire-design.md`

## Finalization Procedure

After all 5 tasks are verified complete, update tracking docs in one batch:

### REVIEW.md updates
For each group (#5, #8, #10, #16, #45):
- Change status to `COMPLETE | Archived`
- Move detailed findings section under a `## Archived Details` heading or remove

### BACKLOG.md updates
- Replace Phase 1 section with: `## Phase 1: COMPLETED (2026-05-03)`
- Renumber Phase 2 → Phase 1, Phase 3 → Phase 2, etc.
- Update summary table

### Archive
Move all plan files from the 5 groups to `docs/plans/archive/`:
```
2026-04-29-pr12-fixes-plan.md
2026-03-29-pr-review-fixes-design.md
2026-04-12-qw-documentation-overhaul-design.md
2026-04-12-qw-docs-phase1-plan.md
2026-04-12-qw-phase1-gapfix-plan.md
2026-04-04-documentation-overhaul-design.md
2026-04-02-quantum-wire-design.md
```
