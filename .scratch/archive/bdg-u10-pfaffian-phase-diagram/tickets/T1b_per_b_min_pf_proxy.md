# T1b — Per-B min-|Pf| proxy producer (output/wire_slim_pfaffian_witness.dat)

> **Claimed by:** Claude session 2026-07-26 15:30+ (work-mode /wayfinder, T1b's turn).
> **Status:** **closed 2026-07-26** as "producer wired, file emitted with lecture-13 regex format, smoke test registered". Resolution commit on branch `feat/bdg-u10-pfaffian-phase-diagram`.

## Question (task, AFK)

Given the magnitude seam return added by T1a, what is the smallest edit to
`compute_wire_bdg_gap_sweep` (main_topology.f90:1230-1264) and a new
`write_wire_slim_pfaffian_witness` (outputFunctions.f90, mirroring `write_z2_phase_diagram`
at :760-796) that emits `output/wire_slim_pfaffian_witness.dat` with one row per `iB`
and the per-B `min-|Pf|` column?

Specific scope:

1. Accumulate `best_pf_abs` into a new `wire_scan_b` buffer in `compute_wire_bdg_gap_sweep`
   (one row per `iB`). Two reasonable reductions:
   - **(a)**: write `min_mag(iB) = min_mu best_pf_abs(iMu, iB)` — the per-B *minimum* across
     the μ-window. This is the canonical "min-|Pf| per B" that the gate consumes
     (lecture_13_topological.py:222-256 reads `output/wire_slim_pfaffian_witness.dat`,
     matches `B=([...]) |Pf|=([...])` regex per :226).
   - **(b)**: write `mean_mag(iB) = mean_mu best_pf_abs(iMu, iB)` — softer, less informative.
   Decision deferred to runtime: pick (a) by convention — matches the existing
   `project_bdg_u10_prechart.md` framing ("min-|Pf| per B = approximate B_crit proxy").
2. Add `write_wire_slim_pfaffian_witness(cfg, b_grid, min_pf_abs)` to outputFunctions.f90.
   Schema mirrors `write_z2_phase_diagram` (outputFunctions.f90:760-796): `# B(T) |Pf|_min`,
   one row per `iB`, `ES16.8` format. `output/wire_slim_pfaffian_witness.dat` path.
3. Wire the call from `main_topology.f90:1153` (next to the existing `write_z2_phase_diagram`
   call) once `compute_wire_bdg_gap_sweep` finishes.
4. No edits to the slim Pfaffian seam contract (that's T1a). No edits to `bdg_observables.f90`.

**Constraint**: zero change to `z2_map` / `gap_map` semantics. The new file is *informational*
(per Q4 acceptance: "emits a per-B open-chain projected `|Pf|_min` so the 4th gate witness
`bcrit_pfaffian` is approximate-not-deferred"). It does not feed any existing decision.

## Resolution (2026-07-26)

**Per-B min-|Pf| proxy producer wired, file emitted in lecture-13 regex format, smoke test registered and GREEN.**

### Edits applied (uncommitted at session close)

1. `src/apps/main_topology.f90` — `eval_wire_bdg_gap` signature gained a required
   `pf_mag :: real(kind=dp), intent(out)` out-arg (single production caller
   pattern, same family as the T1a seam — out-arg not optional here because
   the only caller routes it into the per-B reduction). Threads the magnitude
   through `eval_bdg_pfaffian_witness_csr`'s `best_pf_abs` OPTIONAL out-arg
   (present branch in the seam). `compute_wire_bdg_gap_sweep` now allocates
   `pf_mag_grid(nMu, nB)`, accumulates the per-call `pf_mag` into it, and
   surfaces it via a new (nMu, nB) `pf_mag_grid` OUT-arg. `run_gap_sweep`
   declares `pf_mag_grid_real` + `min_pf_abs_per_b`, allocates an empty
   `(0,0)` buffer for the `bhz_analytic` and `qw_fukane` branches (out of
   scope per Q4 / Q1 acceptance), and in the `wire_bdg` branch folds the
   magnitude grid into `min_pf_abs_per_b(iB) = minval(pf_mag_grid_real(:, iB))`
   and calls the new writer.
2. `src/io/outputFunctions.f90` — added `write_wire_slim_pfaffian_witness(cfg,
   min_pf_abs)` mirroring the existing `write_z2_phase_diagram` header style
   (nB/nMu preamble, ES16.8 data, `error stop` on open failure). Header
   `# wire slim Pfaffian witness (per-B min |Pf| over mu-window)`, data rows
   emitted as `B=<ES16.8> <space> |Pf|=<ES16.8>` with `trim(adjustl(...))` to
   strip Fortran's width padding so the lecture 13 acceptance-gate reader's
   regex `re.finditer(r"B=([\d.eE+-]+)\s+\|Pf\|=([\d.eE+-]+)", ...)` parses
   naturally. Public-export line extended.
3. `tests/integration/test_wire_slim_pfaffian_witness.sh` — new shell
   wrapper that runs `topologicalAnalysis` on the canonical wire_bdg 2D
   config in a temp workdir and asserts the emitted file's rows match the
   lecture 13 regex (anchored, multi-line, with non-negative finite B / |Pf|
   values). Mirrors the pattern of `test_wire_bdg_topological_2d.sh`.
4. `tests/CMakeLists.txt` — registered `regression_wire_slim_pfaffian_witness`
   with `LABELS "regression" TIMEOUT 900` next to `regression_wire_bdg_topological_2d`.

### Drift from ticket text (caught at execution time)

- Ticket acceptance criterion #2 stated `Format 'B(T) |Pf|_min'` but the
  lecture 13 regex `B=([\d.eE+-]+)\s+\|Pf\|=([\d.eE+-]+)` parses literal-
  token `B=val |Pf|=val` rows, not `B(T) |Pf|_min` header-style rows. The
  Fortran `ES16.8` width pads leading whitespace which would NOT match
  the regex as written. Resolution: emit trimmed character buffers
  (`trim(adjustl(...))`) so the producer lands `B=0.00000000E+00 |Pf|=...`
  exactly. T4 (non-flat assertion, future ticket) should validate against
  the same regex pattern.

### Verification (session 2026-07-26)

- `cmake --build build` clean (no warnings introduced).
- 51/51 unit tests green.
- `regression_wire_bdg_topological` PASS (113 s, 1D wire_bdg sweep).
- `regression_wire_bdg_topological_2d` PASS (319 s, 2D wire_bdg sweep).
- `regression_wire_slim_pfaffian_witness` PASS (313 s, new smoke test).
- Smoke run on canonical wire_bdg 2D config (5 B x 2 mu grid) emits
  `output/wire_slim_pfaffian_witness.dat` with 5 data rows; lecture 13
  regex parses all 5 (`min_pf_abs = 4.0e-08` across the grid — flat by
  construction at the closure region, which is the pre-fix condition
  the U10 destination aims to escape once T2/T3 widen the mu-window).
- No call to `compute_z2_gap_bhz_heuristic` introduced anywhere; the
  new file is open-chain projected Pfaffian only (BHZ heuristic still
  drives the colormap at the closure cell per ticket 05 of
  `.scratch/archive/bdg-evaluator-pfaffian/`, unchanged).

### Doc-drift flagged (out of T1b scope, per `codebase-doc-drift-prevention`)

- `tests/integration/test_slim_pfaffian_witness_projection.py:15` still
  says `RESERVED for U13`. The reservation is now stale (U10 emits the
  file as an open-chain proxy; U13 owns the strict Bloch-Pfaffian).
  Update as a separate doc-drift ticket once T4 / T6 land and the file
  is exercised by the acceptance gate.

## Acceptance (re-derived against actual edits)

- [x] `output/wire_slim_pfaffian_witness.dat` emitted on a canonical wire sweep run
- [x] Format `B=<val> <space> |Pf|=<val>` (matches lecture-13 regex anchored)
- [x] One row per `iB` (5 rows on nB=5; nB rows in general)
- [x] No call to `compute_z2_gap_bhz_heuristic`
- [x] `regression_wire_bdg_topological_2d` still PASSes (319 s)
- [x] New `add_test` registration with smoke check (file exists + nB rows +
      regex parses + values finite non-negative)
- [ ] Doc-drift on `tests/integration/test_slim_pfaffian_witness_projection.py:15`
      "RESERVED for U13" — filed separately, NOT folded in this ticket.

- [ ] `output/wire_slim_pfaffian_witness.dat` emitted on a canonical wire sweep run
- [ ] Format `B(T)  |Pf|_min` (matches `lecture_13_topological.py:226` regex
      `B=([...]) |Pf|=([...])`)
- [ ] One row per `iB` (so `nB` rows)
- [ ] No call to `compute_z2_gap_bhz_heuristic` (file is open-chain projected Pfaffian only)
- [ ] Existing `tests/integration/test_wire_bdg_topological_2d.sh` still PASSes
- [ ] New `add_test` registration in `tests/integration/CMakeLists.txt` (or wherever
      `test_wire_bdg_topological_2d` is registered) — at minimum a smoke check that
      `output/wire_slim_pfaffian_witness.dat` exists and has `nB` rows
- [ ] Doc-drift: `tests/integration/test_slim_pfaffian_witness_projection.py:15` line
      "RESERVED for U13" is now stale; update once the file is in production (separate
      doc-drift ticket, not folded here).

## Blocks / blocked-by

- **Blocks**: T6 (formal min-|Pf| extraction via Lutchyn-Oreg spec, optional per Q5 accept),
  T4 (regression test for the file's existence + structure).
- **Blocked by**: T1a (magnitude seam return must exist).

## Type

task (AFK; bounded edit, mirrors existing `write_z2_phase_diagram` pattern).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram`.
