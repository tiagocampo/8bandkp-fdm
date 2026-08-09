# T3 — Drop z2 remap inversion; write Pfaffian-native (-1/+1/0)

## Claimed (2026-07-26, session: T3-grill-via-handoff)

Claimed by Tiago de Campos (tiagocampo@gmail.com) per wayfinder convention. Awaiting grill
resolution on Q3a/Q3b/Q3c before execution. Standing decisions locked (per map §Notes):
BHZ-heuristic fallback removed from `eval_wire_bdg_gap` (wire path); `compute_z2_gap_sweep`
(BHZ-only KTD8) retains heuristic. Native convention `z2 ∈ {-1, +1, 0}` written directly
without remap.

Pre-grill investigation summary (for grilling session context):

**Q3a fact (both paths converge on same file):** `compute_z2_gap_sweep` (BHZ-only,
`main_topology.f90:1130`) writes `z2_map_int` ∈ {0, 1} via `eval_bhz_analytic`
(`topological_analysis.f90:1611-1614`, where `z2=1` ↔ `M_eff < 0`). `compute_wire_bdg_gap_sweep`
(wire BdG, `main_topology.f90:1142`) writes `z2_map_int` ∈ {0, 1} via the post-PR #42 remap at
`main_topology.f90:1420-1423` (where `z2=1` ↔ Pfaffian `s2_sign==-1`). Both paths then convert
integer→real at `main_topology.f90:1156` and write to the SAME `output/z2_phase_diagram.dat`
via `write_z2_phase_diagram` (`main_topology.f90:1162`). T2 added a config-driven
`phase_diagram_file` key so wire_BdG sibling fixture writes to `z2_phase_diagram_phase.dat`
(a separate file) — this means the canonical file is now BHZ-dominant in practice, but the
writer itself is dispatch-blind.

**Q3b fact (downstream consumers already read native):**
- `tests/integration/verify_majorana_polarization.py:74` — reads `z2 == -1` (native).
- `scripts/lecture_13_topological.py:178` — reads parts[0]/parts[3] only (B, gap), never parts[2] (z2).
- `docs/lecture/13-topological-superconductivity.md:872` — text refers to "z2==-1 at mu≈0.6601".
- `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md:1127` — gate precondition reads z2 == -1.
No code path reads `z2 == 0` or `z2 == 1` as "topological"; the BHZ `{0, 1}` convention is
writer-internal only.

**Q3c fact (doc-drift candidates):** All `z2==-1` references in docs/plans/lecture already
match native convention. After T3 ships, the remap block at `main_topology.f90:1420-1434` is
deleted and the BHZ-heuristic comment at lines 1425-1432 (`z2_phase_diagram.dat z2 column at
the B_crit cell should treat it as BHZ-heuristic-decided`) is also deleted. No downstream
doc updates needed; per `codebase-doc-drift-prevention` the comment deletion IS the doc
update (the comment was the inverted-convention stale reference).

## Question (HITL, grill-ticket)

The post-PR #42 wiring at `main_topology.f90:1367-1396` remaps:
```
s2_sign == -1 → z2 = 1   ! topological
s2_sign == +1 → z2 = 0   ! trivial
s2_sign ==  0 → compute_z2_gap_bhz_heuristic(...)  ! closure fallback (U13-stall)
```
…into the `z2` column of `write_z2_phase_diagram` (outputFunctions.f90:760). The remap:
1. Inverts Pfaffian-native semantics (every other witness in `lecture_13_topological.py` uses `z2 == -1` ↔ topological).
2. Forces the gate precondition (`gate_row_colormap_present`, `verify_majorana_polarization.py:74` requires `z2 == -1`) to never match the wire sweep's column.
3. Keeps the BHZ heuristic fallback silently shadowing the decision at closure μ — a U13-stall stub that masquerades as a Pfaffian decision.

Per Q3 acceptance, this ticket removes the remap and writes `s2_sign` directly (`-1/+1/0`) into the `z2` column, and **eliminates the BHZ_heuristic fallback** in the wire path (the fallback is only retained in `compute_z2_gap_sweep` for `bhz_analytic` dispatches, per `compute_z2_gap_sweep` KTD8).

Open questions the grill session must answer:
- **Q3a**: Does `compute_z2_gap_sweep` (which calls `eval_bhz_analytic` directly per `topological_analysis.f90:1522`) write into the same `output/z2_phase_diagram.dat`, or a different file? If same, its column-write logic must use a different schema (e.g. `z2 ∈ {0, 1}` heuristic semantics) — we may need a column comment header to disambiguate.
- **Q3b**: Where in `lecture_13_topological.py` is `z2 == 1` read as "topological" / `z2 == 0` as "trivial"? Update to Pfaffian-native (`z2 == -1`/`z2 == +1`) OR keep both conventions side-by-side?
- **Q3c**: Where else (doc, plans, lecture figure captions) does the inverted `z2 ∈ {0, 1}` convention persist? Update all in the same PR per `codebase-doc-drift-prevention`.

## Acceptance

- [ ] `output/z2_phase_diagram.dat` z2 column is integer `{-1, +1, 0}` for wire BdG dispatches.
- [ ] BHZ_heuristic fallback removed from `eval_wire_bdg_gap` (no call to `compute_z2_gap_bhz_heuristic`).
- [ ] `compute_z2_gap_sweep` (BHZ-only) still calls heuristic — `compute_z2_gap_bhz_heuristic` retained.
- [ ] `verify_majorana_polarization.py:74` `gate_row_colormap_present` reads `z2 == -1` only, `z2 == +1` no-op (existing behavior, just confirmed).
- [ ] `lecture_13_topological.py` updated to read native.
- [ ] Doc-drift: `docs/lecture/13-topological-superconductivity.md` z2 column references updated.

## Resolution (2026-07-27)

Grill resolved: Q3a = (a) header annotation (BHZ path unchanged; preserves 7+ existing unit tests
that pin `compute_z2_gap_bhz_heuristic` to emit `{0,1}`); Q3b = (1) ship native wire convention
directly into z2 column; Q3c = (α) comment deletion is sufficient (`main_topology.f90:1425-1432`
BHZ-heuristic-fallback comment block was the last inverted-convention reference).

**Source changes:**
- `src/apps/main_topology.f90` — `eval_wire_bdg_gap`, lines 1420-1434: collapsed remap+fallback into
  `z2 = s2_sign`. Comment block updated to native-only language (no BHZ-heuristic invocation).
- `src/io/outputFunctions.f90` — `write_z2_phase_diagram`, lines 780-789: added `# z2 semantics:`
  header line that branches on `cfg%topo%sweep_model` (wire_bdg → `{-1,+1,0}` Pfaffian native;
  non-wire_bdg → `{0,1}` heuristic). Header annotates per-file convention; downstream consumers
  read the header to disambiguate.

**Test changes:**
- `tests/integration/test_wire_bdg_topological_phase.sh` — TDD-red pin: asserts z2 column subset
  of `{-1,0,1}` AND at least one row reads -1 (this is the topological cell at low B, pre-B_crit).
  Pre-T3 the remap emitted `{0,1}` so the -1 assertion turned the test red; post-T3 the wire
  path emits native and the test passes.

**Existing-asset impact (no changes needed):**
- `tests/integration/verify_majorana_polarization.py:74` — `z2 = int(float(parts[2]))` then
  `z2 == -1` already reads native.
- `scripts/lecture_13_topological.py:178-196` — reads parts[0] (B), parts[3] (gap); never parts[2]
  (z2). No change needed.
- `docs/lecture/13-topological-superconductivity.md`, `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md`,
  `docs/plans/BACKLOG.md` — every z2 column reference already uses native `z2 == -1` semantics.
  Doc-drift per `codebase-doc-drift-prevention` resolved by the source comment deletion.

**Standing decision added (T3, 2026-07-27):** `wire_bdg` path emits Pfaffian-native `{-1,+1,0}`;
non-wire_bdg paths retain heuristic `{0,1}`; `write_z2_phase_diagram` annotates per-file convention
via header. See `map.md` §Notes.

**Verification (final state):**
- TDD-red: pre-T3 → `FAIL: ... z2 column never reads -1; got values=[1]` (180 error class).
- TDD-green: post-T3 → `PASS: ... native z2 column`.
- Unit: 51/51 green (T3b change is local; no unit-test surface touched).
- Wire-BdG regression: 5/5 green (canonical `_2d` + `wire_slim_pfaffian` + new
  `topological_phase` + `wire_bdg_topological` + `wire_bdg_strain_shift`).
- Plan §U10 + standing decisions need no further amendment — T3 is a clarification, not a
  destination change.

## Blocks / blocked-by

- **Blocks**: T4 (regression test for non-flat) — T4 asserts ≥ 2 distinct `z2` values; post-T3, distinct values are `{-1, +1}` or `{-1, 0}` or `{+1, 0}` — depends on the actual sweep shape.
- **Blocked by**: T2 preferred (sweep window now spans both phase regions; without T2 a single tight window still yields all-zeros from Pfaffian).

## Type
grilling (HITL).

## Branch
`feat/bdg-u10-pfaffian-phase-diagram`.
