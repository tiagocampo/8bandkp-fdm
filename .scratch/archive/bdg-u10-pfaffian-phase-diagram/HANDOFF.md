# U10 — bulk minigap + (B, μ) phase diagram from the Pfaffian — Wayfinder handoff

> **Read this first.** Future `/wayfinder` session entry point for Unit U10 of
> `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`. The seam
> + slim Pfaffian plug-in + heuristic retirement that U10 was waiting on
> shipped on 2026-07-13 via `.scratch/archive/bdg-evaluator-pfaffian/` (7 tickets,
> all resolved, frontier empty). This handoff names the destination, fixes
> the state, marks the loose, and points at the standing preferences so the
> next chart starts grounded.

## Destination

Close U10 end-to-end with the `wire_bdg` / `qw_fukane` sweep paths driving
the slim projected Pfaffian witness (per `eval_bdg_pfaffian_witness_csr`
seam sibling in `bdg_observables.f90`), producing a regenerated, non-flat
**(B, μ) phase diagram** that the existing 4-witness acceptance gate can
verify. The diagram's `z2` column is the slim witness's per-(B, μ) sign;
the `gap_map` column is the (still-windowed) `2·minval(|E|)` near-zero
indicator (see §"Loose" below).

A "shipped" U10 looks like:

- The seam sibling is the canonical per-point invariant for both `wire_bdg`
  and `qw_fukane` sweep paths — `compute_z2_gap` /
  `compute_z2_gap_bhz_heuristic` is only the gap-closure fallback.
- `output/z2_phase_diagram.dat` is regenerated end-to-end by the new path;
  the 4-witness gate (already wired) reads it without modification.
- `regression_wire_bdg_topological_2d` exercises the **slim-witness** phase
  diagram, not just "non-flat"; the L13 lecture script table's `BCRIT`
  values stay consistent across re-runs.
- `tests/unit/test_phase_diagram.pf` (existing) is updated to cover the
  seam-sibling dispatch; new tests pin the gap-vs-invariant caveat.

## State of the tree (what's already shipped)

These are the locked decisions that this map **must not re-litigate**. If a
ticket starts re-deriving one of these, it should be closed and the answer
pointed to instead.

| Locked decision | Ticket / file | Source |
|---|---|---|
| `bdg_observables.f90:eval_bdg_point` is the single BdG per-point invariant SSOT | `.scratch/archive/bdg-evaluator-pfaffian/issues/01-evaluator-api-shape.md` | ticket 01 |
| `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) → s2_sign ∈ {-1, 0, +1}` seam sibling exists | `.scratch/archive/bdg-evaluator-pfaffian/issues/02-slim-pfaffian-location.md` | ticket 02 |
| `bdg_default_pfaffian_floor = 1.0e-12_dp` is the SSOT magic number (replaces the literal at `topological_analysis.f90:1663, 1681, 1766`) | `.scratch/archive/bdg-evaluator-pfaffian/issues/02-slim-pfaffian-location.md` §3 | ticket 02 |
| `compute_z2_gap` / `compute_z2_gap_edge` renamed to `compute_z2_gap_bhz_heuristic` / `compute_z2_gap_edge_bhz_heuristic`; BHZ-only; kept as gap-closure fallback | `.scratch/archive/bdg-evaluator-pfaffian/issues/03-compute-z2-gap-fate.md` | ticket 03 |
| `main_topology.f90:1371` migrated to `eval_bdg_pfaffian_witness_csr`; `:375` / `:1380` renamed to `..._bhz_heuristic`; `compute_z2_gap_sweep` (`:1124`) untouched | `.scratch/archive/bdg-evaluator-pfaffian/issues/04-call-site-audit.md` | ticket 04 |
| 4-witness acceptance gate (1.0 T); slim Pfaffian row read live from `output/z2_phase_diagram.dat` z2 column at mu ≈ 0.6601 ± 0.0001 (z2==-1 first B-row) | `.scratch/archive/bdg-evaluator-pfaffian/issues/05-acceptance-gate-wiring.md` | ticket 05 |
| Strict `@assertTrue` in `test_wire_pfaffian_witness.pf:96` becomes GREEN via non-diagonal synthetic fixture (option a) — decoupled from U13 | `.scratch/archive/bdg-evaluator-pfaffian/issues/06-test-coverage-audit.md` | ticket 06 |
| Lecture 13 §13.0 + §13.7.1 + §13.7.4 + new §13.7.5 disclosure text updated; `scripts/lecture_13_topological.py:section_wire_rung` reads colormap z2 column; gate `PFAFFIAN_DEFERRED` branch stripped | `.scratch/archive/bdg-evaluator-pfaffian/issues/07-plan-backlog-close-out.md` | ticket 07 |

Code / artifacts on disk today:

- `src/physics/bdg_observables.f90` — seam with `eval_bdg_pfaffian_witness_csr` + `eval_bdg_kitaev_majorana` siblings + `bdg_default_pfaffian_floor` SSOT.
- `src/physics/topological_analysis.f90` — `wire_pfaffian_witness` + `wire_pfaffian_witness_sweep` (dense-path, unit-test only); `compute_z2_gap_bhz_heuristic` / `compute_z2_gap_edge_bhz_heuristic` (gap-closure fallback); `compute_z2_gap_sweep` (BHZ-only, calls `eval_bhz_analytic` directly, untouched).
- `src/apps/main_topology.f90`:
  - `:1115` — `gap_threshold = 1.0e-4_dp` (hardcoded; passed verbatim through `:1230` → `:1277`).
  - `:1122` — `bhz_analytic` dispatch (BHZ-only sweep arm).
  - `:1124` — `compute_z2_gap_sweep` call (BHZ sweep, unchanged).
  - `:1371` — `eval_bdg_pfaffian_witness_csr` migration site (was `wire_pfaffian_witness_sweep`).
  - `:1380` — `compute_z2_gap_bhz_heuristic` call (gap-closure fallback AFTER slim Pf returns 0).
- `src/io/outputFunctions.f90` — `write_z2_phase_diagram` (4 inline `open()` blocks consolidated per Phase-24 follow-up Agent B).
- `output/z2_phase_diagram.dat` — produced by `compute_wire_bdg_gap_sweep`; consumed by `scripts/lecture_13_topological.py:section_wire_rung` (z2==-1 first-B lookup at mu ≈ 0.6601) and `tests/regression/test_wire_bdg_topological_2d.sh` (non-flat assertion; **does not yet exercise the slim-witness column**).
- `tests/regression/regression_wire_bdg_topological_2d` — golden-data check; the .dat is on disk and asserted non-flat (per BACKLOG §Phase 23 + `7312e97`).
- `tests/integration/verify_majorana_polarization.py` — `_gate_row_colormap_present` SKIP precondition helper (per ticket 07).
- `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` — durable pattern record.

## Loose (chart these as tickets on the new map)

These are the U10-specific gaps the new Wayfinder session should grill into
tickets. None are large enough to chart without more grilling — each is
named here as the first-pass fog.

1. **Gap-vs-invariant caveat (parent plan §U10 line 425, "Approach" ¶2).**
   The existing `wire_bdg` `gap_map` is `2·minval(|E|)` over the ±5δ₀-windowed
   set at a single kz=0 — it is the *windowed near-zero indicator*, not
   R9's continuous kz-sweep BZ-minimum. Decide: (a) add a kz-sweep sub-step
   for the true bulk minigap, or (b) label the wire `gap_map` as the
   windowed quantity and defer the BZ-minimum to U13. This sub-decision
   was called out in the parent plan as decided "in U10 against cost" —
   the new map needs to land that decision.

2. **Per-(B, μ) recomputation discipline.** The flat-phase-diagram fix
   (recompute the invariant per (B, μ) with a local config copy, not
   broadcast one evaluation) was an established codebase pattern but the
   sweep contract (`z2_map(nMu, nB)`) doesn't currently enforce it. The
   new map should pin the discipline into a test or invariant.

3. **The `output/z2_phase_diagram.dat` schema.** The lecture script's
   parser at `scripts/lecture_13_topological.py:section_wire_rung` expects
   `mu B z2 gap` (whitespace-separated, 4 columns). The regression's
   `.sh` and Python verifier assume a different layout (skip z2 column —
   see `7312e97`). Reconcile the writers + readers; pin a schema in a
   `tests/unit/test_phase_diagram.pf` assertion.

4. **The L13 lecture script's gate row derives the slim Pfaffian value
   from the colormap; U10's sweep must produce a colormap that the gate
   can read.** Verify the dispatch at `main_topology.f90:1122-1124` keeps
   the BHZ-only `compute_z2_gap_sweep` arm separate from the BdG sweep;
   the new map's tickets should test that the wire sweep (which produces
   `z2_phase_diagram.dat`) and the BHZ sweep (which does not) stay on
   disjoint code paths.

5. **The cross-rung comparison.** The dense-QW rung uses
   `compute_qw_fukane_gap_sweep`; the wire rung uses
   `compute_wire_bdg_gap_sweep`. U10 should verify the slim Pfaffian
   reads identically across both rungs at comparable (B, μ) points —
   catching any divergence in the seam sibling's per-rung semantics.

6. **The `bcrit_pfaffian` vs the phase-diagram B_crit boundary.** The 4-witness
   gate's `bcrit_pfaffian` is read from a single mu ≈ 0.6601 row; U10's
   phase diagram contains the B_crit boundary as a 2D region. The
   transition-detection helper (`detect_z2_transitions` /
   `is_z2_transition` in `topological_analysis.f90`) should produce the
   same `bcrit_pfaffian` at mu ≈ 0.6601 as the gate reads — assert this.

7. **Existing `regression_wire_bdg_topological_2d` test asserts non-flat
   only; U10 should add a slim-witness assertion (e.g., that the z2 column
   flips sign at the Lutchyn-Oreg-predicted boundary).** Update or
   supersede.

8. **Whether the BHZ sweep needs any change.** Per ticket 03 the renamed
   `_bhz_heuristic` functions stay as gap-closure fallback; the BHZ
   sweep at `main_topology.f90:1124` calls `compute_z2_gap_sweep` which
   calls `eval_bhz_analytic` directly. The new map should verify the BHZ
   sweep's output (`bhz_phase_diagram.dat`?) doesn't accidentally become
   a second wire-style colormap that confuses the gate.

## Standing preferences (every session on this map should consult)

- **Skills.** `/superpowers:test-driven-development` for each new test
  (slim-witness assertions are TDD-doubled); `/superpowers:receiving-code-review`
  for any ticket that touches the seam; `/superbrains:subagent-driven-development`
  is the right pattern for the per-(B,μ) recomputation sweep — the
  parallel-rung comparison tickets fit a 2–4-agent fan-out.
- **Domain glossary.** `CONTEXT.md` + `UBIQUITOUS_LANGUAGE.md`. BdG Nambu
  layout and pairing: `src/physics/AGENTS.md` §"BdG Nambu structure".
- **ADRs.** 0001 (no polymorphic builder types — dispatch by enum/plain
  record); 0002 (no new TOML fields — U10 must use existing
  `[bdg]` / `[topology]` sections); 0003 (sweep loop in app, pure
  per-point evaluator in physics — KTD3 of the parent plan); 0005 (one
  stable window per sweep via `apply_solver_window`); 0007 (canonical
  hole block `-conjg(H₀(-k))`); 0008 (BdG P1 invariants, including
  `bdg_default_near_zero_frac = 0.001_dp`).
- **Engineering principles.** Repo root `CLAUDE.md` → "Engineering
  Principles". The BdG-evaluator map established a precedent: SSOT
  promotion for any new magic number (cf. `bdg_default_pfaffian_floor`),
  heuristic retirement via scope-narrow + rename rather than deletion,
  live gate witnesses derived from existing Fortran output (no new
  I/O).
- **Convention reminders.** pFUnit `@assertEqual`/`@assertTrue` are
  single-line (no `&` continuation — `feedback_pfunit_macro_single_line`
  memory); no `stop 1` in new code (use `error stop '<descriptive message>'`);
  no Co-Authored-By trailer on commits (`feedback_no_coauthor_trailer`).

## Out of scope (rule-of-scope decisions for this map)

These belong to other destinations. If a ticket starts to drift here,
close it as out of scope and link from this section.

- **U9 — BdG spectral function + BdG LDOS on BdG.** Different observable,
  different seam (green_functions.f90). U10's `z2_phase_diagram.dat` is
  the invariant signature, not the LDOS. Separate scoped PR.
- **U13 — periodic/Bloch BdG construction.** Per CLAUDE.md Known Issues,
  U13 is a separate scoped PR. U10 works with the slim witness on the
  open-chain finite-difference CSR; the slim → Bloch-Pfaffian swap is
  U13's deliverable. Any ticket that depends on a periodic-along-z
  supercell belongs in U13, not here.
- **BLOCKING-EMPIRICAL polarization verifier gap.** Stays as its own
  `@todo U13` per ticket 07. Polarization is a separate observable
  (Sticlet `P_M`, not invariant).
- **U2 re-scoping as a planning act.** Already absorbed by destination
  (d+ii) of the BdG-evaluator map; the seam shape is locked.
- **Topology module split** (`topological_analysis.f90` ≈ 1740 lines).
  Separate concern; BACKLOG Phase 18 architecture-deepening item.
- **U11 lecture 13 full revamp.** Per the BdG-evaluator map's
  out-of-scope; U10 may update the §13.7 disclosure if the phase diagram
  gets new figures, but the lecture revamp is separate.

## First-pass fog (Not yet specified)

The Wayfinder session should populate this section with anything that
graduates from "I can tell it's coming but can't pin it down" into a
sharp ticket. Initial entries to seed:

- Whether the wire `(B, μ)` colormap's `gap_map` column should be relabeled
  to `gap_windowed` (option b in §Loose #1) or whether a kz-sweep sub-step
  is cheap enough to add (option a). Decision trades cost vs spec fidelity.
- Whether `compute_qw_fukane_gap_sweep` needs an analogous seam migration
  to its per-point invariant, or whether it's already Pfaffian-driven
  (per `compute_z2_fukane_qw_result` which evaluates at TRIM).
- Whether the new sweep emits a QW colormap (`output/qw_z2_phase_diagram.dat`)
  for cross-rung comparison, or whether the gate stays wire-only.
- The interaction with the auto-window gate: `validate_semantic` rejects
  explicitly-set Gershgorin-scale BdG solver windows (per Phase-24
  follow-up Agent A); the phase-diagram sweep's energy window needs to
  honor this constraint.

## Files to Read first

When starting a Wayfinder chart session on this map, Read these in order:

1. **This file** (`.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`) — the destination + state + loose.
2. `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U10 (lines 415–429) — the spec.
3. `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §R9 (Requirements → Observable coverage) — what U10 satisfies.
4. `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §KTD8 — the periodic/Bloch invariant caveat (which U13 owns; U10 works with the open-chain witness).
5. `.scratch/archive/bdg-evaluator-pfaffian/map.md` + the 7 tickets — the locked decisions above.
6. `src/physics/AGENTS.md` §"BdG Nambu structure" + §"Dependency DAG" — the seam's place in the module graph.
7. `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` — the seam-SSOT pattern in durable prose.

## Source

- Parent plan §U10: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` lines 415–429.
- Locked-decisions source: `.scratch/archive/bdg-evaluator-pfaffian/` (7 resolved tickets).
- Status footer: parent plan `status_note` ("U2 closed 2026-07-13 ... **U1, U9, U10 still open**; U13 still explicitly deferred").
- BACKLOG summary: row 26 (U2 close-out) + row 5 + row 530 + row 729 (still-open list).
- REVIEW row 79: parent-plan status (mirrors `status_note`).