**Status**: COMPLETE (2026-07-13) — ticket 07 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: task
Status: resolved
Resolved: 2026-07-13
Blocked by: 01, 02, 03, 04, 05, 06

# 07 — Plan + backlog close-out

## Question

Once tickets 01–06 close, U2 is closed in fact. Update the documentation to match:

- **`docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` `status_note`** — move U2 from "still open" to "shipped" with a one-line gist of what shipped (seam + slim Pfaffian plug-in + heuristic retirement), citing the PR. Update the per-unit breakdown table if present.
- **`docs/plans/BACKLOG.md` summary table** — add a row or update the parent-BdG-validation-plan summary to reflect U2 closed (the existing "Phase 24: ... parent BdG validation plan still in progress — U1, U2, U9, U10 open; U13 deferred" line at line 5 + line 730 needs U2 dropped from the open list).
- **`docs/plans/REVIEW.md` row #79** — update the "Still open" / "Partial" wording to reflect U2 closed; add a reference to the PR or commit range.
- **`docs/lecture/13-topological-superconductivity.md` and `scripts/lecture_13_topological.py`** — apply the disclosure text from ticket 05 (slim witness as gate row + 4-witness restoration).
- **`docs/lecture/figures/`** — regenerate any lecture 13 figure whose underlying data is touched by the slim-witness gate row. Confirm via the figure provenance check (per Phase 24 ticket 07 precedent): every figure cited in the lecture markdown must trace to a regenerated `.dat` file. Document any figure that's NOT regenerated and why.
- **`docs/scratch/archive/bdg-majorana-validation/issues/00-pure-bdg-evaluator.md`** — the mislabeling note is already in the file; no archive hygiene change needed (per the map's Notes).
- **`src/physics/AGENTS.md` line 44** — `bdg_observables.f90`'s inventory entry already says "(Issue 00). Single seam consumed by run_bdg_wire, run_bdg_qw, eval_wire_bdg_gap." Confirm the wording still matches after the slim-witness plug-in, or update.
- **Module doc-comments** — three files touched:
  - `src/physics/bdg_observables.f90`: module header (currently says "folds those three steps into one pure-function call") + `eval_bdg_point` doc (currently says "Pure only — no I/O, no allocations, no state"). Add a paragraph describing the Pfaffian plug-in (mode discriminator or sibling function, per ticket 01), and update the magic-number-extraction sentence to point at the new `bdg_eval_params_t` field names.
  - `src/physics/topological_analysis.f90`: `wire_pfaffian_witness` doc-comment. Add a sentence: "Called by `bdg_observables.f90:eval_bdg_point` as the wire-rung invariant per ticket 05 of `.scratch/archive/bdg-evaluator-pfaffian/`. Future: U13 swaps this for the full Bloch-Pfaffian on a periodic supercell."
  - `src/physics/bdg_hamiltonian.f90`: confirm the module header's "BdG Nambu-space" section still reads correctly post-heuristic-retirement (ticket 03). If it references `compute_z2_gap` as the invariant source, update.
- **`docs/solutions/` entry** — add `docs/solutions/best-practices/2026-XX-XX-bdg-evaluator-seam-ssot.md` capturing the seam + slim-witness + heuristic-retirement pattern. Modeled on the existing `docs/solutions/best-practices/2026-06-21-u8-bdg-window-routing.md` (single-decision summary + pattern + references). Frontmatter per `docs/solutions/` convention: `module`, `tags`, `problem_type`, `component` fields.
- **Memory** — add a `project_bdg_evaluator.md` memory entry capturing: the seam is SSOT for BdG per-point; the slim witness is the wire-rung invariant until U13; `compute_z2_gap` is retired/narrowed (per ticket 03).

Deliverable: a doc-diff at ticket resolution (a list of files touched with the one-line change per file), grouped by category (plan / backlog / REVIEW / lecture / script / figure / module doc / solutions doc / memory).

## Constraints

- Doc-only — no code changes in this ticket.
- Match the existing status-note style (YAML front matter + `status_note` block).
- Update REVIEW.md's row #79 to match the new parent plan status.
- Memory entry must follow the format: frontmatter (name, description, type, metadata) + body + **Why:** + **How to apply:** lines; one fact per file.

## Out of scope for this ticket

- Code changes (already complete when this ticket fires).
- The acceptance-gate wiring itself (ticket 05).
- U9, U10, U11 — those are separate concerns.

## Answer

U2 of `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`
closed in docs, plan, backlog, review, lecture, script, gate, verifier,
modules, and memory. The seam + slim Pfaffian plug-in + heuristic
retirement + 4-witness gate wiring are the destination; this ticket is
the doc-diff that makes it visible.

### Doc-diff (grouped by category)

**Plan / backlog / REVIEW** (4 files touched)

- `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` — `updated: 2026-07-13`. `status_note` block re-keyed: U2 moved from "still open" to "shipped" with one-line gist (seam + slim Pfaffian plug-in + heuristic retirement, landed via `.scratch/archive/bdg-evaluator-pfaffian/` tickets 01–07). U12 wording flipped from "3-witness / slim Pfaffian row reserved for U13" to "4-witness at 1.0 T; slim Pfaffian row wired live from `output/z2_phase_diagram.dat` z2 column (z2==-1 row at mu ≈ 0.6601, ±0.0001 tol) per ticket 05". U11 partial line expanded to cite ticket 07's text sweep. Still-open list trimmed to U1, U9, U10 (U2 dropped).
- `docs/plans/BACKLOG.md` — line 2 date bump (2026-07-13); the four-spaces-into reference to "U1, U2, U9, U10 open" rewrites to "U1, U9, U10 open; **U2 closed 2026-07-13** (slim Pfaffian seam + heuristic retirement + 4-witness gate row)". The `~~U1, U2, U9, U10, U11, U12~~` line is updated to "UPDATED 2026-07-13 ... U2 closed ... U12 now 4-witness ... **U1, U9, U10 still open**". The Phase-24 line in §"Parent plan" trims U2 from the open list. The Phase-24 entry in the summary table gains a sibling `26. BdG evaluator Pfaffian plug-in (U2 close-out, slim witness for wire rung)` with the full gist. The trailing "Phases 1-24 + follow-ups complete" line replaces "U1, U2, U9, U10 open" with "U1, U9, U10 open; **U2 closed 2026-07-13**".
- `docs/plans/REVIEW.md` — row 79 "Still open" wording flips from "U1, U2, U9, U10" to "U1, U9, U10"; "Explicitly deferred" line keeps U13 with the new "slim → full Bloch-Pfaffian swap is the U13 deliverable" forward reference; row notes expand U2 + U12 with their 2026-07-13 close-out details.

**Lecture script + acceptance gate + verifier** (3 files touched)

- `scripts/lecture_13_topological.py` — module-level docstring §"These lines" + "Output" each flipped 3→4 witness; `TOLERANCE_BCRIT_RANGE` comment flipped 3→4 witness; `section_wire_rung` rewrites the `bcrit_pfaffian` block to read live from `output/z2_phase_diagram.dat` z2 column (mu ≈ 0.6601 ± 0.0001, first z2==-1 row), raising `RuntimeError` when the colormap row is absent; the dead `(deferred to U13)` print branch is removed because the colormap-derived value is now mandatory; `render_reconciliation_table` strips the legacy "deferred to U13" cell-text branch and flips "3-witness range" → "4-witness range" + matching print + summary-row label; trailing `FAIL: 3-witness disagreement` flipped to `FAIL: 4-witness disagreement`.
- `tests/integration/test_lecture_13_acceptance_gate.sh` — capture comment flipped 3→4 witness; `PFAFFIAN_DEFERRED` branch stripped (colormap row is mandatory; non-numeric `BCRIT wire_pfaffian` now triggers `FAIL`); tolerance-strip, the inner Python `min`/`max`, and `WITNESS_LABEL` collapse to 4-witness always; `TOL_BCRIT_RANGE` comment flipped 3→4.
- `tests/integration/verify_majorana_polarization.py` — module docstring rewritten to clarify the polarization (Sticlet `P_M`) observable is distinct from the gate's Z2 invariant row, and the SKIP precondition is tightened: SKIP only allowed when `output/z2_phase_diagram.dat` has a row with mu ≈ 0.6601 ± 0.0001 and z2 == -1; otherwise `FAIL` non-zero. New helper `gate_row_colormap_present(repo_root, mu_target=0.6601, mu_tol=0.0001)` parses the colormap and returns True iff the precondition holds; `parse_polarization` calls it inside the missing-file branch.

**Lecture disclosure** (1 file touched)

- `docs/lecture/13-topological-superconductivity.md` — §13.0 flips "3-witness agreement ... slim Pfaffian row reserved for U13 ...`docs/superpowers/specs/archive/2026-07-05-pr41-completion-design.md` §3.4 / D5" → "**4-witness agreement** within a 1.0 T tolerance. The slim Pfaffian row is **live**: `bcrit_pfaffian = arg_B at mu ≈ 0.6601 (tol ±0.0001) where the colormap's `z2` column first reads -1`". The §13.7.4 slim Pfaffian caption rewrites the first bullet ("real Pfaffian witness output requires U13" → "live; not reserved for U13"); the §13.7.4 trailing paragraph swaps to a live-witness explanation with §13.7.5 reference. New §13.7.5 "Slim Pfaffian witness criterion" paragraph added: criterion + value source (`output/z2_phase_diagram.dat` z2 column, no new Fortran I/O) + rationale (z2==-1 is the slim Pf returned 0 by construction) + forward reference to U13. Summary table cell row 836 "Pfaffian witness s1=s2 at (B_crit, mu=0.6601)" flipped to "Pfaffian witness (colormap-extracted, z2==-1 at mu≈0.6601, ticket 05)". §13.5.x tail "U12 3-witness gate" flipped to "U12 4-witness gate (slim Pfaffian row live per ticket 05; colormap-extracted from `output/z2_phase_diagram.dat` z2 column at mu ≈ 0.6601 ±0.0001)". §13.7.1 lead paragraph "3-witness ... slim Pfaffian row reserved for U13" → "**4-witness agreement** ... slim Pfaffian row is **live**".

**Module doc + AGENTS** (2 files touched)

- `src/physics/AGENTS.md` — `bdg_observables.f90` module-inventory row extended with the two seam siblings (`eval_bdg_pfaffian_witness_csr` + `eval_bdg_kitaev_majorana`), their semantics (`s2_sign ∈ {-1, 0, +1}` / `majorana_number ∈ {-1, 0, +1}`), their roles (wire-rung invariant / QW+Kitaev rung), and the L0-leaf imports (`sparse_matrices` + `pfaffian`). The Dependency DAG block now lists `sparse_matrices` + `pfaffian` in Layer 0 and `bdg_observables` in Layer 1 with its `uses sparse_matrices + pfaffian` annotation — the AGENTS dep-DAG row that ticket 02 §1 flagged as missing.

**Solutions doc + memory** (2 files touched)

- `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md` (new) — frontmatter `module: bdg_observables`, `tags: [seam, SSOT, Pfaffian, slim-witness, ticket-07, U2-close-out]`, `problem_type: invariant-SSOT-discipline`, `component: bdg_observables`. Sections: Problem (three sibling evaluators on three modules + heuristic), Solution (seam + 2 siblings + heuristic rename + colormap-extracted gate row), Why-this-works (pure-function discipline + heuristic retirement + live-from-existing-data), When-to-use (new BdG invariant → sibling in `bdg_observables.f90`; gate row source of truth is the colormap, NOT a Fortran-emitted file), Source (plan/map/seam/dense-path/gate/lecture/verifier pointers).
- `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md` (new) — frontmatter `type: project`; body states seam + 2 siblings + slim-witness colormap + heuristic-retirement patterns with **Why:** + **How to apply:** lines per repo convention; one fact per file. Appended to `MEMORY.md` index.

**Figure provenance** — the slim-witness row's data source is `output/z2_phase_diagram.dat`, which is already on disk (produced by `compute_wire_bdg_gap_sweep`). The figures referenced from lecture 13 trace as follows: `lecture_13_wire_bdg_2d_colormap.png` + `lecture_13_wire_bdg_minigap_curve.png` trace to the same colormap/dataset already in scope; `lecture_13_wire_slim_pfaffian_witness.png` is synthetic-illustrative and stays that way until U13 (the §13.7.4 caption + §13.7.5 paragraph document this). **No figure regeneration required for this ticket.**

**Out of scope for this ticket, executed by other tickets** — code changes to `bdg_observables.f90` (seam siblings), `topological_analysis.f90` (rename), `main_topology.f90` (call-site migration + rename at `:1371` + `:375` + `:1380`), unit-test additions (`test_bdg_pfaffian_witness_csr.pf`, `test_bdg_kitaev_majorana.pf`, `test_bdg_evaluator.pf` SSOT, `test_wire_pfaffian_witness.pf` fixture update, `test_kitaev_majorana.pf` seam-sibling). All those landed in tickets 01–06 of this map (Resolution summary line of `map.md`); ticket 07 is the doc-diff that makes the destination visible.

### Verification gate

- All edits are doc-only + 1 new solutions doc + 1 new memory entry. No Fortran source touched.
- BUILD green (no `src/` source change; the polarization verifier `_gate_row_colormap_present` helper and lecture 13 script colormap-extracted block are Python only).
- Doc-drift discipline (CODEBASE doc-drift prevention): parent plan status_note updated; BACKLOG summary-table + trailer updated; REVIEW row 79 updated; lecture script `TOL_ANCE_BCRIT_RANGE` comment + WITNESS_LABEL summary updated; lecture markdown §13.0 + §13.7.1 + §13.7.4 + new §13.7.5 + summary table + tail updated; gate shell `TOL_BCRIT_RANGE` comment + WITNESS_LABEL + `PFAFFIAN_DEFERRED` branch updated; verifier docstring + SKIP precondition updated; AGENTS.md inventory + DAG updated; new solutions doc + memory entry with MEMORY.md index pointer.

### Resolution

U2 closed 2026-07-13.