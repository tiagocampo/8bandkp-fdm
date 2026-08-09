# T4 — Strict invariant → 4-witness acceptance gate

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch;
> the destination row of this map. Materialized from the map's T4 entry.

Type: task
Status: pending
Blocked by: 03, 05

## Question

This is the **destination** of the whole effort: lift the 4-witness acceptance
gate in `scripts/lecture_13_topological.py` from *approximate* (the per-B
`min-|Pf|` proxy U10 shipped, flagged `"approximation. Full Bloch-Pfaffian
deferred U13"`) to **strict S1×S2 agreement** at μ≈0.6601 across B. The gate only
flips once T1 (builder), T2 (seam), T3 (knob/dispatch), and T5 (the strict TDD
test) are all green — hence T4 is blocked by **both** T3 and T5.

Contract a session should nail down (drawing on the map's §State and the memory
row `project_bdg_u10_t6_execution`, where the degenerate-proxy fallback was set):

- **Regenerate the (B, μ) phase diagram.** Run `topologicalAnalysis` with T3's
  `nk_par>1` on the canonical wire-BdG fixture
  (`tests/.../wire_inas_gaas_bdg_topological_phase.toml`, μ≈0.6601 eV over B),
  producing a `z2_phase_diagram.dat` whose `z2` column is the strict
  S1×S2 product, not the flat single-point slim Pfaffian. The all-`-1` colormap
  U10 T4 found is expected to **gain a topological `+1` region** (or closure
  `0`s) where the gap closes at B_crit — the whole point of lifting the
  approximation.
- **Flip the label.** In `scripts/lecture_13_topological.py`, the 4-witness row
  labels at **L267, L278, L290, L301, L311** (and the comments at L223, L294,
  L351) read `"approximation. Full Bloch-Pfaffian deferred U13"`. Flip to a
  `"strict S1×S2"` label now that the deferral is over. Resolve the `@todo U13`
  markers. The degenerate-proxy branch U10 T6 added (the
  "file-present-but-all-at-floor" fallback, `lecture_13_topological.py`
  `section_wire_rung()`, `_PFAFFIAN_DEGENERACY_TOL=1e-6`) **stays** — a future
  FEST-envelope saturation still degrades gracefully to the approximation label;
  it just should no longer fire for the canonical fixture once the strict path is
  wired.
- **Acceptance gate.** The 4-witness gate's `bcrit_pfaffian` row (the row U10 T5
  made auto-detect numeric Pfaffian witness at `lecture_13_topological.py:81`
  → 4-witness when emitted) now reads a **numeric, non-degenerate** B_crit from
  the strict product and asserts it against `bcrit_2d` (the gap-min strategy at
  `lecture_13_topological.py:178-196`). Strict agreement — not "approximate."
- **Tests.** The `test_lecture_13_acceptance_gate.sh` (`tests/integration/`)
  comment blocks U10 T5 annotated (49-56, 79-86, 95-111) flip from the
  approximation phrasing to strict; the shell gate logic itself should be
  unchanged (U10 T5 already made it auto-detect). Keep the BHZ-only
  `compute_z2_gap_sweep` path on `{0,1}` per the U10 T3 standing decision
  (7+ BHZ unit tests pin BHZ to `{0,1}`; do not churn).

Output: regenerated strict phase diagram, flipped label(s) in
`lecture_13_topological.py`, resolved `@todo U13`, the acceptance-gate test
green. The destination row of the map flips to "done." TDD
(`superpowers:test-driven-development`), but the **red** here is T5's strict
test — T4 makes it green.

## Notes for the claiming session

- This ticket earns the destination because it's the user-visible flip
  (approximate → strict). T1–T3 are machinery; T5 is the red.
- The `@todo U13` resolution is a doc-drift obligation: the project memory rule
  `codebase-doc-drift-prevention` says every PR touching behavior must update
  spec + plan status footers. When you close T4, refresh the parent plan's U13
  status footer (`docs/plans/2026-06-14-001-...-plan.md`) and archive the
  finished-map directory under `.scratch/archive/`.
