**Status**: COMPLETE (2026-07-13) — ticket 03 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: grilling
Status: resolved
Claimed-by: Tiago (via wayfinder session 2026-07-13)
Blocked by: 01
Resolved: 2026-07-13

# 03 — Z2-gap helper fate

## Question

`compute_z2_gap` and `compute_z2_gap_edge` live at `topological_analysis.f90:275` and `:318`. `compute_z2_gap_sweep` (`:1468-1521`) is restricted to `sweep_model=bhz_analytic` only and uses `compute_z2_gap` internally. After (d+ii), what happens to them?

Three fates on the table:

**(α) Retire outright.** Delete `compute_z2_gap` and `compute_z2_gap_edge`. Refactor `compute_z2_gap_sweep` to take a generic invariant source via the seam (or fold it into a BHZ-only path that calls the slim/Kitaev witness — likely outside this map's scope).

**(β) Scope-narrow to BHZ.** Keep `compute_z2_gap` only as the BHZ analytic-model invariant; rename for clarity (e.g., `compute_z2_gap_bhz_analytic`); restrict `compute_z2_gap_sweep` to use it; delete `compute_z2_gap_edge` (no current callers per `topological_analysis.f90:19-39` public list).

**(γ) Leave as legacy.** Keep all three untouched. The Pfaffian plug-in is additive; old callers keep their old behaviour. (Risk: future readers confuse the two invariant sources.)

The direct call site at `main_topology.f90:1380` (`compute_z2_gap(eigen_res_local%eigenvalues, gap_threshold)`) is one site to migrate regardless — it's the straggler outside the 4 routed through `eval_bdg_point`. Migration of that one site is ticket 04's scope; this ticket decides what `compute_z2_gap` *becomes* (deleted, BHZ-only, or legacy) once its callers are routed elsewhere.

Decide which fate, and document the consequence for `compute_z2_gap_sweep` (does it need a refactor, or does its `bhz_analytic`-only restriction make the answer trivial?).

## Constraints

- Do not regress `compute_z2_gap_sweep`'s BHZ behavior — the analytic BHZ model is a documented test artifact in `tests/`.
- ADR 0001: no polymorphic evaluator types — if `compute_z2_gap_sweep` is refactored, the seam is an enum-tagged record, not a class.
- The straggler call site `main_topology.f90:1380` migrates regardless of this ticket's decision; ticket 04 owns the migration.

## Out of scope for this ticket

- The migration of any specific call site (ticket 04).
- The acceptance-gate wiring (ticket 05).
- The slim witness's noise-floor tolerance (ticket 02).

## Answer

**Fate: (β) — scope-narrow to BHZ.** All three helpers get a rename that pins down "this is the BHZ heuristic"; their bodies stay unchanged; their callers stay unchanged in this map's scope.

### Premise corrections first

1. **`compute_z2_gap_sweep` does NOT call `compute_z2_gap` internally.** It dispatches on `cfg%topo%sweep_model` and calls `eval_bhz_analytic(B_val, mu_val, cfg, z2, gap, status)` directly (topological_analysis.f90:1501-1509, 1585-1602). The ticket's premise that "it uses `compute_z2_gap` internally" is wrong. So scope-narrowing `compute_z2_gap` does not require touching `compute_z2_gap_sweep` — they're already decoupled.
2. **`compute_z2_gap_edge` is NOT call-free.** Production caller at `main_topology.f90:375`, inside the BHZ analytic-model path. The ticket's premise that "no current callers per `topological_analysis.f90:19-39` public list" was the public-export list, but `main_topology.f90` IS a consumer. The caller is also BHZ-only: it's gated by the same `bhz_*` setup block (lines 350-360) and the `bhz_M_eff` M_eff fallback at line 383-393 is its semantic ground truth.

### The renames

| Old name | New name | Rationale |
|---|---|---|
| `compute_z2_gap` (topological_analysis.f90:275) | `compute_z2_gap_bhz_heuristic` | Hardcodes `n_in_gap >= 2 ⇒ z2=1`. Not a topological invariant. The 6 unit-test callers (test_z2_invariant.pf:39,72,265,298,328,339; test_edge_states.pf:42,78) all pass BHZ eigenvalues — the heuristic coincidentally gives the right answer for BHZ fixtures. Rename makes the limitation explicit; the 6 unit tests stay as-is (out of this map's scope; ticket 07's doc-drift sweep may note them). |
| `compute_z2_gap_edge` (topological_analysis.f90:318) | `compute_z2_gap_edge_bhz_heuristic` | Hardcodes `n_bands = 4` (line 338) — only BHZ uses 4-band sites. 8-band k.p uses 8-band sites. The single production caller (`main_topology.f90:375`) is inside the BHZ analytical-model path with M_eff as ground truth. Rename pins the constraint. |
| `compute_z2_gap_sweep` (topological_analysis.f90:1468) | unchanged | Already BHZ-only (`case ('bhz_analytic')` is the only `case` arm; `case default` hard-`error stop`s). Dispatch stays at `main_topology.f90:1122-1136`'s `select case` on `sweep_model`. ADR 0001's "enum-tagged record, not class" satisfied: the dispatch is a string-mode discriminator, not polymorphism. |

### Why not (α) — retire outright

Retiring `compute_z2_gap` outright would break 6 unit tests with no slim-Pfaffian replacement (the tests are BHZ eigenvalue fixtures; the slim Pfaffian doesn't apply — it's a BdG primitive). Retiring `compute_z2_gap_edge` outright would remove the spatial refinement in the BHZ detection path, forcing every BHZ call to fall through to the M_eff ground truth (which is fine, but a bigger PR). Both retirements are out of this map's scope; ticket 07's doc-drift sweep can file them as a follow-up if desired.

### Why not (γ) — leave as legacy

The destination explicitly states "retired outright, or scope-narrowed to `bhz_analytic` — whichever ticket 03 settles". Legacy is the "neither" option and conflicts with the goal of "single invariant source-of-truth for all BdG per-point work" (the destination's first bullet). Future readers would confuse the heuristic with the slim Pfaffian.

### Consequence for the straggler at `main_topology.f90:1380`

The straggler stays on `compute_z2_gap_bhz_heuristic` after the rename. It's the gap-closure fallback in `eval_wire_bdg_gap`: when the slim Pfaffian witness returns 0 (inconclusive), the code falls back to "≥2 states near zero ⇒ Z2=1" to preserve the open→close→reopen colormap pattern at the B_crit point. This is the **only legitimate use of the heuristic**: at gap closure, even a real topological invariant is inconclusive, so the heuristic is the best available signal. Ticket 04 owns the migration; the rename makes the gap-closure fallback's nature explicit.

### What this unblocks

- **Ticket 04**: the straggler migrates to `compute_z2_gap_bhz_heuristic` (rename only — same body, same site).
- **Ticket 05**: the acceptance gate's wire row uses the slim Pfaffian as the invariant discriminator; the gate's "gap-closure fallback" line stays on the renamed heuristic.
- **Ticket 06**: test coverage audit now has a clear "this is BHZ, not BdG" surface to test against.
- **Ticket 07**: the plan status note can record "compute_z2_gap → compute_z2_gap_bhz_heuristic (BHZ scope)" without ambiguity.

### Resolution

- Fate: **(β) scope-narrow to BHZ**.
- Renames: `compute_z2_gap` → `compute_z2_gap_bhz_heuristic`; `compute_z2_gap_edge` → `compute_z2_gap_edge_bhz_heuristic`.
- `compute_z2_gap_sweep`: unchanged. Already BHZ-only; dispatch stays at `main_topology.f90:1122`.
- Straggler at `main_topology.f90:1380`: keeps the heuristic (renamed), as the gap-closure fallback after slim Pfaffian returns 0.
- Out of scope for this map: retiring the heuristic outright (would break 6 unit tests); refactoring `compute_z2_gap_sweep` (already shape-preserving).