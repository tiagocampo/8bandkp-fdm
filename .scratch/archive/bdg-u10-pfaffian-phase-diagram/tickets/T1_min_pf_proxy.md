# T1 — Wire min-|Pf| per-B proxy producer (informational, no schema break)

> **Claimed by:** Claude session 2026-07-26, branch `feat/bdg-u10-pfaffian-phase-diagram`
> **Status:** resolved 2026-07-26 (see resolution comment below)

## Question (AFK, task-ticket)

Does `compute_wire_bdg_gap_sweep` already expose enough per-(B, μ) Pfaffian magnitudes to write a per-B `min-|Pf|` proxy file (`output/wire_slim_pfaffian_witness.dat`), without changing the slim Pfaffian contract or breaking the `eval_bdg_point` seam?

Specifically: which scalar inside the per-point evaluator (already calling `eval_bdg_pfaffian_witness_csr` per `[main_topology.f90:1379]` and computing `pf_abs` at `[topological_analysis.f90:1710]`) is *already* in scope of the sweep host, and what is the smallest set of edits to (a) accumulate it into a `wire_scan_b` buffer in `compute_wire_bdg_gap_sweep` (one row per `iB`), (b) write the buffer to disk via a new `write_wire_slim_pfaffian_witness` mirroring `write_z2_phase_diagram` (outputFunctions.f90:760 schema only), and (c) not change `eval_bdg_pfaffian_witness_csr`'s public signature?

## Resolution (2026-07-26)

**Answer: NO — the `best_pf` magnitude that `wire_pfaffian_witness_sweep` already computes per (B, μ) is *not* exposed to the sweep host.** The seam currently returns only `s2_sign ∈ {-1, 0, +1}`.

### Evidence (post-#42 merge, `main@257b3c4`)

- `src/physics/bdg_observables.f90:195-205` — `eval_bdg_pfaffian_witness_csr` is a thin wrapper:
  ```fortran
  call wire_pfaffian_witness_sweep(H_bdg_csr, Nbdg, params%pfaffian_floor, s2_sign)
  ```
  The wrapped subroutine accumulates `best_pf` (the max-site `|Pf|` for the chosen band-major scan in `topological_analysis.f90:1685-1716`) but discards it; only `s2_sign` crosses the seam boundary.
- `src/apps/main_topology.f90:1356-1363` — `eval_wire_bdg_gap` does not see `best_pf` either; it routes through `eval_bdg_point` (minigap) and `wire_pfaffian_witness_sweep` (sign) with no channel for the magnitude.
- `src/physics/topological_analysis.f90:1710-1716` — `pf_abs` is computed and cached into `best_pf` *locally*; the winning-site 4×4 subblock (`h_proj_best`) and the Pfaffian value (`pf_val_best`) are also cached locally for the sign extraction at :1720-1728. None of these return-values cross the subroutine boundary.

### Implication

The "no schema break" constraint in T1's question is **not satisfiable as posed**: producing a per-B `min-|Pf|` requires *either* (a) a magnitude return channel through `eval_bdg_pfaffian_witness_csr` / `wire_pfaffian_witness_sweep`, *or* (b) a second pass with its own per-(B, μ) magnitude loop. Option (a) is the cheaper edit (one extra `real(kind=dp)` out-arg) and preserves the seam's "sibling" shape; option (b) requires a second eigendecomposition per gridpoint (~nB×nMu extra FEAST calls, ~5-10× cost — the eigensolver pre-empts the cheap-path argument).

### Decision recording

- **T1 resolution**: closed as "answer is no". The act of producing the per-B proxy file is **not** shipped in this ticket.
- **New follow-up tickets to file** (split T1 into T1a + T1b for boundary clarity):
  - **T1a — Pfaffian magnitude seam return** (task, AFK): add `min_pf_abs` (or `best_pf_abs`) out-arg to `eval_bdg_pfaffian_witness_csr` and `wire_pfaffian_witness_sweep`. One-line schema extension, no caller behaviour change. Mirror `bdg_pfaffian_params_t` factory pattern (bdg_observables.f90:148).
  - **T1b — Per-B `min-|Pf|` proxy producer** (task, AFK, blocked by T1a): accumulate `best_pf` (per-`iB`) in `compute_wire_bdg_gap_sweep`; emit `output/wire_slim_pfaffian_witness.dat` via new `write_wire_slim_pfaffian_witness` mirroring `write_z2_phase_diagram` (outputFunctions.f90:760 schema, B(T)  |Pf|_min columns).
- **Reason for split**: T2/T3 both need a clean "current contract" ground; if T1 also touches the seam, three tickets touch the same public surface, and the diff scope drifts. T1a isolates the contract change; T1b is a pure producer that consumes T1a's seam. U13 still owns the strict Bloch-Pfaffian — the open-chain projected `|Pf|_min` here is the wiring milestone per Q4 acceptance.
- **Doc-drift to flag separately** (per `codebase-doc-drift-prevention`): `tests/integration/test_slim_pfaffian_witness_projection.py:15` reserves `output/wire_slim_pfaffian_witness.dat` for U13. That reservation is now stale (Q4 acceptance reverses it). File a doc-drift ticket — not folded into T1, since U13's scope is "strict invariant" which is different from U10's "open-chain proxy".

### Re-derive ground truth (post-merge drift verification)

All four call sites above are on `main@257b3c4` (post-#42 squash). No drift between this resolution and the merged source.

## Acceptance

- [ ] `output/wire_slim_pfaffian_witness.dat` emitted on a canonical wire sweep run.
- [ ] Format `B(T)  |Pf|_min` (or equivalent, matched against `lecture_13_topological.py:226` regex `B=([…]) |Pf|=([…])`).
- [ ] One row per `iB` (so nB rows).
- [ ] No change to `eval_bdg_pfaffian_witness_csr` public signature.
- [ ] No BHZ_heuristic fallback (file is open-chain projected Pfaffian only).

## Blocks / blocked-by

- **Blocks**: T6 (formal min-|Pf| extraction via Lutchyn-Oreg spec, optional per Q5 accept), T4 (regression test for the file's existence + structure).
- **Blocked by**: none — independent of T2/T3; can run in parallel.

## Type
task (HITL-grillable, AFK-drivable).

## Branch
`feat/bdg-u10-pfaffian-phase-diagram` (to be created at session start with `git switch -c`).
