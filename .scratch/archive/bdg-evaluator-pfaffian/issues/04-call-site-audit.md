**Status**: COMPLETE (2026-07-13) — ticket 04 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: task
Status: resolved
Claimed-by: Tiago (via wayfinder session 2026-07-13, work 04)
Blocked by: 01, 02, 03
Resolved: 2026-07-13

# 04 — Call-site audit

## Question

For each BdG per-point call site in `main_topology.f90`, document:

- The file:line of the call.
- The geometry (wire / QW / BHZ / Kitaev).
- The current shape of the call (eigenvalues only, or eigenvalues + eigenvectors, or eigenvalues + eigenpairs).
- Whether eigenvectors are in scope at the call site (i.e., the surrounding code holds them, or they were already discarded post-solve).
- What the call site does with `result%invariant_flag` afterwards — `print`/`write`/`assert`/aggregate — and whether the meaning change from "0/1 heuristic" to "slim Pfaffian sign" or "Kitaev Majorana sign" affects downstream consumers.
- For the straggler `main_topology.f90:1380` `compute_z2_gap` direct call: the geometry, the `gap_threshold` source, and the migration target (likely `eval_bdg_point` once ticket 03 settles the helper's fate).

Sites to audit (from the wayfinder scouting):

- `main_topology.f90:530` — `eval_bdg_point(eigvals_bdg, evp)`
- `main_topology.f90:552` — `eval_bdg_point(eigvals_bdg, wire_eval_p)`
- `main_topology.f90:826` — `eval_bdg_point(eigvals_bdg, evp)`
- `main_topology.f90:1360` — `eval_bdg_point(eigvals_bdg, evp)`
- `main_topology.f90:1380` — `compute_z2_gap(eigen_res_local%eigenvalues, gap_threshold)` (straggler)

Plus any other `compute_z2_gap*` call sites not in `main_topology.f90` — search `src/` for all of them and audit the same five columns.

Deliverable: a table at ticket resolution with one row per call site, the audit columns filled, and a one-line conclusion per site (migrate to slim witness / migrate to Kitaev wrapper / delete / leave).

## Constraints

- Read-only audit; no code changes.
- Cite file:line for every claim.
- The audit's conclusion depends on tickets 01 (API shape), 02 (slim witness location), 03 (helper fate). All three must close before this ticket's conclusions are actionable.

## Out of scope for this ticket

- The actual migration of any call site (downstream execution after this map's tickets close).
- The acceptance-gate wiring (ticket 05).
- New tests (ticket 06).

## Answer

Read-only audit. Five BdG per-point call sites in `src/apps/main_topology.f90` audited (the four `eval_bdg_point` calls + the `wire_pfaffian_witness_sweep` call at `:1371` not in the ticket's scouting list but in scope per the ticket's question "for each BdG per-point call site"). One BdG straggler (`compute_z2_gap` at `:1380`) audited with `gap_threshold` chain traced. Two non-BdG `compute_z2_gap*` call sites (`:375` in `run_qshe_wire` BHZ analytic path; `:1124` BHZ analytic sweep dispatch) audited for completeness per the ticket's "search `src/` for all of them" clause.

### Scope correction — `src/` callers of `compute_z2_gap*`

`grep -rn 'compute_z2_gap' src/` returns only 7 hits, all in `src/apps/main_topology.f90` and `src/physics/topological_analysis.f90`. The 6 src-side call sites are:

- `main_topology.f90:369` — comment (`! First try: spectral edge-state detection with compute_z2_gap_edge.`)
- `main_topology.f90:375` — call (BHZ analytic, not BdG)
- `main_topology.f90:1124` — call (BHZ analytic sweep dispatch)
- `main_topology.f90:1380` — call (BdG sweep straggler)
- `topological_analysis.f90:275,302,318,384,1468,1501-1521` — definitions + comments, not callers

Unit-test callers (8 sites in `tests/unit/test_z2_invariant.pf`, 2 in `tests/unit/test_edge_states.pf`, 2 in `tests/unit/test_phase_diagram.pf`) are in `tests/`, not `src/`. Ticket scope is `src/` only; unit-test callers are noted in §"Unit-test sites (out of audit scope)" below for context.

### Audit table — BdG per-point sites

| # | File:line | Geometry | Call shape | Eigvecs in scope? | `invariant_flag` (or analogue) consumer | Meaning change affects downstream? | Conclusion |
|---|---|---|---|---|---|---|---|
| 1 | `main_topology.f90:530` | wire BdG (`run_bdg_wire`) | eigenvalues only | YES — `eigen_res_local%eigenvectors` held (used by `compute_majorana_profile` at `:575` and `majorana_polarization` at `:595`) | NONE consumed — only `evr%minigap` extracted | NO — minigap field is invariant-witness-agnostic | **leave** — eigenvalues-only seam unchanged; eigenvectors stay on the wire side for Majorana profile + polarization, never cross the seam |
| 2 | `main_topology.f90:552` | wire BdG (`run_bdg_wire`) | eigenvalues only | YES — same `eigen_res_local%eigenvectors` as site 1 | NONE consumed — only `evr%near_zero_count` extracted | NO — `near_zero_count` is `≥2 near-zero modes` count, distinct from invariant; no Pfaffian/Kitaev path here | **leave** — heuristic near-zero count is a separate observable from the Z2 invariant; the wire Pfaffian witness (site 5) provides the invariant, not this site |
| 3 | `main_topology.f90:826` | QW BdG (`run_bdg_qw`) | eigenvalues only | YES — `H_bdg` holds eigenvectors (line 818); used by `compute_majorana_profile` at `:867` | NONE consumed — only `evr%minigap` extracted | NO — same as site 1 | **leave** |
| 4 | `main_topology.f90:1360` | wire BdG sweep (`eval_wire_bdg_gap`) | eigenvalues only | YES — `eigen_res_local%eigenvectors` consumed downstream by site 5 (`:1371`) for Pfaffian | NONE consumed — only `evr%minigap` extracted | NO — same as site 1 | **leave** |
| 5 | `main_topology.f90:1371` | wire BdG sweep (`eval_wire_bdg_gap`) | **CSR matrix** (`type(csr_matrix)` + `n_full`) | N/A — seam consumes CSR, not eigenvectors | `s2_sign == -1 → z2=1; == +1 → z2=0; == 0 → gap-closure fallback (site 6)`. Result `z2` is written to `result%z2_map` (via `compute_wire_bdg_gap_sweep`/`z2_map_int` cast at `:1145`) and printed at `:424` | YES — currently uses the subroutine `wire_pfaffian_witness_sweep`; ticket 01 + 02's resolution locks the seam sibling `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg_local, params)` as the canonical entry point. Mapping is identical (`s2_sign ∈ {-1, 0, +1}` semantics preserved) | **migrate** — replace `call wire_pfaffian_witness_sweep(...)` with `s2_sign = eval_bdg_pfaffian_witness_csr(...)`. The seam sibling is the ticket-02 thin facade with the same signature semantics + the new `bdg_default_pfaffian_floor` SSOT applied via `params`. Subroutine stays in `topological_analysis.f90` for the unit-test dense path (ticket 02 §8) |

### Audit table — BdG straggler

| # | File:line | Geometry | `gap_threshold` source | Migration target |
|---|---|---|---|---|
| 6 | `main_topology.f90:1380` | wire BdG sweep | `run_gap_sweep:1115` — `gap_threshold = 1.0e-4_dp` (hardcoded; passed through `compute_wire_bdg_gap_sweep:1230` → `eval_wire_bdg_gap:1275,1277`; used at `:1380` only inside this routine) | Per ticket 03 (β scope-narrow): **rename to `compute_z2_gap_bhz_heuristic`**, keep as the gap-closure fallback. NOT migrate to `eval_bdg_point` — the heuristic IS the gap-closure fallback exactly because the slim Pfaffian returns 0 at gap closure (where the colormap should still flag Z2=1). Replacing the fallback with the heuristic inside `eval_bdg_point` would lose the "Pf returned 0 ⇒ inconclusive" semantic and confuse the colormap at B_crit. |

### Audit table — non-BdG `compute_z2_gap*` call sites (for completeness per ticket's "search src/" clause)

| # | File:line | Geometry | Call shape | Eigvecs in scope? | Consumer | Conclusion |
|---|---|---|---|---|---|---|
| 7 | `main_topology.f90:375` | wire QSHE / BHZ analytic (`run_qshe_wire` line 254) | eigvals + eigvecs (passed to `compute_z2_gap_edge`) | YES — `eigen_res_local%eigenvectors` passed at `:376` | `result%z2_invariant = compute_z2_gap_edge(...)` — stored in `result%z2_invariant`, printed at `:424`, compared against `bhz_M_eff < 0` fallback at `:383-393` | **leave** — BHZ analytic 4-band path, not BdG. Ticket 03's rename to `compute_z2_gap_edge_bhz_heuristic` makes the constraint explicit; no migration to seam (BHZ is not BdG; the seam is BdG-only per ticket 01) |
| 8 | `main_topology.f90:1124` | sweep dispatch (`run_gap_sweep`) | N/A — passes `cfg_in`, `B_min/max`, `mu_min/max`, `gap_threshold`, outputs `z2_map_int`/`gap_map_real`/`transitions` | N/A | Aggregated into `result%z2_map`, `result%gap_map`, `result%phase_boundary` (lines 1145, 1146, 1163) | **leave** — `compute_z2_gap_sweep` is already BHZ-only per ticket 03 (case `bhz_analytic` is the only arm; `case default` hard-`error stop`s at `:1135`). No refactor needed; dispatch stays |

### Unit-test sites (out of audit scope, noted for context)

`src/` excludes `tests/`. The ticket's "search `src/`" clause ends at the src boundary. The 12 unit-test callers (`tests/unit/test_z2_invariant.pf` ×6, `test_edge_states.pf` ×2, `test_phase_diagram.pf` ×2) all use `compute_z2_gap*` against BHZ eigenvalue fixtures. They are touched by ticket 03's rename (rename is source-only; tests must follow the renamed symbol) but are NOT migrated to the seam — the seam is BdG, and these tests exercise BHZ heuristics. Ticket 04 records no migration target for them; the rename is a mechanical follow-up at execution time (either inside ticket 07's plan close-out, or as a follow-up commit on top of the rename).

### Per-site one-line conclusions (the ticket's required deliverable)

1. **`:530` (wire BdG kz-sweep)** — leave. Eigenvalues-only seam, minigap consumer, eigenvector scope irrelevant.
2. **`:552` (wire BdG kz-sweep)** — leave. Eigenvalues-only seam, near-zero count consumer, separate observable from invariant.
3. **`:826` (QW BdG k-par)** — leave. Eigenvalues-only seam, minigap consumer, eigenvector scope irrelevant.
4. **`:1360` (wire BdG sweep, gap extraction)** — leave. Eigenvalues-only seam, minigap consumer; downstream Pfaffian call lives at site 5.
5. **`:1371` (wire BdG sweep, Pfaffian)** — **migrate** to `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg_local, params)` (ticket 01 + 02 sibling). Drop-in replacement, same semantics, adds `bdg_default_pfaffian_floor` SSOT via `params`. The subroutine `wire_pfaffian_witness_sweep` stays in `topological_analysis.f90` for `test_wire_pfaffian_witness.pf` dense-path unit tests (ticket 02 §8).
6. **`:1380` (wire BdG sweep, gap-closure fallback)** — rename to `compute_z2_gap_bhz_heuristic` per ticket 03 (β). Keep as the gap-closure fallback AFTER `eval_bdg_pfaffian_witness_csr` returns 0. NOT migrate to seam — heuristic IS the fallback's identity.
7. **`:375` (QSHE wire BHZ analytic, edge-state heuristic)** — leave. Rename to `compute_z2_gap_edge_bhz_heuristic` per ticket 03; no seam migration (BHZ is not BdG).
8. **`:1124` (sweep dispatch)** — leave. `compute_z2_gap_sweep` already BHZ-only; dispatch unchanged.

### Migration map (for the executor of ticket 04's downstream work)

```
:530  eval_bdg_point(eigvals_bdg, evp)                          → unchanged
:552  eval_bdg_point(eigvals_bdg, wire_eval_p)                   → unchanged
:826  eval_bdg_point(eigvals_bdg, evp)                          → unchanged
:1360 evr = eval_bdg_point(eigvals_bdg, evp)                    → unchanged
:1371 call wire_pfaffian_witness_sweep(...)                     → s2_sign = eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg_local, params)
:1380 z2 = compute_z2_gap(eigenvalues, gap_threshold)           → z2 = compute_z2_gap_bhz_heuristic(eigenvalues, gap_threshold)
:375  result%z2_invariant = compute_z2_gap_edge(...)            → compute_z2_gap_edge_bhz_heuristic(...)
:1124 call compute_z2_gap_sweep(...)                            → unchanged
```

Net effect on `main_topology.f90`: 2 source changes (one body migration at `:1371`, one rename at `:1380`) + 1 rename at `:375`. 4 sites unchanged. Total: 3 lines touched in `main_topology.f90` plus the new `eval_bdg_pfaffian_witness_csr` import in the `use bdg_observables` block (line 18).

### What this unlocks for the map

- **Ticket 05 (gate wiring)**: the wire-row gate witness is now unambiguously `eval_bdg_pfaffian_witness_csr` (the new site `:1371`); the gap-closure fallback line is `compute_z2_gap_bhz_heuristic` (site `:1380`).
- **Ticket 06 (test coverage)**: the seam-sibling surface (`eval_bdg_pfaffian_witness_csr` + `eval_bdg_kitaev_majorana`) now has exactly one production caller each — `main_topology.f90:1371` for the wire; no QW/Kitaev caller exists yet in production (the Kitaev wrapper is wired only through unit tests at `tests/unit/test_kitaev_*.pf`). QW rung migration to `eval_bdg_kitaev_majorana` is OUT OF THIS MAP'S SCOPE per map "Out of scope: U11 lecture 13 full revamp" (and the QW rung is currently heuristic-only via `eval_bdg_point` at `:826`; no real QW+Kitaev topology path exists).
- **Ticket 07 (plan close-out)**: doc-drift sweep needs to add the `bdg_observables` import dep on `sparse_matrices` + `pfaffian` in `src/physics/AGENTS.md:46-58` (DAG was missing `bdg_observables` entirely; ticket 02 §1, §10); the rename at `:375` + `:1380` is a clean single-line change; the migration at `:1371` is a drop-in.

### Premise checks

- **Ticket's listed sites**: all 5 listed sites exist at the cited lines (`grep -n 'eval_bdg_point\|compute_z2_gap\b' src/apps/main_topology.f90` matches `:530, :552, :826, :1360, :1380` exactly). Ticket 04 scouting is accurate.
- **Ticket's "Plus any other `compute_z2_gap*` call sites"**: src/ has only the 7 hits listed in §"Scope correction". No additional src/ callers.
- **`gap_threshold` source for `:1380`**: confirmed via `grep -n 'gap_threshold' src/apps/main_topology.f90`. Origin at `run_gap_sweep:1115` = `1.0e-4_dp`. Passed verbatim through `compute_wire_bdg_gap_sweep:1228,1230` to `eval_wire_bdg_gap:1275,1277`. Used at `:1380` only. No upstream config knob.
- **Eigenvector scope at sites 1, 2, 3, 4**: confirmed via code reading — `eigen_res_local%eigenvectors` (sites 1, 2, 4) and `H_bdg` (site 3, holding `bdg_result%eigenvectors` from line 818) are in scope but never passed to `eval_bdg_point`.

### Resolution

- Sites 1, 2, 3, 4, 7, 8 — **leave** (renames per ticket 03 for sites 7 and 6; no seam migration).
- Site 5 (`main_topology.f90:1371`) — **migrate** to `eval_bdg_pfaffian_witness_csr` (the ticket-02 seam sibling).
- Site 6 (`main_topology.f90:1380`) — **rename** to `compute_z2_gap_bhz_heuristic` (ticket 03 β) and keep as gap-closure fallback.
- QW+Kitaev rung has no production caller yet; not on this audit. Out of this map's scope per "Out of scope: U11 lecture 13 full revamp".