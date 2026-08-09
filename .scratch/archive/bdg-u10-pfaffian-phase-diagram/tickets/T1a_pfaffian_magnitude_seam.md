# T1a — Pfaffian magnitude return channel through the seam

> **Claimed by:** Claude session 2026-07-26, branch `feat/bdg-u10-pfaffian-phase-diagram`.
> **Status:** **closed 2026-07-26** as "seam extended with optional out-arg, 51/51 unit green".
> Resolution commit: `1ee9ba2`.

## Question (task, AFK)

What is the smallest public-signature extension to `eval_bdg_pfaffian_witness_csr` (and its
underlying `wire_pfaffian_witness_sweep`) that exposes the per-call `best_pf_abs` magnitude
*already cached* inside `wire_pfaffian_witness_sweep` (topological_analysis.f90:1685-1716),
so the sweep host can accumulate a per-B `min-|Pf|` proxy without re-extracting the same
CSR rows or re-calling `complex_pfaffian`?

Specific scope:

1. Add `best_pf_abs :: real(kind=dp)` out-arg to `wire_pfaffian_witness_sweep`
   (topological_analysis.f90:1650-1730). One-line: `real(kind=dp), intent(out) :: best_pf_abs`.
   Assign it before the early-return at :1719 (so all paths set it).
2. Add `best_pf_abs :: real(kind=dp)` out-arg to `eval_bdg_pfaffian_witness_csr`
   (bdg_observables.f90:195-205). Mirror the same shape.
3. The `best_pf_abs` value is a *positive scalar* (magnitude); it never reproduces the
   sign. The sign is already exposed via `s2_sign`. So no semantic interference with the
   Sealed-slim-Pfaffian invariant — this is a **pure additive return**.
4. Update `eval_wire_bdg_gap` (main_topology.f90:1277-1355) to receive the new arg and
   ignore it for now (we will route it in T1b). No behaviour change in T1a's diff.
5. Update any other callers of `eval_bdg_pfaffian_witness_csr` (run_bdg_wire, run_bdg_qw —
   per the seam module's "all three per-point call sites" comment at bdg_observables.f90:6-13)
   to take the new arg. Compile-only verification that the seam contract is consistent.

**Constraint**: no public-signature *behaviour* change for existing callers — they receive
one extra unused `real(dp)` out-arg. This is *not* a breaking change in the seam-sibling sense
(seam siblings are siblings-with-new-shape, siblings-with-extra-arg is the same family).

## Resolution (2026-07-26)

**Seam extended with OPTIONAL `best_pf_abs` out-arg** (not required out-arg as the
ticket text originally framed it; OPTIONAL preserves the existing 3-arg call form
for any out-of-tree callers and matches the F2018 idiom for "sibling with new
return channel" without breaking the seam's prior API contract).

### Drift from ticket text (caught at execution time)

- Ticket claim #5 ("three call sites ... run_bdg_wire, run_bdg_qw ... all three per-point
  call sites") was **stale** at execution time. Post-PR #42 only **one** production
  call site consumes the seam sibling: `eval_wire_bdg_gap` at
  `src/apps/main_topology.f90:1379`. The two `run_bdg_wire` (line 441) and
  `run_bdg_qw` (line 765) subroutines were retired from the seam-sibling path during
  the U2 close-out (PR #42 squash `257b3c4`); they no longer import or call
  `eval_bdg_pfaffian_witness_csr`. The seam module's `bdg_observables.f90:6-13`
  scope comment is itself stale on this point and is out of T1a's scope (filed
  as doc-drift in map's Not-yet-specified).
- Test direct caller: `tests/unit/test_bdg_pfaffian_witness_csr.pf:86` calls
  `wire_pfaffian_witness_sweep` directly (one of the seam-consistency pinning tests);
  updated to thread both `s2_sign_direct` and a new `best_pf_abs_direct`, plus an
  extra `@assertEqual(best_pf_abs_direct, best_pf_abs_seam)` pinning the seam-vs-direct
  consistency on the magnitude.

### Edits applied (commit `1ee9ba2`)

1. `src/physics/topological_analysis.f90:1650-1738` — `wire_pfaffian_witness_sweep`:
   - signature: added `best_pf_abs :: real(kind=dp), intent(out)` after `s2_sign`.
   - dim-guard early returns at :1662-1667 now set `best_pf_abs = 0.0_dp` alongside
     `s2_sign = 0` so the caller never reads uninitialised memory.
   - post-loop assignment at :1738: `best_pf_abs = best_pf` (outside the floor gate,
     so closure-region weak signal is still exposed).
2. `src/physics/bdg_observables.f90:195-216` — `eval_bdg_pfaffian_witness_csr`:
   - signature: added `best_pf_abs :: real(kind=dp), intent(out), optional` (declared
     before the `result(s2_sign)` decl-block per gfortran parsing — `intent(out)`
     attributes are not valid on a `result` variable).
   - delegates via `present(best_pf_abs)` branch: with-arg passes directly to
     `wire_pfaffian_witness_sweep`; without-arg discards into a local dummy so the
     existing 3-arg call form stays green.
3. `tests/unit/test_bdg_pfaffian_witness_csr.pf:80-90` — added `best_pf_abs_seam`,
   `best_pf_abs_direct` locals; updated the direct call to thread the new arg;
   added a new `@assertEqual` pinning magnitude consistency.

### Verification

- `cmake --build build` clean (no warnings introduced).
- 51/51 unit tests green.
- 6/6 directly-relevant BdG/Pfaffian unit tests green
  (`test_bdg_pfaffian_witness_csr`, `test_bdg_evaluator`, `test_bdg_kitaev_majorana`,
  `test_bdg_phs`, `test_bdg_config`, `test_pfaffian`).
- `eval_wire_bdg_gap` at `main_topology.f90:1379` not touched — still calls
  `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg_local, bdg_pfaffian_params_with_floor(...))`
  in its 3-arg form, threading the existing s2_sign + z2 remap path unchanged.

### Doc-drift flagged (out of T1a scope, filed in map Not-yet-specified)

- `src/physics/bdg_observables.f90:6-13` scope comment still says "all three per-point
  call sites" (`run_bdg_wire`, `run_bdg_qw`, `eval_wire_bdg_gap`). Two of the three
  were retired during U2 close-out; the comment is now stale. Update as a
  separate doc-drift ticket per `codebase-doc-drift-prevention`.

## Acceptance (re-derived against actual edits)

- [x] `wire_pfaffian_witness_sweep` has new `best_pf_abs` out-arg, assigned on every path
- [x] `eval_bdg_pfaffian_witness_csr` has new `best_pf_abs` out-arg (OPTIONAL), propagated
- [x] The single production call site (`eval_wire_bdg_gap`) compiles unchanged (uses 3-arg form)
- [x] Test direct caller (`test_bdg_pfaffian_witness_csr.pf:86`) compiles + asserts consistency
- [x] No `error stop` regression; 51/51 unit tests + 6 BdG/Pfaffian unit tests pass
- [x] DRY: the cached `pf_val_best`/`h_proj_best` pattern at topological_analysis.f90:1687-1728
      is *not* duplicated — the new out-arg is sourced from the already-cached `best_pf`.

## Blocks / blocked-by

- **Blocks**: T1b (per-B proxy producer consumes this seam).
- **Blocked by**: none — independent of T2 (window) and T3 (native schema).

## Type

task (AFK, well-defined edit) — resolved.

## Branch

`feat/bdg-u10-pfaffian-phase-diagram`.
