---
module: bdg_observables
tags: [seam, SSOT, Pfaffian, slim-witness, ticket-07, U2-close-out]
problem_type: invariant-SSOT-discipline
component: bdg_observables
---

# BdG evaluator seam SSOT — slim Pfaffian plug-in + heuristic retirement

## Problem

Three sibling evaluators on three different modules — the heuristic
eigenvalues-only invariant and two ad-hoc per-rung wrappers — with no
single seam consumed by all BdG per-point work. The heuristic
`near_zero_count >= 2` discriminator silently masqueraded as a Z2
invariant on the wire rung. The acceptance gate ran as 3-witness; a
4th slim-Pfaffian row (live site-by-site sweep) was attempted then
reverted (ticket 05) — that row remains reserved for the U13 Bloch-
periodic BdG Pfaffian construction, not this U2 scope.

## Solution

Consolidate the per-point BdG invariant on `bdg_observables.f90` as the
single seam with three faces:

- `eval_bdg_point(eigenvalues, params) → bdg_eval_result_t` —
  eigenvalues-only, minigap + heuristic invariant (kept for QW rung).
- `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) → s2_sign ∈
  {-1, 0, +1}` — wire-rung invariant. S2-projected Pfaffian sign
  (bands 7-8, per k.p block-table SSOT), NOT S1⊗S2 — the S2 strategy
  skips eigendecomp and preserves the H(7,8) asymmetry that yields a
  non-zero Pf. `params` is `bdg_pfaffian_params_t` (factory
  `bdg_pfaffian_params_with_floor`, default `bdg_default_pfaffian_floor`).
  Replaces the heuristic on the wire rung (thin wrapper over
  `wire_pfaffian_witness_sweep` per ticket 04).
- `eval_bdg_kitaev_majorana(H_k_array, k_par_values) → majorana_number ∈
  {-1, 0, +1}` — QW+Kitaev rung. Wraps the existing Kitaev helper.

Scope-narrow `compute_z2_gap` / `compute_z2_gap_edge` to
`compute_z2_gap_bhz_heuristic` / `compute_z2_gap_edge_bhz_heuristic`,
signalling BHZ-only at the call site. Retire the orphan dense
`wire_pfaffian_witness` (+ `s1_project`/`s2_construct` helpers) once its
production call site migrates to the seam sibling (ticket 02).

The Pfaffian site-by-site slim sweep stays as an interim diagnostic;
the wire acceptance gate runs 3-witness (`wire_curve`, `wire_2d`,
`wire_pfaffian` witnesses — `WITNESS_LABEL="3-witness"`). A full B-μ
Bloch-Pfaffian sweep is deferred to U13 (Issue 07), so no 4th live row
ships in U2.

## Why this works

- Pure-function discipline: seam stays one module with three faces,
  not three seams. Dispatch by procedure choice (ADR 0001).
- Layering: seam imports L0 leaves (`sparse_matrices` + `pfaffian`) plus
  one symbol from L3 `topological_analysis` (`wire_pfaffian_witness_sweep`,
  per ticket 04 — kept to avoid re-implementing the CSR-aware S2 row-
  extraction in the seam).
- Heuristic retirement via scope-narrow + rename rather than deletion:
  the heuristic is the gap-closure fallback, not the invariant.
- Live gate witnesses derived from existing Fortran output, no new I/O.
- Magic-number SSOT: `bdg_pfaffian_params_t` + factory
  `bdg_pfaffian_params_with_floor` (`error stop` on `pfaffian_floor <= 0`)
  carry the floor; default `bdg_default_pfaffian_floor = 1.0e-12_dp`. The
  dense-subroutine literals at `topological_analysis.f90:1680`/`:1698`
  retired with `wire_pfaffian_witness` (ticket 02); only the CSR sweep
  site `:1781` consumes the floor.

## Latent bug uncovered during execution (dense variant, now retired)

The dense `wire_pfaffian_witness` (S1+S2 outputs) was extracting
`omega(1:4, 1:4)` from a 16×16 omega structured at `(i, i+8)` indices —
yielding all-zeros by construction (always returned `(0, 0)`, hiding the
strict-test failure). Fixed by using a properly-built local 4×4 omega
(same convention as `wire_pfaffian_witness_sweep`). The dense variant was
then **retired** (ticket 02) once the production call site migrated to the
seam sibling; the fix lives on in the CSR sweep.

## Strict S1+S2 sign agreement is structurally unattainable

For a real-symmetric single-particle BdG with imaginary diagonal s-wave
pairing, the S1 strategy's eigendecomposition returns real orthogonal
eigenvectors, and the off-diagonal pairing contribution to the 4×4
Nambu projection cancels by orthonormality. The S2 strategy avoids
eigendecomp by directly projecting onto bands 7-8, preserving the
H(7,8) asymmetry and yielding a non-zero Pf. Strict
`s1 == s2 .and. s1 /= 0` is therefore impossible; the U2 accepted form
is `s2 ∈ {-1, 0, +1}` with `s2 /= 0` on a non-diagonal synthetic fixture
(the gate-relevant witness). Full S1+S2 strict sign-agreement is the
Majorana-basis Pfaffian problem (Issue 05, deferred to U13).

**User Story 5 (verbatim, spec-of-record at
`.scratch/archive/bdg-evaluator-pfaffian/spec.md:95`):** "As a researcher
running cross-builder BdG unit tests, I want the Pfaffian sibling's
witness assertions to be GREEN today (independent of the deferred
Bloch-periodic BdG construction U13), where the witness is the
S2-projected Pfaffian sign `s2 ∈ {-1, 0, +1}` (range contract pinned by
`s2 >= -1 .and. s2 <= 1` and non-zeroness pinned by `s2 /= 0` on a
non-diagonal synthetic fixture). Full S1+S2 strict sign-agreement is the
Majorana-basis Pfaffian problem (Issue 05), deferred to a follow-up scope
— for U2 the S1 path returns `(0, -1)` on the canonical fixtures, which
is a structurally-valid gap-closure signal in the S1 eigenspace, not a
defect." Test pin: `test_pfaffian_witness_spec_user_story_5_contract` in
`tests/unit/test_bdg_pfaffian_witness_csr.pf`.

## When to use

- New BdG invariant → add a sibling to `bdg_observables.f90`, passing
  the appropriate SSOT params type (`bdg_pfaffian_params_t` for the
  Pfaffian witnesses, `bdg_eval_params_t` for the point evaluator).
  Don't grow the seam into a hub.
- Gate row source of truth is the colormap
  (`output/z2_phase_diagram.dat` z2 column), NOT a separately-emitted
  Fortran file.
- Naming: BHZ-only helpers get `_bhz_heuristic` suffix; never pretend a
  heuristic is a generic Z2 helper.
- Construction-matrix indexing: when building projection helpers, build
  a local 4×4 omega at the projected-subblock indices, NOT extract
  from a larger omega whose structure crosses out of the subblock.

## Source

- Live map + 8 tickets: `.scratch/bdg-u2-actual-ship/`
- Spec-of-record: `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (User
  Story 5 verbatim at line 95)
- Plan: `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md`
- Dense-path witness (retired, ticket 02): was `topological_analysis.f90:1641-1707`
  + `s1_project`/`s2_construct` helpers `:1809-1921`; deleted with the
  production call-site migration to the seam sibling.
- CSR sweep variant (consumes the floor): `topological_analysis.f90` `wire_pfaffian_witness_sweep`
- Seam: `src/physics/bdg_observables.f90` (type at `:62`, factory `:148`)
- Gate: `tests/integration/test_lecture_13_acceptance_gate.sh` (3-witness)
- Verifier: `tests/integration/verify_majorana_polarization.py`
- Lecture: `docs/lecture/13-topological-superconductivity.md` §13.7.4 + §13.7.5
- Verification: ctest `-L unit` 51/51 PASS (unit-count = 51; the 52→51 drop
  is the ctest-target count after retiring `test_wire_pfaffian_witness.pf`).
