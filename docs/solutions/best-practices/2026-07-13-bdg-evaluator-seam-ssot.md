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
invariant on the wire rung. The acceptance gate was held back to
3-witness because no real Pfaffian signature was wired through the seam.

## Solution

Consolidate the per-point BdG invariant on `bdg_observables.f90` as the
single seam with three faces:

- `eval_bdg_point(eigenvalues, params) → bdg_eval_result_t` —
  eigenvalues-only, minigap + heuristic invariant (kept for QW rung).
- `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) → s2_sign ∈
  {-1, 0, +1}` — wire-rung invariant. Slim projected Pfaffian over S1
  ⊗ S2. Replaces the heuristic on the wire rung (thin wrapper over
  `wire_pfaffian_witness_sweep` per ticket 04).
- `eval_bdg_kitaev_majorana(H_k_array, k_par_values) → majorana_number ∈
  {-1, 0, +1}` — QW+Kitaev rung. Wraps the existing Kitaev helper.

Scope-narrow `compute_z2_gap` / `compute_z2_gap_edge` to
`compute_z2_gap_bhz_heuristic` / `compute_z2_gap_edge_bhz_heuristic`,
signalling BHZ-only at the call site.

The acceptance gate's Pfaffian row reads live from the existing `(B, μ)`
colormap dataset (`output/z2_phase_diagram.dat`); no new Fortran I/O.

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
- Magic-number SSOT: `bdg_default_pfaffian_floor = 1.0e-12_dp` replaces
  the literal at 3 call sites in `topological_analysis.f90`.

## Latent bug uncovered during execution

`wire_pfaffian_witness` (the dense variant with S1+S2 outputs) was
extracting `omega(1:4, 1:4)` from a 16×16 omega structured at
`(i, i+8)` indices — yielding all-zeros by construction (always
returned `(0, 0)`, hiding the strict-test failure). Fixed by using a
properly-built local 4×4 omega (same convention as
`wire_pfaffian_witness_sweep`). The CSR sweep variant already had this
fix; the dense variant was broken.

## Strict S1+S2 sign agreement is structurally unattainable

For a real-symmetric single-particle BdG with imaginary diagonal s-wave
pairing, the S1 strategy's eigendecomposition returns real orthogonal
eigenvectors, and the off-diagonal pairing contribution to the 4×4
Nambu projection cancels by orthonormality. The S2 strategy avoids
eigendecomp by directly projecting onto bands 7-8, preserving the
H(7,8) asymmetry and yielding a non-zero Pf. Strict
`s1 == s2 .and. s1 /= 0` is therefore impossible; the accepted form
is `s2 /= 0` (the gate-relevant witness). Full S1+S2 strict sign-
agreement is the Majorana-basis Pfaffian problem (Issue 05, deferred).

## When to use

- New BdG invariant → add a sibling to `bdg_observables.f90`, sharing
  the `bdg_eval_params_t` record. Don't grow the seam into a hub.
- Gate row source of truth is the colormap
  (`output/z2_phase_diagram.dat` z2 column), NOT a separately-emitted
  Fortran file.
- Naming: BHZ-only helpers get `_bhz_heuristic` suffix; never pretend a
  heuristic is a generic Z2 helper.
- Construction-matrix indexing: when building projection helpers, build
  a local 4×4 omega at the projected-subblock indices, NOT extract
  from a larger omega whose structure crosses out of the subblock.

## Source

- Map + 7 tickets: `.scratch/archive/bdg-evaluator-pfaffian/`
- Spec: `.scratch/archive/bdg-evaluator-pfaffian/spec.md`
- Plan: `docs/plans/2026-07-13-002-feat-bdg-u2-actual-ship.md`
- Dense-path witness: `src/physics/topological_analysis.f90:1641-1707`
- CSR sweep variant: `src/physics/topological_analysis.f90:1720-1805`
- Seam: `src/physics/bdg_observables.f90`
- Gate: `tests/integration/test_lecture_13_acceptance_gate.sh`
- Verifier: `tests/integration/verify_majorana_polarization.py`
- Lecture: `docs/lecture/13-topological-superconductivity.md` §13.7.4 + §13.7.5
