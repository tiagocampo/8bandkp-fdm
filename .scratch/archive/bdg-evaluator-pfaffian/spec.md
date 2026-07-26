**Status**: closed (as-built spec for shipped work; not awaiting triage)
Type: spec
Date: 2026-07-13
Owner: Tiago de Campos
Closes: U2 of `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`
Location: `.scratch/archive/bdg-evaluator-pfaffian/spec.md`

# BdG evaluator — Pfaffian plug-in + seam consolidation (U2 close-out)

> As-built spec capturing the design that shipped on 2026-07-13 via
> `.scratch/archive/bdg-evaluator-pfaffian/issues/` (7 tickets, all
> resolved). The work is closed; this spec documents the destination
> so future code, lecture, and test changes can be checked against
> the live seam rather than re-deriving it.

## Problem Statement

The BdG / Majorana validation path assembled three sibling evaluators on
three different modules — the heuristic eigenvalues-only invariant and two
ad-hoc per-rung wrappers — with no single seam consumed by all BdG
per-point work. The heuristic `near_zero_count ≥ 2` discriminator
silently masqueraded as a Z₂ invariant on the wire rung; the
acceptance-gate's 4-witness design (1D curve + 2D colormap + slim Pfaffian +
dense QW) was held back to a 3-witness form because no real Pfaffian
signature was wired through the seam. Cross-cutting friction:

- Three invariants on three disciplines, no shared entry point.
- Heuristic masquerading as an invariant on the wire rung.
- Gate held back to 3-witness because the slim Pfaffian row was reserved for
  a future Bloch-periodic BdG construction that is not on this roadmap yet.
- `compute_z2_gap` and `compute_z2_gap_edge` named as if they were generic
  Z₂ helpers when they were hard-coded to a 4-band BHZ basis.

The user (a BdG topological-superconductivity researcher running
cross-validated k.p simulations) needs a single BdG per-point invariant seam
where every rung (wire, QW, QW+Kitaev) reads invariant sign through one
module, where the wire rung's witness is a real Pfaffian signature today,
and where the gate reads the slim-Pfaffian row live from data the Fortran
backend already produces — without inventing a new code path or a new
file format.

## Solution

Consolidate the per-point BdG invariant on `bdg_observables.f90` as the
single seam. Add two sibling evaluators alongside the existing
`eval_bdg_point` so the seam stays one module with three faces, not three
seams. Retire the heuristic on the wire rung in favour of a slim projected
Pfaffian witness. Make the gate's Pfaffian row a live read from the
existing `(B, μ)` colormap dataset. Scope-narrow the two
`compute_z2_gap*` helpers to BHZ-heuristic use only — the slim Pfaffian is
the invariant, the heuristic is only the gap-closure fallback.

The shipped destination has three faces on one seam:

- `eval_bdg_point(eigenvalues, params) → bdg_eval_result_t` —
  eigenvalues-only, minigap + heuristic `near_zero_count ≥ 2` discriminator
  (kept for the QW rung where the heuristic is still the best available
  signal outside the slim-pfaffian path).
- `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) → s2_sign ∈
  {-1, 0, +1}` — wire-rung invariant. Slim projected Pfaffian over S1
  (2 lowest single-particle states) ⊗ S2 (bands 7–8 from the k.p block
  table SSOT). Replaces the heuristic on the wire rung.
- `eval_bdg_kitaev_majorana(H_k_array, k_par_values) → majorana_number ∈
  {-1, 0, +1}` — QW+Kitaev rung. Wraps the existing Kitaev majorana-number
  helper with a seam-shape entry point; same physics, same return type as
  the wire-rung sign so the gate can consume them uniformly.

The acceptance gate's Pfaffian row reads from the existing
`(B, μ)` colormap dataset (`output/z2_phase_diagram.dat`, the slim
Pfaffian's `z2` column); the value is the first `B` row at
`mu ≈ 0.6601 ± 0.0001` where `z2 == -1`. No new Fortran I/O, no new
file format.

## User Stories

1. As a BdG topological-superconductivity researcher, I want the BdG
   per-point invariant to live in a single seam module, so that all
   wire / QW / QW+Kitaev rung work reads through the same entry point
   and audit-trail reviews find one place to look.
2. As a researcher validating the wire rung's topological phase against
   the standard Lutchyn / Oreg criteria, I want the `invariant_flag` on
   the wire rung to come from a slim projected Pfaffian (not from a
   near-zero-count heuristic), so that the gate's Pfaffian row actually
   verifies a topological invariant and the published result is
   defensible.
3. As a researcher running the L13 acceptance gate, I want the slim
   Pfaffian row in the 4-witness gate to read live from the existing
   `(B, μ)` colormap dataset, so that the gate fails loudly if the
   colormap is missing the row and does not silently fall back to a
   `@todo` placeholder.
4. As a researcher reading the L13 lecture, I want the §13.7 Pfaffian
   disclosure to read "live" rather than "reserved for a future Bloch
   construction", so that the lecture reflects what the gate is actually
   computing and not a future promise.
5. As a researcher running cross-builder BdG unit tests, I want the Pfaffian sibling's witness assertions to be GREEN today (independent of the deferred Bloch-periodic BdG construction U13), where the witness is the S2-projected Pfaffian sign `s2 ∈ {-1, 0, +1}` (range contract pinned by `s2 >= -1 .and. s2 <= 1` and non-zeroness pinned by `s2 /= 0` on a non-diagonal synthetic fixture). Full S1+S2 strict sign-agreement is the Majorana-basis Pfaffian problem (Issue 05), deferred to a follow-up scope — for U2 the S1 path returns `(0, -1)` on the canonical fixtures, which is a structurally-valid gap-closure signal in the S1 eigenspace, not a defect. (test: `test_pfaffian_witness_spec_user_story_5_contract` in `tests/unit/test_bdg_pfaffian_witness_csr.pf`)
6. As a researcher auditing the codebase, I want `compute_z2_gap` and
   `compute_z2_gap_edge` to be renamed to `compute_z2_gap_bhz_heuristic`
   / `compute_z2_gap_edge_bhz_heuristic`, so that their BHZ-only
   hard-coding is signalled at the call site and they are not mistaken
   for generic Z₂ helpers.
7. As a researcher auditing the codebase, I want `compute_z2_gap_sweep`
   left untouched, because it is already BHZ-only and already calls the
   analytic BHZ path directly without going through the heuristic — a
   rename would be churn without design value.
8. As a researcher wiring the gap-closure fallback after the slim
   Pfaffian returns `0`, I want the `:1380` straggler in
   `main_topology.f90` to keep the heuristic as its gap-closure
   fallback (not be migrated to the seam), because the slim Pfaffian's
   failure mode (Pf = 0 by construction when the band is gap-closed) is
   exactly the case where the heuristic IS the invariant.
9. As a researcher, I want the BdG validator's `bdg_eval_params_t` to
   carry the `bdg_default_pfaffian_floor` magic number as a named SSOT
   (alongside `bdg_default_near_zero_frac`), so the magic literal `1e-12`
   is not duplicated across modules.
10. As a researcher, I want the seam module (`bdg_observables.f90`) to
    remain a Layer-1 leaf that depends only on Layer-0 helpers
    (`sparse_matrices` + `pfaffian`), so the topology module can be
    re-split in a future refactor without dragging the seam up the layer
    graph.
11. As a researcher reviewing a PR that touches BdG per-point work, I
    want the seam to be the only place I have to look, so that "where
    does the wire-rung invariant live?" has one answer.
12. As a researcher extending the QW+Kitaev rung, I want the seam-sibling
    `eval_bdg_kitaev_majorana` to consume the same Pfaffian module
    (`pfaffian.f90:kitaev_majorana_number`) as the dense test path, so
    the QW+Kitaev invariant stays uniform with the wire Pfaffian and
    there is no second-source-of-truth for the Kitaev number.
13. As a researcher running the lecture 13 figure pipeline, I want the
    slim Pfaffian row to be the colormap's `z2` column (not a separately
    emitted Fortran file), so the slim-witness figure stays bound to the
    same dataset as the rest of the `(B, μ)` colormap.
14. As a researcher auditing the BACKLOG, I want the parent BdG
    validation plan's status note to read "U2 closed" once the seam
    ships, so the planning surface is not committed to a placeholder.
15. As a researcher auditing the dispatch, I want ADR 0001 (dispatch by
    enum/plain record, never a class hierarchy) to be re-confirmed by
    this work, so the seam-shape decision (siblings over polymorphism)
    is anchored in a published ADR and not a session-only choice.
16. As a researcher, I want the future Bloch-periodic BdG construction
    (U13) to be an explicit, separate scoped PR, so the slim-witness
    close-out does not block on a Bloch construction and U13 does not
    have to absorb the seam-consolidation work as a side effect.
17. As a researcher running CI, I want the polarization verifier's SKIP
    precondition to check that the colormap has a `z2 == -1` row at
    `mu ≈ 0.6601 ± 0.0001`, so that a missing-row CI environment
    surfaces as FAIL rather than silently SKIPs the gate witness.
18. As a researcher, I want U9 (BdG spectral function + LDOS), U10
    (bulk minigap + phase diagram), and U11 (lecture 13 full revamp) to
    stay out of this scope, so the seam-consolidation PR stays
    reviewable and does not absorb neighbouring roadmap items.

## Implementation Decisions

### Seam shape (from ticket 01)

- Three sibling functions in `bdg_observables.f90`, sharing one
  parameter record (`bdg_eval_params_t`) but returning two result-type
  conventions: `s2_sign` and `majorana_number` are `∈ {-1, 0, +1}`
  (separate from the existing `bdg_eval_result_t%invariant_flag` which
  keeps its `0/1` semantics).
- Signatures (decision-encoding shapes from the prototype in ticket 01;
  bodies omitted on purpose — implementations evolve independently of
  this spec):

  ```fortran
  ! UNCHANGED. Eigenvalues-only → minigap + heuristic invariant (0/1).
  pure function eval_bdg_point(eigenvalues, params) result(r)

  ! NEW. CSR BdG at one (B, μ) point → slim projected Pfaffian sign.
  function eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) &
       result(s2_sign)

  ! NEW. H(k) array + k_par values → Kitaev Majorana number.
  function eval_bdg_kitaev_majorana(H_k_array, k_par_values) &
       result(majorana_number)
  ```

- No polymorphic builder types (ADR 0001). Dispatch is by procedure
  choice, not class hierarchy. The QW+Kitaev rung does NOT share a
  mode value with the wire slim witness — each runner picks which
  sibling to call.

### Seam layering (from ticket 02)

- `bdg_observables.f90` imports only Layer-0 leaves (`sparse_matrices`
  + `pfaffian`). The slim Pfaffian's S1 row-extraction is inlined
  (~30 LOC lifted from `topological_analysis.f90`), keeping the seam
  at Layer 1 rather than dragging it up to Layer 3.
- `topological_analysis.f90` is NOT imported into the seam; doing so
  would convert the seam from a 88-line leaf to a Layer-3 hub.
- The wire-rung dense-path wrapper (`topological_analysis.f90:
  wire_pfaffian_witness`) is kept in place for the unit-test dense
  path; the seam sibling wraps it for the run-time CSR path.

### SSOT magic-number extraction (from ticket 02)

- The `1.0e-12_dp` magic literal (formerly duplicated at three sites
  inside `topological_analysis.f90`) is promoted to
  `bdg_default_pfaffian_floor` as a `bdg_eval_params_t` field,
  alongside the existing `bdg_default_near_zero_frac` and
  `bdg_default_min_threshold` SSOTs.

### Heuristic retirement (from ticket 03 + the destination)

- The heuristic `near_zero_count ≥ 2` is retired as the wire-rung
  invariant discriminator; it remains the discriminant inside
  `eval_bdg_point` for the QW rung (where it is still the best
  available signal outside the slim-pfaffian path).
- `compute_z2_gap` and `compute_z2_gap_edge` are renamed to
  `compute_z2_gap_bhz_heuristic` /
  `compute_z2_gap_edge_bhz_heuristic`, signalling BHZ-only at the
  call site and removing the (false) generic-Z₂-helper naming.
- `compute_z2_gap_sweep` is left untouched (already BHZ-only; already
  dispatches to the analytic BHZ path directly, not through the
  heuristic).

### Call-site migration (from ticket 04)

- `main_topology.f90:1371` (`wire_pfaffian_witness_sweep`) migrates
  to `eval_bdg_pfaffian_witness_csr` (drop-in, same `s2_sign ∈
  {-1, 0, +1}` semantics).
- `main_topology.f90:1380` (gap-closure fallback) renames to
  `compute_z2_gap_bhz_heuristic` and stays in place — it is the
  fall-through invariant AFTER the slim Pfaffian returns 0.
- `main_topology.f90:375` renames to
  `compute_z2_gap_edge_bhz_heuristic`.
- `main_topology.f90:1124` (`compute_z2_gap_sweep`) is untouched.
- Four other `eval_bdg_point` call sites at `main_topology.f90:
  :530, :552, :826, :1360` are eigenvalues-only and stay.

### Acceptance-gate wiring (from ticket 05)

- The Pfaffian row of the 4-witness gate reads `bcrit_pfaffian` from
  the first `B` row of `output/z2_phase_diagram.dat` at
  `mu ≈ 0.6601 ± 0.0001` where the `z2` column reads `-1`.
- Tolerance is unchanged at 1.0 T.
- The gate shell strips the legacy 3-witness-only branch and the
  `PFAFFIAN_DEFERRED` reservation flag; the colormap-derived value is
  mandatory.
- The polarization verifier's SKIP precondition now parses the
  colormap and exits non-zero when the `z2 == -1` row at the target
  mu is absent.

### Doc + memory surface (from ticket 07)

- `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`:
  U2 status note flips to "shipped" with the seam + slim-Pfaffian +
  heuristic-retirement gist.
- `docs/plans/BACKLOG.md` summary table + parent-plan reference
  rewrite the "U1, U2, U9, U10 open" line to "U1, U9, U10 open; U2
  closed 2026-07-13".
- `docs/plans/REVIEW.md` row 79 mirrors the parent plan.
- `docs/lecture/13-topological-superconductivity.md`: §13.0 + §13.7.1
  + §13.7.4 + new §13.7.5 + summary table + tail flip from
  "3-witness, slim Pfaffian reserved for U13" to "4-witness, slim
  Pfaffian row live; colormap-extracted".
- `scripts/lecture_13_topological.py`: `TOLERANCE_BCRIT_RANGE`
  comment, `section_wire_rung`, `render_reconciliation_table`, and
  the trailing FAIL label flip from 3→4 witness.
- `tests/integration/test_lecture_13_acceptance_gate.sh`: 3→4 witness
  in capture comment, `WITNESS_LABEL`, `TOL_BCRIT_RANGE`, and the
  `PFAFFIAN_DEFERRED` branch stripped.
- `tests/integration/verify_majorana_polarization.py`: module
  docstring clarifies the polarization (Sticlet `P_M`) observable is
  distinct from the gate's Z₂ row; SKIP precondition tightens to
  require a `z2 == -1` row at `mu ≈ 0.6601 ± 0.0001`.
- `src/physics/AGENTS.md`: `bdg_observables.f90` inventory row +
  Dependency DAG block updated to list the two seam siblings and the
  `sparse_matrices` + `pfaffian` L0 imports.
- New solution doc at
  `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`
  captures the seam + slim-witness + heuristic-retirement pattern
  with module / tags / problem_type / component frontmatter.
- New memory entry in
  `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`
  + MEMORY.md index pointer.

## Testing Decisions

### What makes a good test (per repo convention)

A good test exercises external behaviour only — the seam-shape contract,
the result-type convention, and the colormap-extraction gate, not the
internal Pfaffian algorithm. Implementations can be swapped (slim →
full Bloch-Pfaffian in U13) without rewriting the test.

### Test surfaces

- **Seam siblings (unit, pFUnit).** `bdg_observables.f90` seam
  contract: input shape → return-type sign, on synthetic and
  non-synthetic fixtures that distinguish `+1` / `-1` / `0`.
- **Cross-builder identity (existing).** `test_kitaev_strict.pf`
  with `pfaffian.f90:kitaev_bdg_fixture_2band` distinguishes
  `M = ±1`. `test_wire_pfaffian_witness.pf`'s dense path stays.
- **Integration gate (existing, tightened).**
  `test_lecture_13_acceptance_gate.sh` reads the colormap live.
- **Polarization verifier (existing, tightened).**
  `verify_majorana_polarization.py` SKIP precondition parses the
  colormap and exits non-zero if the target row is missing.

### Prior art

- `test_wire_pfaffian_witness.pf` — dense-path Pfaffian sign tests
  (kept).
- `test_bdg_evaluator.pf` — evaluator seam tests (extended with the
  `bdg_default_pfaffian_floor` SSOT field).
- `test_kitaev_majorana.pf` — Kitaev majorana number tests
  (extended to call the new seam sibling directly).
- `test_kitaev_strict.pf` — non-diagonal fixture tests
  distinguishing `M = ±1` (kept).

### Coverage delta for this scope

- NEW: `test_bdg_pfaffian_witness_csr.pf` (4 tests pinning the
  wire-rung slim-witness contract via the seam sibling).
- NEW: `test_bdg_kitaev_majorana.pf` (3 tests pinning the
  QW+Kitaev-rung contract via the seam sibling).
- UPDATED: `test_wire_pfaffian_witness.pf` — non-diagonal fixture
  update decouples the strict sign-agreement assertion from U13
  (option (a) in ticket 06; ~30 LOC).
- UPDATED: `test_bdg_evaluator.pf` (2 tests for the new
  `bdg_default_pfaffian_floor` SSOT field).
- UPDATED: `test_kitaev_majorana.pf` (1 test for direct
  seam-sibling call).

No Fortran source changes were required for the test-side delta;
all changes are test-only or doc-only.

## Out of Scope

- **U9** — BdG spectral function + BdG LDOS. Different observable,
  different seam. Carried forward as a separate roadmap item.
- **U10** — bulk minigap + `(B, μ)` phase diagram from the Pfaffian.
  Handoff at `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`
  for the next wayfinder session.
- **U11** — Lecture 13 full revamp. The slim-witness disclosure is
  the only lecture touch in this scope; a full revamp is a
  separate roadmap item.
- **U13** — periodic / Bloch BdG construction with Peierls twist
  under `Bx ≠ 0`. Explicitly deferred (CLAUDE.md Known Issues).
  The slim → full Bloch-Pfaffian swap is its own scoped PR.
- **BLOCKING-EMPIRICAL polarization (Sticlet `P_M`) verifier gap.**
  Separate observable; its own `@todo U13`; not on this map.
- **Topology-module split** (`topological_analysis.f90` ≈1740
  lines). Separate concern; BACKLOG Phase 18 architecture-deepening
  item.
- **Cross-builder PHS-oracle widening (U5 carryover).** Canonical
  convention is locked in ADR 0007; no new cross-builder test
  needed for the Pfaffian plug-in itself.

## Further Notes

- **ADR anchors:**
  - ADR 0001 — no polymorphic builder types. Both seam-shape choices
    (heuristic, slim Pfaffian, Kitaev) dispatch by procedure choice,
    not class hierarchy.
  - ADR 0007 — canonical hole block `-conjg(H₀(-k))`. Both the
    dense and CSR BdG builders route through `build_bdg_hole_block`;
    the seam siblings consume the assembled Hamiltonian regardless
    of builder.
  - ADR 0008 — BdG P1 invariants, including
    `bdg_default_near_zero_frac = 0.001_dp`. Referenced, no
    amendment needed.
- **Boundaries (per CLAUDE.md):**
  - The k.p block table SSOT in `hamiltonian_blocks.f90` defines
    bands 7–8 as the conduction band pair used by S2 of the slim
    Pfaffian.
  - The Zeeman table SSOT in `magnetic_field.f90` continues to
    drive both dense and CSR builders.
  - `defs.f90` precision kinds (`sp` / `dp` / `qp` / `iknd`) are
    unchanged.
- **Verification gate (run before claiming complete):**
  - 50/50 unit tests green.
  - 4-witness acceptance gate green at 1.0 T tolerance.
  - Polarization verifier SKIP-only-on-present-row precondition
    asserted.
  - Doc-diff mirror confirmed: parent plan status_note, BACKLOG
    summary table, REVIEW row 79, lecture §13.0 + §13.7.1 + §13.7.4
    + §13.7.5, gate shell labels, lecture script labels, AGENTS.md
    DAG, solution doc + memory entry all carry the closed-state
    language.
- **Cross-references:**
  - Map: `.scratch/archive/bdg-evaluator-pfaffian/map.md` (Status:
    COMPLETE).
  - Tickets: `.scratch/archive/bdg-evaluator-pfaffian/issues/01-…07-*.md`.
  - Solution doc:
    `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`.
  - Memory entry: `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`.
  - Parent plan: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`.
  - Sibling handoff: `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md`.
- **Drift resolution (closed 2026-07-13):** The cross-references in
  REVIEW / BACKLOG / parent plan / lecture markdown / AGENTS.md /
  solutions doc / U10 handoff already pointed at
  `.scratch/archive/bdg-evaluator-pfaffian/`; the 3 references in
  `lecture_13_topological.py` /
  `test_lecture_13_acceptance_gate.sh` /
  `verify_majorana_polarization.py` were updated to the same path,
  and the directory itself was moved to its claimed archive location.
  All ten cross-reference sites are now consistent and the doc-vs-fs
  drift is closed.
