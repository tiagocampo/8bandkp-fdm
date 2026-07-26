**Status**: COMPLETE (2026-07-13) — ticket 01 of `.scratch/archive/bdg-evaluator-pfaffian/`.

Type: grilling
Status: resolved
Claimed-by: Tiago (via wayfinder session 2026-07-13)
Blocked by: —
Resolved: 2026-07-13

# 01 — Evaluator API shape

## Question

What is the public API of the Pfaffian-plugged `eval_bdg_point` (and any siblings)?

Two shapes are on the table:

**(A) Single mode-discriminated function.** Extend `bdg_eval_params_t` with a `mode` field (`heuristic | slim_pfaffian | kitaev_majorana`); `eval_bdg_point` branches internally on it. One call site, one struct, one place to look. The `invariant_flag` field's meaning becomes mode-dependent.

**(B) Two sibling functions, one result struct.** Keep `eval_bdg_point(eigenvalues, params)` eigenvalues-only and heuristic-only. Add `eval_bdg_pfaffian_witness(eigenvectors, params)` returning the same `bdg_eval_result_t` with the `invariant_flag` driven by the Pfaffian (slim or Kitaev). The caller picks which to call. Two call sites, two functions, one struct.

ADR 0001 (no polymorphic types) rules out a third shape — a class hierarchy with virtual dispatch. Both A and B dispatch on a plain enum-tagged record or a procedure choice.

Decide:

- Which shape (A or B) — and if A, whether `mode` lives in `bdg_eval_params_t` or in a sibling record.
- Whether `invariant_flag` keeps its current `0/1` semantics or is renamed/repurposed (e.g., `pfaffian_sign ∈ {-1, +1, 0}`, where `0` means inconclusive / near-gap-closure).
- Whether the QW + Kitaev call sites share one entry point (a `eval_bdg_kitaev_majorana(eigenpairs, params)` sibling) or get folded into B's `eval_bdg_pfaffian_witness`.
- Whether the slim witness on the wire gets its own mode value or shares with the QW/Kitaev mode (likely shares — it's the same Pfaffian code path conceptually, only the input shape differs).

The decision here unblocks tickets 02, 03, 04, 05, 06 — every downstream decision depends on the API shape.

## Constraints

- ADR 0001: dispatch by enum-tagged plain record, never a class hierarchy.
- The seam must remain pure (`pure function`): no I/O, no allocations, no state.
- `bdg_eval_params_t`'s existing fields (`delta_0`, `near_zero_frac`) are public; their semantics must remain backward-compatible or be deprecated loudly.
- `main_topology.f90`'s 4 existing call sites must keep their lambda-style one-liners (or be refactored consistently).

## Out of scope for this ticket

- The slim witness's noise-floor tolerance (sub-decision of ticket 02).
- The exact witness criterion in the acceptance gate (ticket 05).
- Which call sites migrate (ticket 04).
- The retirement of `compute_z2_gap` (ticket 03).

## Answer

**Shape B — sibling functions, all in `bdg_observables.f90`.** Three functions, one shared `bdg_eval_params_t`, two shared result-type conventions (`s2_sign` / `majorana_number ∈ {-1, 0, +1}`, separately from the existing `bdg_eval_result_t`).

### The three siblings

```fortran
! UNCHANGED. Eigenvalues-only → minigap + heuristic invariant (0/1). Stays pure.
pure function eval_bdg_point(eigenvalues, params) result(r)

! NEW. CSR BdG matrix at one (B, μ) point → S2 Pfaffian sign.
! Wraps wire_pfaffian_witness_sweep from topological_analysis.f90:1705.
function eval_bdg_pfaffian_witness_csr(H_bdg_csr, n_full, params) result(s2_sign)
  type(csr_matrix), intent(in) :: H_bdg_csr
  integer, intent(in) :: n_full
  type(bdg_eval_params_t), intent(in) :: params
  integer :: s2_sign  ! ∈ {-1, 0, +1}; -1 = topological, +1 = trivial, 0 = inconclusive

! NEW. Stack of BdG matrices across a k-par sweep → Kitaev Majorana number.
! Wraps kitaev_majorana_number from pfaffian.f90:120.
function eval_bdg_kitaev_majorana(H_k_array, k_par_values, params) result(majorana_number)
  complex(kind=dp), intent(in) :: H_k_array(:,:,:)
  real(kind=dp), intent(in) :: k_par_values(:)
  type(bdg_eval_params_t), intent(in) :: params
  integer :: majorana_number  ! ∈ {-1, 0, +1}; same convention as s2_sign
```

### Sub-decisions

1. **Shape (A) vs (B): B.** The three operationally distinct signatures — eigenvalues-only / CSR-matrix / full k-array — are structurally incompatible with a single mode-discriminated function. A mode field would either force every caller to pass dummy matrices (`H_k_array(:,:,:)` for the eigenvalues-only case), or grow `bdg_eval_params_t` into a polymorphic-by-struct shape that ADR 0001 explicitly rejected. Sibling functions match each input shape cleanly and keep the existing seam's purity intact.

2. **`invariant_flag` keeps its current `0/1` semantics in `bdg_eval_result_t`.** Widening it to `{-1, 0, +1}` would break the existing 4 call sites (`main_topology.f90:530, 552, 826, 1360`) which read it as the heuristic "≥2 near-zero modes" flag and feed `result%n_majorana`. The new Pfaffian/Kitaev results return their own integers (`s2_sign`, `majorana_number`) — the caller's job is to map them onto `topological_result%z2_invariant` (with `−1 → 1`, `+1 → 0`, `0 → fallback to gap-closure heuristic`) at the call site, where the SC-minigap fallback already lives (`main_topology.f90:1377-1380`).

3. **QW + Kitaev gets its own sibling, not folded into the wire slim witness.** Different input shape (`H_k_array(:,:,:)` for k-par sweep vs. CSR matrix at one point), different mathematical primitive (Lutchyn-Oreg sign-of-det, ADR 0008 §1, vs. projected Pfaffian on bands 7-8 per k.p block table SSOT). Folding would force the wire caller to pass a degenerate 1-element k-array, which is wasteful and obscures the intent.

4. **Wire slim witness does NOT share a mode value with QW/Kitaev.** Different signatures, different return types, different code paths. The fact that both compute a Pfaffian sign is a math-level observation; at the API level they're distinct entry points. ADR 0001's spirit (no polymorphism, single-source-of-truth data) is satisfied by `bdg_observables.f90` re-exporting the underlying primitives as a thin facade — not by sharing a discriminator.

### `bdg_eval_params_t` is NOT extended

The `mode` field stays out. The existing fields (`delta_0`, `near_zero_frac`) stay exactly as-is. The new siblings take `bdg_eval_params_t` for symmetry (and to leave room for ticket 02's sub-decision on the slim witness's noise-floor tolerance), but only `delta_0` is consumed — `near_zero_frac` is ignored by the Pfaffian/Kitaev code paths. This keeps the 4 existing lambda-style one-liners (`bdg_eval_params_with_delta(cfg%bdg%delta_0)`) untouched and preserves backward compatibility per ticket constraint.

### Purity note

The new siblings cannot be `pure function` because `complex_pfaffian` (`pfaffian.f90:79-114`) allocates an internal `Asym` workspace. This is a soft loss vs. the existing `eval_bdg_point` seam, but the seam itself stays pure and the new siblings remain "no I/O, no allocations, no state" at the API level (internal allocations are encapsulated inside `complex_pfaffian`). ADR 0001's spirit is preserved: dispatch by enum-tagged plain record or sibling procedure choice, not a class hierarchy.

### `pure` placement

`eval_bdg_point`, `bdg_eval_params_with_delta`, `q_zero_tol` keep their `pure function` declarations verbatim. `eval_bdg_pfaffian_witness_csr` and `eval_bdg_kitaev_majorana` are plain `function` (no `pure`). This matches the existing purity of the underlying primitives they wrap (`wire_pfaffian_witness_sweep` is a `subroutine`, `kitaev_majorana_number` is a non-pure `function`).

### What this unlocks

- **Ticket 02** can now place the slim witness as a thin re-implementation (or copy) of `wire_pfaffian_witness_sweep`'s S2 path inside `bdg_observables.f90` — the seam becomes the single landing zone.
- **Ticket 03** has a clear path: retire `compute_z2_gap` / `compute_z2_gap_edge` by removing them from the wire-rung code paths that now route through `eval_bdg_pfaffian_witness_csr`, while keeping `compute_z2_gap_sweep` (BHZ-only) in scope only if its own grilling resolves it.
- **Ticket 04** has a fixed audit target: 4 call sites for `eval_bdg_point` (unchanged) + 1 wire rung site for `eval_bdg_pfaffian_witness_csr` (new) + N QW/Kitaev sites for `eval_bdg_kitaev_majorana` (new, depending on rungs).
- **Ticket 05** can wire the gate row directly: `wire rung invariant discriminator = eval_bdg_pfaffian_witness_csr(...)`.
- **Ticket 06** has a stable API surface to write tests against.

### Resolution

- Shape: **B** (sibling functions).
- `invariant_flag`: keeps `0/1` semantics in `bdg_eval_result_t`.
- QW + Kitaev: own sibling (`eval_bdg_kitaev_majorana`).
- Wire slim + QW/Kitaev: do not share a mode value; distinct siblings.
- Location: all three live in `bdg_observables.f90`. The new siblings are thin facades over the existing primitives (`wire_pfaffian_witness_sweep`, `kitaev_majorana_number`); the primitives stay in their current modules to avoid dep-graph reshuffles.
- `bdg_eval_params_t`: unchanged.