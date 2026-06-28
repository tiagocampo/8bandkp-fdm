# ADR 0007: BdG hole-block canonical convention is the QW-dense form `-conjg(H₀(-k))`

## Status

Accepted (2026-06-27, grill session for `2026-06-14-001-feat-bdg-majorana-validation-plan.md`)

## Context

The two BdG builders — the dense QW path (`build_bdg_hamiltonian_qw`) and the
CSR wire path (`build_bdg_hamiltonian_wire`) — construct the hole block
differently:

- **Wire CSR**: `H_hole = -H₀ᵀ(+k)` (transpose, no conjugation)
- **QW dense**: `H_hole = -conjg(H₀(-k))` (conjugation, evaluated at -k)

Both satisfy the class-D particle-hole-symmetry constraint
`C H_BdG C⁻¹ = −H_BdG` (verified numerically at generic k with Zeeman and
Peierls in `test_bdg_phs_at_finite_bx`, post-PR40) — so the two forms are
demonstrably equivalent. But the divergence is a maintenance hazard: a
change to one builder's hole-block construction can drift from the other
silently, with the divergence invisible to k=0-only tests.

The plan `2026-06-14-001-feat-bdg-majorana-validation-plan.md` identifies
this as U4 (R3): the hole block must be derived from the class-D PHS
constraint and unified across both builders. U1's 2026-06-15 probe found
the divergence is **real but tertiary** — the all-zero minigap headline
bug is μ + window (PR40 fixed), not the hole block. U4 is still required
for correctness/trust, but not on the critical path to killing the
all-zero.

## Decision

Adopt the **QW dense form `-conjg(H₀(-k))` as the canonical hole-block
convention**, in three layers:

1. **Pinning oracle (Layer A)** — `test_bdg_phs_at_finite_bx` (and any
   generic-k extensions) is the authoritative PHS oracle. Both builders
   must satisfy `C H_BdG C⁻¹ = −H_BdG` to numerical tolerance at generic
   k, with and without Zeeman and Peierls. The convention is documented in
   `src/physics/AGENTS.md` and in each builder's module header. No
   Hamilton-construction code is changed by this layer.

2. **Shared wrapper (Layer B)** — extract a single `build_bdg_hole_block`
   procedure in `src/physics/bdg_hamiltonian.f90` that both the dense QW
   builder and the CSR wire builder call. The wrapper implements the QW
   dense form `-conjg(H₀(-k))` for any caller. The wire CSR builder's
   inline `-H₀ᵀ(+k)` is replaced by the shared wrapper; the QW dense
   builder's existing construction is replaced by the shared wrapper.
   No polymorphic types (ADR 0001); this is procedures + data, not a
   class hierarchy.

3. **Canonical form choice (Layer D)** — inside the Layer B wrapper, the
   canonical form is **`H_hole = -conjg(H₀(-k))`**. Rationale: it is the
   k≠0-general form (matches Leijnse-Flensberg Eq. 38), the published
   derivation cold, and the dense-QW rung has fewer existing tests to
   update than the wire CSR rung. The wire CSR tests will be updated to
   the unified convention (Layer D is the "flip the wire convention to
   match the dense" decision; it is the implementation strategy inside
   Layer B, not a separate code change).

The Pfaffian structure matrix ω (U3) and the Sticlet polarization signs
(U6) are derived once for the 16N Nambu layout under this convention
(KTD7 in the plan), and any relative spin-sector sign in `P_M` is pinned
here, not assumed.

## Why

- **Generalization.** `-conjg(H₀(-k))` is correct at k≠0; `-H₀ᵀ(+k)` is its
  transpose-form equivalent only when the Hamiltonian is real-valued at the
  evaluation point (true at k=0 in the absence of Peierls phases, false
  generically). Choosing the k≠0-general form up front avoids future
  surprises when the wire rung is evaluated at non-zero kz (U13 deferred
  but pending).
- **Published derivation.** Leijnse-Flensberg Eq. 38 derives
  `-conjg(H₀(-k))` from the Nambu construction. Choosing the derived
  form is consistent with the codebase's other SSOT derivations (k.p
  block table in `hamiltonian_blocks.f90`, Bir-Pikus in
  `strain_solver.f90`, Zeeman in `magnetic_field.f90`).
- **Test locality.** The dense QW rung (U7 bisection control) has fewer
  BdG tests than the wire CSR rung; updating dense is cheaper than
  updating wire.
- **PHS oracle-driven, not a priori.** The convention is *pinned* by the
  numerical PHS check (Layer A), not chosen and then verified. The
  oracle is authoritative; the convention is whatever passes it.
- **ADR 0001 compliance.** Shared procedures, not polymorphic builder
  types.

## Consequences

- The wire CSR builder's existing hole-block construction is replaced by
  a call to the shared wrapper. The `-H₀ᵀ(+k)` form is removed from the
  codebase; the canonical form is `-conjg(H₀(-k))`.
- Existing convention-pinning pFUnit tests
  (`test_bdg_wire_bx_hole_block_negative_transpose`,
  `test_bdg_qw_particle_hole_nonzero_k`) are revised to the unified
  convention. The convention-agnostic `test_bdg_hermiticity_*` tests are
  unaffected.
- `src/physics/AGENTS.md` BdG Nambu contract line is updated from
  `-H₀ᵀ` to the canonical form. The pairing-embed
  (`pairing_partner`/`pairing_sign`, already module data) and the μ-shift
  are documented as shared operations.
- The Sticlet polarization spin-sector signs (U6) and the Pfaffian
  structure matrix ω (U3) are derived once for the 16N Nambu layout under
  this convention. They are co-pinned by the closed-form Kitaev harness
  (AE1) and the Lutchyn-Oreg sign witness on the wire rung (U13 deferred
  but applies when the wire Pfaffian lands).
- `bdg_hamiltonian.f90` carries six deprecated `stop 1` statements —
  migrate them to `error stop` opportunistically while the file is open,
  or defer explicitly. Approval-gated per `CLAUDE.md` Boundaries
  (Hamiltonian-construction code).

## Implementation Record

- **Date**: 2026-06-28
- **Branch / commit**: `feat/bdg-validation-pass2`
- **Layer B wrapper**: `src/physics/bdg_hamiltonian.f90::build_bdg_hole_block` (private; called by both `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw`)
- **Layer D canonical form**: `H_hole = -conjg(H₀(-k))` (the QW-dense form)
- **PHS oracle**: `tests/unit/test_bdg_phs.pf::test_hole_block_is_time_reversed_*` — 4/4 RED → 4/4 GREEN at generic k under all four field combinations (was RED pre-Issue-03, rel_resid ≈ 0.12578; now GREEN, rel_resid = 0.0)
- **Convention-pinning tests updated**: `test_bdg_wire_bx_hole_block_negative_transpose` (revised to canonical form via independent `ZB8bandGeneralized(-kz)` reference); `test_bdg_qw_particle_hole_nonzero_k` (unchanged, already canonical)
- **Cross-builder spectral identity**: not directly tested; a clean cross-builder identity test would require passing a single shared 8x8 H0 to both builders with the wire reduced to a single spatial point (deferred; was removed in Issue 02 Fix Round 1)
- **Side effects**:
  - All six `stop 1` statements in `bdg_hamiltonian.f90` migrated to `error stop`.
  - `Vz_delta = Vz_opt − Vz_cfg` double-counting guard preserved; explanatory comment added to prevent automated-review false positives.
  - BdG Hermiticity tests (`test_bdg_hermiticity_*`, `test_bdg_qw_hermitian`) GREEN with and without Peierls.
  - BdG regression tests (`regression_wire_bdg_strain_shift`, `regression_wire_bdg_topological`, `regression_wire_insb_gfactor`, `regression_wire_dense_sparse_consistency`) GREEN.
  - Strain validation tests (`strain_validation_wire`, `strain_validation_wire_quantitative`) GREEN.
- **Approval**: Tiago de Campos — "Approved: hole-block unification per ADR 0007 Layers B+D on 2026-06-28"

### Follow-up issues (filed as separate, scoped PRs)

1. **Standard `tau_x K` PHS tests in `test_bdg_phs.pf` (`test_phs_wire_*`)**: these tests verify `C H C⁻¹ = -H` with `C = tau_x K`. The canonical form `-conjg(H₀(-k))` does NOT satisfy this constraint at generic k with Peierls phases. Class-D PHS is `tau_x H^T(-k) tau_x = -H(k)` (requires k-transformation); the simple test was wrong from the start but happened to pass with the old wire convention. Follow-up: rewrite to use the correct class-D PHS operator.
2. **`test_bdg_hamiltonian::test_bdg_phs_at_finite_bx`**: same root cause as (1). Same fix applies.
3. **`test_krylov_snapshots::test_snapshot_wire_peierls`**: reference snapshot was generated with the OLD wire convention. Reference data needs regeneration in a follow-up.

## Alternatives Considered

### A. Pinning oracle only (no wrapper, no convention choice)

Document the divergence in `AGENTS.md`; leave both builders as-is; rely on
the U5 PHS oracle to prove current correctness.

**Pros**: No approval-gated code change; tiny cost.
**Cons**: Does not prevent future divergence (someone changes one builder,
not the other, and the convention drifts); does not satisfy the user's
stated goal of "making sure the BdG is correctly implemented" beyond
verification.

### B. Wire CSR convention wins (flip QW dense)

Adopt `-H₀ᵀ(+k)` as canonical; flip the dense QW builder.

**Pros**: Wire is the primary fixture for the ladder (U7, U8, U9, U10
all run on wire); minimizes disruption to the higher-traffic path.
**Cons**: Not the k≠0-general form; not the published derivation; more
dense-QW tests to update than wire-CSR tests (false premise — see
decision rationale).

### C. Both kept, marked by flag

Add a `hole_block_convention` flag; one builder applies a correction to
match the other; both kept in source.

**Pros**: Preserves both forms for reference.
**Cons**: Confusing; does not reduce surface area; adds complexity
without removing divergence. Rejected.

### D. Sign flip without shared wrapper

Flip the wire builder's convention independently, without extracting a
shared wrapper.

**Pros**: Minimal change.
**Cons**: Belt-and-suspenders redundancy with Layer B; bypasses the
shared contract. Subsumed by Layer B+D combined.

## Related

- `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`
  (U4, R3, KTD2, KTD7)
- `docs/brainstorms/2026-06-14-bdg-majorana-validation-requirements.md`
  (R3, R4)
- `src/physics/bdg_hamiltonian.f90` — both builders' hole-block
  constructions; site of the shared wrapper
- `src/physics/AGENTS.md` — BdG Nambu contract line
- `tests/unit/test_bdg_hamiltonian.pf` — convention-pinning tests
- `tests/unit/test_bdg_phs.pf` (planned) — generic-k PHS oracle
- Leijnse & Flensberg 2012, arXiv:1206.1736 — Eq. 38 derivation