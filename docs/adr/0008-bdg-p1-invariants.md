# 0008 — BdG/Majorana P1 Stabilization

**Date**: 2026-07-01
**Status**: ACCEPTED (locked decisions from brainstorming)
**Supersedes**: (none — extends ADR 0007)
**Related**: `docs/superpowers/specs/2026-07-01-bdg-p1-fix-design.md`

## Context

The 6-dimension adversarial review of PR #41 (2026-07-01) surfaced 18 P1 / 35 P2 / 19 P3 findings across physics correctness, test quality, architecture, build hygiene, documentation drift, and physics edge cases. Four architectural decisions must be locked before any implementation begins.

## Decisions

### 1. Kitaev Majorana number: Lutchyn-Oreg sign-of-det

The current `kitaev_majorana_number` evaluates `real(pf_val)` on a Pfaffian that is purely imaginary after `(B - B^T)/2` antisymmetrization of the anti-Hermitian product B = H·ω. The result is 0 for every input — silently marking all regimes as "gap closure". Lutchyn-Oreg sign-of-det bypasses the Pf sign extraction by:

```
H_bdg(k) at k=0 and k=π
  → sign(H_bdg) via eigendecomposition (Hermitian → real eigvals)
  → project to odd-parity Nambu subspace (rows/cols 1..n_sp/2)
  → det of odd-parity block
  → sign(det) ∈ {-1, +1}
  → product over k ∈ {0, π}
  → M = sign(product)
```

**Rationale**: literature-canonical, handles diagonal-in-(c,c†) BdG forms where Pf vanishes, makes strict AE1 certification possible. Detects non-Hermitian H_bdg via `error stop` (vs silent failure). Gap closure → M=0 remains legitimate.

### 2. pairing_sign: PUBLIC SSOT

The pairing_sign(8) = [+1,+1,-1,-1,+1,-1,+1,-1] table is duplicated in `bdg_hamiltonian.f90:60-72` and `topological_analysis.f90:876-878`. This DRY violation enabled the block (2,1) sign bug at `bdg_hamiltonian.f90:547`. Make `pairing_sign` and `pairing_partner` PUBLIC in `bdg_hamiltonian.f90`. `topological_analysis.f90` imports via `use bdg_hamiltonian, only: pairing_sign, pairing_partner`.

### 3. Symmetric Peierls on hole block — KEPT (reversed 2026-07-05)

**Reversed 2026-07-05**: Symmetric `add_peierls_coo(-B)` on hole block is KEPT (not removed). The `add_peierls_coo` function applies `exp(-iφ)` where `φ = e*Bx*(y_i-y_j)/hbar`. The `-conjg()` transform in `build_bdg_hole_block` does NOT include Peierls phase.

Empirical PHS oracle verification (rel_resid 1.25e-1 → 0 on removing the call) demonstrated the symmetric call is REQUIRED for class-D PHS at generic k with Bx≠0. The original 'double-count' rationale was disproved 2026-07-01; the investigation note lives at `src/physics/bdg_hamiltonian.f90:406-419` ("KEPT (NOT removed per Task 1.10): the symmetric add_peierls_coo(-B) IS required to make test_bdg_phs pass. The 'double-count' claim in ADR 0008 §3 was wrong").

Companion fixes: `docs/superpowers/specs/2026-07-01-bdg-p1-fix-design.md` line 23 corrected in this PR per spec §6.1; `docs/adr/0007-bdg-hole-block-canonical-convention.md` Layer D gets a Peierls-symmetry footnote per spec §6.2.

### 4. Minigap-pattern regression replaces z2-golden test

`tests/regression/data/wire_inas_gaas_bdg_topological_2d_phase.dat` has `z2=0` at all 10 grid points — same failure mode as the historical Rashba all-zero bug. Replace with minigap-pattern assertion: small minigap near B_crit, large minigap far from it. The Z2 transition is captured by the new `acceptance gate absolute window guard` (BCRIT_MIN ≥ 0.5 T, BCRIT_MAX ≤ 6.0 T for InAs wire).

## Consequences

- Phase 1.4 replaces the broken `real(Pf)` with Lutchyn-Oreg. Existing synthetic-BdG tests still pass; new strict tests pass on real fixtures.
- Phase 1.8 fixes dense-QW block (2,1) sign bug (`pairing_sign(pairing_partner(ib))` → `pairing_sign(ib)`).
- Phase 1.10 removes symmetric Peierls from wire builder; verify PHS oracle stays 4/4 green.
- Phase 2.1 replaces z2 golden with minigap-pattern regression.
- Phase 2.4 adds absolute window guard to the 4-witness acceptance gate.

## Risks

- Lutchyn-Oreg eigendecomposition requires LAPACK `zheev` (complex Hermitian). Already in `linalg.f90`. Cross-check in `pfaffian.f90` for n>12 fallback (Laplace reduction).
- Dense-QW B_crit resolution: B-grid step is 0.5 T → ±0.25 T resolution. Sufficient for the 4-witness gate at 2.0 T range; coarsening this in any future update is a breaking change.
