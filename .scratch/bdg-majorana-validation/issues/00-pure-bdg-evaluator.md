# Issue 02 — PHS oracle catches hole-block divergence at generic k (U5, Layer A of U4)

> **File-numbering note**: This file is `00-pure-bdg-evaluator.md` in `.scratch/bdg-majorana-validation/issues/`. The content is **Issue 02** (PHS oracle) — sourced from `.superpowers/sdd/issue-00-brief.md`. The brief filenames were off by 2 from the issue numbers; this corrected numbering aligns the file's content with its issue number in the dependency graph.

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U5, Layer A of U4)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-00-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md`
**Phase**: PR-A, parallel with Issue 01 (Kitaev harness) — fan-out

## What to build

A numerical particle-hole-symmetry oracle that pins the BdG hole-block convention. Two pure checks, both running at generic `k` (e.g., `k=π/2a`), not just `k=0` — because the wire/QW hole-block divergence is invisible at `k=0`.

1. **PHS check:** `‖C H_BdG C⁻¹ + H_BdG‖ < tol` over the square intersection, reusing the clamped-loop structure of `csr_hermitian_error` with `conjg → C-conjugation`. The particle-hole operator `C = U_C·K` is built from the existing Kramers pairing. Runs with and without Zeeman; with and without Peierls. The four-field combinations are the canonical test matrix.
2. **Cross-builder consistency:** CSR and dense hole blocks yield identical spectra for the same input. The two spectra are derived from independent constructions (not the same constants — the tautology anti-pattern).

The teeth-first discipline: this oracle must **fail on the current divergent convention** at generic k. Writing it to fail before the hole-block unification proves it would have caught the bug; writing it to pass after the unification proves the new convention.

## Acceptance criteria

- [ ] `tests/unit/test_bdg_phs.pf` (new) green on the unified convention (after Issue 03 lands) at generic k under all four field combinations (no Zeeman / Zeeman; no Peierls / Peierls).
- [ ] **Teeth demonstrated:** BEFORE Issue 03 lands, the PHS check fails on the current divergent convention at generic k. Capture the failing diff in the issue report (record the failure numerically).
- [ ] Cross-builder spectra identical for the same `H₀` input (independent constructions, not tautological).
- [ ] PHS oracle reports a tolerance, not a pass/fail boolean — the report shows `‖C H_BdG C⁻¹ + H_BdG‖` numerically so future regressions are debuggable.
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

### `src/physics/bdg_hamiltonian.f90` (433 lines)

- Public: `build_bdg_hamiltonian_1d` (lines 90–308; **hole block at lines 264–283 — uses `-H₀ᵀ(+k)` form**, the wire-CSR convention); `build_bdg_hamiltonian_qw` (lines 330–431; **hole block at lines 420–424 — uses `-conjg(H₀(-k))` form**, the dense-QW convention).
- Module data: `pairing_partner(8) = [4,3,2,1,6,5,8,7]`, `pairing_sign(8) = [+1,+1,-1,-1,+1,-6,5,8,7,-1,+1,-1]`.
- Peierls phase on electron block (line 174–179).
- `Vz_delta` double-counting guard at lines 182–224 (preserve; do NOT fix).

### Existing tests at `tests/unit/test_bdg_hamiltonian.pf`

- `test_bdg_wire_bx_hole_block_negative_transpose` (lines 446–507) — current wire convention: `H(n8+j, n8+i) + H(i,j) < 1e-10`.
- `test_bdg_qw_particle_hole_nonzero_k` (lines 601–649) — current QW convention: `H(n8+row, n8+col) + conjg(H_minus(row, col)) < 1e-12`.
- `test_bdg_phs_at_finite_bx` (lines 930–1046) — existing PHS oracle (not currently gated as authoritative; Issue 02's job to elevate it).
- Convention-agnostic: `test_bdg_hermiticity_*`, `test_bdg_dimension_doubling`, etc.

### Kramers pairing (BdG particle-hole operator C)

The C operator that satisfies `C H_BdG C⁻¹ = -H_BdG` is built from the `pairing_partner`/`pairing_sign` module data in `bdg_hamiltonian.f90`. The 16N×16N BdG matrix in Nambu space has C represented as a permutation-with-signs matrix that maps the electron sector index `(band-1, site)` to the hole sector index `(pairing_partner(band)-1, site)` with sign `pairing_sign(band)`.

## Constraints from CLAUDE.md + ADRs

- **ADR 0001 (fat derived type)**: oracle is a pure procedure inside the pFUnit test file — no polymorphic types, no new module.
- **ADR 0007 (hole-block canonical convention)**: the oracle pins the canonical form chosen by Issue 03 (`-conjg(H₀(-k))`).
- **CLAUDE.md Boundaries**: Layer A (pinning oracle) is documentation/witness only — **does NOT modify Hamilton-construction code**. The convention is verified numerically here; the change is Issue 03's job.
- **CLAUDE.md Code Conventions**: F2018; `pure function` for the oracle core; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.
- **CLAUDE.md Conventions re hermiticity**: reuse the clamped-loop structure of `csr_hermitian_error` with `conjg → C-conjugation` (analogous pattern).

## File ownership (exhaustive)

### New files (you create)
- `tests/unit/test_bdg_phs.pf` — the PHS oracle test. The pure PHS-check procedure lives inside this file (a `module test_phs_oracle` or a `contains`-block procedure).

### Modified files
- `src/physics/bdg_hamiltonian.f90` — module-header doc-comment update ONLY (no Hamilton-construction code change). Add a line stating the canonical convention: "Hole block = -conjg(H₀(-k)) per ADR 0007; convention pinned by tests/unit/test_bdg_phs.pf".
- `tests/CMakeLists.txt` — `add_pfunit_ctest(test_bdg_phs ...)` with LABEL "unit".

### NOT modified (out of scope)
- The two builders' hole-block construction code (Issue 03 owns that change).
- `src/physics/topological_analysis.f90` (Issue 04/07 own).
- `src/apps/main_topology.f90` (Issue 03/05/06/07 own).
- The convention-pinning tests `test_bdg_wire_bx_hole_block_negative_transpose` and `test_bdg_qw_particle_hole_nonzero_k` (Issue 03 will update these).

## Suggested test layout (illustrative)

```fortran
module test_phs_oracle
  use pFUnit
  use kind_module, only: dp
  use hamiltonian_wire, only: bdg_hamiltonian_wire_t  ! if applicable
  use bdg_hamiltonian, only: pairing_partner, pairing_sign, build_bdg_hamiltonian_1d, &
                             build_bdg_hamiltonian_qw
  implicit none
contains

  ! Build the C operator as a 16N x 16N matrix from pairing_partner/pairing_sign
  pure function build_phs_operator(N) result(C)
    integer, intent(in) :: N
    complex(dp) :: C(16*N, 16*N)
  end function

  ! PHS residual: ||C H_BdG C^{-1} + H_BdG||
  pure function phs_residual(H_bdg, C) result(res)
    complex(dp), intent(in) :: H_bdg(:,:), C(:,:)
    real(dp) :: res
  end function

  @test
  subroutine test_phs_at_finite_k_no_zeeman_no_peierls(...)
    ! generic k = pi/(2*a) [or whatever finite-k your wire/QW fixture uses]
    ! build H_bdg with build_bdg_hamiltonian_1d
    ! assert phs_residual < 1e-10_dp
    ! (NOTE: this WILL FAIL on the current divergent convention — record the failure as evidence)
  end subroutine

  ! ... three more tests for the four-field matrix: {no_Z, Z} x {no_P, P}

  @test
  subroutine test_cross_builder_spectral_identity(...)
    ! Build same H_0 input to both builders independently
    ! assert sort(eigvals_1d) == sort(eigvals_qw) within tolerance
  end subroutine
end module
```

## TDD discipline (mandatory)

1. Write `tests/unit/test_bdg_phs.pf` FIRST.
2. Run the test on the **current** (pre-Issue-03) code — confirm at least ONE of the four-field tests FAILS at generic k (the "teeth demonstrated" requirement). Capture the failure numerically.
3. Document the failure in the report file as evidence (this is the regression baseline — Issue 03's fix must turn it GREEN).
4. Document in the report that the test will pass under the unified convention (after Issue 03 lands). The current code has the divergent convention; this test set is intentionally teeth-first.

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
```

## Out of scope (must NOT do)

- Do NOT modify the hole-block construction code (Issue 03 owns).
- Do NOT modify `src/physics/bdg_hamiltonian.f90`'s Hamiltonian-builder code (only the module-header doc comment).
- Do NOT modify `src/apps/main_topology.f90`.
- Do NOT modify `tests/unit/test_bdg_hamiltonian.pf`.
- Do NOT add `sweep_model` enum value.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-00-report.md`

The report must include:
- Status (DONE_WITH_CONCERNS is acceptable here because the tests will RED until Issue 03 lands)
- Files created / modified
- **Teeth-demonstrated evidence**: the actual numerical failure on the current divergent convention. Capture this BEFORE marking done.
- The four-field test matrix result (each combination: residual norm numerically).
- Cross-builder spectral identity test result.
- Self-review findings
- Commit SHA
- Explicit note: "These tests will be GREEN after Issue 03 lands per ADR 0007; current RED is the regression baseline."

## Risk notes

- The PHS operator C depends on the Nambu ordering. If you build C from `pairing_partner`/`pairing_sign`, verify the sign convention matches ADR 0007's canonical form (`-conjg(H₀(-k))`).
- "Generic k" must be non-zero. `k = π/2a` (or `0.1/Å` for arbitrary lattices) is a good choice — it ensures the divergence between `-H₀ᵀ(+k)` and `-conjg(H₀(-k))` is visible (they agree at k=0 when H₀ is real).
- The cross-builder spectral identity test must use INDEPENDENT constructions — not the same constants. If you accidentally route both builders through a shared helper, the test is tautological and the teeth are lost.

---

## Outcome (as executed)

- **RED pre-Issue-03**: 4 canonical hole-block tests RED with rel_resid ≈ 0.12578.
- **GREEN post-Issue-03**: 4/4 GREEN with rel_resid = 0.0.
- Cross-builder identity test: deferred (removed in Issue 02 Fix Round 1); would require shared H₀ → both builders with wire reduced to single spatial point.
- Standard `tau_x K` PHS tests (`test_phs_wire_*`): rewrite required per Issue 03 fix2 (root cause: canonical form `-conjg(H₀(-k))` does NOT satisfy `tau_x K PHS` at generic k with Peierls; class-D PHS is `C H(k, B) C⁻¹ = -H(-k, -B)` — needs k-transformation; the simple test was wrong from the start).
- Test krylov snapshot regenerated in Issue 03 fix2 (the old reference was generated with the OLD wire convention).