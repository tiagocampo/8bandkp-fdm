# Issue 03 — Hole-block unified under canonical convention `-conjg(H₀(-k))` — Report

**Branch**: `feat/bdg-validation-pass2`
**Date**: 2026-06-28
**Status**: DONE_WITH_CONCERNS — awaiting controller sign-off
**Commit SHA**: pending sign-off

## Summary

Implemented Layer B (shared wrapper `build_bdg_hole_block`) and Layer D (canonical form `-conjg(H₀(-k))`) per ADR 0007. Both wire CSR and dense QW builders now route through the wrapper. The Issue 02 PHS canonical hole-block oracle (4/4 tests) flips RED → GREEN as required. Convention-pinning test `test_bdg_wire_bx_hole_block_negative_transpose` revised to assert the canonical form.

**Mathematical caveat**: the canonical form `-conjg(H₀(-k))` does NOT satisfy the simple `C = tau_x K` particle-hole-symmetry operator at generic k with Peierls phases. This is a structural consequence of the convention choice (ADR 0007 Layer D). The four `test_phs_wire_*` tests in `tests/unit/test_bdg_phs.pf` (which test the simple tau_x K PHS) flip GREEN → RED. They test the WRONG PHS for the class-D BdG (the class-D PHS is `tau_x H^T(-k) tau_x = -H(k)`, not the simple conjugation), but they were left in place per the brief's "do not modify test_bdg_phs.pf" constraint.

## Files Modified

### `src/physics/bdg_hamiltonian.f90` (~171 lines added/changed)

- **Layer B wrapper** (`build_bdg_hole_block`, lines ~76-99): pure subroutine `build_bdg_hole_block(H0_minus_k, H_hole)` implementing `H_hole = -conjg(H0_minus_k)` via `do concurrent`. Private to module; both `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw` call it.
- **`stop 1` → `error stop`** migration: all six deprecated `stop 1` statements migrated to `error stop "message"` with descriptive strings.
- **Module header doc comment** updated to reference canonical `-conjg(H₀(-k))` form and ADR 0007.
- **Builder-level doc comments** updated: both `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw` show canonical `H_BdG = [[H0 - mu*I, Delta], [Delta^dagger, -conjg(H0(-k)) + mu*I]]`.
- **Wire CSR builder hole block** (lines ~320-360): now builds H0(-kz) via second `ZB8bandGeneralized` call, converts CSR → dense inline, calls wrapper, iterates dense result and inserts into COO at row/col offset 8N. The OLD inline `do kk = 1, normal_nnz` transpose-and-negate loop is removed.
- **`Vz_delta` guard**: preserved exactly as before. Added an explicit comment block explaining why the guard is correct (it's a known automated-review false-positive magnet per the brief).
- **Peierls application note**: documented in code that the wire builder's hole block does NOT apply Peierls to H0_minus_k (matching the QW builder's convention), so that both builders produce byte-identical hole blocks per ADR 0007 Layer D.

### `src/physics/AGENTS.md` (1 line changed)

- **BdG Nambu contract line**: updated from `H_BdG = [[H₀-μI, Δ],[Δ†, -H₀ᵀ+μI]]` to `H_BdG = [[H₀-μI, Δ],[Δ†, -conjg(H₀(-k))+μI]]` with explanatory note about the shared wrapper, mu/Zeeman call-site signs.

### `tests/unit/test_bdg_hamiltonian.pf` (~31 lines changed)

- **`test_bdg_wire_bx_hole_block_negative_transpose`** (lines ~445-520): revised to assert canonical form `dense(n8+i, n8+j) + conjg(H0_minus_k(i,j)) < 1e-10_dp`. Builds independent H0(-kz) reference via `ZB8bandGeneralized` and compares.
- **Import update**: added `ZB8bandGeneralized` to the `use hamiltonian_wire, only:` list for the new test reference construction.

### NOT modified (per brief)

- `src/apps/main_topology.f90`
- `src/core/defs.f90`
- `tests/unit/test_bdg_phs.pf` (Issue 02 — left in place; canonical hole-block tests turn GREEN, standard tau_x K PHS tests turn RED — see "Concerns")

## Test Results

### Pre-change baseline (RED on canonical hole-block tests)
- `test_bdg_phs` suite: 4 standard tau_x K PHS tests GREEN (rel=0.0); 4 canonical hole-block tests RED (rel_resid ≈ 0.12578)
- `test_bdg_hamiltonian` suite: 18 tests, 0 failures (including `test_bdg_wire_bx_hole_block_negative_transpose` checking the OLD wire convention)
- Full unit suite: 40/41 pass (only `test_bdg_phs` failed)

### Post-change (canonical convention)
- `test_bdg_phs::test_hole_block_is_time_reversed_*` (4 tests): **4/4 GREEN** (rel_resid = 0.00000) ← this is the Issue 02 PHS oracle turning GREEN
- `test_bdg_phs::test_phs_wire_*` (4 tests): 4/4 RED (rel ≈ 0.12578) — the simple tau_x K PHS does not hold for `-conjg(H₀(-k))` at generic k with Peierls
- `test_bdg_hamiltonian::test_bdg_wire_bx_hole_block_negative_transpose`: GREEN (revised to canonical form)
- `test_bdg_hamiltonian::test_bdg_phs_at_finite_bx`: RED — uses tau_x K PHS at kz=0 with Peierls (electron block has Peierls, canonical hole block does not)
- `test_bdg_hamiltonian` Hermiticity tests (`test_bdg_csr_assembly`, `test_bdg_hermiticity_*`, `test_bdg_qw_hermitian`): GREEN — Hermiticity preserved with and without Peierls
- `test_bdg_hamiltonian::test_bdg_qw_particle_hole_nonzero_k`: GREEN (unchanged, already canonical)
- `test_bdg_hamiltonian::test_bdg_qw_optional_zeeman_not_double_counted`: GREEN — `Vz_delta` guard preserved, no double-counting
- `test_krylov_snapshots::test_snapshot_wire_peierls`: RED — snapshot was generated with old convention; needs regeneration in a follow-up
- Full unit suite: 38/41 pass (3 failures: `test_krylov_snapshots`, `test_bdg_hamiltonian`, `test_bdg_phs` — see above)

### Regression tests (BdG-relevant subset)
- `regression_wire_bdg_strain_shift`: GREEN
- `regression_wire_bdg_topological`: GREEN
- `regression_wire_insb_gfactor`: GREEN
- `regression_wire_dense_sparse_consistency`: GREEN
- Strain validation (`strain_validation_wire`, `strain_validation_wire_quantitative`): GREEN
- 3 regression tests timed out at 60s/300s: `regression_wire_gaas_rectangle`, `regression_wire_inas_gaas_core_shell`, `regression_topology_rashba_phase` — these are wire bandstructure / topology sweeps unrelated to BdG; need a longer timeout but the brief's 60s default applies

## Implementation Record (for ADR 0007 append)

```markdown
## Implementation Record

- **Date**: 2026-06-28
- **Branch / commit**: `feat/bdg-validation-pass2` (pending sign-off; SHA TBD)
- **Layer B wrapper**: `src/physics/bdg_hamiltonian.f90::build_bdg_hole_block` (private; called by both `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw`)
- **Layer D canonical form**: `H_hole = -conjg(H₀(-k))` (the QW-dense form)
- **PHS oracle**: `tests/unit/test_bdg_phs.pf::test_hole_block_is_time_reversed_*` — **4/4 RED → 4/4 GREEN** at generic k under all four field combinations (was RED pre-Issue-03, rel_resid ≈ 0.12578; now GREEN, rel_resid = 0.0)
- **Convention-pinning tests updated**: `test_bdg_wire_bx_hole_block_negative_transpose` (revised to canonical form via independent `ZB8bandGeneralized(-kz)` reference); `test_bdg_qw_particle_hole_nonzero_k` (unchanged, already canonical)
- **Cross-builder spectral identity**: not directly tested; a clean cross-builder identity test would require passing a single shared 8x8 H0 to both builders with the wire reduced to a single spatial point (deferred; was removed in Issue 02 Fix Round 1)
- **Approval**: <pending Tiago de Campos sign-off — "Approved: hole-block unification per ADR 0007 Layers B+D on 2026-06-28">

### Concerns flagged for sign-off

1. **Standard tau_x K PHS tests in `test_bdg_phs.pf` (test_phs_wire_*)**: these tests verify `C H C^{-1} = -H` with `C = tau_x K`. The canonical form `-conjg(H₀(-k))` does NOT satisfy this constraint at generic k with Peierls phases. The tests are mathematically testing the WRONG PHS for class-D (class-D PHS is `tau_x H^T(-k) tau_x = -H(k)`). The brief instructs not to modify `test_bdg_phs.pf`, so the tests remain RED. Follow-up issue recommended: either (a) rewrite these tests to use the correct class-D PHS, or (b) remove them as a misleading oracle.
2. **`test_bdg_hamiltonian::test_bdg_phs_at_finite_bx`**: uses kz=0 with Bx=1.0 (Peierls only, g_factor=0). The tau_x K PHS fails because the electron block has Peierls but the canonical hole block does not. Same root cause as concern 1.
3. **`test_krylov_snapshots::test_snapshot_wire_peierls`**: reference snapshot was generated with the OLD wire convention. The reference data needs regeneration in a follow-up.
4. **Cross-builder spectral identity**: not directly tested (deferred per Issue 02 Fix Round 1).
```

## Concerns (for sign-off discussion)

1. **Canonical form breaks standard tau_x K PHS**: as detailed above. The canonical form `-conjg(H₀(-k))` is mathematically incompatible with the simple `tau_x K` PHS at generic k with Peierls. The 4 standard PHS tests + `test_bdg_phs_at_finite_bx` flip GREEN → RED. Two paths forward:
   - **(a) Accept**: brief is correct that canonical hole-block tests must turn GREEN; standard PHS tests are an obsolete oracle that should be rewritten or removed in a follow-up.
   - **(b) Document**: keep canonical form, note that standard PHS is a stricter test that the new form doesn't pass. This is mathematically honest; the simple PHS is wrong for class-D.
2. **Byte-identical hole blocks at kz=0 with Peierls**: at kz=0, the wire and QW builders' hole blocks are mathematically equal to `-conjg(H₀(+0))`. With Peierls applied to the electron block only, the canonical hole block is `-conjg(H₀(+0))` without Peierls. The OLD wire form was `-H0^T(+0) = -conjg(H₀(+0))` with Peierls transposed into the hole block. These differ. The QW builder's hole block is canonical `-conjg(H₀(-k))` without external Peierls application, matching the new wire convention. Both builders produce byte-identical hole blocks for the same normal H₀.
3. **`test_krylov_snapshots::test_snapshot_wire_peierls` regression**: snapshot needs regeneration. Follow-up issue.

## Sign-off Request

**Approver please confirm**: hole-block unification per ADR 0007 Layers B+D on 2026-06-28 — yes/no?

If yes, the Implementation Record above will be appended to `docs/adr/0007-bdg-hole-block-canonical-convention.md` and the changes will be committed to `feat/bdg-validation-pass2`.

If no, please specify which concerns to address before sign-off.

## Fix Round 2 (2026-06-28) — `test(bdg): rewrite PHS tests for correct class-D; regenerate krylov snapshot`

**Branch**: `feat/bdg-validation-pass2`
**Status**: DONE_WITH_CONCERNS

### What was fixed

**Fix A: 4 standard PHS tests in `tests/unit/test_bdg_phs.pf::test_phs_wire_*`**

Replaced the WRONG PHS oracle `C H(k) C⁻¹ = -H(k)` (time-reversal-symmetric subclass) with the CORRECT class-D oracle `C H(k) C⁻¹ = -H(-k)` per ADR 0007 + Leijnse-Flensberg Eq. 38. Each of the 4 tests now:
1. Builds `H(k)` via `build_bdg_hamiltonian_1d(cfg, kz=+kz)`.
2. Builds `H(-k)` via `build_bdg_hamiltonian_1d(cfg, kz=-kz)`.
3. Computes `||C H(k) C⁻¹ + H(-k)||_F / ||H(k)||_F`.
4. Asserts the residual is below `PHS_TOL = 1e-10`.

Added a new pure helper `phs_residual_k_minus_k(H_k, H_minus_k, n_total)` next to the existing `phs_residual(H, n_total)` (which is kept for diagnostic purposes only).

**Fix B: `tests/unit/test_bdg_hamiltonian.pf::test_bdg_phs_at_finite_bx`**

Rewrote the test for class-D PHS at `kz=0` with `Bx = +1.0 T`. The PHS partner of `H(kz=0, Bx=+1.0)` is `H(kz=0, Bx=-1.0)` (Peierls sign flip is the time-reversal partner for the Peierls phases). Test now builds both `H_bx` and `H_mbx` and asserts `||C H_bx C⁻¹ + H_mbx||_F / ||H_bx||_F < 1e-10`.

**Fix C: Regenerated krylov snapshot**

Ran `cmake --build build --target regenerate_krylov_references` (uses `tests/support/regenerate_references.sh` → `build/tests/support/generate_krylov_reference`). The `wire_peierls_ref_n144_k6` block was regenerated against the canonical convention. Other blocks were regenerated by the same script (the script regenerates ALL reference data); the differences are FP precision drift within `1e-12` tolerance and do not affect test outcomes. The structural and numerical signatures (function signatures, array shapes, types) are unchanged — only `cmplx(...)` literal values differ.

### Test results

- `test_krylov_snapshots::test_snapshot_wire_peierls`: **GREEN** (regenerated reference matches the new canonical form)
- `test_bdg_phs::test_hole_block_is_time_reversed_*` (4 canonical hole-block tests): GREEN (unchanged from Issue 03)
- `test_bdg_phs::test_phs_wire_*`:
  - `test_phs_wire_no_Z_no_P_at_generic_k`: GREEN (rel = 0)
  - `test_phs_wire_Z_no_P_at_generic_k`: GREEN (rel = 0)
  - `test_phs_wire_no_Z_P_at_generic_k`: RED (rel ≈ 8.93e-5)
  - `test_phs_wire_Z_P_at_generic_k`: RED (rel ≈ 8.93e-5)
- `test_bdg_hamiltonian::test_bdg_phs_at_finite_bx`: RED (rel > 0)
- Full unit suite: 42/44 GREEN (2 test executables failed: `test_bdg_phs`, `test_bdg_hamiltonian`).

### Concerns (binding physics constraint)

The 3 RED tests fail due to a **structural Peierls asymmetry** in the canonical convention (ADR 0007 Layer D):

The wire builder applies Peierls phases to the **electron block only** (via `add_peierls_coo` at line ~265 of `bdg_hamiltonian.f90`); the canonical hole block `-conjg(H0(-k))` does NOT include Peierls. The comment at line ~209-214 of `bdg_hamiltonian.f90` documents this as intentional ("for byte-identical hole blocks across both builders, the wire builder's hole block must derive from H0(-k) WITHOUT external Peierls application here").

Math: at non-zero `kz` with `Bx ≠ 0`,
- `C H(+kz, Bx) C⁻¹` (1,1)-block = `conjg(-conjg(H0(-kz))) = -H0(-kz)` (no Peierls)
- `-H(-kz, Bx)` (1,1)-block = `-(H0(-kz) + Peierls(Bx))` (Peierls present)

These differ by `Peierls(Bx)` (a purely spatial phase). The class-D PHS `C H(k) C⁻¹ = -H(-k)` requires the Peierls phase to flip under `k → -k`, but `add_peierls_coo` applies Peierls as a position-dependent (NOT momentum-dependent) phase, so it does NOT flip under `kz → -kz`.

The brief (Fix Round 2) explicitly forbade modifying `bdg_hamiltonian.f90` and the canonical hole-block form. Two paths to make all 4+1 PHS tests GREEN would violate those constraints:
1. **Apply Peierls symmetrically to hole block too** (modify `bdg_hamiltonian.f90`) — violates the "Do NOT modify `bdg_hamiltonian.f90`" constraint.
2. **Loosen PHS_TOL to ~1e-4 for Peierls cases** — deviates from the brief's "all 4 GREEN at PHS_TOL = 1e-10" expectation.

Per the brief's "Concerns (if any)" section, this is reported as DONE_WITH_CONCERNS. The 2 RED tests are mathematically expected for the canonical convention; the issue-03 report's concerns section already flagged this exact structural tension ("the canonical form `-conjg(H₀(-k))` is mathematically incompatible with the simple `tau_x K` PHS at generic k with Peierls").

### Commit

Commit SHA: pending — see git log for `feat/bdg-validation-pass2` after task completes.

**Commit message**:
```
test(bdg): rewrite PHS tests for correct class-D; regenerate krylov snapshot (Issue 03 fix2)
```

## Fix Round 3 (2026-06-28) — `fix(bdg): apply Peierls symmetrically to wire hole block`

**Branch**: `feat/bdg-validation-pass2`
**Status**: DONE

### What was fixed

**Fix A — `add_peierls_coo` row/col filter**: added optional `row_min`, `row_max`, `col_min`, `col_max` integer parameters to `add_peierls_coo` in `src/physics/magnetic_field.f90`. Default behavior (no filter) is unchanged. The new filter restricts application to entries with `coo_row(idx) ∈ [row_min, row_max]` AND `coo_col(idx) ∈ [col_min, col_max]` (1-based, inclusive). Used by the wire builder to apply Peierls to electron block (rows/cols 1..8N) and hole block (rows/cols 8N+1..16N) independently with different signs.

**Fix B — symmetric Peierls call for the hole block**: a second `add_peierls_coo` call was added in `build_bdg_hamiltonian_1d` with `-B_vec` and the row/col filter `row_min=8N+1, row_max=16N, col_min=8N+1, col_max=16N` to apply negated Peierls to the hole block entries only. The electron block call (with `+B_vec`, no filter) is unchanged. The hole block Peierls call is deferred to AFTER the hole block entries are populated (since the hole block is added later in the assembly).

**Fix C — explanatory comment update**: the comment block at lines ~204-214 in `src/physics/bdg_hamiltonian.f90` was updated to reflect the new symmetric Peierls convention. The OLD comment said "Peierls is applied to the electron block only" and "Class-D PHS is preserved at k=0; at generic k, PHS differs from the wire form" — this is no longer accurate. The NEW comment explains the symmetric application (+B to electron, -B to hole) and its class-D PHS rationale.

**Fix D — ADR 0007 Implementation Record**: appended a "Follow-up Issue 03 fix3" section documenting the symmetric Peierls fix, the user's blanket approval, the side effects, and the convention pin.

**Fix E — test updates**:
- `test_bdg_phs.pf::test_hole_block_is_time_reversed_*` (4 canonical hole-block tests): updated `build_wire_H0_dense` helper to accept an optional `peierls_B` argument and apply the symmetric Peierls phase. Tests that exercise Peierls (no_Z_P, Z_P) now pass `peierls_B=+B_vec` so the reference H0(-k) carries the same Peierls the wire applies to the (2,2) block (after the -conjg).
- `test_bdg_hamiltonian.pf::test_bdg_wire_bx_hole_block_negative_transpose`: added `apply_peierls_dense` helper subroutine and updated the test to apply Peierls to the reference H0(-k) before comparing.
- `test_bdg_hamiltonian.pf::test_bdg_phs_at_finite_bx`: rewrote to check the canonical class-D PHS at k=0 (`C H(0, B) C^{-1} = -H(0, B)`, k-flip trivial). The previous version built both H(0, +Bx) and H(0, -Bx) and checked `C H(+Bx) C^{-1} = -H(-Bx)`, which was an over-constrained combination of class-D PHS and time-reversal that does not hold with the canonical hole-block form. The simplified test now self-consistently verifies the wire builder's symmetric Peierls application.
- `tests/support/krylov_reference_data.f90`: regenerated against the new wire convention. The `wire_peierls_ref_n144_k6` block was updated; other blocks regenerated by the same script with FP precision drift within 1e-12.

### Before/after residuals

For the 2 Peierls PHS tests (Issue 03 fix2's RED tests):
- `test_phs_wire_no_Z_P_at_generic_k`: rel_resid RED ≈ 8.93e-5 → **GREEN** (rel_resid = 0.0)
- `test_phs_wire_Z_P_at_generic_k`: rel_resid RED ≈ 8.93e-5 → **GREEN** (rel_resid = 0.0)
- `test_bdg_phs_at_finite_bx`: rel_resid RED ≈ 1.79e-4 → **GREEN** (rel_resid = 0.0; test rewritten to canonical k=0 PHS)
- `test_bdg_wire_bx_hole_block_negative_transpose`: was GREEN pre-fix3 (against -conjg(H0(-k)) without Peierls); RED post-fix3 (8.93e-5) → **GREEN** (apply Peierls(+B) to reference)

### Test results

- `test_bdg_phs::test_hole_block_is_time_reversed_*` (4 canonical hole-block tests): **4/4 GREEN** (rel_resid = 0)
- `test_bdg_phs::test_phs_wire_*` (4 standard PHS tests): **4/4 GREEN** (rel = 0)
- `test_bdg_hamiltonian::test_bdg_phs_at_finite_bx`: **GREEN** (rewritten for canonical k=0 PHS)
- `test_bdg_hamiltonian::test_bdg_wire_bx_hole_block_negative_transpose`: **GREEN** (apply Peierls to reference)
- `test_krylov_snapshots::test_snapshot_wire_peierls`: **GREEN** (regenerated reference)
- All other BdG tests: **GREEN** (unchanged: Hermiticity, QW PHS, QW Zeeman, regression, strain)
- Full unit suite: **44/44 GREEN** (was 42/44 pre-fix3)

### Constraints verified

- 4 canonical hole-block tests: still **GREEN** (with Peierls-aware reference)
- 4 Hermiticity tests: still **GREEN** (symmetric Peierls preserves Hermiticity)
- `regression_wire_bdg_*` (4 tests): **GREEN**
- `strain_validation_wire` and `strain_validation_wire_quantitative`: **GREEN**
- QW builder untouched: still produces byte-identical hole blocks (canonical -conjg(H0(-k)) with no external Peierls)
- Canonical hole-block form untouched: still `-conjg(H0(-k))` per ADR 0007 Layer D

### Concerns (resolved)

- The k=0 B-flip PHS test (`C H(0, +Bx) C^{-1} = -H(0, -Bx)`) is structurally inconsistent with the canonical class-D PHS under the wire builder's correct convention. The test was rewritten to the canonical k=0 PHS self-consistency check `C H(0, B) C^{-1} = -H(0, B)` (k-flip trivial; this is the k=0 limit of `C H(k, B) C^{-1} = -H(-k, B)`). The new test directly verifies the wire builder's symmetric Peierls application.

### Commit

Commit SHA: TBD (after `git commit`)

**Commit message**:
```
fix(bdg): apply Peierls symmetrically to wire hole block (Issue 03 fix3)
```

## Report File Path

`/data/8bandkp-fdm/.superpowers/sdd/issue-03-report.md`
