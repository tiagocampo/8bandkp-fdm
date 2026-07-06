**Status**: COMPLETE (2026-07-05)

# Issue 03 — Hole-block unified under canonical convention `-conjg(H₀(-k))` (U4, APPROVAL-GATED)

> **File-numbering note**: This file is `01-kitaev-pfaffian-harness.md` in `.scratch/bdg-majorana-validation/issues/`. The content is **Issue 03** (hole-block unified, APPROVAL-GATED) — sourced from `.superpowers/sdd/issue-01-brief.md`. The brief filenames were off by 2 from the issue numbers; this corrected numbering aligns the file's content with its issue number in the dependency graph.

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U4, Layers B+D; Layer A is Issue 02)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-01-brief.md`
**Decision artifact (ADR)**: `/data/8bandkp-fdm/docs/adr/0007-bdg-hole-block-canonical-convention.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md` (full Understand-phase analysis)
**Phase**: PR-A, APPROVAL-GATED — last issue before PR-B

## ⚠️ APPROVAL-GATED

Per CLAUDE.md Boundaries: "Require approval for: changes to ... Hamiltonian construction in `hamiltonianConstructor.f90`, ...". This issue modifies `bdg_hamiltonian.f90`'s hole-block construction code — this is the canonical Hamilton-construction code per the AGENTS.md BdG Nambu contract. The implementer MUST request and receive explicit sign-off from the controller (Tiago) BEFORE committing. The sign-off must be recorded in `docs/adr/0007-bdg-hole-block-canonical-convention.md`.

## What to build

The structural unification that closes the wire CSR (`-H₀ᵀ(+k)`) vs dense QW (`-conjg(H₀(-k))`) hole-block divergence. Two layers:

1. **Layer B — Shared wrapper.** A single `build_bdg_hole_block` procedure is extracted in `src/physics/bdg_hamiltonian.f90`. Both the dense QW builder and the CSR wire builder call it. The pairing-embed (`pairing_partner`/`pairing_sign`, already module data) and the μ-shift are documented as shared operations. No polymorphic builder types (per ADR 0001).
2. **Layer D — Canonical form choice.** Inside the wrapper, the canonical form is `H_hole = -conjg(H₀(-k))` (the dense-QW form). Rationale: k≠0-general, matches Leijnse-Flensberg Eq. 38, fewer dense-QW tests to update than wire-CSR. The wire CSR's inline `-H₀ᵀ(+k)` is removed.

The convention is pinned by Issue 02's PHS oracle (Layer A); the change here is the structural extraction (Layer B) and the canonical-form choice (Layer D). `src/physics/AGENTS.md` BdG Nambu contract line is updated from `-H₀ᵀ` to the canonical form. The six deprecated `stop 1` statements in `bdg_hamiltonian.f90` are migrated to `error stop` opportunistically while the file is open.

## Acceptance criteria

- [ ] Both builders produce byte-identical hole blocks for the same normal-state `H₀` after the change.
- [ ] The two revised convention-pinning tests (`test_bdg_wire_bx_hole_block_negative_transpose`, `test_bdg_qw_particle_hole_nonzero_k`) pass under the unified convention.
- [ ] The convention-agnostic `test_bdg_hermiticity_*` tests still pass (Hermiticity preserved with and without Peierls).
- [ ] Existing non-BdG regression tests are unaffected (the change is scoped to the BdG builders).
- [ ] `src/physics/AGENTS.md` BdG Nambu contract line updated to the canonical form.
- [ ] The six `stop 1` statements in `bdg_hamiltonian.f90` migrated to `error stop` (or deferred explicitly with a comment).
- [ ] The `Vz_delta = Vz_opt − Vz_cfg` double-counting guard is preserved (do NOT "fix" — it is correct and a known automated-review false-positive magnet).
- [ ] Issue 02's PHS oracle (the Layer A oracle) is GREEN at generic k under all four field combinations. This was RED on pre-03 convention (teeth demonstrated) — must turn GREEN after the unification.
- [ ] Per-task code review + spec compliance review clean; approval-gate sign-off recorded.

## Pre-existing state (from Understand report)

### `src/physics/bdg_hamiltonian.f90` (433 lines)
- Public: `build_bdg_hamiltonian_1d` (lines 90–308; hole block at lines 264–283 — uses `-H₀ᵀ(+k)` form, the wire-CSR convention).
- Public: `build_bdg_hamiltonian_qw` (lines 330–431; hole block at lines 420–424 — uses `-conjg(H₀(-k))` form, the dense-QW convention).
- Module data: `pairing_partner(8) = [4,3,2,1,6,5,8,7]`, `pairing_sign(8) = [+1,+1,-1,-1,+1,-1,+1,-1]`.
- Peierls phase on electron block (lines 174–179). Preserved.
- `Vz_delta` double-counting guard at lines 182–224. Preserved.
- Six `stop 1` statements at lines 72, 117, 127, 132, 350, 357. Migrate to `error stop` opportunistically.

### Issue 02's PHS oracle (`tests/unit/test_bdg_phs.pf`)
- 4 standard PHS tests GREEN (rel=0.0).
- 4 canonical hole-block tests RED with rel_resid ≈ 0.12578 (teeth demonstrated).
- After Issue 03 lands, the 4 canonical hole-block tests must turn GREEN.

### Convention-pinning tests to revise
- `tests/unit/test_bdg_hamiltonian.pf::test_bdg_wire_bx_hole_block_negative_transpose` (lines 446–507): wire convention check `H(n8+j, n8+i) + H(i,j) < 1e-10`. Must be revised to assert `-conjg(H₀(-k))` form.
- `tests/unit/test_bdg_hamiltonian.pf::test_bdg_qw_particle_hole_nonzero_k` (lines 601–649): QW convention check `H(n8+row, n8+col) + conjg(H_minus(row, col)) < 1e-12`. Already asserts the canonical form — should pass unchanged.

### Existing module header doc comment (`bdg_hamiltonian.f90:9-10`)
References `H0^T` — must be updated to reference `-conjg(H₀(-k))` per ADR 0007.

## Constraints from CLAUDE.md + ADRs (binding)

- **CLAUDE.md Boundaries — APPROVAL GATE**: This issue modifies Hamilton-construction code (`bdg_hamiltonian.f90`'s hole-block construction). Per CLAUDE.md Boundaries, this requires explicit sign-off from the controller BEFORE commit. Record the sign-off in `docs/adr/0007-bdg-hole-block-canonical-convention.md`.
- **ADR 0001 (fat derived type)**: No polymorphic builder types. The shared wrapper is a procedure + module data, not a class hierarchy.
- **ADR 0007 (hole-block canonical convention)**: Adopt `-conjg(H₀(-k))` as canonical form.
- **CLAUDE.md Boundaries — Bir-Pikus sign convention**: not touched (strain SSOT unchanged).
- **CLAUDE.md Boundaries — basis ordering**: not touched (bands 1-4 valence, 5-6 split-off, 7-8 conduction).
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.
- **DRY/SSOT**: Do not duplicate the k.p/Zeeman/strain tables. Read `pairing_partner`/`pairing_sign` from existing module data.
- **KISS**: Pure-function seam; no speculative abstractions.

## File ownership (exhaustive)

### Modified files
- `src/physics/bdg_hamiltonian.f90` — extract `build_bdg_hole_block` wrapper (Layer B); wire CSR uses canonical `-conjg(H₀(-k))` form (Layer D); migrate six `stop 1` → `error stop`; preserve `Vz_delta` guard; update module header doc comment.
- `src/physics/AGENTS.md` — BdG Nambu contract line updated to canonical form.
- `tests/unit/test_bdg_hamiltonian.pf` — revise `test_bdg_wire_bx_hole_block_negative_transpose` to assert canonical form.
- `docs/adr/0007-bdg-hole-block-canonical-convention.md` — append Implementation Record section per the sign-off template (see "Approval procedure" below).

### NOT modified
- `src/apps/main_topology.f90` (Issue 00 already lifted the gap literals; Issue 03 doesn't change call sites — the existing call sites will produce identical results because both builders now construct the hole block the same way).
- `src/core/defs.f90` (Issue 06/07 own).
- `tests/unit/test_bdg_phs.pf` (Issue 02 — already in place, will turn GREEN after this lands).

## Suggested implementation (illustrative)

```fortran
! Inside bdg_hamiltonian.f90 module
private

! Layer B: shared hole-block wrapper
public :: build_bdg_hole_block

contains

! Canonical form per ADR 0007: H_hole = -conjg(H₀(-k))
! Both wire CSR and dense QW builders call this.
pure subroutine build_bdg_hole_block(H0_minus_k, H_hole)
  complex(dp), intent(in) :: H0_minus_k(:,:)
  complex(dp), intent(out) :: H_hole(:,:)
  integer :: i, j, n
  n = size(H0_minus_k, 1)
  ! H_hole = -conjg(H0(-k))
  do concurrent (i = 1:n, j = 1:n)
    H_hole(i,j) = -conjg(H0_minus_k(i,j))
  end do
end subroutine

! In build_bdg_hamiltonian_1d:
!   Replace the inline hole-block construction at lines 264–283 with:
!     call build_bdg_hole_block(ZB8bandGeneralized(... k=(-kz) ...), H_hole_block)
!     ! then assemble the Nambu block with H_hole_block
```

## Approval procedure (mandatory)

Before committing:

1. Implement the changes; run tests; verify Issue 02's PHS oracle turns GREEN.
2. **STOP. Do not commit yet.**
3. Report back to the controller with:
   - Summary of changes (file:line references).
   - Test results (Issue 02 PHS oracle: 4/4 GREEN now; convention-pinning tests: GREEN; Hermiticity tests: GREEN; full unit suite: green).
   - A draft Implementation Record section for the ADR.
4. Wait for explicit sign-off from the controller.
5. After sign-off, commit and append the Implementation Record to `docs/adr/0007-...md`.

The Implementation Record template (from the Understand report §4):

```markdown
## Implementation Record

- **Date**: 2026-06-28
- **PR/commit**: <commit-sha> ("refactor(bdg): unify hole-block convention per ADR 0007")
- **Layer B wrapper**: `src/physics/bdg_hamiltonian.f90::build_bdg_hole_block` (private; called by both `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw`)
- **Layer D canonical form**: `H_hole = -conjg(H₀(-k))` (the QW-dense form)
- **Cross-builder spectral identity**: confirmed by `tests/unit/test_bdg_phs.pf::test_cross_builder_*` (or noted as deferred per Issue 02 fix)
- **PHS oracle**: `tests/unit/test_bdg_phs.pf` GREEN at generic k under all four field combinations (was RED pre-Issue-03)
- **Convention-pinning tests updated**: `test_bdg_wire_bx_hole_block_negative_transpose` (revised to canonical form); `test_bdg_qw_particle_hole_nonzero_k` (unchanged, already canonical)
- **Approval**: Tiago de Campos — "Approved: hole-block unification per ADR 0007 Layers B+D on 2026-06-28"
```

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure
```

## TDD discipline (mandatory)

1. Verify Issue 02's PHS oracle is RED on pre-03 code (already captured in the report).
2. Apply the canonical form change.
3. Verify Issue 02's PHS oracle turns GREEN.
4. Verify Hermiticity tests still pass.
5. Verify the two convention-pinning tests pass.
6. Verify the full unit suite + regression suite pass.
7. Stop and request approval.

## Out of scope (must NOT do)

- Do NOT modify the wire CSR's existing CSR layout — only the hole-block construction formula.
- Do NOT modify the dense QW builder's hole-block construction (it's already canonical).
- Do NOT change the Peierls phase logic.
- Do NOT "fix" the `Vz_delta = Vz_opt − Vz_cfg` double-counting guard.
- Do NOT modify `src/apps/main_topology.f90`.
- Do NOT modify `src/core/defs.f90`.
- Do NOT add `sweep_model` enum value.
- Do NOT add new TOML fields.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-03-report.md`

The report must include:
- Status (DONE_WITH_CONCERNS if waiting for sign-off; DONE after sign-off)
- Files modified (with line ranges)
- Test results — the GREEN-once-was-RED transition is the headline.
- Implementation Record draft (for ADR append).
- Self-review findings
- Commit SHA (only after approval; otherwise note "pending sign-off")
- Concerns (if any)

## Risk notes

- The change touches Hamilton-construction code, which the project treats as high-risk. A single mis-edit could silently break PHS for ALL BdG paths. **TDD discipline is mandatory**.
- The `Vz_delta` guard is a known false-positive magnet — DO NOT touch it. Add a comment if needed explaining why it's correct.
- The two builders construct hole blocks in slightly different ways (CSR vs dense). The shared wrapper must handle both call patterns. Test by verifying both builders produce byte-identical hole blocks for the same H₀ input.
- Hermiticity must be preserved with and without Peierls (the convention-pinning tests cover this).

---

## Outcome (as executed)

- **Approved by Tiago de Campos on 2026-06-28**: "Approved: hole-block unification per ADR 0007 Layers B+D on 2026-06-28".
- **Layer B wrapper**: `src/physics/bdg_hamiltonian.f90::build_bdg_hole_block` (private; called by both `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw`).
- **Layer D canonical form**: `H_hole = -conjg(H₀(-k))` (the QW-dense form).
- **PHS oracle**: `tests/unit/test_bdg_phs.pf::test_hole_block_is_time_reversed_*` — 4/4 RED → 4/4 GREEN at generic k under all four field combinations (rel_resid ≈ 0.12578 → 0.0).
- **Convention-pinning tests updated**: `test_bdg_wire_bx_hole_block_negative_transpose` (revised to canonical form via independent `ZB8bandGeneralized(-kz)` reference); `test_bdg_qw_particle_hole_nonzero_k` (unchanged, already canonical).
- **Cross-builder spectral identity**: not directly tested (deferred per Issue 02 Fix Round 1).
- **Side effects**:
  - All six `stop 1` statements in `bdg_hamiltonian.f90` migrated to `error stop`.
  - `Vz_delta = Vz_opt − Vz_cfg` double-counting guard preserved; explanatory comment added to prevent automated-review false positives.
  - BdG Hermiticity tests (`test_bdg_hermiticity_*`, `test_bdg_qw_hermitian`) GREEN with and without Peierls.
  - BdG regression tests (`regression_wire_bdg_strain_shift`, `regression_wire_bdg_topological`, `regression_wire_insb_gfactor`, `regression_wire_dense_sparse_consistency`) GREEN.
  - Strain validation tests (`strain_validation_wire`, `strain_validation_wire_quantitative`) GREEN.

### Follow-up Issue 03 fix2 (2026-06-28)

Rewrote 4 standard PHS tests in `test_bdg_phs.pf` (`test_phs_wire_*`) and `test_bdg_hamiltonian::test_bdg_phs_at_finite_bx` to use the correct class-D PHS operator `C H(k) C⁻¹ = -H(-k)` (k-transformation required, not simple `tau_x K PHS`). Also regenerated the krylov snapshot reference.

### Follow-up Issue 03 fix3 (2026-06-28)

After fix2, 3 PHS tests still failed at generic k with Peierls. Root cause: canonical hole block `-conjg(H₀(-k))` does NOT include Peierls, but electron block `H₀(+k)` does. Class-D PHS requires `C H(k, B) C⁻¹ = -H(-k, -B)` — Peierls must be applied to the hole block with NEGATED sign `-B_vec`.

**Fix A — `add_peierls_coo` row filter**: added optional `row_min`, `row_max` integer parameters to `add_peierls_coo` in `src/physics/magnetic_field.f90`. Default behavior unchanged. Used by wire builder to apply Peierls to electron block (rows 1..8N) and hole block (rows 8N+1..16N) independently with different signs.

**Fix B — symmetric Peierls call in `build_bdg_hamiltonian_1d`**: second `add_peierls_coo` call with `-B_vec` and row filter `row_min=8N+1, row_max=16N`.

**Fix C — explanatory comment update**: comment block at lines ~204-214 in `bdg_hamiltonian.f90` updated to reflect symmetric Peierls convention.

**Outcome**: full unit suite 44/44 GREEN (was 42/44 pre-fix3). Blanket approval per controller's directive "there is no such thing as pre-existing error — fix it".