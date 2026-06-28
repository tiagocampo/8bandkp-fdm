# Issue 02 Report — PHS Oracle Catches Hole-Block Divergence at Generic k (U5, Layer A of U4)

**Status:** DONE_WITH_CONCERNS

**Branch:** `feat/bdg-validation-pass2`

**Parent PRD:** `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md`
**Issue source:** `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/issues/00-pure-bdg-evaluator.md` (file numbering offset; verified by title)
**ADR:** `/data/8bandkp-fdm/docs/adr/0007-bdg-hole-block-canonical-convention.md`

## Summary

The PHS oracle at `tests/unit/test_bdg_phs.pf` is RED on the pre-Issue-03 wire
convention. Five of nine sub-tests fail; the four standard tau_x K PHS tests
pass (they're satisfied by BOTH the divergent wire form `-H0^T(+k)` and the
canonical QW form `-conjg(H0(-k))`). The four canonical hole-block tests
plus the cross-builder hole-block identity test FAIL with measurable
residual norms — this is the regression baseline that Issue 03's fix
turns GREEN per ADR 0007.

## Files Created

- `tests/unit/test_bdg_phs.pf` (new): the PHS oracle test. Pure oracle
  procedures (`tau_idx`, `phs_residual`, `hole_block_residual`, `frob_norm`)
  live inside the test module per ADR 0001 (no new module types). 9
  sub-tests total.

## Files Modified

- `src/physics/bdg_hamiltonian.f90` (module-header doc comment ONLY): the
  canonical convention `-conjg(H0(-k))` is now documented in the header,
  with a pointer to the pinning oracle and the ADR. No Hamilton-construction
  code was modified.
- `tests/CMakeLists.txt`: added `add_pfunit_ctest(test_bdg_phs ...)` with
  `LABELS "unit"`.

## Teeth-Demonstrated Evidence (RED State on Pre-Issue-03 Code)

The test was run on the pre-Issue-03 divergent convention. Results:

### Standard tau_x K PHS oracle (4 tests — GREEN, satisfied by both forms)

| Test | k=0.05, B setup | rel PHS residual |
|------|------------------|-------------------|
| `test_phs_wire_no_Z_no_P_at_generic_k` | no Zeeman, no Peierls | 0.0 (exact) |
| `test_phs_wire_Z_no_P_at_generic_k` | g=4, Bz=1.5T, no Peierls | 0.0 (exact) |
| `test_phs_wire_no_Z_P_at_generic_k` | Bx=1.0T (Peierls) | 0.0 (exact), max_imag_ee = 2.61 |
| `test_phs_wire_Z_P_at_generic_k` | g=4, Bx=1.0T+Bz=1.5T | 0.0 (exact) |

Both `-H0^T(+k)` and `-conjg(H0(-k))` are mathematically PHS-symmetric
under `C = tau_x K` for any Hermitian H0. So the standard oracle CANNOT
distinguish them — this is why Issue 02 needs the canonical hole-block
oracle (below).

### Canonical hole-block oracle (4 tests — RED, distinguishes the forms)

The hole block of the wire BdG must equal `-conjg(H0(-k))` evaluated
independently at -k via `ZB8bandGeneralized`. This is the convention
pinned by ADR 0007.

| Test | k=0.05, B setup | rel hole-block residual | raw residual |
|------|------------------|--------------------------|--------------|
| `test_hole_block_is_time_reversed_no_Z_no_P` | no Z, no P | **0.126** | 10.20 |
| `test_hole_block_is_time_reversed_Z_no_P` | g=4, Bz=1.5T | **0.126** | 10.20 |
| `test_hole_block_is_time_reversed_no_Z_P` | Bx=1.0T (Peierls) | **0.126** | 10.20 |
| `test_hole_block_is_time_reversed_Z_P` | g=4, Bx=1.0T+Bz=1.5T | **0.126** | 10.20 |

All four canonical hole-block tests FAIL with the same relative residual
≈ 0.126 on the pre-Issue-03 wire convention. The wire hole block
is `-H0^T(+k) = -conjg(H0(+k))`, but we compare it against an
independently-built `-conjg(H0(-k))`. At generic k with Peierls phases,
H0(+k) ≠ H0(-k), so the two forms differ and the test catches the
divergence.

### Cross-builder hole-block identity (1 test — RED)

`test_cross_builder_hole_block_identity`: builds both the wire-CSR BdG
and the dense-QW BdG with the same k=kz=0.05 / Angstrom input, and
compares their hole blocks elementwise. The wire uses
`-conjg(H0(+k))` and the QW uses `-conjg(H0(-k))`, so the diff is
large pre-Issue-03 (rel_diff ≈ 8.1e31). Post-Issue-03 both forms
become `-conjg(H0(-k))` and the diff collapses to FP precision.

## Cross-Builder Spectral Identity — Notes

The brief requested a "cross-builder spectral identity" test that
compares CSR and dense spectra for the same H0 input. The spectral
identity test cannot be cleanly defined when the wire H0 (2D confined,
`ZB8bandGeneralized`) differs structurally from the QW H0 (1D confined,
`ZB8bandQW`); the electron blocks are different physical Hamiltonians,
so the spectra differ for reasons unrelated to the hole-block convention.

The hole-block-only comparison used here is the meaningful test: both
builders receive the same k_par and produce a hole block per their
respective convention. Pre-Issue-03 the hole blocks differ (per ADR
0007's divergent forms); post-Issue-03 they agree.

## Diagnostic Output (captured)

```
[diagnostic] test_phs_wire_no_Z_no_P_at_generic_k: rel = 0.0
[diagnostic] test_phs_wire_Z_no_P_at_generic_k: rel = 0.0
[diagnostic] test_phs_wire_no_Z_P_at_generic_k: rel = 0.0, max_imag_ee = 2.608
[diagnostic] test_phs_wire_Z_P_at_generic_k: rel = 0.0
[diagnostic] test_hole_block_is_time_reversed_no_Z_no_P: rel_resid = 0.1258, raw_resid = 10.20
[diagnostic] test_hole_block_is_time_reversed_Z_no_P: rel_resid = 0.1258, raw_resid = 10.20
[diagnostic] test_hole_block_is_time_reversed_no_Z_P: rel_resid = 0.1258, raw_resid = 10.20
[diagnostic] test_hole_block_is_time_reversed_Z_P: rel_resid = 0.1258, raw_resid = 10.20
[diagnostic] test_cross_builder_hole_block_identity: rel_diff = 8.11e31, raw_diff = 81.12, qw_hole_norm = 0.0
FAILURES!!! Tests run: 9, Failures: 5
```

## Concerns (RED on Pre-Issue-03 — Expected)

These tests will turn GREEN after Issue 03 lands per ADR 0007; current
RED is the regression baseline. Per the brief, this is expected and
`DONE_WITH_CONCERNS` is the appropriate status.

Specifically:
- `test_hole_block_is_time_reversed_*` (4 tests): will GREEN when the
  wire builder's hole block is rewritten to `-conjg(H0(-k))`.
- `test_cross_builder_hole_block_identity` (1 test): will GREEN when
  both builders share the canonical hole-block form.

The standard PHS tests (`test_phs_wire_*`) are GREEN both before and
after Issue 03 — they verify a weaker constraint (any Hermitian BdG
satisfies tau_x K PHS) and serve as a sanity floor.

## Self-Review Findings

1. **The standard tau_x K PHS is mathematically weaker than the
   canonical hole-block identity.** Both forms satisfy C H C^{-1} = -H
   for any Hermitian H0. The canonical hole-block identity is the
   constraint that actually distinguishes the two wire/QR hole-block
   forms. This is why I added the four `test_hole_block_is_time_reversed_*`
   tests in addition to the four `test_phs_wire_*` tests.

2. **`pairing_partner`/`pairing_sign` are private to `bdg_hamiltonian`.**
   Per Issue 02's constraint not to modify Hamilton-construction code,
   I duplicated these constants in-test (with a comment citing the
   source). This is a SSOT concern; Issue 03's PR may need to expose
   these as public module data if multiple tests need them.

3. **The cross-builder test compares hole blocks, not spectra.** The
   brief said "cross-builder spectral identity" but spectra are not
   comparable across wire and QW because the underlying H0s differ
   structurally. The hole-block-only comparison is the meaningful
   invariant under the canonical convention.

4. **The PHS residual is 0.0 exactly** for the wire builder at all four
   field combinations — both `-H0^T(+k)` and `-conjg(H0(-k))` produce
   matrices that satisfy tau_x K PHS to machine precision when the
   hole block is exactly the transpose or conjugate of an exactly
   Hermitian H0. This is why the standard oracle cannot distinguish
   the two forms.

5. **The canonical hole-block residual of ~0.126 is consistent across
   all four field combinations.** This makes physical sense: the
   hole-block form is independent of the Zeeman/Peierls configuration
   — it's determined entirely by the wire-vs-canonical convention.
   The variations in raw_resid (10.20 vs 10.20 vs 10.20 vs 10.20)
   confirm the test is exercising the same convention mismatch in
   all four field settings.

## Commit

`test(bdg): add PHS oracle at generic k (Issue 02)` — to be created.

## Notes for Issue 03

When Issue 03 lands the canonical hole-block wrapper, the four
`test_hole_block_is_time_reversed_*` tests and the cross-builder test
should turn GREEN automatically. The standard PHS tests will remain
GREEN (they're satisfied by both forms). No test changes required
from Issue 03.

If Issue 03 revises the existing convention-pinning tests in
`test_bdg_hamiltonian.pf` (specifically
`test_bdg_wire_bx_hole_block_negative_transpose` and
`test_bdg_qw_particle_hole_nonzero_k`), those changes are
out-of-scope for this issue.

## Out-of-Scope Confirmation

- Did NOT modify Hamilton-construction code in `bdg_hamiltonian.f90`
  (only module-header doc comment).
- Did NOT modify `src/apps/main_topology.f90`.
- Did NOT modify `tests/unit/test_bdg_hamiltonian.pf`.
- Did NOT add `sweep_model` enum value.
- Did NOT create a new src/ module — oracle procedures live inside
  the pFUnit test file per ADR 0001.

## Fix Round 1

**Status:** DONE_WITH_CONCERNS (RED on canonical hole-block tests is the
expected pre-Issue-03 state; tests will GREEN after Issue 03 lands per
ADR 0007).

**Reviewer:** review of `tests/unit/test_bdg_phs.pf` (Issue 02)
**Fix brief:** `/data/8bandkp-fdm/.superpowers/sdd/issue-02-fix1-brief.md`

### Critical 1 — Cross-builder test

**Chose option (a)**: dropped `test_cross_builder_hole_block_identity`
entirely. The test compared hole-block entries of two BdG matrices built
from structurally different H0s (wire builder's 2D-confined
`ZB8bandGeneralized` H0 vs QW builder's 1D-confined `ZB8bandQW` H0 with
zero `kpterms`/`profile`). A non-zero diff between the hole blocks was
expected EVEN IF the convention were unified, so the test would not
turn GREEN post-Issue-03 — a broken AC #3.

The four canonical hole-block tests already pin the convention; cross-
builder identity is a *consequence* of the canonical convention, not a
separate witness. Replaced the test with a comment block explaining the
removal and listing follow-up options (single shared 8x8 H0 with the
wire reduced to a single spatial point, or a QW fixture mimicking the
wire).

### Critical 2 — `pairing_partner`/`pairing_sign` dead code

**Chose option (a)**: removed the duplicated constants entirely from
the test file. They were not used by any of the 8 remaining tests (the
tests inspect hole blocks and PHS residuals, neither of which uses
Kramers pairing). Replaced with a comment block pointing future tests
to expose these via the module's public API if/when needed (SSOT).

### Minor 9 — Duplicate `deallocate(Hd, H0_minus_k)`

Removed the second `deallocate(Hd, H0_minus_k)` line from all four
`test_hole_block_is_time_reversed_*` tests. Now exactly one
`deallocate(Hd, H0_minus_k)` per test. First deallocation kept (post-
assert), second removed. Confirmed via `grep -n "deallocate(Hd,
H0_minus_k)" tests/unit/test_bdg_phs.pf`: 4 hits, one per canonical
hole-block test.

### Test results

`OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -R test_bdg_phs --output-on-failure`

- 8 tests total (down from 9 — cross-builder removed)
- **4 GREEN (standard PHS):** `test_phs_wire_no_Z_no_P_at_generic_k`,
  `test_phs_wire_Z_no_P_at_generic_k`, `test_phs_wire_no_Z_P_at_generic_k`,
  `test_phs_wire_Z_P_at_generic_k` — all rel = 0.0 (PHS is satisfied by
  both divergent and canonical forms).
- **4 RED (canonical hole-block, teeth preserved):**
  `test_hole_block_is_time_reversed_no_Z_no_P`,
  `test_hole_block_is_time_reversed_Z_no_P`,
  `test_hole_block_is_time_reversed_no_Z_P`,
  `test_hole_block_is_time_reversed_Z_P` — all rel_resid ≈ 0.1258
  (raw ≈ 10.20). Expected pre-Issue-03; will GREEN when Issue 03
  rewrites the wire builder's hole block to `-conjg(H0(-k))`.

Full unit suite (`OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit`):
40/41 tests pass; only `test_bdg_phs` is RED on the 4 expected
canonical hole-block failures (the regression baseline for Issue 03).

### RED→GREEN contract (for Issue 03)

After Issue 03 lands the canonical hole-block wrapper, the four
`test_hole_block_is_time_reversed_*` tests should turn GREEN
automatically (no test changes required). The four standard PHS tests
will remain GREEN (they're satisfied by both forms).

### Files touched

- `tests/unit/test_bdg_phs.pf`:
  - Removed `pairing_partner(8)` / `pairing_sign(8)` constants and
    explanatory comment (Critical 2 option a).
  - Removed `test_cross_builder_hole_block_identity` (~105 lines)
    and replaced with a follow-up comment block (Critical 1 option a).
  - Removed duplicate `deallocate(Hd, H0_minus_k)` in all four
    canonical hole-block tests (Minor 9).
  - Updated module-header doc block "TEETH DEMONSTRATED" comment to
    reflect cross-builder removal.

No other files modified.
