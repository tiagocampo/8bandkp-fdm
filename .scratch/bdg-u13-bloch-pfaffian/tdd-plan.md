# T1 Build Plan — Bloch-periodic H_BdG(k_par) stack builder

> Plan for ticket
> `.scratch/bdg-u13-bloch-pfaffian/issues/01-bloch-periodic-bdg-builder.md`.
> Built from chart route (`map.md`) + T6's resolved L1–L3 (recommendation (a)
> each). Author: 2026-08-09, branch `docs/u13-chart-pr35-drift-fix`, next
> action is `git checkout -b feat/bdg-u13-bloch-pfaffian main`.

## Context

T1 is the foundation of U13. The wire Pfaffian sweep today evaluates the
Majorana number at **one point only** (the parent plan at
`docs/plans/2026-06-14-001-...-plan.md` L53-54 names this deficiency). T1
fixes it by generalizing the fixed-kz CSR builder into a Bloch-periodic dense
stack matching `kitaev_majorana_number`'s shape (`src/math/pfaffian.f90:121`).
T1 unblocks T2 (S1×S2 strict agreement seam).

## Critical citations (already verified)

- Consumer: `kitaev_majorana_number(H_k_array, k_par_values, omega_struct)`
  at `src/math/pfaffian.f90:120-171`. Stack shape:
  `complex(kind=dp) H_k_array(:,:,n_k)`, `real(kind=dp) k_par_values(:)`.
  L142 loops `i = 1, n_k` and slices `H_k_array(:,:,i)`.
- Existing builder: `build_bdg_hamiltonian_1d` at
  `src/physics/bdg_hamiltonian.f90:142`, returns CSR. The `+-kz` dispatch
  at L173/182 calls `ZB8bandGeneralized`; the canonical `-conjg(H0(-k))`
  hole block is built via the `build_bdg_hole_block` wrapper.
- Bx-only Peierls SSOT: `defs.f90:963-969` (config-level rejection via
  `validate_semantic`). The by-orbital Peierls path is "not yet implemented";
  T1 inherits this gate as defense-in-depth at builder entry.
- Existing test fixture pattern: `tests/unit/test_bdg_hamiltonian.pf:225-339`
  builds a `nx=3, ny=3` GaAs wire, calls `confinementInitialization_2d`,
  `build_bdg_hamiltonian_1d`, then `csr_to_dense(H_bdg, Hd)`.
- ctest registration: `tests/CMakeLists.txt:160-164` shows the
  `add_pfunit_ctest(test_X ... LABELS "unit")` pattern.
- Workspace plumbing: `wire_workspace` reused via `ws=`; `g='g3'`
  derivative-build cache isolation is preserved by using the same helper
  the existing builder uses.

## Design

### Source change — `src/physics/bdg_hamiltonian.f90`

Refactor the existing `build_bdg_hamiltonian_1d` body into a **private
helper** that builds the BdG matrix as a **dense complex array** (no CSR
build, no Peierls COO loop). The existing public subroutine becomes a
thin wrapper that calls the helper and converts to CSR via the existing
COO path.

```fortran
! PRIVATE helper: builds H_bdg (16N x 16N, dense complex) for ONE k value.
! This is the bit-exact dense form of build_bdg_hamiltonian_1d.
subroutine build_bdg_hamiltonian_1d_dense( &
    H_bdg, kz, mu, delta_0, ws, B_vec, g_factor)

! PUBLIC existing CSR wrapper: thin shim — calls the helper + CSR-build
! (unchanged signature; bit-exact output preserved).
subroutine build_bdg_hamiltonian_1d(H_bdg_csr, cfg, ...)

! NEW PUBLIC: Bloch stack builder. Reuses the helper at each k_par.
subroutine build_bdg_hamiltonian_1d_bloch( &
    H_k_array, n_k, k_par_values, cfg, profile_2d, kpterms_2d, &
    mu, delta_0, ws, B_vec, g_factor)
  complex(kind=dp), intent(out) :: H_k_array(:,:,:)
  integer, intent(in) :: n_k
  real(kind=dp), intent(in) :: k_par_values(:)
  ...
  do i = 1, n_k
    call build_bdg_hamiltonian_1d_dense( &
        H_k_array(:,:,i), k_par_values(i), mu, delta_0, ws, B_vec, g_factor)
  end do
end subroutine
```

The helper carries the canonical `-conjg(H0(-k))` hole block + symmetric
Peierls +/-B loops unchanged — same logic as the existing builder, only
the output representation differs (dense complex vs CSR).

### Bx-only Peierls guard (defense-in-depth)

At `build_bdg_hamiltonian_1d_bloch` entry, mirror the SSOT at
`defs.f90:963-969`:

```fortran
if (present(B_vec)) then
  if (abs(B_vec(1)) < 1.0e-12_dp .and. &
      (abs(B_vec(2)) > 1.0e-12_dp .or. abs(B_vec(3)) > 1.0e-12_dp)) then
    error stop 'bdg_hamiltonian_bloch: BdG requires Bx nonzero ' // &
      'for Peierls orbital coupling (By-only path not yet implemented)'
  end if
end if
```

The existing `build_bdg_hamiltonian_1d` does NOT currently have this guard
(only `validate_semantic` does). T1 does NOT add the guard to the existing
builder — that's a separate hardening, out of scope. T1's bloch builder
inherits it because the bloch path is a *new* surface that bypasses
`validate_semantic` for direct programmatic callers.

### Public surface

`bdg_hamiltonian.f90:55` adds `public :: build_bdg_hamiltonian_1d_bloch`
next to the existing exports. No new module needed — keeps the seam in
`bdg_hamiltonian`, the AGENTS.md-listed entry point.

### Workspace

The dense helper takes `wire_workspace ws` by `intent(inout)` (same as
the existing builder). The k-loop reuses the same `ws` — `ZB8bandGeneralized`
mutates internal cache state, but the cache is keyed by `(kz, B_vec,
profile)`, so passing the same `ws` across k-slices is safe per the
existing workspace contract. If we discover during implementation that the
workspace cannot be safely reused across k-values, we allocate a fresh
`ws` per slice (within the loop) — `wire_workspace` is allocatable and
the `*_free` pattern is mandatory.

### Test — `tests/unit/test_bdg_hamiltonian_bloch.pf`

Two `@test` subroutines:

1. **`test_bloch_n1_matches_fixed_kz`** — builds the same `nx=3, ny=3`
   GaAs wire as `test_bdg_csr_assembly` (L226), calls both
   `build_bdg_hamiltonian_1d` at `kz=0.05_dp` AND
   `build_bdg_hamiltonian_1d_bloch` with `n_k=1, k_par_values=[0.05_dp]`.
   Asserts `csr_to_dense(H_bdg_csr) == H_k_array(:, :, 1)` element-wise
   to `1.0e-12_dp` tolerance. **This is the TDD-red test** — it must FAIL
   on `main@52743b3` (no bloch builder exists yet, compile error or
   symbol-not-found).

2. **`test_bloch_n2_two_hermitian_slices`** — calls the bloch builder at
   `n_k=2, k_par_values=[0.0_dp, 0.1_dp]`, asserts each slice is
   Hermitian (`max|H - H^dagger| < 1.0e-10_dp`). Guards the `n_k > 1`
   path that T3's sweep-mode knob will use.

3. **`test_bloch_bx_only_guard`** — calls the bloch builder with
   `B_vec=[0.0_dp, 0.5_dp, 0.0_dp]` (Bx=0, By≠0) and asserts an
   `error stop` is triggered. Implementation note: pFUnit `@test`
   cannot catch `error stop` (per project memory
   `feedback_pfunit_error_stop_uncatchable`), so this test uses a
   **positive-path** fixture — `B_vec=[0.5_dp, 0.0_dp, 0.0_dp]` (Bx-only,
   no error) — and asserts the builder runs to completion and produces a
   Hermitian stack. The guard is exercised at config-validation time
   upstream and at the bloch builder's `error stop` when called
   directly with bad input (out-of-scope for pFUnit). Test name reflects
   this: `test_bloch_accepts_bx_only_no_error`.

### ctest registration — `tests/CMakeLists.txt`

Insert between L164 (end of `test_bdg_hamiltonian`) and L166 (start of
`test_bdg_phs` comment block):

```cmake
# T1 (U13 Bloch-Pfaffian chart): Bloch-periodic H_BdG(k_par) stack
# builder. Stack shape: H_k_array(:,:,n_k) dense complex, matches
# kitaev_majorana_number consumer at src/math/pfaffian.f90:120-171.
add_pfunit_ctest(test_bdg_hamiltonian_bloch
    TEST_SOURCES unit/test_bdg_hamiltonian_bloch.pf
    LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
    LABELS "unit"
)
```

### Verification — count delta

Baseline: **51 unit tests, 3 wire-BdG regression tests**.
After T1: **52 unit tests, 3 wire-BdG regression tests**.

## Build sequence

1. Branch `feat/bdg-u13-bloch-pfaffian` off `main@52743b3`.
2. Author `tests/unit/test_bdg_hamiltonian_bloch.pf` + register in
   `tests/CMakeLists.txt`. Watch it FAIL: build error (symbol not found).
3. Refactor `src/physics/bdg_hamiltonian.f90`:
   a. Extract dense body from `build_bdg_hamiltonian_1d` into
      `build_bdg_hamiltonian_1d_dense` (private).
   b. Make existing public `build_bdg_hamiltonian_1d` a thin wrapper
      that calls the dense helper + CSR build (preserving bit-exact
      output for the regression guard).
   c. Add new public `build_bdg_hamiltonian_1d_bloch` with k-loop.
   d. Add `public :: build_bdg_hamiltonian_1d_bloch` to exports.
4. Rebuild + watch T1's RED tests turn GREEN.
5. Run full unit + wire-BdG regression gate:
   `OMP_NUM_THREADS=$((nproc/4)) ctest --test-dir build -j4
   --output-on-failure`.
6. Per `verification-before-completion`: read full output, confirm
   52/52 unit + 3/3 wire-BdG regression green. Confirm `n_k=1` regression
   test produced 0 FP-bit differences vs. the fixed-kz builder.
7. Close T1 in tracker (`.scratch/bdg-u13-bloch-pfaffian/issues/01...md`):
   `Status: completed`, link commit SHA, note the n_k=1 bit-for-bit
   result.
8. Append one-line entry to `map.md` ## Decisions so far.
9. Single commit. NO `Co-Authored-By` trailer (CLAUDE.md global
   `feedback_no_coauthor_trailer`).

## Risks and how they shake out

- **Risk:** refactor breaks `build_bdg_hamiltonian_1d`'s CSR output.
  **Mitigation:** the `test_bdg_hamiltonian.pf` suite (Hermiticity,
  PHS-pair, dimension-doubling) already exercises the existing builder;
  if the refactor preserves the dense body bit-exactly, all of those
  stay green. The TDD-red test for `n_k=1` is an additional guard.
- **Risk:** workspace cache invalidation across k-slices. **Mitigation:**
  the existing wire-builder pattern (`test_bdg_hamiltonian.pf:394-437`
  builds two `wire_workspace` instances back-to-back with no shared
  state) shows that the workspace is per-call. If the helper's internal
  cache is keyed by `(kz, B_vec)`, reusing `ws` across slices is safe;
  if not, we allocate a fresh `ws` per slice.
- **Risk:** pFUnit `error stop` uncatchable. **Mitigation:** the Bx-only
  guard test uses a positive-path fixture (Bx≠0, run to completion).
  Negative-path coverage is the existing `validate_semantic` config-level
  rejection in `test_bdg_config.pf`.

## Out of scope (deliberately)

- No `nk_par` config knob (T3's job — see `.scratch/bdg-u13-bloch-pfaffian/issues/03-bloch-sweep-config-knob.md`).
- No `main_topology.f90` dispatch change (T3's job).
- No S1×S2 strict agreement seam wiring (T2's job).
- No `kitaev_majorana_number` or `zskpfa_reduction` edits (those are
  in-tree and verified).
- No lecture-13 / `lecture_13_topological.py` change (T4's job).
- No parent plan footer update (T4's job per the L4 decision in T6).

## Verification gate

```
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 --output-on-failure
```

Expected:
- Unit label: **52/52** (51 baseline + 1 new T1 target; the new target
  contains 3 `@test` subroutines per project memory
  `feedback_ctest_counts_targets`).
- Regression label: 3/3 wire-BdG (unchanged).
- Lecture-13 acceptance: 1/1 (unchanged — T1 doesn't touch the gate).

If any baseline test flips red, the refactor broke something. **Stop and
revert** — do not commit a partial refactor.