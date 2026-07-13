# PR #39 Review Fixes — Design

- **Date:** 2026-06-21
- **Branch:** `fix/pr39-review-fixes` (off `refactor/architecture-deepening`, PR #39)
- **Status:** Approved (brainstorm 2026-06-21)
- **Supersedes / companions:** builds on `2026-06-13-eigensolver-review-fixes-design.md` (the dispatch cleanup that already landed in PR #39)

## Context

A max-effort code review of PR #39 ("FEAST parity everywhere except bulk + tested
foundation under the k.p block table") surfaced 15 findings. PR #39's default
(AUTO/DENSE) paths are bit-identical and all 123 ctests pass; the risk is
concentrated in the **opt-in FEAST paths**, where the per-sweep reconciliation
logic is unguarded against partial FEAST returns and inconsistent across sites,
plus two unenforced spec invariants.

This spec covers the corrective fixes. Two larger refactors and one dedup are
**deferred** to `docs/plans/BACKLOG.md` ("PR #39 Review — Deferred Refactors").

## Scope

**In this pass (Groups A + B + safe C):**

| # | Finding | Fix group |
|---|---------|-----------|
| 1, 2, 4, 10 | FEAST↔DENSE-INDEX reconciliation duplicated, divergent, OOB | A1 |
| 3 | `info=3` populates result; wire g-factor + SC loop use truncated spectrum | A2 |
| 5 | `bulk + FEAST` not permanently rejected | A3 |
| 7 | `apply_solver_window` OR vs `validate()` AND → degenerate partial window | A4 |
| 13 | `kp_scalar_block` silent-0 vs sibling `error stop` | B2 |
| 8 | Cached dense-LAPACK buffer non-conformance on shrinking N (latent) | B3 |
| 14, 11 | YAGNI dead code (`adopt_precomputed`, `asw_evals`); duplicated margin constants | C1 |

**Deferred to BACKLOG:** #9 (`resolve_kp_term` truly drives the CSR formula),
#12 (decompose `compute_strain_wire`), #15 (dedup the sweep-window derivation
block across 5 sites — not currently buggy; needs a build-callback abstraction).

## Locked design decisions

1. **#3 — root-cause in `solve_feast`, not at call sites.** The result is
   populated only on genuine convergence. Every caller already does the right
   thing with `nev_found = 0`; the wire SC loop yields zero contribution for an
   info=3 k-point (accepted — safer than charging ahead with a truncated
   spectrum).
2. **Reconciliation — one shared helper, strict.** A single
   `eigensolver_result_band_slice` replaces the three divergent hand-rolled
   blocks; a partial spectrum is a hard `error stop` (no silent zero-fill). The
   QW band-structure path switches to it (correctness fix on the opt-in FEAST
   path; AUTO/DENSE golden outputs untouched).
3. **#7 — reject partial windows in `validate()`.** "Both-or-neither" becomes an
   explicit config contract; `apply_solver_window`'s OR override is then
   provably safe and is left untouched.

## Components

### A1 — Reconciliation helper + OOB guard  (`src/math/eigensolver.f90`)

New public free subroutine:

```fortran
subroutine eigensolver_result_band_slice(result, il, iu, evals_out, evecs_out)
  type(eigensolver_result), intent(in) :: result
  integer, intent(in) :: il, iu
  real(kind=dp), intent(out), contiguous :: evals_out(:)        ! (nev_target)
  complex(kind=dp), intent(out), contiguous :: evecs_out(:,:)   ! (N, nev_target)
  integer :: nev_target, idx_lo
  nev_target = iu - il + 1
  if (result%nev_found >= iu) then
    idx_lo = il                          ! FEAST full in-window spectrum
  else if (result%nev_found == nev_target) then
    idx_lo = 1                           ! DENSE+INDEX pre-sliced (zheevx range 'I')
  else
    error stop 'eigensolver_result_band_slice: FEAST window truncated — got ' // &
      '<M> eigenvalues, need <iu> to extract bands [<il>,<iu>]; widen the energy window'
  end if
  evals_out(1:nev_target)   = result%eigenvalues(idx_lo : idx_lo+nev_target-1)
  evecs_out(:,1:nev_target) = result%eigenvectors(:, idx_lo : idx_lo+nev_target-1)
end subroutine
```

**Call sites refactored** to use it (removing the inline `idx_lo` / lowest-M
logic):
- `src/apps/main_optics.f90` — QW optics k_par sweep
- `src/apps/main.f90` — Landau k-sweep
- `src/apps/main.f90` — QW band-structure CSR sweep (drops `eig(1:M)=eigenvalues(1:M)`)

The `M < nev_target` guards already present at each site stay (they fire first,
with the existing "insufficient eigenvalues" message); the helper's partial-case
`error stop` covers the previously-OOB `nev_target < M < iu` window.

Bulk (8×8) is **not** routed through this — trivial, unchanged.

### A2 — info=3 root-cause  (`src/math/eigensolver.f90:424`)

Change the result-population guard in `solve_feast`:

```fortran
! was: if (M > 0 .and. M <= M0 .and. info >= 0) then
if (M > 0 .and. M <= M0 .and. (info == 0 .or. info == 2)) then
```

Now `converged` and `nev_found` are consistent: a non-converged solve yields
`nev_found = 0`. By the time this guard runs, `info = 3` can only mean the retry
loop exhausted with `M0 < N` (genuinely incomplete), so returning zero is honest.
No call-site changes. (`test_eigensolver.pf` is all converged cases; the
empty-matrix test already covers the `converged=.false., nev_found=0` shape.)

### A3 — bulk + FEAST permanent rejection  (`src/core/defs.f90`, `validate()`)

Next to check I15:

```fortran
if (trim(cfg%confinement) == 'bulk' .and. trim(cfg%solver%method) == 'FEAST') then
  error stop 'validate_simulation_config: bulk is always 8x8 dense; FEAST is not ' // &
    'supported for bulk (use method = AUTO or DENSE)'
end if
```

Enforces PRD invariant 1 ("bulk never offers FEAST; rejection is permanent").
AUTO→DENSE for bulk is unchanged.

### A4 — partial-window rejection  (`src/core/defs.f90`, `validate()`)

```fortran
if ((cfg%solver%emin == 0.0_dp) .neqv. (cfg%solver%emax == 0.0_dp)) then
  error stop 'validate_simulation_config: set both solver emin and emax, or neither ' // &
    '(0 = auto). A partial energy window is ambiguous.'
end if
```

`apply_solver_window` (OR override) is left as-is — now provably safe because
validate guarantees both-or-neither. The existing I11 (`emin < emax` when both
nonzero) stays.

### B2 — `kp_scalar_block` failure altitude  (`src/physics/hamiltonianConstructor.f90`)

Drop `pure` from `kp_scalar_block` and make its `case default` an `error stop`,
matching `kp_block_dense:727` exactly. (`error stop` is illegal in a `pure`
function; purity is valueless on an 8×8 / 52-entry lookup.)

### B3 — shrinking-N cache guard  (`src/math/eigensolver.f90:846`)

`if (N > self%cached_n)` → `if (N /= self%cached_n)`. Same-N (the hot sweep
path) still reuses buffers; any N change reallocates to exact size — correct for
both growth and shrinkage, eliminating the latent non-conformable assignment.

### C1 — remove dead code (fixes #14 **and** #11) 

- Delete `wire_setup_adopt_precomputed` + its `public` export
  (`src/physics/wire_setup.f90:127`).
- Delete `asw_evals` + its `module procedure` line in the `apply_solver_window`
  interface + the private `asw_apply_margin` helper
  (`src/math/eigensolver.f90:574,503`). **This also resolves #11**: the
  duplicated `margin_frac` / `margin_floor` constants lived only in
  `asw_apply_margin`, so they vanish with it. DRY-by-deletion over extracting
  constants with a single remaining user.

## Error-handling philosophy

All new guards are `error stop` with descriptive messages — matches the codebase
convention (ADR 0002: no silent corrections). The two `validate()` additions
(A3, A4) live where all config validation already lives.

## Testing

**Unit (pFUnit) — `tests/unit/test_eigensolver.pf`:**
- `eigensolver_result_band_slice`: full-spectrum case (`idx_lo = il`) and
  pre-sliced case (`idx_lo = 1`) return the expected contiguous slice.
- A2 (info=3 → `nev_found = 0`): best-effort — hard to reliably force `info=3`
  in a unit test; the empty-matrix test already covers the
  `converged=.false., nev_found=0` shape, so this is regression-protected by the
  caller behavior.

**Rejection tests (shell, exit-code) — `tests/integration/test_validate_rejects_bad_configs.sh`:**
- `bulk + method=FEAST` rejected (A3).
- `[solver] emin` set without `emax` (and vice-versa) rejected (A4).
- `eigensolver_result_band_slice` partial-spectrum `error stop` (A1) — a small
  driver program or an extension of the rejection harness.

**Integration (Python verify) — `tests/integration/`:**
- **New:** `verify_qw_bandstructure_dense_feast.py` — DENSE vs FEAST QW
  band-structure agree within tolerance (locks the A1 QW-CSR behavior change;
  mirrors `verify_qw_optics_dense_feast.py`, incl. the dispatch-diagnostic gate).
- **Regression (must stay green):** `verify_qw_optics_dense_feast.py`,
  `verify_qw_gfactor_dense_feast.py`, `verify_landau_dense_feast.py`,
  `verify_dense_sparse_all_geometries.py`. Golden regression outputs unchanged
  (AUTO/DENSE default paths untouched).

**Gate:** `123 + new` ctests, `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest -j4 --output-on-failure`.

## Acceptance criteria

1. `eigensolver_result_band_slice` exists, is used at all 3 sweep sites, and the
   QW band-structure FEAST path returns the gap-straddling `[il, iu]` bands.
2. No OOB read is reachable on a partial FEAST result (covered by the strict
   `error stop` + a rejection test).
3. `solve_feast` returns `nev_found = 0` on `info=3`; wire g-factor and SC loop
   no longer use a truncated spectrum.
4. `bulk + FEAST` and partial `[solver]` windows are rejected by `validate()`
   with clear messages.
5. `kp_scalar_block` error-stops on an unknown tag; dense solver cache is correct
   under shrinking N; `adopt_precomputed`, `asw_evals`, `asw_apply_margin` and
   the duplicated margin constants are gone.
6. All existing ctests + golden outputs unchanged; new verify/rejection tests
   pass.

## Out of scope (deferred — see BACKLOG)

- #9 — make `resolve_kp_term` truly drive the CSR formula.
- #12 — decompose `compute_strain_wire` (~520 lines).
- #15 — dedup the sweep-window derivation block across 5 sites.

---

## Post-Execution Outcomes (2026-06-21)

Executed via subagent-driven development + systematic debugging. **Final: 124/124 ctest green** (was 123; +1 new `verify_qw_bandstructure_dense_feast`); `tests/regression/data/` empty diff (golden AUTO/DENSE outputs unchanged). Commits `1f851c2`→`303d6a6` on `refactor/architecture-deepening`.

**Findings fixed:** #7, #5, #4, #1/#2/#10 (the `reconcile_band_slice` helper), #13, #8, #14.

**#3 (A2) — REVERTED as a non-bug** (`b045e7e`). The locked premise was false: callers requiring full convergence already gate on `result%converged` (`main_gfactor.f90:120`, `sc_loop.f90:215,732`), and FEAST `info=3`'s truncated spectrum holds the correct LOW bands (the wide auto/Gershgorin window spans the full spectrum). Zeroing `nev_found` on `info=3` broke the wire band-structure sweep and the BdG solve. Confirmed by single-variable experiments (both pass with the guard reverted). Restored `info >= 0` with an inline comment documenting the design decision. **#3 is dropped from this pass.** Design lesson: FEAST `info=3` (subspace too small) is common here; the truncated low-M spectrum is usable; do not zero `nev_found`.

**#11 — NOT resolved.** The duplicated `margin_frac`/`margin_floor` constants live in `asw_apply_margin`, which is LIVE — two `@test`s exercise it via the generic `apply_solver_window` interface (the plan's safety-gate grep missed generic-interface callers). Deleting it would break the build. Needs a constant-extraction refactor → BACKLOG. Consequently **C1's scope reduced**: only `wire_setup_adopt_precomputed` (genuinely dead) was deleted; `asw_evals`/`asw_apply_margin` kept.

**Task 7 (A1d) test follow-up** (`303d6a6`). `reconcile_band_slice` requires a full-spectrum window (AUTO envelope, `emin=emax=0`) + `m0` large enough to converge; narrow user ENERGY windows can't honor global-index `[il,iu]` extraction. Updated `verify_qw_sparse_solver` and `verify_dense_sparse_all_geometries` QW cases to the AUTO envelope, and removed `verify_qw_feast_truncation_warning` (it guarded the "only the lowest will be kept" warning Task 7 replaced with exact `[il,iu]` extraction). The QW DENSE path retains lowest-M (correct for DENSE+INDEX pre-slice); DENSE+ENERGY vs FEAST+ENERGY would now differ but DENSE+ENERGY is not a production path.

**Deferred remains:** #9, #12, #15 (BACKLOG), plus #11 (above).
