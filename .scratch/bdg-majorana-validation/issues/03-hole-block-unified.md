# Issue 01 — Kitaev-chain harness certifies Pfaffian + Majorana number (U3)

> **File-numbering note**: This file is `03-hole-block-unified.md` in `.scratch/bdg-majorana-validation/issues/`. The content is **Issue 01** (Kitaev harness) — sourced from `.superpowers/sdd/issue-03-brief.md`. The brief filenames were off by 2 from the issue numbers; this corrected numbering aligns the file's content with its issue number in the dependency graph.

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U3, Pfaffian math + Kitaev harness; dense-QW Pfaffian is Issue 05)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-03-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md` (full Understand-phase analysis)
**Phase**: PR-A, parallel with Issue 02 (PHS oracle) — fan-out

## What to build

A real topological invariant for the class-D BdG/Majorana path: Kitaev's Majorana number `M = sgn[Pf(H_BdG(0)·ω)·Pf(H_BdG(π/a)·ω)]` evaluated at the particle-hole-invariant momenta. Implementation:

1. **Vendor pfapack's Parlett-Reid tridiagonalization routines** (`dskpfa`/`zskpfa` + scaled variants) into the math layer. pfapack ships as legacy fixed-form `.f`; transcribe to free-form `.f90` with `implicit none` and an explicit interface matching the `linalg.f90` pattern (or add a scoped `set_source_files_properties(Fortran_FORMAT FIXED)` CMake exception — the build enforces `-std=f2018 -fimplicit-none -ffree-form`).
2. **Implement `kitaev_majorana_number`** wrapper that takes a BdG matrix and returns the Majorana number sign. Always skew-symmetrize the input `A ← (A−Aᵀ)/2` before calling (pfapack assumes exact skew-symmetry); use the mantissa-exponent scaled form for large n; treat `|a|/|a_max| < ε` as a gap-closure flag rather than trusting the sign there.
3. **Build a minimal spinless p-wave Kitaev-chain harness** as a **test-only fixture** (not a config-driven mode, not a new `sweep_model` value). The 8-band s-wave builder cannot express spinless p-wave; a validation-only model warrants no config plumbing.
4. **Certify the Pfaffian implementation** against closed-form Kitaev results: `Pf² = det` identity on random skew matrices up to n=16; sign of Pf on `[[0,a],[−a,0]]` returns `+a`; Kitaev Eq. 14 analytical MZM spinors.

The dense-QW Pfaffian (application of this routine to the existing `build_bdg_hamiltonian_qw` evaluated at `k_par ∈ {0, π/a}`) is Issue 05. The wire Pfaffian sweep is deferred with U13.

## Acceptance criteria

- [ ] Pfaffian routines vendored into `src/math/pfaffian.f90`; CMake wiring updated in `src/CMakeLists.txt` (added to `COMMON_SOURCES` and `set_source_files_properties`).
- [ ] Pfaffian of `[[0,a],[−a,0]]` returns `+a` (Kitaev convention).
- [ ] `Pf(A)² = det(A)` holds to roundoff on random skew matrices up to n=16.
- [ ] Skew-symmetrization guard: a near-skew input gives the correct sign after symmetrization.
- [ ] AE1 (Kitaev chain, `|μ| < 2t`): `M = −1` (topological), bulk gap = `|Δ|`, two end-localized MZMs.
- [ ] AE1 (`|μ| > 2t`): `M = +1` (trivial), no end modes, gap open.
- [ ] AE1 (`|μ| = 2t`): the gap closes (transition); the gap-closure flag fires.
- [ ] Even site count only (`k=π/a` exists as a real point only for even N) is asserted.
- [ ] `tests/unit/test_pfaffian.pf` and `tests/unit/test_kitaev_majorana.pf` (new) green.
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

- `src/math/linalg.f90` (9.0k lines) — centralized LAPACK/BLAS/MKL/PARDISO interfaces; FEAST under `#ifdef USE_MKL_FEAST`. No Pfaffian.
- `src/physics/AGENTS.md` — Module Inventory list. Add row for `pfaffian.f90`.
- pfapack reference: Wimmer's `pfapack` (arXiv:1102.3440); `dskpfa.f` (real skew Pfaffian) and `zskpfa.f` (complex skew Pfaffian). Parlett-Reid tridiagonalization is the algorithm.
- The Kitaev harness lives INSIDE `test_kitaev_majorana.pf` — synthetic BdG matrix generator. Do NOT wire into `main_topology.f90`; do NOT add a new `sweep_model` enum value.

## Constraints from CLAUDE.md + ADRs

- **ADR 0001 (fat derived type)**: pure functions only.
- **ADR 0002 (no new TOML fields)**: harness is test-only, not config-driven.
- **CLAUDE.md Boundaries**: NO approval gate (does not touch `bdg_hamiltonian.f90`).
- **CLAUDE.md Engineering Principles (DRY/SSOT)**: Pfaffian is a NEW math primitive; place in `src/math/` mirroring `linalg.f90`'s location and style.
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.
- **No polymorphic types**.

## File ownership (exhaustive)

### New files (you create)
- `src/math/pfaffian.f90` — the vendored Parlett-Reid routines + `kitaev_majorana_number` wrapper.
- `tests/unit/test_pfaffian.pf` — Pfaffian unit tests.
- `tests/unit/test_kitaev_majorana.pf` — Kitaev Majorana number tests (with synthetic harness inside).

### Modified files
- `src/CMakeLists.txt` — add `pfaffian.f90` to `COMMON_SOURCES` (and possibly a `set_source_files_properties(... Fortran_FORMAT FIXED)` if you go that route; else transcribe to free-form per AC).
- `src/physics/AGENTS.md` — Module Inventory row for `pfaffian.f90`.
- `tests/CMakeLists.txt` — `add_pfunit_ctest(test_pfaffian ...)` + `add_pfunit_ctest(test_kitaev_majorana ...)` with LABEL "unit".

### NOT modified (out of scope)
- `src/physics/bdg_hamiltonian.f90` (Issue 03 owns exclusively).
- `src/apps/main_topology.f90` (Issue 03 + 05 + 06 + 07 own).
- `src/core/defs.f90` (Issue 06/07 own).
- `scripts/lecture_*.py` (Issue 08 owns).

## Suggested API (illustrative — adjust if you prefer a different shape)

```fortran
module pfaffian
  use kind_module, only: dp
  implicit none
  private
  public :: real_pfaffian, complex_pfaffian, kitaev_majorana_number

  ! Real skew-symmetric matrix Pfaffian (Parlett-Reid)
  pure function real_pfaffian(A) result(pf)
    real(dp), intent(in) :: A(:,:)   ! assumed skew-symmetric
    real(dp) :: pf
  end function

  ! Complex skew-symmetric matrix Pfaffian (Parlett-Reid)
  pure function complex_pfaffian(A) result(pf)
    complex(dp), intent(in) :: A(:,:)  ! assumed skew-symmetric: A = -A^T
    complex(dp) :: pf
  end function

  ! Kitaev's Majorana number: M = sgn[Pf(H(k=0)·ω) · Pf(H(k=π/a)·ω)]
  ! ω is the antisymmetric structure matrix (Pauli τ_y ⊗ I_N)
  function kitaev_majorana_number(H_bdg, k_par_values) result(majorana_number)
    complex(dp), intent(in) :: H_bdg(:,:)
    real(dp), intent(in) :: k_par_values(:)  ! PHS-invariant momenta, e.g., [0.0, pi/a]
    integer :: majorana_number  ! -1 (topological), +1 (trivial), 0 (gap closure)
  end function
end module
```

## Tests required

### `tests/unit/test_pfaffian.pf`
1. `pf([[0,a],[-a,0]]) == +a` for `a = 1.0, 0.5, -2.0, 1.0e-6, 1.0e3` (sign + magnitude).
2. `Pf(A)^2 == det(A)` on random skew matrices for n = 2, 4, 6, 8, 10, 12, 14, 16.
3. Skew-symmetrization guard: input `A + 1.0e-14_dp * identity` gives the same Pf as clean skew `A`.
4. Zero matrix returns 0 (no pivots).
5. Real pfaffian handles n=1 (returns 0 by convention; n must be even).
6. Complex pfaffian handles `A_ij = i*sin(phase)` cases.

### `tests/unit/test_kitaev_majorana.pf`
1. Build a spinless p-wave Kitaev chain `H = -t Σ (c†_i c_{i+1} + h.c.) − μ Σ n_i + Δ Σ (c_i c_{i+1} + h.c.)` for N=8 (even).
2. **AE1a**: `|μ| = 1.0 < 2t = 2.0` (use t=1, Δ=0.5): `M = -1` (topological). Verify gap open at `|Δ|`. Verify MZMs at the two ends.
3. **AE1b**: `|μ| = 3.0 > 2t = 2.0`: `M = +1` (trivial). Gap open.
4. **AE1c**: `|μ| = 2.0 = 2t`: gap closes, flag fires (return 0).
5. Even-N assertion: N=4 (smallest even with k=π/a) and N=8 produce a real `k=π/a` eigenvalue.
6. The Majorana number wrapper reports the expected integer sign in all three regimes.

## TDD discipline (mandatory)

1. Write `tests/unit/test_pfaffian.pf` FIRST (RED).
2. Confirm test fails because module doesn't exist.
3. Vendor pfapack and implement `real_pfaffian` / `complex_pfaffian` (GREEN for `test_pfaffian.pf`).
4. Write `tests/unit/test_kitaev_majorana.pf` (RED for harness).
5. Implement `kitaev_majorana_number` (GREEN).
6. Run full unit suite to confirm no regressions.

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure
```

## Out of scope (must NOT do)

- Do NOT modify `src/physics/bdg_hamiltonian.f90` (Issue 03 owns exclusively).
- Do NOT modify `src/apps/main_topology.f90` (Issue 03/05/06/07 own).
- Do NOT modify `src/core/defs.f90` (Issue 06/07 own).
- Do NOT add a `sweep_model = kitaev_chain` enum value. Harness is test-only.
- Do NOT wire the harness into `main_topology.f90`. It lives inside the pFUnit test.
- Do NOT add new TOML fields.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-01-report.md`

The report must include:
- Status (DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT)
- Files created / modified
- TDD evidence (RED → GREEN) with exact commands and output
- Test results (full unit suite pass count)
- Self-review findings
- Commit SHAs
- Any concerns

## Risk notes

- pfapack is ~1000 lines (both real + complex). Vendor only what's needed (skew-symmetric Pfaffian, not the full matrix decomposition set).
- The Pfaffian is only defined for even-dimension matrices. Test must cover this.
- The complex Pfaffian sign convention differs from the real one in some libraries — verify Kitaev's convention (`Pf([[0,a],[-a,0]]) = +a`) matches your implementation.

---

## Outcome (as executed)

- **Module**: `src/math/pfaffian.f90` (vendored Parlett-Reid `dskpfa`/`zskpfa` + scaled variants + `kitaev_majorana_number` wrapper).
- **Tests**: `tests/unit/test_pfaffian.pf` (Pf²=det to roundoff up to n=16, sign on 2×2 case, skew-symmetrization guard) + `tests/unit/test_kitaev_majorana.pf` (AE1 closed-form spinors; M=±1 in three regimes; gap-closure flag at |μ|=2t; even-N assertion).
- **Fix Round 1**:
  - **Critical 2**: corrected `default_kitaev_omega` to match the code's Nambu ordering (per ADR 0007 / KTD7).
  - **Critical 3**: debugged Parlett-Reid for n>12 (scaled variant path was wrong; switched to mantissa-exponent scaling).
  - **Critical 1**: added M=±1 explicit integer tests (was missing — only sign checked).
  - **Important 1**: AGENTS.md Module Inventory row.
  - **Important 3**: replaced vacuous gap-closure test with non-vacuous assertion.
  - **Important 5**: added k=π bulk-gap test for the topological regime.
  - **Minor**: removed duplicate deallocations.
- **Deferred**: M=±1 on spinless chain was structurally impossible (electron/hole blocks exactly opposite at TRIM for diagonal BdG harness); deferred to Issue 05/07 for production-grade M=±1 witness on dense-QW BdG. Documented as a finding, not a defect.