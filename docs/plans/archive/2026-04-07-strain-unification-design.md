# Strain Code Unification Design

**Date:** 2026-04-07
**Branch:** feature/docs-overhaul
**Scope:** Refactor Bir-Pikus strain to eliminate formula duplication, block topology duplication, and parameter threading.

## Problem

The Bir-Pikus strain code has three copies of the same physics:

1. **`ZB8bandBulk`** (hamiltonianConstructor.f90:1086-1178): Computes P_eps, Q_eps, R_eps, S_eps inline for a single grid point, inserts into 8x8 dense matrix. ~90 lines.
2. **`ZB8bandQW`** (hamiltonianConstructor.f90:808-875): Receives precomputed `bir_pikus_blocks`, inserts into 8Nx8N dense matrix. ~68 lines.
3. **`insert_strain_coo`** (hamiltonianConstructor.f90:2085-2293): Same as QW but into COO sparse format. ~209 lines.

Additionally, `bir_pikus_blocks` is threaded as an optional parameter through 4+ subroutines (ZB8bandQW, ZB8bandGeneralized, self_consistent_loop, self_consistent_loop_wire).

## Design

### 1. Formula unification: `bp_scalar` + `compute_bp_scalar`

Add a scalar struct and a pure computation function:

```fortran
! defs.f90 -- pure data, no deps
type :: bp_scalar
  real(dp) :: delta_Ec, delta_EHH, delta_ELH, delta_ESO
  complex(dp) :: R_eps, S_eps
  real(dp) :: QT2_eps
end type
```

```fortran
! strain_solver.f90 -- single source of truth
pure function compute_bp_scalar(params, eps_xx, eps_yy, eps_zz, &
    eps_xy, eps_xz, eps_yz) result(s)
  type(paramStruct), intent(in) :: params
  real(dp), intent(in) :: eps_xx, eps_yy, eps_zz, eps_xy, eps_xz, eps_yz
  type(bp_scalar) :: s
  ! P_eps, Q_eps, T_eps, R_eps, S_eps, diagonal shifts
end function
```

- `compute_bir_pikus_blocks` calls `compute_bp_scalar` per grid point, distributes fields to the 7 arrays.
- `ZB8bandBulk` calls `compute_bp_scalar` once instead of inline formulas.

### 2. Block topology helper: `add_bp_strain_dense`

Key insight: bulk (8x8) and QW (8Nx8N) use the same indexing formula. Band `b` at grid point `ii` is at row `b*N + ii`. For bulk, N=1 and ii=1.

```fortran
! hamiltonianConstructor.f90
subroutine add_bp_strain_dense(HT, ii, N, bp)
  complex(dp), intent(inout) :: HT(:,:)
  integer, intent(in) :: ii, N
  type(bir_pikus_blocks), intent(in) :: bp
  ! Inserts all 32 entries for grid point ii
end subroutine
```

Call sites:
- `ZB8bandBulk`: creates a 1-element `bir_pikus_blocks`, calls `add_bp_strain_dense(HT, 1, 1, bp)`.
- `ZB8bandQW`: `do ii = 1, N; call add_bp_strain_dense(HT, ii, N, bp); end do`.
- `insert_strain_coo` stays separate (COO format is fundamentally different).

### 3. Move `bir_pikus_blocks` into `simulation_config`

Move the type definition to `defs.f90` (no circular dep -- `strain_solver` uses `defs`, not vice versa). Add as a field of `simulation_config`:

```fortran
! defs.f90 -- in simulation_config
type(bir_pikus_blocks) :: strain_blocks
```

Remove `strain_blocks` optional parameter from:
- `ZB8bandQW` (signature simplified)
- `ZB8bandGeneralized` (signature simplified)
- `self_consistent_loop` (no more threading)
- `self_consistent_loop_wire` (no more threading)

Hamiltonian constructors access `cfg%strain_blocks` directly (they already receive `cfg`).

`main.f90` / `main_gfactor.f90` populate `cfg%strain_blocks` instead of a local variable.

### 4. Remove deprecated `apply_pikus_bir`

The old diagonal-only profile-modification approach is superseded by `compute_bir_pikus_blocks` + Hamiltonian insertion. It has a sign convention mismatch and is dead code in production. Remove it and update tests.

## Files changed

| File | Change |
|---|---|
| `src/core/defs.f90` | Add `bp_scalar` type, move `bir_pikus_blocks` from strain_solver, add field to `simulation_config` |
| `src/physics/strain_solver.f90` | Add `compute_bp_scalar`, refactor `compute_bir_pikus_blocks` to use it, remove `bir_pikus_blocks` type (now in defs), remove `apply_pikus_bir` |
| `src/physics/hamiltonianConstructor.f90` | Add `add_bp_strain_dense`, refactor bulk/QW strain insertion to use it, access `cfg%strain_blocks` instead of parameter |
| `src/physics/sc_loop.f90` | Remove `strain_blocks` parameter, access via `cfg` |
| `src/apps/main.f90` | Populate `cfg%strain_blocks` instead of local `bp`, remove from call sites |
| `src/apps/main_gfactor.f90` | Same as main.f90 |
| `tests/unit/test_strain_solver.pf` | Update tests for removed `apply_pikus_bir`, test `compute_bp_scalar` |

## Implementation order

1. **Add `bp_scalar` to defs.f90** -- no existing code changes, just new type
2. **Add `compute_bp_scalar` to strain_solver.f90** -- new pure function
3. **Move `bir_pikus_blocks` to defs.f90** -- move type, add to simulation_config, update strain_solver imports
4. **Refactor `compute_bir_pikus_blocks`** -- call `compute_bp_scalar` internally
5. **Add `add_bp_strain_dense` helper** -- new subroutine in hamiltonianConstructor
6. **Refactor `ZB8bandBulk`** -- use `compute_bp_scalar` + `add_bp_strain_dense`
7. **Refactor `ZB8bandQW`** -- use `add_bp_strain_dense`, access `cfg%strain_blocks`
8. **Refactor `ZB8bandGeneralized`** -- access `cfg%strain_blocks`
9. **Refactor sc_loop** -- remove parameter, access via cfg
10. **Update main.f90 / main_gfactor.f90** -- populate cfg%strain_blocks
11. **Remove `apply_pikus_bir`** -- delete function, update tests
12. **Build and test** -- verify all 24 tests pass
