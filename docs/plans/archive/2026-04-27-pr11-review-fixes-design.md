# PR #11 Review Fixes Design

Date: 2026-04-27
Branch: `feature/hamiltonian-performance-refactor`

## Problem

PR #11 introduces `wire_workspace` fast-path COO assembly that eliminates per-kz CSR allocation overhead. The fast path for terms S, PP, and R uses element-wise array arithmetic on kpterms with **different sparsity patterns**, triggering gfortran F2003 `realloc-lhs` that silently corrupts CSR structures.

Empirical confirmation: on an 11x11 grid, `blk_S%values` shrinks from 480 to 242 entries after the first fast-path call. `nnz` and `rowptr` still reference 480 entries.

Additional issues: uninitialized `diag_pos`, silent COO overflow, `prev_wire_eval` UB on zero eigenvalues, dead workspace fields, no fast-path test coverage.

## Fix: Scatter Maps

### Core idea

During workspace initialization (slow path), pre-compute integer index arrays that map each kpterm's nonzero positions to the correct positions in the union-structure workspace block. In the fast path, scatter values through these maps instead of doing element-wise arithmetic on differently-shaped arrays.

### Scatter map construction

`build_scatter_map(src, dst, map)` iterates through `src`'s CSR. For each nonzero at `(row, col)`, binary-searches `dst`'s sorted columns in that row. Complexity: O(NNZ_src * log nnz_per_row), done once at init. Fails with `stop 1` if entry not found (dst is a superset, so this should never happen).

### Workspace type changes

Add 6 integer arrays to `wire_workspace`:

```fortran
integer, allocatable :: scatter_S_14(:)   ! kp14 -> blk_S
integer, allocatable :: scatter_S_15(:)   ! kp15 -> blk_S
integer, allocatable :: scatter_PP_12(:)  ! kp12 -> blk_PP
integer, allocatable :: scatter_PP_13(:)  ! kp13 -> blk_PP
integer, allocatable :: scatter_R_16(:)   ! kp16 -> blk_R
integer, allocatable :: scatter_R_11(:)   ! kp11 -> blk_R
```

### Fast-path rewrite for S (pattern for PP, R)

```fortran
blk%values = 0.0_dp
blk%values(ws%scatter_S_14) = coeff1 * kpterms_2d(14)%values
blk%values(ws%scatter_S_15) = blk%values(ws%scatter_S_15) + coeff2 * kpterms_2d(15)%values
```

### Conjugate-transpose terms: eliminate blk_temp dependency

SC, RC, PM now read from their forward block directly instead of using `blk_temp` as scratch:

- `build_kp_term_SC`: conjugate-transpose `ws%blk_S` -> `blk`
- `build_kp_term_RC`: conjugate-transpose `ws%blk_R` -> `blk`
- `build_kp_term_PM`: conjugate-transpose `ws%blk_PP` -> `blk`

This means `build_kp_term_SC` no longer calls `build_kp_term_S` internally. Fast-path call order becomes: build Q -> T -> S -> R -> PZ -> PP -> A, then SC = conj(S), RC = conj(R), PM = conj(PP).

`blk_temp` is only written once per kz-point for `0.5*(Q+T)` where sparsity matches.

## Additional Fixes

### C2: COO overflow -> `stop 1`
Replace `print *, "WARNING"` + `return` with `stop 1` in `insert_csr_block`, `insert_csr_block_scaled`, `insert_profile_diagonal`, and `insert_strain_coo`. The nnz estimate is exact + 20% buffer; overflow means something is fundamentally wrong.

### C3: Initialize diag_pos
```fortran
allocate(ws%diag_pos(N))
ws%diag_pos = 0
do i = 1, N
  ...
end do
! Assert all found
if (any(ws%diag_pos == 0)) then
  print *, 'ERROR: diag_pos: missing diagonal entry'
  stop 1
end if
```

### I2: Guard prev_wire_eval assignment
Move lines 375-376 inside the `if (eigen_res%nev_found > 0)` block.

### I3: Remove dead coo_to_csr / coo_nnz_in fields
Remove declarations from type, deallocation from free.

### I4/I5: feast_workspace hardening
Add `N` field and check `fw%N == N` alongside `fw%M0 == M0`.

### Suggestions
- `case default; error stop` in `insert_strain_coo` select-case
- Cache `build_strain_table()` as module-level save variable
- Guard workspace double-init: `if (ws%initialized) error stop`
- Remove dead `nnz_est` from ZB8bandGeneralized outer scope

## Test Plan

1. **Fast-path vs slow-path comparison**: Build wire Hamiltonian at kz=0.05 twice (slow, then fast). Assert identical CSR element-by-element.

2. **Multi-kz consistency**: Build at kz1 (slow, init), then kz2 (fast). Verify fast-path matches fresh slow-path at kz2.

3. **All existing tests pass**: 15 unit + 22 regression. Wire golden data needs regeneration since the old fast-path produced wrong results.

4. **Strain table coverage**: Extend test to spot-check one entry per group (diagonal, S_eps, R_eps, VB-SO coupling).

## Files Modified

- `src/physics/hamiltonian_wire.f90`: scatter maps, fast-path rewrites, diag_pos init, COO overflow, dead field removal
- `src/math/eigensolver.f90`: feast_workspace N validation
- `src/apps/main.f90`: prev_wire_eval guard
- `tests/unit/test_hamiltonian_2d.pf`: fast-path comparison tests
- `tests/unit/test_strain_solver.pf`: extended strain table coverage
