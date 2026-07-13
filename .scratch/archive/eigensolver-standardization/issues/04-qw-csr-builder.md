# QW CSR builder — FEAST on QW end-to-end

**Type:** AFK

## Parent

PRD: `.scratch/eigensolver-standardization/PRD.md`

## What to build

Complete vertical slice enabling sparse FEAST solver on QW geometry: build a CSR-form QW Hamiltonian from 1D kpterms, integrate it into simulation_setup's QW path, and verify eigenvalues match the trusted dense builder.

**New module `hamiltonian_qw`** — Builds the QW Hamiltonian in CSR format using the same 52-entry kp block table as the dense and wire builders. The key insight: QW has only 1D spatial coupling (z-direction FD stencils), making R, RC, PP, PM purely diagonal (no scatter maps needed). The existing COO insertion layer from `hamiltonian_wire.f90` is fully grid-dimensionality-agnostic and reused directly.

- **`qw_workspace` type**: 12 pre-allocated CSR blocks (blk_Q through blk_A, blk_diff, blk_temp), COO buffers, COO-to-CSR cache, diagonal position maps for Q/T and A. No scatter maps (unlike wire). Fast path: updates CSR block values in-place on subsequent k-points, O(NNZ) rebuild.

- **10 kp-term assembly functions** (`build_kp_term_Q_1d` through `build_kp_term_A_1d`): each with slow path (first call, allocates) and fast path (cached, updates values). QW formulas derived from existing dense `ZB8bandQW`:
  - Q/T: tridiagonal FD stencil + diagonal k² term
  - S/SC/PZ: scalar × off-diagonal tridiagonal
  - R/RC/PP/PM: purely diagonal (kx/ky scalar multipliers, no FD stencil)
  - A: tridiagonal FD stencil + diagonal k² term

- **1D kpterms→CSR conversion**: Converts existing dense `kpterms(nz,nz,10)` tridiagonal/diagonal matrices to CSR by extracting 3 diagonals. O(N) per term.

- **`ZB8bandQW_csr` subroutine**: Calls the 10 assembly functions, then reuses COO insertion helpers from `hamiltonian_wire.f90` (insert_main_blocks, insert_profile_diagonal, insert_strain_coo, insert_zeeman_coo, finalize_coo_to_csr). Slow path on first call, fast path on subsequent k-points.

**Expose COO helpers** — Make `insert_main_blocks`, `insert_profile_diagonal`, `insert_strain_coo`, `insert_zeeman_coo`, `finalize_coo_to_csr`, and COO cache types public in `hamiltonian_wire.f90`.

**QW sparse path in simulation_setup** — When method is FEAST for QW (user explicitly sets `method = "FEAST"` for QW), use `ZB8bandQW_csr` instead of dense `ZB8bandQW`, build CSR, and dispatch through `solver%solve_sparse`.

**Build integration** — Add `physics/hamiltonian_qw.f90` to `COMMON_SOURCES` in `src/CMakeLists.txt`.

## Acceptance criteria

- [ ] `src/physics/hamiltonian_qw.f90` exists, compiles, and is linked into all executables
- [ ] `qw_workspace` type with init/free/cache; fast path works (workspace%initialized = .true. on second call)
- [ ] All 10 `build_kp_term_*_1d` functions produce correct CSR blocks
- [ ] `ZB8bandQW_csr` produces eigenvalues matching `ZB8bandQW` (dense) to 1e-10 tolerance at k=0 and k≠0
- [ ] COO insertion helpers are public in `hamiltonian_wire.f90`
- [ ] QW with `method = "FEAST"` in `[solver]` section works end-to-end (simulation_setup dispatches to sparse path)
- [ ] Unit test: `test_qw_csr_vs_dense` passes (eigenvalue comparison)
- [ ] Integration test: `verify_qw_sparse_solver.py` passes (QW with FEAST matches QW with DENSE)
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #02 (unified interface must exist so QW CSR can be consumed by FEAST via solve_sparse)
