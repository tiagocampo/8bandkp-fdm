# Delete ARPACK + unguard LDOS/PARDISO

**Type:** AFK

## Parent

PRD: `.scratch/eigensolver-standardization/PRD.md`

## What to build

Remove all ARPACK infrastructure. ARPACK (~280 lines) is guarded by `#ifdef USE_ARPACK` and never compiled in CI — it requires `find_package(ARPACK)` which is never satisfied. Delete the `arpack_solver_t` type, `solve_arpack` subroutine, `solve_arpack_dispatch`, ARPACK cases in `make_eigensolver` and `solve_sparse_evp`, ARPACK interfaces (`znaupd`, `zneupd`) in linalg.f90, CMake detection (`find_package(ARPACK)`, `FindARPACK.cmake`), and ARPACK unit tests. The LDOS computation (`compute_ldos_csr`) uses only PARDISO, not ARPACK — unguard `pardiso_c` and `compute_ldos_csr` from `#ifdef USE_ARPACK` so they are always available.

This is pure subtraction with one unguard. Zero physics impact. All existing tests continue to pass.

## Acceptance criteria

- [ ] No `#ifdef USE_ARPACK` remains in any source file
- [ ] No `arpack_solver_t`, `solve_arpack`, or `solve_arpack_dispatch` exists
- [ ] No `'ARPACK'` case in `make_eigensolver` or `solve_sparse_evp`
- [ ] `pardiso_c` interface is public in `linalg.f90` without preprocessor guard
- [ ] `compute_ldos_csr` is public in `green_functions.f90` without preprocessor guard
- [ ] No `find_package(ARPACK)` in any CMakeLists.txt; `cmake/FindARPACK.cmake` deleted
- [ ] No test references ARPACK (`test_arpack_fallback` deleted, `cfg%method = 'ARPACK'` replaced with `'DENSE'`)
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #00 (branch must exist)
