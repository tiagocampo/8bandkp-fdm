## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Create a new `hamiltonian_blocks.f90` module in `src/physics/` that defines the 8-band zinc-blende k.p block structure as data. This is the foundation for Phase C: both the dense builder and the COO builder will import this table.

Define named constants for k.p terms: `KP_Q`, `KP_R`, `KP_S`, `KP_T`, `KP_P`, `KP_PZ`, `KP_PM`, `KP_PP`, `KP_A`, `KP_DIAGONAL`. Define a `kp_entry` type with fields: `row_band`, `col_band` (0-indexed, 0-7), `kp_term` (named constant), `prefactor` (complex), `use_conjg` (logical). Implement `get_kp_block_table()` returning the full 52-entry table covering all 8x8 band-pair interactions. Include 8 band-offset entries as `KP_DIAGONAL` entries.

The module depends only on `defs.f90` (for `dp` kind). No physics module dependencies.

Write pFUnit tests verifying: entry count is 52, all band pairs are covered, no duplicate (row, col) pairs for non-Hermitian entries, Hermitian symmetry holds (entry (i,j) has conjugate prefactor to (j,i)), and every `kp_term` constant appears at least once.

## Acceptance criteria

- [ ] New file `src/physics/hamiltonian_blocks.f90` with `private` default and explicit `public` exports
- [ ] Named constants for all k.p terms (KP_Q, KP_R, KP_S, KP_T, KP_P, KP_PZ, KP_PM, KP_PP, KP_A, KP_DIAGONAL)
- [ ] `kp_entry` type with row_band, col_band, kp_term, prefactor, use_conjg
- [ ] `get_kp_block_table()` returns 52 entries (cached)
- [ ] Module depends only on `defs.f90`
- [ ] CMakeLists.txt updated to compile the new module
- [ ] New pFUnit test file with tests for entry count, band-pair coverage, Hermitian symmetry, no duplicates, kp_term coverage
- [ ] All existing 33 unit tests still pass (no import chain breakage)

## Blocked by

None - can start immediately (new module with no existing consumers)
