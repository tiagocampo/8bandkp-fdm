## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Replace all three Bir-Pikus strain implementations in the dense Hamiltonian builder with a single table-driven subroutine that reads `get_strain_table()` from `strain_solver`. Currently the dense module has three variants: (1) an inline copy in `ZB8bandQW` done to avoid compiler temporaries with -O3+OpenMP, (2) a standalone `add_bp_strain_dense` subroutine, and (3) a scalar `apply_bp_strain_inline` for bulk. All three encode the same 32 Bir-Pikus entries as manual HT() assignments.

Write a new `apply_strain_table_dense` subroutine that iterates over the 32-entry strain table, looks up the field value from the `bir_pikus_blocks` type, applies the prefactor and conjugation, and writes to HT(). Use scalar field lookups (not array temporaries) to avoid the original -O3+OpenMP stack corruption. Replace all three call sites. Delete the old `add_bp_strain_dense` and `apply_bp_strain_inline` subroutines.

The `hamiltonianConstructor` module already imports `compute_bp_scalar` and `bir_pikus_blocks_free` from `strain_solver`. Add imports for `strain_entry`, `get_strain_table`.

## Acceptance criteria

- [ ] `apply_strain_table_dense` subroutine reads `get_strain_table()` and applies all 32 entries to the dense matrix
- [ ] The inline strain block in `ZB8bandQW` is removed and replaced by a call to `apply_strain_table_dense`
- [ ] `add_bp_strain_dense` subroutine is removed
- [ ] `apply_bp_strain_inline` (bulk scalar variant) is removed or replaced
- [ ] All 33 unit tests pass (`ctest -L unit`)
- [ ] QW and bulk regression tests pass (no physics regression)
- [ ] OpenMP test: `OMP_NUM_THREADS=4 ctest -L unit` passes (no stack corruption from temporaries)
- [ ] A targeted unit test compares old vs new strain output element-wise to machine precision

## Blocked by

None - can start immediately
