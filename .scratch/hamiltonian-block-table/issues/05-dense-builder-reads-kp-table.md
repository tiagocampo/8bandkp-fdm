## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Rewrite the k.p block construction in `ZB8bandQW` to read the k.p block table from `hamiltonian_blocks` instead of using 52 hardcoded HT() assignments. The dense builder iterates over the 52-entry table, looks up each k.p term value from the precomputed `kpterms` arrays (populated by `confinementInitialization`), applies the prefactor and conjugation, and writes to HT(row_band*N + ii, col_band*N + ii).

The builder handles derivative modes (dH/dk for g-factor velocity operators): when in derivative mode, apply the power rule to the k.p term before emitting. Only entries with non-zero k-dependence contribute in derivative mode.

Also rewrite `ZB8bandBulk` to use the same table for the 8x8 bulk Hamiltonian (scalar variant, no spatial grid).

Verify eigenvalues are unchanged: run QW and bulk regression tests and the convergence test suite (U4, U5). Richardson-extrapolated values must agree within 1e-12 eV.

## Acceptance criteria

- [ ] `ZB8bandQW` imports `get_kp_block_table` from `hamiltonian_blocks` and uses it for k.p block construction
- [ ] The 52 hardcoded HT() assignments in `ZB8bandQW` are replaced by a data-driven loop over the table
- [ ] `ZB8bandBulk` also uses the table for bulk Hamiltonian construction
- [ ] Derivative modes (g='g', g='g3') work correctly via builder-side power-rule handling
- [ ] All 33 unit tests pass
- [ ] QW regression tests pass with unchanged eigenvalues
- [ ] Convergence test suite (U4, U5) Richardson values unchanged within 1e-12 eV

## Blocked by

- `.scratch/hamiltonian-block-table/issues/04-create-hamiltonian-blocks-module.md`
