## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Rewrite the k.p block insertion in `hamiltonian_wire` to read the k.p block table from `hamiltonian_blocks` instead of using 52 hardcoded `insert_csr_block*` calls. The COO builder iterates over the 52-entry table, gets the CSR block for each k.p term from the precomputed blocks, and inserts via the existing COO insertion mechanism.

Also handle the derivative mode (`insert_g3_blocks`, 20 entries for dH/dkz) using the same table with builder-side derivative rules.

Verify wire eigenvalues are unchanged: run wire regression tests and the wire convergence test (U6). Richardson-extrapolated values must agree within 1e-12 eV.

## Acceptance criteria

- [ ] `hamiltonian_wire` imports `get_kp_block_table` from `hamiltonian_blocks`
- [ ] The 52 `insert_csr_block*` calls in `insert_main_blocks` are replaced by a data-driven loop over the table
- [ ] `insert_g3_blocks` (derivative mode, dH/dkz) uses the same table with builder-side power-rule handling
- [ ] All 33 unit tests pass (including `test_hamiltonian_2d`)
- [ ] Wire regression tests pass with unchanged eigenvalues
- [ ] Wire convergence test (U6) Richardson values unchanged within 1e-12 eV

## Blocked by

- `.scratch/hamiltonian-block-table/issues/04-create-hamiltonian-blocks-module.md`
