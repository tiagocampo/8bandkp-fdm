## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Extract the 8-entry Zeeman diagonal into a data-driven table in `strain_solver`, and make both the dense and COO paths use it. Currently the Zeeman splitting values (HH=-1.5, LH=+0.5, SO=-0.5, CB=+1.0) are encoded inside `compute_zeeman_vz` as an implicit mapping. The wire module already has `insert_zeeman_coo` as a data-driven loop over 8 bands. The dense module has three inline copies: in `ZB8bandQW`, `ZB8bandBulk`, and `ZB8bandLandau`.

Add a `zeeman_entry` type to `strain_solver` (alongside `strain_entry`) with fields: `band_index` (0-7), `g_multiplier` (the spin-dependent factor). Add `get_zeeman_table()` returning the 8-entry table. Both dense and COO Zeeman insertion routines read this table. `compute_zeeman_vz` remains as the computation routine but its band-dependent values are now visible as data.

## Acceptance criteria

- [ ] `zeeman_entry` type defined in `strain_solver` with `band_index` and `g_multiplier` fields
- [ ] `get_zeeman_table()` returns the 8-entry table (cached, same pattern as strain table)
- [ ] Dense Zeeman insertion in `ZB8bandQW` uses the table instead of inline values
- [ ] Dense Zeeman insertion in `ZB8bandBulk` uses the table
- [ ] COO Zeeman insertion in `insert_zeeman_coo` uses the table (replaces hardcoded Vz loop)
- [ ] All 33 unit tests pass
- [ ] g-factor regression tests pass (Zeeman is critical for g-factor calculation)

## Blocked by

None - can start immediately
