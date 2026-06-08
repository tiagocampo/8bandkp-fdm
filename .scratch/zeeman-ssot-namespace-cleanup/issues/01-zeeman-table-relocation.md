# Fix Zeeman single-source-of-truth violation (C3)

**Type:** AFK

## What to build

The Zeeman diagonal coefficients are defined in two places: `get_zeeman_table()` in `strain_solver.f90` (the declared SSOT) and hard-coded in `compute_zeeman_vz()` in `magnetic_field.f90`. The BdG module uses the hard-coded copy, creating a sync hazard.

Move the entire Zeeman table infrastructure (type, builder, cached accessor, init subroutine) from `strain_solver.f90` to `magnetic_field.f90` — its domain-correct home (Zeeman is magnetic, not strain). Rewrite `compute_zeeman_vz` to read from the table instead of hard-coding. Update all import sites: `hamiltonianConstructor.f90`, `hamiltonian_wire.f90`, and `main.f90` must now import Zeeman symbols from `magnetic_field` instead of `strain_solver`.

Additionally, the Bir-Pikus `field_id` dispatch (7-way `select case`) is character-for-character duplicated between `apply_strain_table_dense` in `hamiltonianConstructor.f90` and `insert_strain_coo` in `hamiltonian_wire.f90`. Extract a shared `lookup_bp_field` pure function into `strain_solver.f90` (the module that defines the `field_id` values) and have both builders call it.

Key design decisions:
- `compute_zeeman_vz` loses its `pure` attribute (the table accessor reads a SAVE cache). Acceptable — its only caller (`bdg_hamiltonian.f90`) is non-pure, and the cache is pre-warmed before any OpenMP fork.
- No circular dependency risk: `strain_solver` does not import `magnetic_field`.
- `init_strain_cache` stays in `strain_solver.f90` (it's strain-related, not Zeeman).

Commit as: `refactor: move Zeeman table to magnetic_field and extract lookup_bp_field`

## Acceptance criteria

- [ ] `zeeman_entry` type, `build_zeeman_table`, `get_zeeman_table`, `init_zeeman_cache` all live in `magnetic_field.f90`, not `strain_solver.f90`
- [ ] `compute_zeeman_vz` reads from `get_zeeman_table()` instead of hard-coding coefficients
- [ ] `compute_zeeman_vz` no longer has `pure` attribute
- [ ] `hamiltonianConstructor.f90` and `hamiltonian_wire.f90` import Zeeman symbols from `magnetic_field`, not `strain_solver`
- [ ] `main.f90` imports `init_zeeman_cache` from `magnetic_field`, not `strain_solver`
- [ ] `lookup_bp_field` is a public pure function in `strain_solver.f90`, used by both dense and COO strain table appliers
- [ ] No duplicated `select case (field_id)` blocks remain in `hamiltonianConstructor.f90` or `hamiltonian_wire.f90`
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #00 (branch must exist)
