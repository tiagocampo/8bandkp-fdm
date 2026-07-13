# Update documentation for Zeeman SSOT relocation and namespace changes

**Type:** AFK

## What to build

Update all project documentation to reflect the code changes from issues #01 and #02.

**CLAUDE.md changes:**
- Global Invariants: change Zeeman SSOT from `strain_solver.f90` → `magnetic_field.f90`
- Boundaries: update "NEVER change the strain or Zeeman block tables" to reference `magnetic_field.f90` for Zeeman
- Module dependency graph: `magnetic_field.f90` now exports `zeeman_entry`, `get_zeeman_table`, `init_zeeman_cache`; `strain_solver.f90` no longer exports Zeeman

**src/physics/AGENTS.md changes:**
- `magnetic_field.f90` entry: add Zeeman table to public API, note `compute_zeeman_vz` reads from table (not pure)
- `strain_solver.f90` entry: remove Zeeman from public API, add `lookup_bp_field`
- `finitedifferences.f90` entry: note 10 symbols now PRIVATE, list the 6 public ones
- `utils.f90` entry: add `get_unit` and `ensure_output_dir` to public API
- `outputFunctions.f90` entry: remove `get_unit` and `ensure_output_dir`, note they moved to `utils.f90`

Commit as: `docs: update CLAUDE.md and AGENTS.md for Zeeman SSOT and namespace changes`

## Acceptance criteria

- [ ] CLAUDE.md Zeeman SSOT invariant points to `magnetic_field.f90`
- [ ] CLAUDE.md Boundaries section references `magnetic_field.f90` for Zeeman
- [ ] CLAUDE.md dependency graph reflects the relocation
- [ ] `src/physics/AGENTS.md` reflects Zeeman in `magnetic_field` and not in `strain_solver`
- [ ] `src/physics/AGENTS.md` reflects `lookup_bp_field` in `strain_solver`
- [ ] No stale references to Zeeman table in `strain_solver.f90` remain in any documentation file

## Blocked by

- #01 (Zeeman relocation)
- #02 (namespace cleanup)
