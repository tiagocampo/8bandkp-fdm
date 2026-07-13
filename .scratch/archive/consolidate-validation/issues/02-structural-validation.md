# Issue 2: Consolidate structural validation into validate()

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 1, 2, 3, 4, 5, 6, 7, 8, 21
**GitHub issue**: #21

## What to build

Add 8 config-agnostic checks to the existing `simulation_config_validate` type-bound procedure in `defs.f90`. These checks require no app context — they apply regardless of which executable runs. All use `error stop` with contextual messages including the actual value and what was expected.

Checks to add (V1–V8 from ADR 0002):

| ID | Check | Error message pattern |
|----|-------|-----------------------|
| V1 | `evnum > 8` for bulk confinement | `evnum (=X) exceeds 8 for bulk confinement (8x8 Hamiltonian)` |
| V2 | `num_cb > NUM_CB_STATES * grid%npoints()` for QW | `num_cb (=X) exceeds maximum (Y) for QW confinement` |
| V3 | `num_vb > NUM_VB_STATES * grid%npoints()` for QW | `num_vb (=X) exceeds maximum (Y) for QW confinement` |
| V4 | Wave vector mode is one of: kx, ky, kz, kxky, kxkz, kykz, k0 | `wave_vector mode 'X' not recognized (expected one of: kx, ky, kz, kxky, kxkz, kykz, k0)` |
| V5 | B-sweep step > 0 when B-sweep is present | `B_sweep step must be positive, got X` |
| V6 | `z_min(i) < z_max(i)` for all material layers | `material layer I: z_min (X) >= z_max (Y)` |
| V7 | SC + bulk confinement is invalid | `SC loop requires confinement='qw', got 'bulk'` |
| V8 | Electric field requires `z(1) /= 0` (guard with `allocated(cfg%z)`) | `electric field requires z(1) /= 0` |

These checks replace the old behavior of silently capping values (V1–V3) or printing warnings (V7). All are now hard errors.

Key implementation notes:
- V1: Only fires when `conf_direction(cfg%confinement) == 'n'` (bulk). The current code silently caps `evnum` to 8 and forces `num_cb=2, num_vb=6` — this silent correction is eliminated.
- V2–V3: `NUM_CB_STATES` and `NUM_VB_STATES` are named constants from `defs.f90` (values 2 and 6).
- V4: The wave vector mode check is currently duplicated in `main.f90` and `main_gfactor.f90` with different valid sets. Consolidate to the full set of 7 modes.
- V5: Only check when B-sweep is present (array allocated and size >= 3).
- V6: New check — `z_min >= z_max` for any material layer is a crash-level risk. Loop over `cfg%num_layers`.
- V7: Currently a warning in `simulation_setup.f90` — upgrade to hard error.
- V8: Guard with `allocated(cfg%z)` — bulk mode may not allocate the z array. Only check when external_field is electric type.

## Acceptance criteria

- [ ] All 8 checks added to `simulation_config_validate` in `defs.f90`
- [ ] Each check uses `error stop` with contextual message including actual value
- [ ] Invalid configs that previously passed now fail at parse time with clear messages
- [ ] Existing 19 checks in `validate()` remain unchanged
- [ ] All 34 unit tests pass (existing configs are valid, so new checks do not fire)

## Blocked by

None — can start immediately.
