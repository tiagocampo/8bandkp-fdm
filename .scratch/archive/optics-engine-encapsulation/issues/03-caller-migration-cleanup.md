Status: ready-for-agent

## Parent

PRD: `.scratch/optics-engine-encapsulation/PRD.md`

## What to build

Migrate `main_optics.f90` from the old interface to the new `optics_engine` type. Then remove all old code: `save` variables, old routine signatures, old public exports.

Caller changes in `main_optics.f90` (all three confinement branches — bulk, QW, wire):
- Add `type(optics_engine) :: oe` declaration
- Replace `call optics_init(cfg%optics)` with `call optics_init(oe, cfg%optics)`
- Replace `call optics_accumulate(cfg%optics, ...)` with `call optics_accumulate(oe, ...)` (and same for spontaneous, gain, ISBT)
- Replace `call optics_finalize(cfg%optics)` with `call optics_apply_prefactor(oe)`
- For QW exciton block: replace `E_grid` with `oe%E_grid`, `alpha_te` with `oe%alpha_te`, `alpha_tm` with `oe%alpha_tm`, `nE` with `oe%nE`
- After exciton block (or directly if no exciton): add `call optics_write_output(oe)`
- Replace `call optics_cleanup()` with `call optics_free(oe)`

After caller is migrated, remove from `optical_spectra.f90`:
- All 15 `save` variable declarations
- Old `optics_init(optcfg)` (single-arg version)
- Old `optics_cleanup()`
- Old `optics_finalize(optcfg)`
- Old accumulate/gain/ISBT wrappers
- `public :: E_grid, alpha_te, alpha_tm, nE` exports
- `public :: optics_cleanup` and `public :: optics_finalize` exports

Verify `exciton.f90` needs no changes — it takes arrays as arguments.

Run the full regression suite and verify the 96/98 baseline still holds.

## Acceptance criteria

- [ ] `main_optics.f90` uses `type(optics_engine) :: oe` throughout all three confinement branches
- [ ] Pipeline in each branch: `optics_init(oe, cfg%optics)` → accumulate loop → `optics_apply_prefactor(oe)` → optional exciton → `optics_write_output(oe)` → `optics_free(oe)`
- [ ] QW exciton block accesses arrays via `oe%E_grid`, `oe%alpha_te`, `oe%alpha_tm`, `oe%nE`
- [ ] All 15 `save` variable declarations removed from `optical_spectra.f90`
- [ ] Old single-arg `optics_init`, `optics_cleanup`, `optics_finalize`, and old wrapper routines removed
- [ ] Old public exports (`E_grid`, `alpha_te`, `alpha_tm`, `nE`) removed from public statement
- [ ] `exciton.f90` has zero changes
- [ ] Full regression suite passes (96/98 baseline)
- [ ] Build succeeds with `-std=f2018` enforcement

## Blocked by

- Issue 02 (accumulate routines + finalize split must exist)
