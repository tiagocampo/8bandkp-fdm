Status: ready-for-agent

## Parent

PRD: `.scratch/optics-engine-encapsulation/PRD.md`

## What to build

Define the `optics_engine` derived type in `optical_spectra.f90` and implement its creation and destruction lifecycle. The type has public allocatable components for all accumulation arrays, a stored `optics_config`, and scalar state for gain quasi-Fermi levels.

`optics_init(oe, optcfg)` stores the config, allocates the energy grid and accumulation arrays (unconditionally: `E_grid`, `alpha_te`, `alpha_tm`, `alpha_isbt`; conditionally based on config flags: gain, spontaneous, spin-resolved), and zeros everything. `optics_free(oe)` deallocates all arrays and resets scalars. A `final optics_engine_finalize` subroutine delegates to `optics_free`.

The old module-level `save` variables and old public routines remain untouched — both interfaces coexist at this stage. The new type is additive only.

Write pFUnit tests in a new test file that verify: init allocates and zeros arrays, free deallocates, double-free is safe, finalizer fires on scope exit, conditional arrays are allocated/omitted correctly based on config flags.

## Acceptance criteria

- [ ] `optics_engine` type defined with all components from the PRD (public allocatable arrays, stored config, scalar gain state)
- [ ] `optics_init(oe, optcfg)` stores config, allocates correct arrays based on config flags, zeros all arrays
- [ ] `optics_free(oe)` deallocates all arrays, resets `nE=0`, `gain_fermi_computed=.false.`, `mu_e=0`, `mu_h=0`
- [ ] Finalizer on `optics_engine` delegates to `optics_free`
- [ ] Double-free is safe (calling `optics_free` on already-freed type does nothing)
- [ ] pFUnit tests pass: init allocates, free deallocates, finalizer works, conditional allocation correct
- [ ] Old `save` variables and old public routines (`optics_init(optcfg)`, `optics_cleanup()`, etc.) remain untouched and functional
- [ ] Existing test suite (96/98 baseline) still passes

## Blocked by

None - can start immediately.
