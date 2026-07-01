Status: ready-for-agent

## Parent

PRD: `.scratch/optics-engine-encapsulation/PRD.md`

## What to build

Add new module procedures that take `type(optics_engine), intent(inout) :: oe` as their first argument, replacing all module-level `save` variable access with `oe%` component access. Each new routine reads config from `oe%optcfg` instead of taking `optcfg` as a separate argument.

New routines to add:
- `optics_accumulate(oe, eigvals, eigvecs, k_weight, vel, numcb, numvb, fermi_level)` — reads `oe%nE`, `oe%E_grid`, `oe%optcfg%confinement`, etc.; writes `oe%alpha_te`, `oe%alpha_tm`, spin-resolved arrays
- `optics_accumulate_spontaneous(oe, eigvals, eigvecs, k_weight, vel, numcb, numvb, fermi_level)` — same pattern with `oe%spont_te/tm`
- `compute_gain_qw(oe, eigvals, eigvecs, k_weight, vel, numcb, numvb, carrier_density)` — reads/writes `oe%gain_fermi_computed`, `oe%mu_e`, `oe%mu_h`, `oe%alpha_gain_te/tm`
- `compute_isbt_absorption(oe, eigvals, eigvecs, vel, numcb, numvb, k_weight, fermi_level)` — writes `oe%alpha_isbt`
- `gain_reset(oe)` — resets `oe%gain_fermi_computed`, `oe%mu_e`, `oe%mu_h`

Split `optics_finalize` into two new routines:
- `optics_apply_prefactor(oe)` — applies the physical prefactor to all arrays. No file I/O. Uses `oe%optcfg%refractive_index` instead of passed argument.
- `optics_write_output(oe)` — writes all spectra files from `oe%` components.

The old routines become thin wrappers that delegate to the new ones (old `optics_accumulate(optcfg, ...)` calls new `optics_accumulate(oe, ...)` by reading module-level `save` vars into a temporary type, etc.). This keeps the old caller (`main_optics.f90`) functional until the migration in issue 03.

Write TDD tests with minimal fixtures: small eigenvalue/eigenvector arrays (e.g., 2 CB + 2 VB states for a bulk 8x8 system), CSR velocity matrices with known values. Verify: accumulate produces non-zero arrays, apply_prefactor scales values correctly, write_output produces files with correct content.

## Acceptance criteria

- [ ] All five accumulate/gain/ISBT routines have new overloads taking `oe` as first arg, reading config from `oe%optcfg`
- [ ] `optics_apply_prefactor(oe)` applies prefactor to all arrays without file I/O
- [ ] `optics_write_output(oe)` writes all spectra files from type components
- [ ] Old routines remain functional as thin wrappers delegating to new routines
- [ ] No module-level `save` variable access in new routines — all via `oe%` components
- [ ] pFUnit tests pass: accumulate with fixture data produces expected non-zero values, apply_prefactor scales correctly
- [ ] Existing test suite (96/98 baseline) still passes

## Blocked by

- Issue 01 (type skeleton + init/free must exist)
