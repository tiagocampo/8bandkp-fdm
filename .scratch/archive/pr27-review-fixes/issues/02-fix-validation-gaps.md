# Issue 02: Fix band_idx bulk gap + fermi_mode silent fallback + rejection tests

**Priority:** Critical (C2 + C3)
**Commit:** `fix(validation): extend band_idx check to bulk, reject unknown fermi_mode`
**Files:** `src/core/defs.f90`, `src/io/input_parser.f90`, `tests/integration/test_validate_rejects_bad_configs.sh`

## Problem

**C2:** S1 check in `validate_semantic` (defs.f90:814) guards `band_idx` range check with `confinement == 'wire' .or. confinement == 'qw'`, excluding bulk. Bulk gfactor accesses `cb_state(:, bandIdx+1)` — `band_idx=2, num_cb=2` passes validation but crashes at runtime.

**C3:** `fermi_mode` `case default` in `parse_sc` (input_parser.f90:501-502) silently falls back to `charge_neutrality` instead of `error stop`. Violates ADR 0002 "no silent corrections".

## Fix

- [ ] **C2:** In `defs.f90:814`, remove the confinement guard so `band_idx` check applies to all modes:
  ```fortran
  ! Remove this line:
  if (trim(cfg%confinement) == 'wire' .or. trim(cfg%confinement) == 'qw') then
  ```
  Keep the range check `if (cfg%band_idx < 1 .or. cfg%band_idx + 1 > cfg%bands%num_cb)` without the guard.
  Remove the corresponding `end if`.

- [ ] **C3:** In `input_parser.f90:501-502`, replace `case default; cfg%sc%fermi_mode = 0` with:
  ```fortran
  case default
    error stop 'parse_sc: unknown fermi_mode ''' // trim(fermi_mode_str) // &
      '''. Valid values: charge_neutrality, fixed'
  ```

- [ ] **Test C2:** Add rejection test case — TOML with `confinement='bulk'`, `which_band=0`, `band_idx=2`, `num_cb=2`, run `gfactorCalculation`, expect non-zero exit + "bandIdx" in output

- [ ] **Test C3:** Add rejection test case — TOML with `[sc]` section containing `fermi_mode = 'typo'`, run `bandStructure`, expect non-zero exit + "fermi_mode" in output

## Verification

- [ ] New rejection tests pass (non-zero exit + expected message)
- [ ] Existing validation tests still pass
- [ ] Existing gfactor regression tests still pass
