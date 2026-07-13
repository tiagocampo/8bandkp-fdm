# Issue 6: Add type-mismatch warning helper for toml-f optional lookups

## Parent

PRD: Input System Architecture Cleanup (`.scratch/input-system-cleanup/PRD.md`)

## What to build

Add a private helper subroutine `check_optional_stat(stat, key, section)` to input_parser.f90 that prints a warning when `stat /= 0` on optional field lookups. Currently, when a user writes a wrong-type value (e.g., `linewidth_lorentzian = "thick"`), toml-f returns `stat /= 0` and the code silently falls through to the default. The helper provides visibility:

```fortran
subroutine check_optional_stat(stat, key, section)
  integer, intent(in) :: stat
  character(len=*), intent(in) :: key, section
  if (stat /= 0) then
    print *, 'Warning: [', trim(section), '] key ''', trim(key), &
      ''' has wrong type, using default'
  end if
end subroutine
```

Apply this helper to optional `get_value` calls (those with default values and `stat=stat`). Do NOT apply to required-field lookups (those use `require_string`/`require_int`/`require_real` which already error-stop). Do NOT apply to `requested=.false.` table/array lookups (those check `associated()` instead of stat).

This is an incremental improvement. Apply to the optional lookups within parse routines. The helper is small and the application is mechanical.

## Acceptance criteria

- [ ] `check_optional_stat` helper exists in input_parser.f90
- [ ] Warning message includes section name and key name
- [ ] Optional `get_value` calls with defaults call `check_optional_stat` after `get_value`
- [ ] Required-field lookuses (`require_*` helpers) are unchanged
- [ ] `requested=.false.` table/array lookups are unchanged
- [ ] A config with a type-mismatched optional key produces a warning on stderr/stdout
- [ ] A config with correct types produces no warnings
- [ ] All 34 unit tests pass
- [ ] All ~29 regression tests pass

## Blocked by

None — can start immediately. Lowest priority of the 6 issues.
