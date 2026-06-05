# Issue 1: Add `conf_direction()` function and remove `conf_dir` field

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 5, 6

## What to build

Add a `conf_direction()` pure function to `defs.f90` that maps the confinement string ('bulk', 'qw', 'wire', 'landau') to a direction character ('n', 'z', 'z', 'x'). Then remove the `conf_dir` field from `simulation_config` entirely — replace all `cfg%conf_dir` references across the four executables and physics modules with calls to `conf_direction(cfg%confinement)`. Remove the four mode-branch assignments in `input_parser.f90` that set `cfg%conf_dir`.

The function signature:
```
pure function conf_direction(conf) result(d)
  character(len=*), intent(in) :: conf
  character(len=1) :: d
  select case(trim(conf))
  case('bulk');   d = 'n'
  case('qw','wire'); d = 'z'
  case('landau'); d = 'x'
  end select
end function
```

This is the first of three normalizations eliminating mode-dependent poly-fields from `simulation_config`.

## Acceptance criteria

- [ ] `conf_direction()` pure function added to `defs.f90`, exported as `public`
- [ ] `conf_dir` field removed from `simulation_config` type
- [ ] All `cfg%conf_dir` references replaced with `conf_direction(cfg%confinement)` across `main.f90`, `main_gfactor.f90`, `main_optics.f90`, `main_topology.f90`, `simulation_setup.f90`, `confinement_init.f90`, and any other consumers
- [ ] Parser assignments `cfg%conf_dir = ...` removed from `input_parser.f90`
- [ ] All 34+ unit tests pass
- [ ] All regression tests pass (`ctest --test-dir build`)
- [ ] No remaining references to `conf_dir` in `src/` (grep confirms zero hits)

## Blocked by

None — can start immediately.
