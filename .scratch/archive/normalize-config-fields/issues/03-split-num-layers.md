# Issue 3: Split `num_layers` — wire uses explicit `wire%num_regions`

**Type**: AFK
**Blocked by**: Issue 2 (same files: `defs.f90`, `input_parser.f90`)
**User stories**: 1, 2, 4, 6, 7

## What to build

Narrow `cfg%num_layers` to mean material layers only. Currently wire mode silently sets it to the region count (`nreg`), conflating two concepts. After this change:

- Bulk/landau: `num_layers = 1` (unchanged)
- QW: `num_layers = nmat` (count of `[[material]]` entries — unchanged)
- Wire: `num_layers = 1` (was `nreg`; regions stay in `cfg%wire%num_regions`)

Wire consumers that currently branch on `num_layers` must switch to explicit confinement-mode checks:
- `sc_loop.f90`: wire doping iteration uses `cfg%wire%num_regions`
- `topological_analysis.f90`: wire midpoint computed from `cfg%wire%num_regions`
- `simulation_setup.f90`: eigensolver choice (`num_layers == 1` → zheev) replaced with `confinement == 'bulk'` check
- `defs.f90` validation: wire array-size checks use `cfg%wire%num_regions`

## Acceptance criteria

- [ ] Wire parser branch sets `cfg%num_layers = 1` instead of `nreg`
- [ ] `sc_loop.f90` wire path iterates over `cfg%wire%num_regions` with explicit confinement guard
- [ ] `topological_analysis.f90` wire path computes midpoint from `cfg%wire%num_regions`
- [ ] `simulation_setup.f90` eigensolver choice branches on `confinement` instead of `num_layers == 1`
- [ ] `defs.f90` validation: wire array-size checks use `cfg%wire%num_regions`
- [ ] Unit test added for `validate_semantic` with wire config verifying region-based validation
- [ ] All 34+ unit tests pass
- [ ] All regression tests pass (`ctest --test-dir build`)
- [ ] No remaining wire-path code that uses `cfg%num_layers` to mean region count

## Blocked by

- Issue 2 (same files: `defs.f90`, `input_parser.f90`)
