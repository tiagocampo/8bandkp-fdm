# ISBT InGaAs/InAlAs Material Sweep — 3 Widths × QCL-Relevant System

**Type:** AFK
**Blocked by:** #04 (ISBT GaAs/AlGaAs width sweep)
**GitHub:** #34

## What to Build

Extend the ISBT benchmark to cover a second material system: Ga₀.₄₇In₀.₅₃AsW / Al₀.₄₇In₀.₅₃AsW (InP lattice-matched, Winkler parameters). This is a real quantum cascade laser material system with a narrower gap than GaAs/AlGaAs, producing stronger non-parabolicity and larger deviation from the infinite-well formula.

Create 3 new TOML configs (50, 100, 200 Å) using the existing `qw_ingaas_algaas_strained_optics.toml` as template, flipping ISBT=true and adjusting E_min/E_max for the narrower gap. Include `[strain]` section with reference to InP substrate (lattice-matched, so strain is zero, but the section activates strain calculation for the optics pipeline).

Extend `verify_isbt_benchmark.py` to handle the second material. Update the deviation figure to show both materials on the same axes — the InGaAs/InAlAs curve should show larger deviation from infinite-well than GaAs/AlGaAs at every width. Add COVERAGE annotations for the InGaAs/InAlAs ISBT cells.

## Acceptance Criteria

- [ ] 3 new TOML configs created: 50, 100, 200 Å Ga₀.₄₇In₀.₅₃AsW / Al₀.₄₇In₀.₅₃AsW QWs with ISBT enabled and [strain] ref=InP
- [ ] `verify_isbt_benchmark.py` extended with InGaAs/InAlAs material sweep
- [ ] All 6 data points (2 materials × 3 widths) pass quantitative assertions
- [ ] Figure updated: both materials on same axes, showing material-dependent deviation
- [ ] InGaAs/InAlAs ISBT cells added to `validation_universe.yml` with COVERAGE annotations
- [ ] ctest registration covers all 6 configs
- [ ] All existing tests still pass

## Blocked by

- #04 (ISBT GaAs/AlGaAs width sweep must exist first)

## User Stories Covered

- #4: Material dependence of ISBT deviation (InGaAs/InAlAs narrow gap)
- #6: InGaAs/InAlAs uses existing Ga₀.₄₇In₀.₅₃AsW / Al₀.₄₇In₀.₅₃AsW parameter sets
- #7: Extended benchmark catches regressions
