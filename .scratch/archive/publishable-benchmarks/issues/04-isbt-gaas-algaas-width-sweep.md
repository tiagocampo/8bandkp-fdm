# ISBT GaAs/AlGaAs Width Sweep — 50Å and 200Å Configs + Deviation Figure

**Type:** AFK
**Blocked by:** #01 (ISBT GaAs/AlGaAs tracer bullet)
**GitHub:** #32

## What to Build

Extend the ISBT tracer bullet to cover 50 Å and 200 Å GaAs/Al₀.₃Ga₀.₇As QWs. Create two new TOML configs (the 100 Å config already exists from the tracer bullet). Extend `verify_isbt_benchmark.py` to loop over all three widths. Update the figure to show `|z₁₂|_8band / |z₁₂|_inf` vs well width — the key publishable figure showing how the 8-band deviation from the infinite-well formula depends on confinement.

The deviation should increase for narrower wells (stronger band mixing, more wavefunction penetration into barriers).

## Acceptance Criteria

- [ ] Two new TOML configs created: 50 Å and 200 Å GaAs/Al₀.₃Ga₀.₇As QWs with ISBT enabled
- [ ] `verify_isbt_benchmark.py` extended to loop over all three widths (50, 100, 200 Å)
- [ ] All three width points pass the quantitative assertion
- [ ] Figure updated: ratio vs well width plot (x-axis: width in Å, y-axis: |z₁₂|_8band / |z₁₂|_inf)
- [ ] ctest registration updated to cover all three configs
- [ ] All existing tests still pass

## Blocked by

- #01 (ISBT GaAs/AlGaAs tracer bullet must exist first)

## User Stories Covered

- #3: Benchmark covers multiple well widths (50, 100, 200 Å)
- #5: Benchmark table shows ratio vs well width for GaAs/AlGaAs
