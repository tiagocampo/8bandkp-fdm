# ISBT GaAs/AlGaAs Tracer Bullet — Single-Point z-Dipole vs Infinite-Well

**Type:** AFK
**Blocked by:** None — can start immediately
**GitHub:** #29

## What to Build

A thin end-to-end tracer bullet for the ISBT oscillator strength benchmark. Create a single TOML config for a 100 Å GaAs / Al₀.₃Ga₀.₇As QW with `[optics] ISBT = true`. Run the `opticalProperties` executable. Parse `output/isbt_transitions.dat` to extract the e1→e2 z-dipole `|⟨1|z|2⟩|` and oscillator strength `f₁₂`. Compare against the infinite square well analytical formulas:

- z-dipole: `|⟨1|z|2⟩|_inf = 16L / (9π²)` where L is the well width
- oscillator strength: `f₁₂ = (2m₀E₁₂ / ℏ²) |⟨1|z|2⟩|²`

Assert that the 8-band result is within a physically justified tolerance of the analytical value (expected ratio < 1 due to finite barrier penetration and band mixing). Generate a figure showing the ISBT absorption spectrum with the analytical transition energy marked. Register as ctest target `verification_isbt_benchmark` under labels `verification;standard-star`. Promote the ISBT cell from `aspirational` to `required` in `validation_universe.yml` (GaAs/AlGaAs QW). Add COVERAGE annotations.

This is the tracer bullet — once this works, the width sweep and material sweep are incremental extensions.

## Acceptance Criteria

- [ ] TOML config created for 100 Å GaAs/Al₀.₃Ga₀.₇As QW with ISBT enabled
- [ ] Verification script `verify_isbt_benchmark.py` runs `opticalProperties`, parses ISBT output, compares z-dipole and f₁₂ against infinite-well analytical
- [ ] Quantitative assertion passes: |z₁₂|_8band / |z₁₂|_inf < 1 with documented tolerance
- [ ] Figure generated: ISBT absorption spectrum with analytical transition energy marker
- [ ] ctest target `verification_isbt_benchmark` registered under labels `verification;standard-star`
- [ ] ISBT cell promoted from `aspirational` to `required` in `validation_universe.yml`
- [ ] COVERAGE annotation added to the verification script
- [ ] All existing 111 tests still pass

## Blocked by

None — can start immediately.

## User Stories Covered

- #1: Compare ISBT z-dipole against infinite-well analytical
- #2: Compare oscillator strength against analytical formula
- #7: ISBT benchmark catches regressions in optical pipeline
- #8: ISBT cell promoted to required in coverage matrix
