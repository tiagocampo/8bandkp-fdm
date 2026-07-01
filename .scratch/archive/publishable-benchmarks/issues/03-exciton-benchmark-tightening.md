# Exciton Benchmark Tightening — Miller ≤15% + Bastard 1982 + Universe Registration

**Type:** AFK
**Blocked by:** None — can start immediately
**GitHub:** #31

## What to Build

Tighten the existing exciton binding energy benchmark and register it in the validation coverage matrix. Three changes:

1. **Tighten Miller tolerance:** The existing `test_exciton_convergence.py` compares Richardson-extrapolated Eb against Miller et al. 1985 (`Eb_ref = 9.0 meV`) at 30% tolerance. The convergence data suggests the actual agreement is much better. Tighten the `MILLER_TOLERANCE` constant from 0.30 to 0.15.

2. **Add Bastard 1982 comparison:** Add a second quantitative check comparing Eb against Bastard PRB 1982 tabulated values for GaAs QWs with Al₀.₃Ga₀.₇As barriers. Bastard provides Eb at specific well widths — use the 100 Å value as the primary reference. Expected tolerance: ≤20% (the Bastard model uses a simpler effective-mass approach).

3. **Register in universe:** Add `exciton_Eb` as a new observable in `validation_universe.yml` metadata.observables. Add a required cell: `exciton_Eb / QW / GaAs/AlGaAs / required` with references to both Miller 1985 and Bastard 1982.

## Acceptance Criteria

- [ ] `test_exciton_convergence.py` Miller tolerance tightened from 0.30 to 0.15
- [ ] Bastard 1982 quantitative reference added with ≤20% tolerance assertion
- [ ] Both Miller and Bastard checks pass at 100 Å GaAs/Al₀.₃Ga₀.₇As QW
- [ ] `exciton_Eb` observable added to `validation_universe.yml` metadata.observables
- [ ] Required cell added: `exciton_Eb / QW / GaAs/AlGaAs / required`
- [ ] COVERAGE annotation added or verified in the convergence test
- [ ] All existing 111 tests still pass (the tightened tolerance must not break CI)

## Blocked by

None — can start immediately.

## User Stories Covered

- #16: Miller tolerance tightened to ≤15%
- #17: Bastard 1982 quantitative comparison added
- #18: Well-width trend covered (leverages lecture_14 infrastructure)
- #20: Tightened tolerance catches regressions
- #21: exciton_Eb registered in coverage matrix
