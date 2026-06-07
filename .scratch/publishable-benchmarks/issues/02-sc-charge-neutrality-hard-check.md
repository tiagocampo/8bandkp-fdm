# SC Charge Neutrality Hard Check — Conservation Law Benchmark

**Type:** AFK
**Blocked by:** None — can start immediately
**GitHub:** #30

## What to Build

A thin end-to-end tracer bullet for the SC Schrödinger-Poisson benchmark suite. Reuse the existing GaAs/AlAs QW config from `test_sc_convergence.py` (FDstep=51, ND=5e18 cm⁻³ in the GaAs well, SC tolerance=1e-8). Run `bandStructure` with SC enabled. Parse `output/sc_charge.dat` to extract n_e(z). Integrate via trapezoidal rule. Assert charge neutrality as a hard pass/fail check:

```
|∫n_e dz - ∫N_D dz| / |∫N_D dz| < 1%
```

where `∫N_D dz = N_D × L_well`. This is an exact conservation law — the SC solver must satisfy it. Register as ctest target `verification_sc_benchmark` under labels `verification;standard-star`. Add `charge_neutrality` observable to `validation_universe.yml` as a required cell (GaAs/AlAs QW). Add COVERAGE annotations.

## Acceptance Criteria

- [ ] Verification script `verify_sc_benchmark.py` runs `bandStructure` with SC enabled, parses charge density, integrates n_e(z)
- [ ] Charge neutrality assertion passes: `|∫n_e - ∫N_D| / |∫N_D| < 1%`
- [ ] ctest target `verification_sc_benchmark` registered under labels `verification;standard-star`
- [ ] `charge_neutrality` observable cell added to `validation_universe.yml` as required (QW, GaAs/AlAs)
- [ ] COVERAGE annotation added to the verification script
- [ ] All existing 111 tests still pass

## Blocked by

None — can start immediately.

## User Stories Covered

- #9: Charge neutrality enforced as hard pass/fail (< 1%)
- #14: SC benchmark catches regressions in charge conservation
- #15: charge_neutrality registered in coverage matrix
