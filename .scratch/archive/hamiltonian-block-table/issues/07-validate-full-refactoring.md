## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Run the full convergence test suite (U4: QW grid, U5: QW order, U6: wire, U7: SC, U8: exciton) after all phases are complete. Compare Richardson-extrapolated values against the pre-refactoring baseline to verify no physics regression. Document any changes, particularly the known S5 (InAs/GaSb broken-gap) negative convergence rates and SC CB1_shift 40% GCI failure as calibration gaps (not caused by C4).

Run the full regression test suite and all 33 unit tests. Run with OpenMP (`OMP_NUM_THREADS=4`) to verify no stack corruption from the table-driven approach.

## Acceptance criteria

- [ ] All convergence tests pass: `ctest -L convergence` green
- [ ] All 33 unit tests pass: `ctest -L unit` green
- [ ] All regression tests pass: `ctest -L regression` green
- [ ] OpenMP test: `OMP_NUM_THREADS=4 ctest -L unit` green
- [ ] S4 (GaAs/AlGaAs) Richardson values unchanged within 1e-12 eV for CB1_energy, subband_spacing, gz
- [ ] S6 (InAs/GaAs strained) Richardson values unchanged within 1e-12 eV
- [ ] Wire (S7) Richardson values unchanged within 1e-12 eV
- [ ] SC and exciton Richardson values unchanged within 1e-12 eV
- [ ] S5 negative convergence rates documented as known calibration gap (pre-existing)
- [ ] SC CB1_shift 40% GCI documented as known calibration gap (pre-existing)

## Blocked by

- `.scratch/hamiltonian-block-table/issues/02-dense-strain-adopts-strain-table.md`
- `.scratch/hamiltonian-block-table/issues/03-dense-zeeman-adopts-zeeman-table.md`
- `.scratch/hamiltonian-block-table/issues/05-dense-builder-reads-kp-table.md`
- `.scratch/hamiltonian-block-table/issues/06-coo-builder-reads-kp-table.md`
