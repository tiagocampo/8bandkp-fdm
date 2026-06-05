## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Fix the wire convergence test (`test_wire_convergence.py`) so it works with both absolute and relative build directory paths. Currently the test crashes with `FileNotFoundError` when invoked as `python3 test_wire_convergence.py build-nolto .` because `star_helpers.run_exe` constructs the executable path as `build_dir/src/bandStructure`, which fails for relative paths. Other convergence tests (U4, U5, U7, U8) handle this correctly — follow their pattern.

After fixing, run the test directly and verify it produces non-empty JSON output with CB1_energy and gz Richardson extrapolation results for S7 (InAs wire at 4 grid levels: 16x16, 21x21, 26x26, 31x31).

## Acceptance criteria

- [ ] `python3 tests/integration/test_wire_convergence.py build-nolto .` completes without FileNotFoundError
- [ ] `python3 tests/integration/test_wire_convergence.py $(pwd)/build-nolto $(pwd)` also works (no regression)
- [ ] `tests/integration/convergence_results/wire_convergence.json` contains non-empty CB1_energy results with Richardson extrapolation, GCI, and max_observed_rate
- [ ] `tests/integration/convergence_results/wire_convergence.json` contains non-empty gz results
- [ ] `ctest -L convergence -R wire` passes (no regression)

## Blocked by

None - can start immediately
