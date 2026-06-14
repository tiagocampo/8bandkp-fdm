# Clean namespace pollution — narrow imports, hide unused FD exports, relocate file utilities (C4)

**Type:** AFK

## What to build

Three coordinated mechanical changes that reduce compilation coupling and clean up the module interface layer. Zero behavioral change. Do this after issue #01 so the `only:` lists reflect the updated Zeeman import locations.

**1. Make 10 unused FD symbols PRIVATE.** `finitedifferences.f90` exports 16 symbols but only 6 have external callers. Remove the 10 unused ones from `public ::` declarations. The module already uses `implicit private`, so simply removing them from the public list makes them private. The 6 that stay public: `Identity`, `buildFD2ndDerivMatrix`, `buildFD1stDerivMatrix`, `buildStaggeredD1Inner`, `buildStaggeredD1Outer`, `interpolateToHalfPoints`.

**2. Move `get_unit` and `ensure_output_dir` from `outputFunctions.f90` to `utils.f90`.** These are file-level utilities with no dependency on output formatting. Moving them breaks the backwards dependency where `input_parser.f90` and 4 physics modules imported `outputFunctions` solely for these two helpers. After the move, `outputFunctions.f90` imports them from `utils.f90`. All external importers (`optical_spectra.f90`, `scattering.f90`, `exciton.f90`, `main_optics.f90`, `input_parser.f90`) change their import target. The 4 broad `use outputFunctions` importers (`main.f90`, `main_gfactor.f90`, `main_topology.f90`) should no longer rely on getting `get_unit`/`ensure_output_dir` from there.

**3. Add `only:` to all 22 broad `use definitions` imports.** The heaviest user references only 16 of 88 public symbols — no pragmatic cutoff needed. For each file, determine the actual symbols used by grepping for `definitions` public symbols within the file. Short symbol names (`c`, `e`, `m0`) require careful verification to avoid false positives from local variables. Do this step last, after parts 1 and 2, so the `only:` lists incorporate all prior import changes.

Commit as: `refactor: namespace cleanup — only imports, private FD exports, file utils to utils`

## Acceptance criteria

- [ ] No module imports `definitions` without an `only:` clause (`grep -rn 'use definitions$' src/` returns nothing)
- [ ] `finitedifferences.f90` exports exactly 6 public symbols
- [ ] `get_unit` and `ensure_output_dir` live in `utils.f90`, not `outputFunctions.f90`
- [ ] No module imports `get_unit` or `ensure_output_dir` from `outputFunctions`
- [ ] `input_parser.f90` does not import `outputFunctions` (it imports from `utils` instead)
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #00 (branch must exist)
- #01 should be committed first (Zeeman import changes reflected in the `only:` lists)
