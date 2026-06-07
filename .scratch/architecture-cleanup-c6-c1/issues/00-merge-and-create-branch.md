# Merge Phase 21 and Create Architecture Cleanup Branch

**Type:** AFK
**Blocked by:** None — can start immediately

## What to Build

Merge the current `feat/publishable-benchmarks-phase21` branch (8 commits: validation + review fixes, PR #35) into `main`. Then fork a new branch `refactor/architecture-cleanup` from the merged `main`. Verify all 113 tests pass on the new branch.

## Acceptance Criteria

- [ ] `feat/publishable-benchmarks-phase21` merged into `main` via PR #35 (or direct merge if already approved)
- [ ] New branch `refactor/architecture-cleanup` created from `main`
- [ ] `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl` succeeds on the new branch
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all 113 tests

## Blocked by

None — can start immediately.

## User Stories Covered

- #10: simulation_setup handles all 4 modes (prerequisite)
