# Create feature branch and verify clean build

**Type:** AFK

## What to build

Create the feature branch `refactor/zeeman-ssot-namespace-cleanup` from `main`. All subsequent issues in this campaign will be committed to this branch. After all issues are merged into the branch, a PR will be opened against `main`.

Configure with tests enabled, build, and run the full test suite to establish a clean baseline.

## Acceptance criteria

- [ ] Branch `refactor/zeeman-ssot-namespace-cleanup` exists off `main`
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

None — can start immediately.
