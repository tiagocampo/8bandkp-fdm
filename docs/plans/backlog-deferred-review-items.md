---
title: Backlog — Deferred Review Items from PR #27
date: 2026-05-30
status: backlog
source: /tmp/8bandkp-fdm-pr27-review-handoff.md
---

# Deferred Review Items

Items identified during the PR #27 review that were deferred to future PRs.

## Architectural Debt

| ID | Description | Why Deferred | Effort |
|----|-------------|-------------|--------|
| I1 | `main.f90` bypasses `simulation_setup` for Landau + k-sweep | Works correctly; refactor is architectural debt, not a bug | Medium |
| I2 | `defs.f90` at ~1,039 lines — 3x the 300-line guideline | Natural split point unclear; don't split without a clear decomposition | Medium |
| I3 | Duplicate confinement string checks in multiple modules | Low priority — works correctly, just repetitive | Low |
| I4 | `grid_ngrid`/`grid%npoints()` redundancy | Legacy wrapper still called from 8 source files; migration is gradual | Low |

## Testing Gaps

| ID | Description | Why Deferred | Effort |
|----|-------------|-------------|--------|
| I17 | No wire optics path tested in any unit test | Requires full wire optics setup; regression tests cover integrated path | High |
| I18 | `hamiltonian_blocks` has no dedicated unit test | Block table consumed by builders; integration tests exercise it indirectly | Medium |
| I19 | Convergence suite doesn't cover all physics modules | Expanding coverage incrementally per feature | Medium |
| I20 | ISBT test assertions were weak — improved by C8 but could be stronger | C8 added proper physics tests; further hardening is incremental | Low |

## Code Cleanup (Backlog)

| ID | Description | Why Deferred | Effort |
|----|-------------|-------------|--------|
| I13-ext | ~60 `stop 1` remaining in physics/math/io modules | Non-executable files; bulk sed replacement. Done for `simulation_setup.f90` + 4 executables in PR #27. | Low |

## Completed in PR #27

See git log for commits addressing: C1–C12, I5–I12, I15–I16, I21, I23–I27, I6, I13, I14 (verified already handled), I22.
