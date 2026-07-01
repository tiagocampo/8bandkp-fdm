# 05 — Confinement-assembly unit tests (C7 gate)

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 3, C7). This issue is the
**hard gate** for #06.

## What to build

The finite-difference k.p term assembly and the block-tag → block-formula
interpretation are only tested transitively through downstream eigenvalue
checks. A subtly wrong but still-Hermitian assembly could shift eigenvalues past
tolerance. Add **direct unit tests** that establish a safety net before any
assembly refactor:

- Known FD stencils; diagonal material-parameter placement; profile/params
  construction for single-material and heterostructure configs.
- Hermiticity of assembled terms.
- FD-order 2/4/6 consistency.
- **Direct block-tag → formula tests** asserting the derived identities
  (difference `Q − T`, half-sum `0.5·(Q + T)`) at the source — not via
  eigenvalues.

These tests must pass against the **current, un-refactored** code, establishing
the net before #06 touches anything.

## Acceptance criteria

- [ ] New confinement-assembly unit-test file covering FD stencils, material
      placement, Hermiticity, and FD-order consistency.
- [ ] Direct block-tag → formula tests assert the difference and half-sum
      identities at the source.
- [ ] All new tests pass against current un-refactored code.

## Blocked by

None — can start immediately.
