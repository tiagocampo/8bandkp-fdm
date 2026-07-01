# Issue 1: Parser guard: reject fd_step < 2 before division-by-zero

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 19, 22
**GitHub issue**: #20

## What to build

Add a fail-fast guard in the QW material parser that rejects `fd_step < 2` immediately before the division `cfg%delta = cfg%totalSize / real(cfg%fd_step - 1, dp)` executes. Currently, `fd_step = 1` causes a division-by-zero crash. The existing `validate()` check (`fd_step >= 3` for QW) runs too late — after the division has already happened.

The guard uses `error stop` with a contextual message: `QW fd_step must be >= 2, got <value>`. This is in addition to the existing `validate()` check (belt-and-suspenders).

This is check P1 from ADR 0002 (`docs/adr/0002-consolidate-validation.md`).

## Acceptance criteria

- [ ] Parser guard added before the `cfg%delta` division in `parse_materials_qw`
- [ ] Guard uses `error stop` with contextual message including the actual value
- [ ] Existing `validate()` check for `fd_step >= 3` remains unchanged
- [ ] Config with `fd_step = 1` and `confinement = "qw"` triggers the guard at parse time
- [ ] Config with `fd_step >= 3` and `confinement = "qw"` passes the guard
- [ ] All 34 unit tests pass

## Blocked by

None — can start immediately.
