# Issue 6: Update docs for consolidated validation

**Type**: AFK
**Blocked by**: Issue 4, Issue 5
**User stories**: 20, 21
**GitHub issue**: #25

## What to build

Mechanical documentation updates to reflect the completed validation consolidation. No judgment calls — all changes are precisely specified below.

### 1. ADR 0002 status flip

In `docs/adr/0002-consolidate-validation.md`, change:
- `Status: Proposed` → `Status: Accepted`

### 2. CLAUDE.md — add validation bullet to Architecture section

In the "Key design concepts" section, after the bullet that describes `validate()` and `validate_semantic()` (or after the `confinementInitialization` bullet if no validation bullet exists), add:

```markdown
- **Config validation** consolidated in `validate()` (19 structural + 8 new checks) and `validate_semantic(cfg, app_name)` (5 existing + 10 new app-specific checks) in `defs.f90`. Executables only handle runtime errors (LAPACK info codes, FEAST convergence, file I/O). All config-level checks use `error stop` with contextual messages. No silent corrections. See ADR 0002 (`docs/adr/0002-consolidate-validation.md`) for the full checklist.
```

### 3. CLAUDE.md — remove stale "validation split" language

Search CLAUDE.md for any references to "validation split across" executables, or silent correction/warning behavior for `evnum`, `num_cb`, `num_vb`, or `fd_step`. If found, remove or replace with "rejected by `validate()`." If not found, skip this step.

### 4. input-reference.md — check for stale silent-correction docs

In `docs/reference/input-reference.md`, search for the terms "cap", "warn", "silently", or "correct" in the context of `evnum`, `num_cb`, `num_vb`, `fd_step`, or `B_sweep`. If any such language exists describing silent correction behavior, replace with "rejected at parse time by `validate()`." If not found, skip this step.

## Acceptance criteria

- [ ] ADR 0002 status changed from `Proposed` to `Accepted`
- [ ] One new bullet added to CLAUDE.md Architecture section describing consolidated validation
- [ ] No stale references to silent corrections remain in CLAUDE.md or input-reference.md
- [ ] All 34+ unit tests pass (documentation changes should not affect tests)

## Blocked by

- Issue 4 (remove old checks from executables)
- Issue 5 (pFUnit tests for new validators)
