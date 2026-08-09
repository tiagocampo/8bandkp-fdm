# P7 — TOML validation for new writer-path keys (IMPORTANT, defense-in-depth)

## Question (task, AFK)

T2 added two config-driven writer paths: `[topology] phase_diagram_file` and `slim_pfaffian_witness_file` (`src/core/defs.f90:339-340`, parser at `src/io/input_parser.f90:641-654`). Both default to canonical values; the sibling fixture overrides them.

The fields are `character(len=64)`. No validation guards against:

1. **Empty string** `phase_diagram_file = ""`: writer attempts `output/' // trim(adjustl('')) = 'output/'`. The `open()` call fails loud (`error stop 'cannot open phase_diagram.dat'`) but the error message references the wrong file name — confusing for users.

2. **Whitespace string** `phase_diagram_file = "   "`: same as empty.

3. **Absolute path** `phase_diagram_file = "/tmp/foo.dat"`: writer concatenates `output/' // '/tmp/foo.dat'` = `output//tmp/foo.dat` — broken path. The user's intent is silently ignored.

4. **Path traversal** `phase_diagram_file = "../escape"`: writer emits to `output/../escape` — silently escapes output dir.

5. **Length truncation** > 64 chars: `trim(...)` silently truncates.

Per `CLAUDE.md` "no silent corrections" rule, these need explicit `error stop` in `validate_semantic(cfg, app_name)`.

## Acceptance

- [ ] Add `validate_semantic` checks for `cfg%topo%phase_diagram_file` and `cfg%topo%slim_pfaffian_witness_file`:
  - Reject empty / whitespace-only strings.
  - Reject absolute paths (leading `/`).
  - Reject path traversal (`..`).
  - Reject paths exceeding 64 chars (or extend the field length if the truncation risk is too aggressive).
  - Each check produces `error stop 'phase_diagram_file: <reason>'` with descriptive message.
- [ ] Add a rejection-class shell test `tests/integration/test_topology_writer_paths_validation.sh` that drives `topologicalAnalysis` with each malformed input and asserts exit code != 0 (per `feedback_pfunit_error_stop_uncatchable`).
- [ ] Update `tests/unit/test_topology_parser.pf` with two new `@test` cases pinning the rejection behavior.

## Cross-references

- T2 ticket: `tickets/T2_mu_window.md` (the field addition)
- ADR 0002: `docs/adr/0002-consolidate-validation.md` (consolidated validation pattern)
- CLAUDE.md: "All config validation via `validate()` + `validate_semantic()` in `defs.f90` — no silent corrections"

## Blocks / blocked-by

- **Blocks**: nothing.
- **Blocked by**: none.

## Type

task (AFK).

## Branch

`feat/bdg-u10-pfaffian-phase-diagram` (PR #43 review follow-up).

## Claimed by

Unclaimed (2026-08-08 review pass).
