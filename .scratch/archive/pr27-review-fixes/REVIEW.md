# REVIEW — PR #27 Review Fixes (Batch 2)

**Status: COMPLETE** — audited 2026-06-22 against `feat/bdg-u8-window-routing`; final open item (Issue 03/I2 `b_field%components` silent fallback) closed by PR40 C1 follow-up commit on 2026-06-27. Archived to `.scratch/archive/pr27-review-fixes/` on 2026-06-27.
**Overall:** 4 of 4 issues fully done.

> This PRD was reviewed against the current working tree. The completed parts are
> listed for traceability; only the open item is actionable.

## What is DONE (verified in HEAD)

- **Issue 01 — thread-safety (C1):** `init_kp_block_cache` / `init_strain_cache` /
  `init_zeeman_cache` are called before all four OpenMP forks in `main.f90`
  (Landau FEAST envelope `:461-463`, Landau k-sweep `:508-510`, QW FEAST `:841-843`,
  QW DENSE `:890-892`). Init subroutines defined in `hamiltonian_blocks.f90:149`,
  `strain_solver.f90:142`, `magnetic_field.f90:74`.
- **Issue 02 — validation gaps:** bulk `band_idx` guard (`defs.f90:880-900`, branches
  on `which_band`, applies to all confinement modes); `fermi_mode` typo now
  `error stop` (`input_parser.f90:503-504`); rejection tests V9/V10 added
  (`test_validate_rejects_bad_configs.sh:139-174`).
- **Issue 03 — partial (see open item):**
  - I1 done — 0× `stop 1` remain in `input_parser.f90` (24× `error stop`).
  - I2 polygon done — `check_optional_stat` on `x`/`y` (`input_parser.f90:301,303`).
  - I4 done — `type(kp_entry) :: table(52)` stack-allocated (`hamiltonian_wire.f90:1029,1151`).
- **Issue 04 — Codex P1 responses:** inline review comments posted on PR #27
  ("False positive — investigated and confirmed safe" on `b_field` copy and
  `strainSubstrate`), per `gh api`.

## What is LEFT (open item)

### Issue 03 / I2 — b_field components must error on type mismatch

The PRD (`issues/03-parser-cleanup.md`, "I2 b_field", lines 25-30) requires the three
B-field component extractions to use `check_optional_stat` so a type mismatch errors
instead of silently defaulting to 0.0 (User Story #6).

**Current state — still the silent fallback.** `src/io/input_parser.f90:451-456`:

```fortran
call get_value(comp_arr, 1, cfg%b_field%components(1), stat=stat)
if (stat /= 0) cfg%b_field%components(1) = 0.0_dp
call get_value(comp_arr, 2, cfg%b_field%components(2), stat=stat)
if (stat /= 0) cfg%b_field%components(2) = 0.0_dp
call get_value(comp_arr, 3, cfg%b_field%components(3), stat=stat)
if (stat /= 0) cfg%b_field%components(3) = 0.0_dp
```

Note the inconsistency: the very next line, `g_factor` (`:459`), already uses
`check_optional_stat`. The sibling deliverable I2 polygon is also done. Only the
three `components` lines were missed.

**Adaptation note:** the PRD's illustrative snippet assumes separate `Bx`/`By`/`Bz`
keys, but the actual parser reads a `components` **array** (indices 1/2/3). Keep the
array access; only swap the error handling. Suggested fix:

```fortran
call get_value(comp_arr, 1, cfg%b_field%components(1), 0.0_dp, stat=stat)
call check_optional_stat(stat, 'components[1]', 'b_field')
call get_value(comp_arr, 2, cfg%b_field%components(2), 0.0_dp, stat=stat)
call check_optional_stat(stat, 'components[2]', 'b_field')
call get_value(comp_arr, 3, cfg%b_field%components(3), 0.0_dp, stat=stat)
call check_optional_stat(stat, 'components[3]', 'b_field')
```

(Check whether a default-value form of `get_value` for array elements exists; if not,
keep the `0.0_dp` default and let `check_optional_stat` gate it, matching the
`g_factor` precedent on `:458-459`.)

**Verification to add:**
- All existing tests still pass (`ctest --test-dir build`).
- Optionally add a rejection case (e.g. V12) in
  `tests/integration/test_validate_rejects_bad_configs.sh` with a non-numeric
  `components` entry to lock in the new error behavior (mirrors V9/V10).
