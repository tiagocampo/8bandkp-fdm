# T3 — Config-driven `nk_par` / `bloch_sweep` knob + dispatch

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch;
> unblocks T4. Materialized from the map's T3 entry.

Type: task
Status: resolved
Owner: tiago
Blocked by: 02
Resolved: 2026-08-09

## Question

Surface the Bloch-periodic sweep the builder (T1) and seam (T2) now support as a
**config knob** so the topology executable can opt into the strict S1×S2 product
without disturbing the fixed-kz default path. Scope isolation is the non-negotiable
invariant: `nk_par=1` (or the knob absent) must reproduce the fixed-kz
build & sweep **bit-for-bit** before and after U13.

Contract a session should nail down (drawing on the map's §State and T6's L2
decision):

- **TOML surface.** Add `nk_par` (integer, default `1`) to the `[topology]`
  section (and/or a `bloch_sweep` boolean switch — whichever T6's L2 decision
  resolves as the cleaner lever). Optional sections are enabled by presence —
  no separate enable flag (`CLAUDE.md` §Input file). The parser entries go in
  `src/io/input_parser.f90`; the derived-type field in `defs.f90`'s
  `topology_config` (read `src/core/defs.f90` for the exact sub-type name before
  editing — `defs.f90` derived types are a **Require-approval** boundary per
  `CLAUDE.md` Boundaries).
- **Validation.** Add the structural + semantic checks to `validate()` /
  `validate_semantic(cfg, app_name)` in `defs.f90` per the consolidated-validation
  SSOT (ADR 0002, `docs/adr/0002-consolidate-validation.md`): `nk_par >= 1`,
  integer, and (if L2 ties it to a BZ period) bounded by the wire supercell
  period. `error stop` with a contextual message — no silent corrections.
- **Dispatch.** `main_topology.f90` branches on `nk_par`: `== 1` → existing
  fixed-kz path (the three call sites `:480, :718, :1346` the map §State
  records); `> 1` → T1's Bloch stack → T2's strict seam. The dispatch is where
  the scope-isolation regression guard lives: a test asserting the `nk_par=1`
  output matches the pre-U13 fixed-kz output byte-identical.
- **Scope-isolation regression guard (the real deliverable here).** A
  regression test (`tests/regression/` or `tests/integration/`, per the AGENTS
  there) that runs the canonical wire-BdG fixture at `nk_par=1` and diffs
  against a golden pre-U13 output. This guard is what lets U13 ship without
  re-validating every existing wire-BdG result. **It must fail** if the build
  path diverges, which is exactly the regression we want to catch.

Output: the knob, its validation, the dispatch, and the `nk_par=1`
bit-for-bit regression test. TDD (`superpowers:test-driven-development`): write
the regression test **red** against current output, wire the knob, confirm
green, then add a `nk_par>1` smoke. Read `src/physics/AGENTS.md` and
`tests/integration/AGENTS.md` before modifying.

## Answer

- **TOML surface:** `nk_par` (int, default 1), `k_par_min` (real, default 0.0), `k_par_max` (real, default 0.0) in `[bdg]`. Optional keys; defaults preserve the U10 bit-exact fixed-kz path. Uniform grid: when `nk_par > 1`, k_par_values span `[k_par_min, k_par_max]` in `nk_par` equal steps.
- **Validation:** `validate_semantic` rejects `nk_par < 1` with `error stop`; rejects `nk_par > 1` with `k_par_max <= k_par_min`. Both rejections use the existing `error stop` + buf-string pattern at `defs.f90:991-1005`.
- **Dispatch:** `compute_wire_bdg_gap_sweep` (`main_topology.f90:1283-1291`) branches on `cfg_in%bdg%nk_par`: `<= 1` → existing `eval_wire_bdg_gap` (zero behavior change); `> 1` → new `eval_wire_bdg_gap_bloch` (Bloch stack + T2 strict seam).
- **Regression guard:** `tests/regression/test_bdg_u13_nk_par_equiv.sh` runs the canonical wire-bdg 2D fixture twice (nk_par absent vs nk_par=1) and diffs all emitted outputs byte-for-byte. Registered in ctest as `regression_bdg_u13_nk_par_equiv`.
- **Smoke:** `tests/regression/test_bdg_u13_nk_par_sweep.sh` runs the canonical wire-bdg fixture with nk_par=4 over `[-0.001, +0.001]` 1/A (collapsed to Gamma so FEST's tight window converges). Full physics-correct sweep coverage is T4's destination.

## Files & commits

- `e17fe11` `test(bdg): U13 T3 — nk_par=1 vs absent byte-equivalence regression guard` — `tests/regression/test_bdg_u13_nk_par_equiv.sh` + ctest entry
- `5896455` `fix(test): populate cfg%params in wire spectral tests via paramDatabase` — prerequisite UB fix (uncovered during Task 2 dispatch)
- `9336dac` `feat(bdg): U13 T3 — nk_par/k_par_min/k_par_max config knob + validation` — `src/core/defs.f90` (bdg_config fields + validation) + `src/io/input_parser.f90` (parse_bdg keys)
- `f4af568` `refactor(sparse_matrices): extract public dense_to_csr helper` — DRY refactor per user approval (Task 3 Step 2)
- `3db44cc` `feat(bdg): U13 T3 — nk_par dispatch (fixed-kz vs Bloch stack)` — `src/apps/main_topology.f90` (dispatch + eval_wire_bdg_gap_bloch) + `tests/regression/test_bdg_u13_nk_par_sweep.sh` + ctest entry
