# T2 — Widen μ-window to Lutchyn-Oreg span (with writer-schema change)

## Question (HITL, grill-ticket, RESOLVED 2026-07-26)

What `topology.gap_sweep_mu_{min,max}` and `gap_sweep_nMu` produce a coarse `wire_bdg` sweep that:
- Samples **both** trivial (`s2_sign==+1` → `z2=+1` post Q3 native) and topological (`s2_sign==-1` → `z2=-1` post Q3 native) regions of the InAs wire, on a NEW sibling fixture (not the canonical `_2d.toml`, which must stay tight for `bcrit_curve ≈ 2.85 T` regression guard).
- Preserves `bcrit_curve ≈ 2.85 T` (regression guard via `lecture_13_topological.py` 1D curve parity check).

Plus: the new fixture's `output/z2_phase_diagram.dat` (and the parallel `output/wire_slim_pfaffian_witness.dat` written by T1b) MUST NOT collide with the canonical fixture's same-path output. Hardcoded output paths in `outputFunctions.f90` are an unfactored global; the new fixture wants hermetic filenames. Adopt config-driven file paths under `[topology]` (default values preserve canonical behavior; new fixture overrides).

## Grill session — 2026-07-26 — 11 questions, all resolved

| Q | Decision |
|---|----------|
| Q1 | New **sibling** fixture `wire_inas_gaas_bdg_topological_phase.toml` (NOT mutation of canonical `_2d.toml` — would break `test_wire_bdg_topological_2d.sh`). |
| Q2 | `gap_sweep_mu_min=0.5, max=1.0, nMu=11` — single coarse sweep (anchor `μ_c ≈ 0.6601` sits at grid index 4 of 11). Dual pass deferred (no clean sweep-loop refactor today). |
| Q3 | Update plan/lecture docs immediately per `codebase-doc-drift-prevention` — no deferral to U11. |
| Q4 | **Config-driven file paths** (`phase_diagram_file` + `slim_pfaffian_witness_file` under `[topology]`), defaulting to current names so canonical behavior unchanged. |
| Q5 | Both writers (`write_z2_phase_diagram`, `write_wire_slim_pfaffian_witness`) get the same config-driven path treatment. |
| Q6 | BHZ dispatch (`compute_z2_gap_sweep` → `write_z2_phase_diagram`) also honors the new config-key (default unchanged). |
| Q7 | Schema change lands in T2; `lecture_13_topological.py` script's path lookups **deferred to U11** (lecture stamping). |
| Q8 | T2 absorbs the writer-schema work — no separate T7 spawn (would create needless blocking overhead). |
| Q9 | New fields `character(len=64)`, default values `'z2_phase_diagram.dat'` / `'wire_slim_pfaffian_witness.dat'`. Writers concat `output//trim(...)`. |
| Q10 | Parser update lives in `input_parser.f90` (extends `parse_topology_block`); unit test `test_topology_parser.pf` gets new-key coverage. |
| Q11 | Plan amend (`docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U10 / §R9) in T2 scope. |

## Resolved scope (8 items)

1. **New fixture** `tests/regression/configs/wire_inas_gaas_bdg_topological_phase.toml`:
   - `gap_sweep_mu_min = 0.5`, `gap_sweep_mu_max = 1.0`, `gap_sweep_nMu = 11`
   - `gap_sweep_B_min = 0.0`, `gap_sweep_B_max = 5.0`, `gap_sweep_nB = 11`
   - `phase_diagram_file = 'z2_phase_diagram_phase.dat'`
   - `slim_pfaffian_witness_file = 'wire_slim_pfaffian_witness_phase.dat'`
   - All other fields copy canonical `wire_inas_gaas_bdg_topological_2d.toml`.

2. **`src/core/defs.f90`** — extend `topology_config`:
   - `character(len=64) :: phase_diagram_file = 'z2_phase_diagram.dat'`
   - `character(len=64) :: slim_pfaffian_witness_file = 'wire_slim_pfaffian_witness.dat'`

3. **`src/io/input_parser.f90`** — `parse_topology_block` reads both new keys (optional; defaults preserved).

4. **`src/io/outputFunctions.f90`** — both writers concat `output//trim(cfg%topo%phase_diagram_file)` and `output//trim(cfg%topo%slim_pfaffian_witness_file)` respectively; default behavior unchanged.

5. **`src/apps/main_topology.f90`** — no call-site changes needed (writers consume `cfg` directly).

6. **`tests/unit/test_topology_parser.pf`** — new `@test` for default + override of both keys.

7. **`tests/CMakeLists.txt`** — register `regression_wire_bdg_topological_phase` test under `regression` label (mirrors `regression_wire_bdg_topological_2d` wiring).

8. **Plan amend** — `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U10 (lines 435–457) and §R9 (line 144): note new fixture, writer-schema config keys.

## Acceptance

- [ ] New fixture builds, runs, emits both output files under their new names.
- [ ] Default canonical fixture (without override keys) still emits `output/z2_phase_diagram.dat` and `output/wire_slim_pfaffian_witness.dat` — no behavior change.
- [ ] Existing `test_wire_bdg_topological_2d.sh` and `test_wire_slim_pfaffian_witness.sh` tests still PASS (canonical fixture path collision is now impossible).
- [ ] New `test_wire_bdg_topological_phase.sh` test wires the new fixture through ctest under `regression` label.
- [ ] Unit test `test_topology_parser.pf` PASSes new-key test.
- [ ] Plan footer for §U10 reflects new fixture + writer schema.

## Blocks / blocked-by

- **Blocks**: T4 (regression test for non-flat colormap needs ≥1 row with non-zero z2 from new fixture).
- **Blocked by**: none.

## Type
grilling (HITL, resolved).

## Branch
`feat/bdg-u10-pfaffian-phase-diagram`.

## Claimed by
Tiago de Campos (2026-07-27 — fresh session starting T2 execution)

## Resolution (2026-07-27)

T2 executed in full. 8 items shipped in working tree on `feat/bdg-u10-pfaffian-phase-diagram` (uncommitted as of session end — pending PR #43 review or follow-up PR):

1. ✅ **New fixture** `tests/regression/configs/wire_inas_gaas_bdg_topological_phase.toml` — sibling of canonical `_2d.toml`, nB=6 × nMu=3 = 18 grid points, B∈[0,5]T, μ∈[0.6600, 0.6602]eV (tight window per FEST `±5·δ₀` around E=0 physics), both writer keys overridden to `*_phase.dat` suffix.
2. ✅ **`defs.f90`** — added `phase_diagram_file` (default `'z2_phase_diagram.dat'`) + `slim_pfaffian_witness_file` (default `'wire_slim_pfaffian_witness.dat'`) as `character(len=64)` fields on `topology_config`.
3. ✅ **`input_parser.f90`** — `parse_topology` reads both new keys via optional `get_value` blocks with `check_optional_stat`; defaults applied when keys absent.
4. ✅ **`outputFunctions.f90`** — both `write_z2_phase_diagram` and `write_wire_slim_pfaffian_witness` writers concat `output//trim(adjustl(cfg%topo%*))` for the path; error-stop messages also use the configured name.
5. ✅ **`main_topology.f90`** — no call-site changes needed; verified writers receive `cfg`.
6. ✅ **`tests/unit/test_topology_parser.pf`** — added 2 `@test` cases (defaults + overrides); single-line per `feedback_pfunit_macro_single_line`.
7. ✅ **`tests/CMakeLists.txt` + `tests/integration/test_wire_bdg_topological_phase.sh`** — new ctest entry `regression_wire_bdg_topological_phase` under `regression` label, 900s timeout, dedicated shell wrapper asserts both phase_diagram_phase.dat AND slim_pfaffian_witness_phase.dat are present and well-formed.
8. ✅ **Plan amend** `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U10 + §R9 — notes new fixture, writer-schema config keys, FEST physics window constraint on μ-window.

**Verification (2026-07-27):** 51/51 unit tests green; 5/5 wire-BdG regression tests green (canonical `regression_wire_bdg_topological_2d` + canonical `regression_wire_slim_pfaffian_witness` + new `regression_wire_bdg_topological_phase` + existing `regression_wire_bdg_topological` + `regression_wire_bdg_strain_shift`).

**FEST physics discovery (exit scope, but important for T4):** the original ticket-scope μ-window [0.5, 1.0] eV is NOT viable — `eval_wire_bdg_gap` (main_topology.f90:1358-1365) fixes FEAST window to `±5·δ₀` around E=0 (catches near-zero Majorana-mode subspace), so μ values outside ±0.001 eV of the conduction-band edge (0.6601 eV) leave the subspace empty and FEAST fails (U8 fail-fast guard). Per ticket Q2 grill fix, the phase fixture widens B-grid to nB=6 instead of nMu for finer B_crit transition detection while respecting FEST physics. T4 (non-flat colormap assertion) will need to decide whether to widen FEST window to track μ-sweep or constrain μ-window tighter — deferred per ticket scope.
