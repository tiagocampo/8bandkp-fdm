# Issue 008: Topology + BdG parser

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 12

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add `[topology]` section with nested sub-sections and `[bdg]` section. This replaces the most complex part of the old parser — 220 lines of deeply nested peek/backspace logic.

### `[topology]` section
- Base fields: `mode` ("qhe"/"qshe"/"bdg"), `compute_chern`, `compute_hall`, `qwz_u`, `compute_z2`, `z2_method`, `bhz_M`, `bhz_d`, `extract_edge_states`, `edge_E_window`, `compute_ldos`, `ldos_eta`, `ldos_E_range`, `ldos_num_E`
- Nested optional sub-sections (presence = enabled):
  - `[topology.gap_sweep]`: `B_range` (3-element array), `mu_range` (3-element array), `sweep_model`
  - `[topology.conductance]`: `method`, `berry_nk`, `landauer_energy`
  - `[topology.spectral]`: `k_range` (3-element), `E_range` (4-element: min, max, npts, eta)

### `[bdg]` section
- `mu`, `delta_0`, `gauge`, `B_sweep` (3-element), `kz`

Rewrite `test_topology_parser.pf` and `test_bdg_config.pf` with TOML test cases.

## Acceptance criteria

- [ ] `[topology]` base fields parse correctly for all 3 modes
- [ ] `[topology.gap_sweep]`, `[topology.conductance]`, `[topology.spectral]` parse as optional sub-sections
- [ ] `[bdg]` section parses correctly
- [ ] `test_topology_parser.pf` rewritten with TOML test cases (all 5 existing tests ported)
- [ ] `test_bdg_config.pf` rewritten with TOML test case
- [ ] All topology/BdG regression configs converted to TOML
- [ ] Topology simulation outputs match pre-migration golden data

## Blocked by

- Issue 002 (bulk tracer bullet)
