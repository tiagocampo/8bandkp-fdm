# Issue 5: pFUnit tests for consolidated validate() and validate_semantic() checks

**Type**: AFK
**Blocked by**: Issue 2, Issue 3
**User stories**: 22
**GitHub issue**: #24

## What to build

Write pFUnit passing-case tests for the new checks added in issues 2 and 3. Follow the existing pattern: construct a minimal valid config, call the validator, assert it does not crash. pFUnit 4.x cannot catch `error stop`, so only passing cases are tested.

**Tests for `validate()` checks (add to `test_defs.pf`):**

One test per new check. Each constructs a valid config that exercises the acceptance path:

| Test | Config setup | Validates |
|------|-------------|-----------|
| `test_validate_bulk_evnum_ok` | bulk, evnum=8 | V1: bulk accepts evnum=8 |
| `test_validate_qw_num_cb_ok` | QW with grid, num_cb within max | V2: QW accepts valid num_cb |
| `test_validate_qw_num_vb_ok` | QW with grid, num_vb within max | V3: QW accepts valid num_vb |
| `test_validate_wave_vector_modes` | Parameterized or per-mode tests for all 7 modes | V4: all 7 wave vector modes accepted |
| `test_validate_bsweep_positive` | B-sweep with positive step | V5: positive B-sweep step accepted |
| `test_validate_zmin_lt_zmax` | Material with z_min < z_max | V6: valid material layers accepted |
| `test_validate_sc_qw_ok` | QW with SC enabled | V7: SC with QW accepted |
| `test_validate_efield_z_nonzero` | QW with electric field, z(1) /= 0 | V8: valid electric field accepted |

**Tests for `validate_semantic()` checks (add to `test_parameters.pf`):**

One test per new check. Each constructs a valid config and calls `validate_semantic`:

| Test | App name | Config setup | Validates |
|------|----------|-------------|-----------|
| `test_validate_semantic_gfactor_bandidx_wire_ok` | gfactor | wire, valid bandIdx | S1 |
| `test_validate_semantic_topology_qshe_qw_ok` | topologicalAnalysis | QW, QSHE Z2 | S2 |
| `test_validate_semantic_topology_bdg_ok` | topologicalAnalysis | wire, BdG enabled | S3, S4 |
| `test_validate_semantic_topology_spectral_ok` | topologicalAnalysis | wire, spectral params valid | S5, S6, S7 |
| `test_validate_semantic_topology_sweep_ok` | topologicalAnalysis | QW, valid sweep_model | S8 |
| `test_validate_semantic_topology_conductance_ok` | topologicalAnalysis | valid conductance_method | S9 |
| `test_validate_semantic_topology_mode_ok` | topologicalAnalysis | valid topology mode | S10 |

Follow the existing pattern from `test_parameters.pf:228-269` and `test_defs.pf:174-210`. Each test constructs a `type(simulation_config)`, sets the minimum required fields, calls the validator, and asserts `.true.` on success.

## Acceptance criteria

- [ ] ~8 new passing-case tests in `test_defs.pf` for `validate()` checks
- [ ] ~7 new passing-case tests in `test_parameters.pf` for `validate_semantic()` checks
- [ ] Each test follows the existing pattern (construct valid config, call validator, assert .true.)
- [ ] All 34+ unit tests pass (original 34 + new tests)

## Blocked by

- Issue 2 (structural validation in validate())
- Issue 3 (app-semantic validation in validate_semantic())
