Status: ready-for-agent

## Parent

PRD: `.scratch/simulation-setup-orchestration/PRD.md`

## What to build

Final integration validation and cleanup. This slice has two deliverables: a pFUnit integration test proving the setup module's full pipeline works, and an ADR documenting the architectural decision.

**pFUnit integration test:** Add a test to `test_simulation_setup.pf` that exercises the full pipeline: `read_config` → `simulation_setup_init` → `setup_build_H` → `setup_solve_kpoint` for a simple bulk GaAs config. Assert that the eigenvalues at k=0 match the known 8-band GaAs bandgap (Eg ≈ 1.519 eV at Γ). This test proves the module's external API is correct end-to-end. The existing ~91 regression tests already prove the refactored apps produce identical output — no new integration script is needed.

**Remove backward-compat wrapper:** Remove the `read_and_setup` wrapper from `input_parser.f90`. Rename all remaining callers to use `read_config` directly. This includes `main_topology.f90` and the unit test files (`test_bdg_config.pf`, `test_topology_parser.pf`). Parser tests that need `profile`/`kpterms` call `read_config` then do their own minimal confinement init inline — they do not use the setup module.

**ADR:** Create `docs/adr/0001-simulation-setup-fat-type.md` documenting the architectural decision: fat derived type with allocatable components vs. polymorphic subtypes/strategy pattern. The PRD section "Type design: fat type with allocatable components" is the spec for the ADR content. The ADR covers: the trade-off (simplicity vs. type safety), the dispatch-seam for future Candidate 4 (unified Hamiltonian block structure), and the `read_config` split rationale.

## Acceptance criteria

- [ ] pFUnit integration test: full pipeline (read_config → setup_init → build_H → solve_kpoint) produces correct GaAs bulk eigenvalues at k=0
- [ ] `read_and_setup` wrapper removed from `input_parser.f90`
- [ ] `main_topology.f90` calls `read_config` directly
- [ ] `test_bdg_config.pf` and `test_topology_parser.pf` call `read_config` directly with inline confinement init where needed
- [ ] No remaining references to `read_and_setup` anywhere in the codebase
- [ ] `docs/adr/0001-simulation-setup-fat-type.md` exists with the fat-type trade-off decision
- [ ] All existing tests pass

## Blocked by

- Issue 02 (velocity matrices + optics migration)
- Issue 03 (wire path + bandStructure migration)
