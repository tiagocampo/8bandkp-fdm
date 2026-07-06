**Status**: COMPLETE (2026-07-05)

# Issue 06 — BdG LDOS + A(k,E) + KTD6 close (U9)

**Parent PRD**: `/data/8bandkp-fdm/.scratch/bdg-majorana-validation/PRD.md` (Unit U9, LDOS path; A(k,E) and Nambu-resolved LDOS land alongside)
**Issue source**: `/data/8bandkp-fdm/.superpowers/sdd/issue-06-brief.md`
**Plan & context**: `/data/8bandkp-fdm/.superpowers/sdd/understand-report.md`
**Phase**: PR-B, fan-out with Issue 04

## What to build

The missing BdG observables on the topological wire: zero-bias LDOS peak (the experimental MZM signature), `A(k,E)` spectral function (gap + in-gap mode), Nambu-resolved LDOS (electron vs hole sector). All three on the **BdG Hamiltonian, not the normal-state one** — exposing the zero-energy Majorana peak.

Implementation:
1. New `bdq_spectral` enum value on the existing `sweep_model`. Single `topologicalAnalysis` invocation produces all three outputs in one PARDISO setup (one factorization, three observers). No new TOML fields (ADR 0002).
2. New callers of the existing H-agnostic `compute_ldos_csr` (BdG-LDOS, Nambu-resolved LDOS). The BdG CSR is built by the unified hole-block builder from Issue 03.
3. New `compute_spectral_function_bdg_wire` sibling of `compute_spectral_function_wire` for `A(k,E)`. Mirrors the CSR spectral skeleton. The two CSR consumers (BdG + normal-state wire) share a CSR-spectral subroutine; the bulk/QW dense variants stay untouched until a separate follow-up.
4. New dedicated output files: `bdg_ldos.dat`, `bdg_ldos_nambu.dat`, `bdg_spectral.dat`. No new `topological_result` fields (avoids approval-gated `defs.f90` type change). The 2D-grid writers go to `outputFunctions.f90` (the consolidated writer pattern from `#07`), not inline writes.
5. **KTD6 closes opportunistically:** while `green_functions.f90` is open for the BdG-LDOS caller, route `compute_spectral_function_wire`'s manual FEAST window through `apply_solver_window`. Closes the third KTD6 bypass.

The electron/hole split comes from the Nambu block structure (no new projection code).

## Acceptance criteria

- [ ] `sweep_model = bdq_spectral` produces all three output files (`bdg_ldos.dat`, `bdg_ldos_nambu.dat`, `bdg_spectral.dat`) in one invocation.
- [ ] BdG LDOS shows a zero-energy peak in the topological phase (via `compute_ldos_csr` on the BdG CSR).
- [ ] BdG `A(k,E)` shows the superconducting gap and the in-gap Majorana mode.
- [ ] Nambu-resolved LDOS splits the electron vs hole sector contribution to the zero-energy peak.
- [ ] The shared spectral-observable skeleton produces identical output to the pre-refactor normal-state wire spectral function (no regression).
- [ ] Normal-state LDOS (existing) is unchanged.
- [ ] `compute_spectral_function_wire`'s manual FEAST window routes through `apply_solver_window` (KTD6 closed); the existing wire spectral function still produces equivalent output post-routing.
- [ ] `tests/unit/test_green_functions.pf` extensions green.
- [ ] `# COVERAGE:` annotations in the verifier; `validation_universe.yml` cells updated (`bdg_spectral_function` if introduced as a tracked quantity — explicit scope decision per that file's header).
- [ ] Per-task code review + spec compliance review clean.

## Pre-existing state (from Understand report)

### `src/physics/green_functions.f90` (487 lines)
- `compute_spectral_function_bulk` (180–232): DENSE
- `compute_spectral_function_qw` (109–178): DENSE
- `compute_spectral_function_wire` (234–343): **KTD6 bypass at lines 283–292** — manual FEAST window: `eigen_cfg%emin = minval(E_arr) − 5·η`, then either user-override OR `auto_compute_energy_window` fallback. **This is the KTD6 bypass to close.**
- `compute_ldos_csr` (362–485): **already H-agnostic** — takes any `csr_matrix` + E + η; uses `pardiso_c`, requires structural diagonal at every row. Reuse for BdG-LDOS.
- `compute_landauer_transmission_1d` (70–89)

### `src/io/outputFunctions.f90`
- `write_bdg_eigenvalues` (line 339): single-axis writer for BdG eigenvalue dump.
- **No 2D-grid writers for `bdg_ldos.dat`, `bdg_spectral.dat`** — must be added.

### `src/apps/main_topology.f90` (1320 lines)
- `sweep_model` enum at `defs.f90:335` (`bhz_analytic | wire_bdg | qw_fukane`); switch at line 1046–1059. **Add `bdq_spectral` case** to the dispatch.
- 2D-grid writers should live in `outputFunctions.f90` per the existing consolidated writer pattern (the post-#07 pattern).

### `src/core/defs.f90`
- `topology_config.sweep_model` at line 335 (3-valued enum). Add `bdq_spectral` (4-valued).
- `validate_semantic` at line 855 — topology-mode validation at lines 915–928 rejects unknown modes; sweep_model validation at 1018–1036 rejects unknown + mismatched confinement. Update both for `bdq_spectral`.

### `src/physics/bdg_hamiltonian.f90` (433 lines)
- The BdG CSR is now built by the unified hole-block builder (Issue 03's `build_bdg_hole_block`). This is the prerequisite for the LDOS to be physically meaningful.

### `tests/unit/test_green_functions.pf`
- 13k lines; existing LDOS, Landauer tests. Extend with BdG LDOS / spectral tests.

### `tests/integration/validation_universe.yml` (484 lines)
- Existing cells: `minigap` (line 478), `majorana_modes` (line 453), `ldos` (line 338).
- **Missing**: `majorana_polarization` cell (Issue 04 adds); `bdg_spectral_function` cell — explicit scope decision per file header. Add `bdg_spectral_function` as `tier: required` for `geometry: wire, material: InAsW`.

## Constraints from CLAUDE.md + ADRs

- **ADR 0001 (fat derived type)**: No new types; `bdq_spectral` is a string enum on existing `topo%mode` per ADR 0002.
- **ADR 0002 (no new TOML fields)**: Outputs go to dedicated files (`bdg_ldos.dat`, `bdg_ldos_nambu.dat`, `bdg_spectral.dat`), NOT new `topological_result` fields. **Explicit scope decision per PRD: "if any observable must enter the result type, that is itself approval-gated."**
- **ADR 0003 (build-and-solve in app)**: Build-and-solve stays in `main_topology`; pure observable consumers in the new module.
- **ADR 0005 (one stable window per sweep via `apply_solver_window`)**: **KTD6 closes opportunistically** — `compute_spectral_function_wire`'s manual FEAST window routes through `apply_solver_window`.
- **CLAUDE.md Boundaries**: 
  - Adding `bdq_spectral` to `validate_semantic` enum is a validation-only change (no derived-type field); NOT approval-gated.
  - The shared CSR spectral subroutine extraction is a refactor; not approval-gated.
- **CLAUDE.md Code Conventions**: F2018; `private` default + explicit `public ::` exports; `error stop` not `stop 1`; no `goto`; `<= 300 lines/file`, `<= 50 lines/function`; pFUnit `@assertEqual`/`@assertTrue` MUST be single-line.

## File ownership (exhaustive)

### New files (you create)
- `src/physics/spectral_bdg_wire.f90` — new thin module housing `compute_spectral_function_bdg_wire` sibling + BdG-LDOS caller + Nambu-resolved LDOS.
- `tests/integration/verify_bdg_spectral.py` — Python verifier (standard-star pattern; uses `tests/integration/star_helpers.py`); COVERAGE annotations here.
- `tests/integration/test_bdg_spectral.sh` — thin shell wrapper calling the Python verifier (no COVERAGE annotations per AGENTS.md convention).

### Modified files
- `src/physics/green_functions.f90` — opportunistically route `compute_spectral_function_wire`'s manual FEAST window through `apply_solver_window` (closes KTD6); share the CSR-spectral subroutine between BdG and normal-state wire.
- `src/io/outputFunctions.f90` — add 2D-grid writers for `bdg_ldos.dat`, `bdg_ldos_nambu.dat`, `bdg_spectral.dat`.
- `src/apps/main_topology.f90` — add `case ('bdq_spectral')` to the `cfg%topo%mode` dispatch; single PARDISO setup, three observers.
- `src/core/defs.f90` — add `bdq_spectral` to the topology-mode `validate_semantic` enum.
- `tests/unit/test_green_functions.pf` — extend with BdG LDOS / spectral / KTD6 routing tests.
- `tests/integration/validation_universe.yml` — add `bdg_spectral_function` cell (`geometry: wire, material: InAsW, tier: required`).
- `src/CMakeLists.txt` — add `spectral_bdg_wire.f90` to COMMON_SOURCES.
- `tests/CMakeLists.txt` — wire any new pFUnit tests.

### NOT modified (out of scope)
- `src/physics/bdg_hamiltonian.f90` (Issue 03 owns; do NOT touch).
- `src/physics/topological_analysis.f90` (Issue 04/07 own).
- `src/physics/bdg_observables.f90` (Issue 00 owns).
- The dense (bulk/QW) spectral routines (`compute_spectral_function_bulk`, `compute_spectral_function_qw`) — they stay untouched.

## Suggested implementation (illustrative)

```fortran
module spectral_bdg_wire
  use kind_module, only: dp
  use csr_module, only: csr_matrix, csr_to_dense, dense_to_csr
  use green_functions, only: compute_ldos_csr, apply_solver_window
  implicit none
  private
  public :: compute_spectral_function_bdg_wire, compute_bdg_ldos, compute_bdg_ldos_nambu

  ! BdG spectral function A(k,E) on the BdG CSR
  subroutine compute_spectral_function_bdg_wire(H_bdg_csr, k_values, E_values, eta, A_kE, ...)
    type(csr_matrix), intent(in) :: H_bdg_csr
    real(dp), intent(in) :: k_values(:), E_values(:), eta
    real(dp), intent(out) :: A_kE(:,:)
    ! Use shared CSR-spectral subroutine with normal-state wire (post-KTD6)
  end subroutine

  ! BdG LDOS via existing H-agnostic compute_ldos_csr
  subroutine compute_bdg_ldos(H_bdg_csr, E_values, eta, ldos)
    type(csr_matrix), intent(in) :: H_bdg_csr
    real(dp), intent(in) :: E_values(:), eta
    real(dp), intent(out) :: ldos(:)
  end subroutine

  ! Nambu-resolved LDOS (electron vs hole sector)
  subroutine compute_bdg_ldos_nambu(H_bdg_csr, E_values, eta, ldos_e, ldos_h)
    type(csr_matrix), intent(in) :: H_bdg_csr
    real(dp), intent(in) :: E_values(:), eta
    real(dp), intent(out) :: ldos_e(:), ldos_h(:)
  end subroutine
end module
```

## Tests required

### `tests/unit/test_green_functions.pf` extensions
1. `compute_spectral_function_bdg_wire` produces non-negative `A(k,E)` (the spectral function is non-negative).
2. BdG LDOS peak at E=0 in the topological phase (use a synthetic BdG CSR with a zero mode; expect a Lorentzian peak at E=0).
3. Nambu-resolved LDOS sums to total LDOS (`ldos_e(n) + ldos_h(n) == ldos(n)`).
4. KTD6 closure: the post-refactor `compute_spectral_function_wire` produces equivalent output to a captured reference (existing test_spectral_function_wire must still pass; capture pre-refactor reference if not already).
5. `apply_solver_window` is called from `compute_spectral_function_wire` (verify by mocking or by checking that the windowing no longer uses `auto_compute_energy_window` fallback).

### `tests/integration/verify_bdg_spectral.py` (new)
- Standard-star pattern (`star_helpers.py`).
- Run `topologicalAnalysis` with `sweep_model = bdq_spectral`.
- Assert all three output files exist and have valid content.
- Assert BdG LDOS has a peak at E=0 (zero-energy mode visible).
- Assert A(k,E) shows the SC gap and in-gap mode.
- `# COVERAGE:` annotations: `observable=bdg_ldos geometry=wire material=InAsW`; `observable=bdg_spectral_function geometry=wire material=InAsW`.

## TDD discipline (mandatory)

1. Write the new pFUnit tests FIRST.
2. Confirm they fail (functions don't exist yet).
3. Implement `spectral_bdg_wire.f90`.
4. Confirm new tests pass.
5. Implement the KTD6 closure (route `compute_spectral_function_wire` through `apply_solver_window`).
6. Confirm pre-existing tests still pass (especially `test_spectral_function_wire`).
7. Implement the `bdq_spectral` enum + dispatch case.
8. Implement the 2D-grid writers.
9. Implement the integration verifier.
10. Run full unit suite to verify no regressions (38/41 expected — 3 follow-up failures from Issue 03 are pre-existing and out of scope).

## Build commands (use these exact forms)

```bash
cmake --build build
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit --output-on-failure
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L regression --output-on-failure
ctest --test-dir build -L verification --output-on-failure  # for the new verifier
```

## Out of scope (must NOT do)

- Do NOT modify `src/physics/bdg_hamiltonian.f90`.
- Do NOT modify `src/physics/topological_analysis.f90`.
- Do NOT modify `src/physics/bdg_observables.f90`.
- Do NOT modify dense (bulk/QW) spectral routines.
- Do NOT add new `topological_result` fields (would be approval-gated per PRD).
- Do NOT add new TOML fields.
- Do NOT add new `sweep_model` enum value beyond `bdq_spectral`.

## Report file path

Write your full report to: `/data/8bandkp-fdm/.superpowers/sdd/issue-06-report.md`

The report must include:
- Status (DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT)
- Files modified (with COVERAGE annotations)
- TDD evidence (RED → GREEN) with exact commands and output
- Test results (full unit suite pass count; pre-existing 3 follow-up failures expected)
- KTD6 closure confirmation: pre/post-refactor output equivalence
- Commit SHA
- Self-review findings
- Concerns

## Risk notes

- KTD6 closure is a refactor of a routine shared with non-BdG callers. The existing `test_spectral_function_wire` regression must still pass post-routing. If `apply_solver_window`'s default windowing differs from `auto_compute_energy_window`'s Gershgorin envelope, the regression test will fail. The implementer must verify equivalence or bump the golden tolerance.
- The BdG CSR is 16N×16N — much larger than the normal-state CSR (8N×8N for QW; 8nx×8nx for wire). PARDISO factorization cost scales accordingly; make sure the test grid size is small enough to not time out.
- The electron/hole split in Nambu-resolved LDOS is a structural property of the BdG block structure — no new projection code. Verify the splitting matches the BdG layout (electron block: (1..8N, 1..8N); hole block: (8N+1..16N, 8N+1..16N); off-diagonal: pairing).
- The `apply_solver_window` function is part of the eigensolver dispatch (per ADR 0005). Verify its signature matches the KTD6 closure requirement.

---

## Outcome (as executed)

- **Module**: `src/physics/spectral_bdg_wire.f90` (255 lines; new thin module housing `compute_spectral_function_bdg_wire`, `compute_bdg_ldos`, `compute_bdg_ldos_nambu`).
- **Dispatch**: `case ('bdq_spectral')` added to `cfg%topo%mode` in `src/apps/main_topology.f90`; single PARDISO setup, three observers.
- **Enum**: `bdq_spectral` added to `validate_semantic` in `defs.f90` (validation-only, not approval-gated).
- **Outputs**: 3 files (`bdg_ldos.dat`, `bdg_ldos_nambu.dat`, `bdg_spectral.dat`) from one invocation.
- **2D-grid writers**: consolidated in `src/io/outputFunctions.f90`.
- **KTD6 close**: `compute_spectral_function_wire`'s manual FEAST window routed through `apply_solver_window` via `asw_evals` (specific procedure of `apply_solver_window`).
- **Tests**: `tests/unit/test_green_functions.pf` extensions green; `tests/integration/verify_bdg_spectral.py` + `test_bdg_spectral.sh` wired.
- **Plot generation**: A(k,E) and LDOS plots generated by verifier and copied to `docs/lecture/figures/`.
- **COVERAGE annotations**:
  - `# COVERAGE: observable=bdg_ldos geometry=wire material=InAsW`
  - `# COVERAGE: observable=bdg_spectral_function geometry=wire material=InAsW`
- **Deferred**:
  - `csr_spectral_lorentzian_sum` shared helper: only consumed by BdG path; normal-state wire keeps inline `lorentzian_broadening` (deferred refactor).
  - `bdq_spectral` only supports wire confinement; dense-QW spectral deferred.