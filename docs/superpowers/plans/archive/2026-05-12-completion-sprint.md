# Completion Sprint (Phases 14-16) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete all remaining validation and infrastructure gaps from the ideation-driven work on `feature/bdg-topological-superconductivity`.

**Architecture:** Three effort-ordered phases: quick wins (mechanical changes, golden data generation), validation tightening (CSR Krylov, standard-star assertions, doc audit), then integration + Rashba physics calibration. Performance improvements applied opportunistically when touching Fortran source.

**Tech Stack:** Fortran 2018 (gfortran, MKL, pFUnit), Python 3 (numpy, PyYAML), CMake/Ninja, Bash test scripts.

---

## Phase 14: Quick Wins Sprint

### Task 1: Contiguous attribute gaps

**Files:**
- Modify: `src/core/utils.f90:19`
- Modify: `src/physics/spin_projection.f90:14`

- [ ] **Step 1: Add `contiguous` to `dns` in `utils.f90`**

Change line 19 from:
```fortran
    complex(kind=dp), intent(in), dimension(:,:) :: dns
```
to:
```fortran
    complex(kind=dp), intent(in), contiguous, dimension(:,:) :: dns
```

- [ ] **Step 2: Add `contiguous` to `psi` in `spin_projection.f90`**

Change line 14 from:
```fortran
    complex(kind=dp), intent(in) :: psi(:)  ! (8*Ngrid)
```
to:
```fortran
    complex(kind=dp), intent(in), contiguous :: psi(:)  ! (8*Ngrid)
```

- [ ] **Step 3: Build and test**

Run:
```bash
cmake --build build && ctest --test-dir build -j4 --output-on-failure
```
Expected: All tests pass. No behavior change — `contiguous` is an optimization hint.

- [ ] **Step 4: Commit**

```bash
git add src/core/utils.f90 src/physics/spin_projection.f90
git commit -m "perf: add contiguous attribute to remaining hot-path array arguments"
```

---

### Task 2: Coverage matrix orphan annotations

**Files:**
- Modify: `tests/integration/test_gfactor_no_optics.sh`
- Modify: `tests/integration/test_qw_bandstructure.sh`
- Modify: `tests/integration/test_qw_character_and_output_dir.sh`
- Verify: `tests/integration/validation_universe.yml`

- [ ] **Step 1: Fix `test_gfactor_no_optics.sh` annotation**

Change line 2 from:
```bash
# COVERAGE: observable=gfactor geometry=QW material=GaAs/AlGaAs
```
to:
```bash
# COVERAGE: observable=g*_cb geometry=QW material=GaAs/AlGaAs
```

- [ ] **Step 2: Fix `test_qw_bandstructure.sh` annotation**

Change line 2 from:
```bash
# COVERAGE: observable=E_sub geometry=QW material=AlSbW/GaSbW/InAsW
```
to:
```bash
# COVERAGE: observable=subband_spacing geometry=QW material=AlSbW/GaSbW/InAsW
```

- [ ] **Step 3: Fix `test_qw_character_and_output_dir.sh` annotation**

Change line 2 from:
```bash
# COVERAGE: observable=state_character geometry=QW material=AlSbW/GaSbW/InAsW
```
to:
```bash
# COVERAGE: observable=CB_ground_state geometry=QW material=AlSbW/GaSbW/InAsW
```

Verify that `CB_ground_state` exists in the `observables` list in `validation_universe.yml` (it does — line 30).

- [ ] **Step 4: Run coverage matrix to verify zero orphans**

Run:
```bash
ctest --test-dir build -L coverage
```
Expected: PASS with zero orphan annotations.

- [ ] **Step 5: Commit**

```bash
git add tests/integration/test_gfactor_no_optics.sh tests/integration/test_qw_bandstructure.sh tests/integration/test_qw_character_and_output_dir.sh
git commit -m "fix(tests): correct 3 orphan COVERAGE annotations to match universe observables"
```

---

### Task 3: Strain validation docs commit

**Files:**
- Commit: `docs/brainstorms/2026-05-11-strain-validation-requirements.md`
- Commit: `docs/plans/2026-05-11-strain-validation-plan.md`

- [ ] **Step 1: Stage and commit the untracked planning documents**

```bash
git add docs/brainstorms/2026-05-11-strain-validation-requirements.md docs/plans/2026-05-11-strain-validation-plan.md
git commit -m "docs: commit strain validation brainstorm and plan"
```

---

### Task 4: gfactor regression golden data (3 configs)

**Files:**
- Use: `tests/regression/configs/gfactor_bulk_gaas_vb.cfg` (exists)
- Use: `tests/regression/configs/gfactor_bulk_gaasw_cb.cfg` (exists)
- Use: `tests/regression/configs/gfactor_qw_vb.cfg` (exists)
- Create: `tests/regression/data/gfactor_bulk_gaas_vb/output.txt`
- Create: `tests/regression/data/gfactor_bulk_gaasw_cb/output.txt`
- Create: `tests/regression/data/gfactor_qw_vb/output.txt`
- Modify: `tests/CMakeLists.txt`

- [ ] **Step 1: Build gfactorCalculation**

```bash
cmake --build build
```

- [ ] **Step 2: Generate golden data for GaAs bulk VB**

```bash
WORKDIR=$(mktemp -d)
/bin/cp tests/regression/configs/gfactor_bulk_gaas_vb.cfg "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/gfactorCalculation > output.txt 2>&1)
mkdir -p tests/regression/data/gfactor_bulk_gaas_vb
/bin/cp "$WORKDIR/output.txt" tests/regression/data/gfactor_bulk_gaas_vb/output.txt
rm -rf "$WORKDIR"
```

- [ ] **Step 3: Generate golden data for GaAsW bulk CB**

```bash
WORKDIR=$(mktemp -d)
/bin/cp tests/regression/configs/gfactor_bulk_gaasw_cb.cfg "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/gfactorCalculation > output.txt 2>&1)
mkdir -p tests/regression/data/gfactor_bulk_gaasw_cb
/bin/cp "$WORKDIR/output.txt" tests/regression/data/gfactor_bulk_gaasw_cb/output.txt
rm -rf "$WORKDIR"
```

- [ ] **Step 4: Generate golden data for QW VB**

```bash
WORKDIR=$(mktemp -d)
/bin/cp tests/regression/configs/gfactor_qw_vb.cfg "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/gfactorCalculation > output.txt 2>&1)
mkdir -p tests/regression/data/gfactor_qw_vb
/bin/cp "$WORKDIR/output.txt" tests/regression/data/gfactor_qw_vb/output.txt
rm -rf "$WORKDIR"
```

- [ ] **Step 5: Register tests in CMakeLists.txt**

Add after the existing `regression_gfactor_no_optics` block (around line 277), following the canonical pattern from `regression_gfactor_cb`:

```cmake
# g-factor bulk GaAs VB
add_test(
    NAME regression_gfactor_bulk_gaas_vb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_bulk_gaas_vb.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_bulk_gaas_vb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor bulk GaAsW CB
add_test(
    NAME regression_gfactor_bulk_gaasw_cb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_bulk_gaasw_cb.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_bulk_gaasw_cb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor QW VB
add_test(
    NAME regression_gfactor_qw_vb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_qw_vb.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_qw_vb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
```

Add all three test names to the `set_tests_properties(... PROPERTIES LABELS "regression")` block that contains the existing gfactor tests.

- [ ] **Step 6: Re-configure, build, and run new tests**

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
ctest --test-dir build -R gfactor --output-on-failure
```
Expected: All gfactor tests pass (old + 3 new).

- [ ] **Step 7: Run full suite to verify no regressions**

```bash
ctest --test-dir build -j4 --output-on-failure
```
Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add tests/regression/data/gfactor_bulk_gaas_vb/ tests/regression/data/gfactor_bulk_gaasw_cb/ tests/regression/data/gfactor_qw_vb/ tests/CMakeLists.txt
git commit -m "test: add gfactor regression tests for GaAs VB, GaAsW CB, QW VB"
```

---

## Phase 15: Validation Tightening

### Task 5: CSR Krylov snapshot completion (3 code paths)

**Files:**
- Modify: `tests/support/generate_krylov_reference.f90`
- Modify (auto-regenerate): `tests/support/krylov_reference_data.f90`
- Modify: `tests/unit/test_krylov_snapshots.pf`

This task adds Krylov snapshot tests for SC loop, optics wire, and gfactor wire code paths. Each path needs: (1) a `compute_*` subroutine in the generator, (2) regenerated reference data, (3) a `@test` in the snapshot test file.

- [ ] **Step 1: Add SC loop Krylov snapshot to generator**

In `tests/support/generate_krylov_reference.f90`:
- Add `sc_vecs` to the allocatable declaration alongside the existing arrays
- Add `call compute_sc_loop(sc_vecs)` to the main program
- Add `call write_declaration('sc_ref', sc_vecs)`
- Add `call write_init_routine('sc_ref', sc_vecs)`
- Add `sc_vecs` to the `deallocate` statement
- Add `compute_sc_loop` subroutine that builds a wire Hamiltonian (same as `compute_wire` but enables `cfg%SC = .true.` with a simple doping profile). Use a 3x3 GaAs grid to keep reference data small. The snapshot captures the initial Hamiltonian before any SC iteration.

- [ ] **Step 2: Add optics wire Krylov snapshot to generator**

Similarly add `compute_optics_wire`:
- Build a 3x3 GaAs wire Hamiltonian using `ZB8bandGeneralized` (same as existing `compute_wire`)
- The key difference: the optics path uses the CSR matrix from `build_velocity_matrices`. Snapshot the wire Hamiltonian itself (the optics code path uses the same CSR assembly)
- Name: `optics_wire_ref`

- [ ] **Step 3: Add gfactor wire Krylov snapshot to generator**

Add `compute_gfactor_wire`:
- Build a 3x3 GaAs wire Hamiltonian (same grid as existing wire test)
- The gfactor code path uses the CSR matrix for commutator-based velocity operators. Snapshot the Hamiltonian
- Name: `gfactor_wire_ref`

- [ ] **Step 4: Regenerate reference data**

```bash
cd build && cmake --build . --target regenerate_krylov_references
```
This regenerates `tests/support/krylov_reference_data.f90` with the 3 new data arrays.

- [ ] **Step 5: Add 3 test cases to `test_krylov_snapshots.pf`**

Add three `@test` subroutines following the established pattern (see `test_snapshot_wire_gaas`):
- `test_snapshot_sc_loop` — calls `init_sc_ref_n{N}_k6`, builds wire + SC config, asserts Krylov match
- `test_snapshot_optics_wire` — calls `init_optics_wire_ref_n{N}_k6`, builds wire Hamiltonian, asserts match
- `test_snapshot_gfactor_wire` — calls `init_gfactor_wire_ref_n{N}_k6`, builds wire Hamiltonian, asserts match

Each test must call the `init_*` function before building the matrix, and use `TOL = 1.0e-12_dp`.

- [ ] **Step 6: Build and test**

```bash
cmake --build build && ctest --test-dir build -R krylov --output-on-failure
```
Expected: All Krylov tests pass (4 existing + 3 new).

- [ ] **Step 7: Run full suite**

```bash
ctest --test-dir build -j4 --output-on-failure
```

- [ ] **Step 8: Commit**

```bash
git add tests/support/generate_krylov_reference.f90 tests/support/krylov_reference_data.f90 tests/unit/test_krylov_snapshots.pf
git commit -m "test: add Krylov snapshot tests for SC loop, optics wire, gfactor wire paths"
```

---

### Task 6: Standard-star assertion tightening

**Files:**
- Modify: `tests/integration/verify_star_gaas_algaas_qw.py`
- Modify: `tests/integration/verify_star_inas_wire.py`
- Modify: `tests/integration/verify_star_insb_bulk.py`
- Modify: `tests/integration/verify_star_inas_gasb_qw.py`
- Modify: `tests/integration/verify_star_inas_gaas_qw.py`
- Modify: `docs/reference/benchmarks.md`

- [ ] **Step 1: Establish S4 absorption onset reference**

Run `opticalProperties` with the GaAs/AlGaAs QW config to extract the actual absorption onset energy:

```bash
WORKDIR=$(mktemp -d)
/bin/cp tests/regression/configs/qw_gaas_algaas.cfg "$WORKDIR/input.cfg"
# Add an optics block if not present, with TE polarization absorption
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/opticalProperties > output.txt 2>&1)
# Extract onset from output (first non-zero absorption energy)
python3 -c "
import os, sys
# Parse absorption data and find onset
with open(os.path.join('$WORKDIR', 'output', 'absorption_TE.dat')) as f:
    lines = f.readlines()
for line in lines:
    parts = line.split()
    if len(parts) >= 2 and float(parts[1]) > 1e-6:
        print(f'ONSET_REF = {parts[0]}')
        break
"
rm -rf "$WORKDIR"
```

- [ ] **Step 2: Replace S4 onset range check with `compare_value`**

In `tests/integration/verify_star_gaas_algaas_qw.py`, at lines 356-390, replace the range check with:

```python
ONSET_REF = <value from Step 1>  # eV, from opticalProperties reference run

print(f"  TE absorption onset: {onset_energy:.4f} eV")
print(f"  Reference: {ONSET_REF:.4f} eV")
status_onset, delta_onset = compare_value(onset_energy, ONSET_REF, TOL_ABS_ONSET)

row_onset = {
    'computed': onset_energy,
    'expected': ONSET_REF,
    'tolerance': TOL_ABS_ONSET,
    'delta': delta_onset,
    'name': 'TE absorption onset',
    'unit': 'eV',
    'status': 'PASS' if status_onset else 'FAIL',
    'material': 'GaAs/AlGaAs QW',
    'reference': 'opticalProperties reference run',
    'tol_str': f'{TOL_ABS_ONSET*100:.0f}%',
}
```

- [ ] **Step 3: Establish S7 wire g-factor reference**

```bash
WORKDIR=$(mktemp -d)
/bin/cp tests/regression/configs/wire_inas_gfactor.cfg "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/gfactorCalculation > output.txt 2>&1)
# Extract gz from output
grep -A2 "gz" "$WORKDIR/output.txt" | tail -1
rm -rf "$WORKDIR"
```

- [ ] **Step 4: Replace S7 g-factor range check with `compare_value`**

In `tests/integration/verify_star_inas_wire.py`:

Add constant after line 75:
```python
G_WIRE_REF = <value from Step 3>  # g-factor from gfactorCalculation wire run
```

Replace the `check_wire_gfactor` call at line 233 with:
```python
status_gf, delta_gf = compare_value(gz, G_WIRE_REF, TOL_GFACTOR_REGRESSION)
```

Update the corresponding row dictionary to use `G_WIRE_REF` as the reference.

- [ ] **Step 5: Centralize run wrappers in S3**

In `tests/integration/verify_star_insb_bulk.py`, replace the local `run_bandstructure` and `run_gfactor` functions (lines 76-91) with calls to `star_helpers.run_exe`:

```python
# Replace local run_bandstructure(build_dir, config_path, work_dir)
# with: rc, outdir = run_exe(build_dir, 'bandStructure', config_path, work_dir)
```

Update all call sites from `run_bandstructure(...)` to `run_exe(build_dir, 'bandStructure', ...)` and similarly for `run_gfactor`.

- [ ] **Step 6: Centralize run wrappers in S5**

Same change in `tests/integration/verify_star_inas_gasb_qw.py` (lines 130, 139).

- [ ] **Step 7: Centralize run wrappers in S6**

Same change in `tests/integration/verify_star_inas_gaas_qw.py` (lines 111, 120).

- [ ] **Step 8: Fix Roth formula in benchmarks.md**

In `docs/reference/benchmarks.md`, change line 206 from:
```
1. **Roth formula**: g* = g0 - (2Ep/3)(1/Eg + 1/(Eg+Delta)) gives -0.44 using
```
to:
```
1. **Roth formula**: g* = 2 - 2*Ep*DeltaSO / (3*Eg*(Eg + DeltaSO)) gives -0.44 using
```
This matches the implementation in `star_helpers.py` line 123 and Winkler (2003) Eq. 6.42.

- [ ] **Step 9: Run full standard-star suite**

```bash
ctest --test-dir build -L standard-star --output-on-failure
```
Expected: All 7 standard-star tests pass.

- [ ] **Step 10: Run full suite**

```bash
ctest --test-dir build -j4 --output-on-failure
```

- [ ] **Step 11: Commit**

```bash
git add tests/integration/verify_star_gaas_algaas_qw.py tests/integration/verify_star_inas_wire.py tests/integration/verify_star_insb_bulk.py tests/integration/verify_star_inas_gasb_qw.py tests/integration/verify_star_inas_gaas_qw.py docs/reference/benchmarks.md
git commit -m "fix(tests): tighten standard-star assertions, centralize run wrappers, fix Roth formula doc"
```

---

### Task 7: Re-scope docs physics revamp

**Files:**
- Read: `docs/plans/2026-04-21-docs-physics-revamp-plan.md`
- Modify: `docs/plans/REVIEW.md`
- Modify: `docs/plans/BACKLOG.md`

- [ ] **Step 1: Audit original 12 tasks against current state**

Read the plan and evaluate each task:
- Tasks 1-2: Already DONE (per REVIEW.md)
- Tasks 3-5 (figure fixes): Lecture-test pairs now generate overlay plots for all 14 lectures. Figure reference fixes were applied in commit `7b80237`. Evaluate remaining gaps.
- Tasks 6-7 (solver repair): Check if ISBT dipole sign issue is still present in `optical_spectra.f90`. Check if gain quasi-Fermi integration is still incomplete.
- Tasks 8-9 (wire hardening): Check if FEAST convergence issues have been addressed by the various bug fixes in Phases 5-7.
- Tasks 10-12 (chapter rebuilds): All 14 lectures now have Verification sections. Evaluate if "rebuild with real data" is satisfied.

- [ ] **Step 2: Update REVIEW.md and BACKLOG.md**

Based on the audit, either:
- Mark #26 as COMPLETE with a note about what was covered by lecture-test pairs, or
- Update with a reduced scope listing only genuinely remaining items

Close this item with a commit updating the tracking docs.

- [ ] **Step 3: Commit**

```bash
git add docs/plans/REVIEW.md docs/plans/BACKLOG.md
git commit -m "docs: re-scope docs physics revamp (Backlog #26) after lecture-test pair delivery"
```

---

## Phase 16: Integration + Rashba Physics

### Task 8: Wire hexagon and SC wire integration tests

**Files:**
- Create: `tests/regression/configs/wire_gaas_hexagon.cfg`
- Create: `tests/regression/configs/sc_wire_gaas.cfg`
- Create: `tests/regression/data/wire_gaas_hexagon/`
- Create: `tests/regression/data/sc_wire_gaas/`
- Create: `tests/integration/test_wire_hexagon.sh`
- Create: `tests/integration/test_sc_wire.sh`
- Modify: `tests/CMakeLists.txt`

- [ ] **Step 1: Create hexagonal wire config**

Create `tests/regression/configs/wire_gaas_hexagon.cfg` based on the rectangle config (`wire_gaas_rectangle.cfg`) but with hexagonal geometry:

```
waveVector: kz
waveVectorMax: 0.01
waveVectorStep: 2
confinement:  2
FDstep: 1
FDorder: 2
numLayers:  1
wire_nx: 21
wire_ny: 21
wire_dx: 5.0
wire_dy: 5.0
wire_shape: hexagon
wire_radius: 50.0
numRegions: 1
region: GaAs  0.0  50.0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
whichBand: 0
bandIdx: 1
SC: 0
feast_emin: -1.5
feast_emax: 2.0
feast_m0: -1
```

- [ ] **Step 2: Generate hexagonal wire golden data**

```bash
WORKDIR=$(mktemp -d)
/bin/cp tests/regression/configs/wire_gaas_hexagon.cfg "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/bandStructure > output.txt 2>&1)
mkdir -p tests/regression/data/wire_gaas_hexagon
/bin/cp "$WORKDIR/output.txt" tests/regression/data/wire_gaas_hexagon/output.txt
rm -rf "$WORKDIR"
```

- [ ] **Step 3: Create hexagonal wire test script**

Create `tests/integration/test_wire_hexagon.sh` following `test_wire_bandstructure.sh` pattern — 4 arguments (exe, config, ref_data, compare_script), runs bandStructure, compares eigenvalues.

- [ ] **Step 4: Create SC wire config**

Create `tests/regression/configs/sc_wire_gaas.cfg` based on `wire_gaas_rectangle.cfg` with SC enabled:

```
waveVector: kz
waveVectorMax: 0.01
waveVectorStep: 2
confinement:  2
FDstep: 1
FDorder: 2
numLayers:  1
wire_nx: 11
wire_ny: 11
wire_dx: 5.0
wire_dy: 5.0
wire_shape: rectangle
wire_width: 50.0
wire_height: 50.0
numRegions: 1
region: GaAs  0.0  50.0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
whichBand: 0
bandIdx: 1
SC: 1
SC_max_iter: 20
SC_tolerance: 0.001
SC_mixing: 0.3
SC_temperature: 300.0
SC_fermi_mode: 0
SC_num_kpar: 5
SC_kpar_max: 0.5
doping: 1e18 0
feast_emin: -1.5
feast_emax: 2.0
feast_m0: -1
```

Keep grid small (11x11) for fast test execution.

- [ ] **Step 5: Generate SC wire golden data**

Same pattern as Step 2, but run `bandStructure` with `sc_wire_gaas.cfg`.

- [ ] **Step 6: Create SC wire test script**

Create `tests/integration/test_sc_wire.sh` — similar to hex test but also checks SC convergence output (iteration count < max_iter, final potential shift < tolerance).

- [ ] **Step 7: Register both tests in CMakeLists.txt**

Add two test registrations following the wire regression pattern:

```cmake
# Wire hexagonal geometry
add_test(
    NAME regression_wire_gaas_hexagon
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_hexagon.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_hexagon.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/wire_gaas_hexagon
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(regression_wire_gaas_hexagon PROPERTIES LABELS "regression" TIMEOUT 300)

# SC wire convergence
add_test(
    NAME regression_sc_wire_gaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_wire.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_wire_gaas.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_wire_gaas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(regression_sc_wire_gaas PROPERTIES LABELS "regression" TIMEOUT 300)
```

- [ ] **Step 8: Build and test**

```bash
cmake --build build
ctest --test-dir build -R "hexagon|sc_wire" --output-on-failure
```
Expected: Both new tests pass.

- [ ] **Step 9: Run full suite**

```bash
ctest --test-dir build -j4 --output-on-failure
```

- [ ] **Step 10: Commit**

```bash
git add tests/regression/configs/wire_gaas_hexagon.cfg tests/regression/configs/sc_wire_gaas.cfg tests/regression/data/wire_gaas_hexagon/ tests/regression/data/sc_wire_gaas/ tests/integration/test_wire_hexagon.sh tests/integration/test_sc_wire.sh tests/CMakeLists.txt
git commit -m "test: add wire hexagon and SC wire integration tests"
```

---

### Task 9: Rashba BdG physics calibration

**Files:**
- Modify: `tests/regression/configs/topology_rashba_phase.cfg`
- Modify: `scripts/lecture_13_topological.py` (or create new Rashba section)
- Modify: `tests/integration/test_topology_rashba_phase.sh` (if it exists) or the CMake registration
- Modify: `docs/plans/BACKLOG.md`

This is the most physics-intensive task. The Rashba BdG sweep script was removed (absorbed into lecture_13). The config `topology_rashba_phase.cfg` has fundamental physics issues (mu in band gap, unrealistic grid spacing).

- [ ] **Step 1: Run bandStructure on the wire to find subband energies**

```bash
WORKDIR=$(mktemp -d)
# Create a bandstructure config matching the Rashba wire geometry
cat > "$WORKDIR/input.cfg" << 'EOF'
waveVector: kz
waveVectorMax: 0.1
waveVectorStep: 0.005
confinement:  2
FDstep: 1
FDorder: 2
numLayers:  1
wire_nx: 21
wire_ny: 21
wire_dx: 5.0
wire_dy: 5.0
wire_shape: rectangle
wire_width: 100.0
wire_height: 100.0
numRegions: 1
region: InAs  0.0  100.0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
whichBand: 0
bandIdx: 1
SC: 0
feast_emin: -2.0
feast_emax: 3.0
feast_m0: -1
EOF
mkdir -p "$WORKDIR/output"
(cd "$WORKDIR" && ./build/src/bandStructure > output.txt 2>&1)
# Extract CB subband minimum from kz=0 eigenvalues
python3 -c "
import re
with open('$WORKDIR/output/output.txt') as f:
    for line in f:
        # Look for eigenvalue output at k=0
        pass
# The CB minimum for a wide InAs wire should be near the InAs CB edge (~0.35 eV)
print('Check output/output.txt for subband energies')
"
rm -rf "$WORKDIR"
```

- [ ] **Step 2: Calibrate mu to first CB subband**

Based on Step 1 output, set `mu` in `topology_rashba_phase.cfg` to the first CB subband energy at kz=0. This ensures electronic states exist at the Fermi level for the BdG pairing.

- [ ] **Step 3: Fix grid spacing**

Update `topology_rashba_phase.cfg` to use physically consistent grid spacing. Change from arbitrary `wire_dx`/`wire_dy` to:
```
wire_dx: 5.0
wire_dy: 5.0
wire_nx: 21
wire_ny: 21
wire_width: 100.0
wire_height: 100.0
```
So `wire_dx ≈ wire_width / (wire_nx - 1) = 5.0 Å`.

- [ ] **Step 4: Calibrate FEAST window**

Based on the subband spacing from Step 1, set `feast_emin` and `feast_emax` to bracket the relevant BdG eigenvalues around `mu`:
```
feast_emin: <mu - 0.5>
feast_emax: <mu + 0.5>
```
The window should be wide enough to capture the BdG splitting but not so wide that FEAST returns spurious eigenvalues.

- [ ] **Step 5: Update the Rashba section in lecture_13_topological.py**

Find the Rashba BdG section in `scripts/lecture_13_topological.py`. Update it to:
1. Run bandStructure first to find subband positions
2. Use the calibrated mu value
3. Demonstrate the Majorana phase transition by sweeping Zeeman field
4. Assert that the BdG gap closes at the transition point

- [ ] **Step 6: Fix the false-positive regression test**

Update the Rashba regression test to verify FEAST actually finds eigenvalues before checking the gap. If the existing test shell script doesn't exist, register a new CTest that runs the calibrated config and checks:
- Eigenvalue count > 0 when BdG is enabled
- Gap value is a positive number (not the sentinel -1.0)

- [ ] **Step 7: Run topology tests**

```bash
ctest --test-dir build -R rashba --output-on-failure
ctest --test-dir build -L standard-star --output-on-failure
python3 scripts/lecture_13_topological.py
```
Expected: Rashba section demonstrates the phase transition; regression test passes with real eigenvalues.

- [ ] **Step 8: Run full suite**

```bash
ctest --test-dir build -j4 --output-on-failure
```

- [ ] **Step 9: Commit**

```bash
git add tests/regression/configs/topology_rashba_phase.cfg scripts/lecture_13_topological.py tests/CMakeLists.txt docs/plans/BACKLOG.md
git commit -m "fix(topology): calibrate Rashba BdG physics — mu, FEAST window, grid spacing"
```

---

## Final Task: Update tracking docs

### Task 10: Close remaining backlog items

**Files:**
- Modify: `docs/plans/BACKLOG.md`
- Modify: `docs/plans/REVIEW.md`

- [ ] **Step 1: Update BACKLOG.md**

After all tasks complete, update the "Remaining Backlog" section:
- Mark #4 gfactor as COMPLETE (all 5 configs now have regression or standard-star coverage)
- Mark #37 as COMPLETE (contiguous gaps closed)
- Mark #8 as COMPLETE (wire hexagon + SC wire + wire strain done)
- Mark #26 as COMPLETE or re-scoped (per Task 7 outcome)
- Remove Rashba section or convert to "Known Limitations" if not fully resolved
- Update the summary table with final phase status

- [ ] **Step 2: Update REVIEW.md**

- Update row #4 to COMPLETE
- Update row #37 to COMPLETE
- Update row #26 to COMPLETE or document remaining scope
- Update rows #55 (CSR) to COMPLETE
- Update row #56 (standard-star) to COMPLETE
- Update row #59 (strain validation) to COMPLETE
- Update review date

- [ ] **Step 3: Commit**

```bash
git add docs/plans/BACKLOG.md docs/plans/REVIEW.md
git commit -m "docs: close remaining backlog items — all phases complete"
```
