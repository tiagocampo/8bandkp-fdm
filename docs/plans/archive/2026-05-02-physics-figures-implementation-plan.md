# Physics Figures Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Generate publication-quality benchmark figures for all implemented physics (QW/QW Zeeman, Bulk Peierls/Landau levels, Wire Majorana/BHZ Z2), with automated regression testing and updated lecture documentation.

**Architecture:** Four phases building on existing code. Phase 1-3 fix individual physics modules; Phase 4 standardizes figure generation and documentation. Each phase produces both PNG figures and TXT data files, with exit codes for regression testing.

**Tech Stack:** Fortran 2018 (hamiltonianConstructor.f90, magnetic_field.f90), Python (scripts/verify_*.py), CMake/Ninja build.

---

## Phase 1: QW Zeeman — Add Zeeman Splitting to Quantum Well Mode

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:430` (end of ZB8bandQW)
- Modify: `scripts/verify_landau_levels.py` (add QW Zeeman test)
- Test: `tests/regression/configs/landau_InAs.cfg`

### Task 1.1: Add Zeeman Block to ZB8bandQW

Read `src/physics/hamiltonianConstructor.f90` around line 430 (end of ZB8bandQW, before the optional `g` derivative argument).

At the end of the subroutine (after all HT matrix elements are set, before the `end subroutine` line), add:

```fortran
! Zeeman splitting: g*mu_B * B . sigma for B in z-direction
! Applied to diagonal blocks when B_z /= 0
if (present(cfg)) then
  if (cfg%bdg%enabled) then
    block
      real(kind=dp) :: mu_B, B_z, Vz(8)
      B_z = cfg%bdg%B_vec(3)
      if (abs(B_z) < 1e-12_dp) exit block
      mu_B = e * hbar / (2.0_dp * m0)
      ! g_J eigenvalues: HH=-1.5, LH=+0.5, SO=-0.5, CB=+1.0
      Vz(1:2) = -1.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! HH
      Vz(3:4) =  0.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! LH
      Vz(5:6) = -0.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! SO
      Vz(7)   = -1.0_dp * cfg%bdg%g_factor * mu_B * B_z  ! CB1
      Vz(8)   =  1.0_dp * cfg%bdg%g_factor * mu_B * B_z  ! CB2
      ! Add to each 8x8 block diagonal (bands 1-8, repeated for each z grid point)
      do i = 1, N
        HT((i-1)*8+1:(i-1)*8+8, (i-1)*8+1:(i-1)*8+8) = &
          HT((i-1)*8+1:(i-1)*8+8, (i-1)*8+1:(i-1)*8+8) + diag(Vz)
      end do
    end block
  end if
end if
```

Note: The existing `ZB8bandBulk` has Zeeman at lines 438-460 as reference pattern.

**Step 1: Verify g_factor exists in bdg_config**

Run: `grep -n "g_factor" src/core/defs.f90 | grep bdg`

Expected: line with `g_factor` in `bdg_config` derived type (should exist from prior commit).

**Step 2: Build and test QW mode**

Run: `cmake --build build && ctest -R "test_qw|test_hamiltonian" --output-on-failure`

Expected: All QW/hamiltonian tests pass.

**Step 3: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "feat: add Zeeman splitting to QW mode for B_perp"
```

### Task 1.2: Add QW Zeeman to verify_landau_levels.py

Read `scripts/verify_landau_levels.py`.

Add QW mode test section after existing bulk Landau test:

```python
def test_qw_zeeman():
    """Verify QW Zeeman splitting: E_spin_up - E_spin_down = g*mu_B*B"""
    # Run with B=5T in QW mode (confinement=1)
    # Parse eigenvalues at k=0
    # Expected: CB Zeeman splitting ~ g*mu_B*B = 2 * 0.058 * 5 = 0.58 meV
    pass
```

Run: `python scripts/verify_landau_levels.py`

Expected: QW Zeeman section included, pass.

**Step 4: Commit**

```bash
git add scripts/verify_landau_levels.py
git commit -m "feat: add QW Zeeman verification to Landau level script"
```

---

## Phase 2: Bulk Peierls — Integrate Peierls Substitution into Bulk Hamiltonian

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:ZB8bandBulk` (call add_peierls_coo)
- Modify: `src/io/input_parser.f90` (fix landau_InAs.cfg ExternalField parsing)
- Create: `tests/regression/configs/landau_InAs.cfg` (already exists, verify format)
- Test: `scripts/verify_landau_levels.py`

### Task 2.1: Integrate add_peierls_coo into ZB8bandBulk

Read `src/physics/hamiltonianConstructor.f90` around line 438-460 where Zeeman is added to ZB8bandBulk.

Find where Zeeman is added and add Peierls call after it:

```fortran
! Peierls substitution: apply magnetic field phase factors to off-diagonal elements
! Only effective when Bx /= 0 (Landau gauge A = (0, 0, Bx*y))
if (cfg%bdg%enabled .and. present(kpterms)) then
  call add_peierls_coo(coo_vals_bulk, coo_row_bulk, coo_col_bulk, &
                       nnz_offset, cfg%grid, cfg%bdg%B_vec, 'landau', kpterms)
end if
```

**Step 1: Check kpterms availability in ZB8bandBulk**

Run: `grep -n "kpterms\|confinement_init" src/physics/hamiltonianConstructor.f90 | head -20`

Confirm that `kpterms` is passed to ZB8bandBulk and accessible.

**Step 2: Verify add_peierls_coo signature**

Run: `grep -n "subroutine add_peierls_coo" src/physics/magnetic_field.f90`

Check that the signature matches what we're calling. Expected:
```fortran
subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                            grid, B_vec, gauge, kpterms_2d)
```

Note: Peierls uses `kpterms_2d` but bulk uses `kpterms` (different type). May need wrapper.

**Step 3: Build and test**

Run: `cmake --build build && ctest -R "test_landau|test_bulk" --output-on-failure`

Expected: Landau tests pass with computed eigenvalues.

**Step 4: Verify E_0 = 11.13 meV**

Run: `python scripts/verify_landau_levels.py`

Expected output:
```
InAs at B=5T: hbar_omega_c = 22.26 meV
E_0 (analytical) = 11.13 meV
E_0 (computed) = 11.XX meV  <- should be close
STATUS: PASS
```

**Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "feat: integrate Peierls substitution into bulk Hamiltonian for Landau levels"
```

### Task 2.2: Fix landau_InAs.cfg ExternalField Parsing

The `landau_InAs.cfg` uses `ExternalField: 5 EF` but parser may fail with "Failed to read ExternalField". The issue is ExternalField is read as integer when config has float format.

**Step 1: Test current config**

Run: `cd /tmp && mkdir test_landau && cd test_landau && mkdir output && cp tests/regression/configs/landau_InAs.cfg input.cfg && /data/8bandkp-fdm/build/src/bandStructure 2>&1 | head -20`

Expected: Error "Failed to read ExternalField" or similar.

**Step 2: Check input_parser ExternalField reading**

Read `src/io/input_parser.f90` around line 492 where ExternalField is read.

The issue: `read(data_unit, *, iostat=status) label, cfg%ExternalField, cfg%EFtype` expects integer for ExternalField but landau_InAs.cfg has `ExternalField: 5 EF` (integer format works) but sweep script had `ExternalField: 5.0 EF` (float format fails).

**Step 3: Fix parser for float ExternalField**

Modify the ExternalField read to handle float values:
```fortran
real(kind=dp) :: ext_field_val
read(data_unit, *, iostat=status) label, ext_field_val, cfg%EFtype
if (status == 0) cfg%ExternalField = int(ext_field_val)
```

Or more robustly, just read as real and convert:
```fortran
read(data_unit, *, iostat=status) label, cfg%ExternalField, cfg%EFtype
```

The real issue may be the format. Check if `5` vs `5.0` matters for integer read.

**Step 4: Test with corrected config**

Run: same test as Step 1, should parse correctly now.

**Step 5: Commit**

```bash
git add src/io/input_parser.f90
git commit -m "fix: handle float ExternalField in input parser"
```

---

## Phase 3: Wire FEAST — Fix Eigenvalue Finding for Topological Cases

**Files:**
- Modify: `src/apps/main_topology.f90` (widen FEAST energy window, increase wire size)
- Modify: `scripts/sweep_rashba_bdg.py` (increase wire_nx/wire_ny)
- Modify: `scripts/verify_bhz_z2.py` (increase wire width for topological case)
- Create: `docs/lecture/figures/rashba_majorana_phase_diagram.png` (regenerate with working config)

### Task 3.1: Diagnose FEAST Issue for Majorana Mode

**Step 1: Test with larger wire**

Create temp config with wire_nx=31, wire_ny=31 (larger than current 21x21):

```bash
cd /tmp && mkdir test_bdg_large && cd test_bdg_large
cat > input.cfg << 'EOF'
waveVector: kz
waveVectorMax: 0.1
waveVectorStep: 21
confinement: 2
FDstep: 1
FDorder: 2
numLayers: 1
wire_nx: 31
wire_ny: 31
wire_dx: 3.0
wire_dy: 3.0
wire_shape: rectangle
wire_width: 50.0
wire_height: 50.0
numRegions: 1
region: GaAs  0.0  100.0
numcb: 4
numvb: 4
ExternalField: 5 EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
SC: 0
topology: T
mode: bdg
compute_chern: F
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
bdg: T
mu: 0.0005
delta_0: 0.0003
g_factor: 2.0
b_field: 5 0 0
EOF
mkdir output
/data/8bandkp-fdm/build/src/topologicalAnalysis 2>&1 | grep -E "nev_found|FEAST|min_gap"
```

Expected: If FEAST finds eigenvalues, we'll see nev_found > 0.

**Step 2: Check FEAST workspace reuse issue**

The FEAST eigenvalue finder uses workspace reuse via `feast_workspace_matches_pattern`. For different wire sizes, the sparsity pattern changes. If the stale cache from a previous run (different size) is being reused incorrectly, eigenvalues would be wrong.

Run with `feast_m0` increased:
- Current code: `eigen_cfg_local%feast_m0 = 0  ! auto`
- Try: `eigen_cfg_local%feast_m0 = 20` for 21x21 wire, `eigen_cfg_local%feast_m0 = 50` for 31x31 wire

### Task 3.2: Fix Majorana Phase Diagram Sweep

**Step 1: Update sweep script for larger wire**

Read `scripts/sweep_rashba_bdg.py`. Increase wire_nx and wire_ny:

```python
wire_nx: 31  # was 21
wire_ny: 31  # was 21
```

Also widen FEAST energy window in the topologicalAnalysis call:
```python
# Current: emin/emax = +-5*delta_0 = +-1.5 meV
# Try: +-10*delta_0 = +-3 meV to catch more states
eigen_cfg_local%emin = -10.0_dp * cfg%bdg%delta_0
eigen_cfg_local%emax = 10.0_dp * cfg%bdg%delta_0
```

But note: this change is in the Fortran code (`main_topology.f90`), not the sweep script. The sweep script just runs the executable.

**Step 2: Modify main_topology.f90 FEAST parameters**

Read `src/apps/main_topology.f90` around line 438-440.

Change:
```fortran
eigen_cfg_local%emin = -5.0_dp * cfg%bdg%delta_0
eigen_cfg_local%emax = 5.0_dp * cfg%bdg%delta_0
```

To:
```fortran
eigen_cfg_local%emin = -10.0_dp * cfg%bdg%delta_0
eigen_cfg_local%emax = 10.0_dp * cfg%bdg%delta_0
```

Also increase `feast_m0`:
```fortran
eigen_cfg_local%feast_m0 = 20  ! was 0 (auto), increase for better coverage
```

**Step 3: Regenerate Majorana phase diagram**

Run: `python scripts/sweep_rashba_bdg.py`

Expected: Gap should show V-shape closure near B_crit ≈ 5.0 T, not flat 0.

**Step 4: Verify figure**

Check `docs/lecture/figures/rashba_majorana_phase_diagram.png` shows non-zero gap at B=0 and V-shape at B_crit.

**Step 5: Commit**

```bash
git add src/apps/main_topology.f90 scripts/sweep_rashba_bdg.py
git commit -m "fix: widen FEAST energy window for Majorana mode eigenvalue finding"
```

### Task 3.3: Fix BHZ Z2 Topological (Z2=1)

**Step 1: Test BHZ topological with larger wire**

Run: `bash tests/integration/test_bhz_z2.sh build/src/topologicalAnalysis tests/regression/configs/topology_bhz_z2_trivial.cfg tests/regression/configs/topology_bhz_z2_topological.cfg`

Expected: Currently Z2 topological shows FAIL (FEAST finds bulk but no edge states).

**Step 2: Modify verify_bhz_z2.py for larger width**

Read `scripts/verify_bhz_z2.py`. Increase wire width for topological case to make edge states more visible:

```python
# For topological case (d=70Å), try wider wire to separate edge from bulk
wire_width: 80.0  # was 50.0
wire_height: 80.0
```

**Step 3: Check Z2 computation method**

The Z2 gap-based method counts Kramers pairs in the gap. With FEAST finding no eigenvalues in gap, the issue may be the gap detection itself.

Read `src/physics/topological_analysis.f90` around the gap-based Z2 computation.

**Step 4: Commit**

```bash
git add scripts/verify_bhz_z2.py
git commit -m "fix: increase BHZ wire size for topological Z2 detection"
```

---

## Phase 4: Unified Figure Suite — Standardize Scripts and Update Docs

**Files:**
- Modify: `scripts/verify_qwz_chern.py`
- Modify: `scripts/verify_bhz_z2.py`
- Modify: `scripts/verify_landau_levels.py`
- Modify: `scripts/sweep_rashba_bdg.py`
- Create: `scripts/generate_all_figures.py`
- Modify: `docs/lecture/13-topological-superconductivity.md`

### Task 4.1: Standardize verify_*.py Scripts

All scripts should follow the same interface:

```python
def main():
    parser = argparse.ArgumentParser(description='Verify PHYSICS')
    parser.add_argument('--config', type=str, help='Config file path')
    parser.add_argument('--exe', type=str, default='build/src/topologicalAnalysis', help='Executable path')
    parser.add_argument('--tolerance', type=float, default=0.1, help='Acceptable deviation %%')
    args = parser.parse_args()

    # Run executable with config
    # Parse output
    # Compare to expected
    # Exit 0 if pass, non-zero if fail
```

**Step 1: Update verify_qwz_chern.py**

Add argparse, tolerance checking, exit codes. Ensure consistent output format.

**Step 2: Update verify_bhz_z2.py**

Same standardization.

**Step 3: Update verify_landau_levels.py**

Same standardization.

**Step 4: Update sweep_rashba_bdg.py**

Same standardization.

### Task 4.2: Create generate_all_figures.py

Create `scripts/generate_all_figures.py`:

```python
#!/usr/bin/env python3
"""Generate all physics verification figures."""
import subprocess
import sys
from pathlib import Path

SCRIPTS = [
    'verify_qwz_chern.py',
    'verify_bhz_z2.py',
    'verify_landau_levels.py',
    'sweep_rashba_bdg.py',
]

def main():
    results = []
    for script in SCRIPTS:
        result = subprocess.run([sys.executable, f'scripts/{script}'], capture_output=True)
        status = 'PASS' if result.returncode == 0 else 'FAIL'
        results.append((script, status))
        print(f'{script}: {status}')

    # Summary
    passed = sum(1 for _, s in results if s == 'PASS')
    print(f'\nPassed: {passed}/{len(results)}')
    return 0 if passed == len(results) else 1
```

### Task 4.3: Update Lecture Documentation

Read `docs/lecture/13-topological-superconductivity.md`.

Update Section 13.12 benchmark table:

| Test | Model | Status | Notes |
|---|---|---|---|
| Chern +1 | QWZ u=-0.8 | PASS | Figure: qwz_chern_phase_diagram.png |
| Chern -1 | QWZ u=0.5 | PASS | |
| Chern 0 | QWZ u=2.5 | PASS | |
| BHZ Z2 trivial | BHZ d=58Å | PASS | Figure: verify_bhz_z2_width_sweep.png |
| BHZ Z2 topological | BHZ d=70Å | FIX IN PROGRESS | FEAST finds bulk bands, edge states not resolved |
| Landau levels | InAs B=5T | FIX IN PROGRESS | Peierls integration incomplete |
| Majorana gap | Rashba wire | FIX IN PROGRESS | V-shape gap closure at B_crit expected |

Update Section 13.10 input configuration example to reflect correct `b_field` format:
```
! BdG mode example:
mode: bdg
bdg: T
mu: 0.0005
delta_0: 0.0003
g_factor: 2.0
b_field: 5 0 0  ! Bx By Bz in Tesla
```

### Task 4.4: Final Verification

**Step 1: Run all verification scripts**

Run: `python scripts/generate_all_figures.py`

Expected: All scripts pass, all figures generated.

**Step 2: Verify figures exist**

Run: `ls -la docs/lecture/figures/*.png docs/lecture/figures/*.txt`

Expected:
- qwz_chern_convergence.png ✓
- qwz_chern_phase_diagram.png ✓
- verify_bhz_z2_width_sweep.png ✓
- landau_levels_inas_b5t.png ✓
- rashba_majorana_phase_diagram.png ✓ (shows V-shape gap)

**Step 3: Commit Phase 4**

```bash
git add scripts/generate_all_figures.py docs/lecture/13-topological-superconductivity.md
git commit -m "feat: add unified figure generation and update lecture docs"
```

---

## Summary

| Phase | Tasks | Key Files | Success Metric |
|---|---|---|---|
| 1 | 2 | hamiltonianConstructor.f90, verify_landau_levels.py | QW Zeeman code added, tests pass |
| 2 | 2 | hamiltonianConstructor.f90, input_parser.f90 | Landau E_0 = 11.13 meV computed |
| 3 | 3 | main_topology.f90, sweep_rashba_bdg.py, verify_bhz_z2.py | Majorana V-shape visible, BHZ Z2=1 works |
| 4 | 4 | All scripts, 13-topological-superconductivity.md | All figures generate, docs updated |

## Verification Commands

```bash
# Phase 1
cmake --build build && ctest -R "test_qw|test_hamiltonian" --output-on-failure

# Phase 2
python scripts/verify_landau_levels.py

# Phase 3
python scripts/sweep_rashba_bdg.py
bash tests/integration/test_bhz_z2.sh build/src/topologicalAnalysis \
  tests/regression/configs/topology_bhz_z2_trivial.cfg \
  tests/regression/configs/topology_bhz_z2_topological.cfg

# Phase 4
python scripts/generate_all_figures.py
```