# Magnetic Field Integration — Landau Levels Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Enable magnetic field effects (Landau quantization + Zeeman splitting) across all three confinement modes: Wire (2D grid), QW (1D grid), and Bulk (k-space).

**Architecture:** Phase 1 verifies existing wire Peierls via regression test. Phase 2 adds Zeeman to QW mode. Phase 3 implements material-generic analytical Landau markers for bulk. Phase 4 executes the BdG sweep script for Majorana phase diagram. Phase 5 is deferred.

**Tech Stack:** Fortran 2018 (hamiltonianConstructor.f90), Python (verify_landau_levels.py, sweep_rashra_bdg.py), CMake/Ninja build.

---

## Task 1: Verify Wire Peierls (Phase 1)

**Files:**
- Run: `bash tests/integration/test_bhz_z2.sh build/src/topologicalAnalysis tests/regression/configs/topology_bhz_z2_trivial.cfg tests/regression/configs/topology_bhz_z2_topological.cfg`

**Step 1: Run BHZ Z2 integration test**

Run: `bash tests/integration/test_bhz_z2.sh build/src/topologicalAnalysis tests/regression/configs/topology_bhz_z2_trivial.cfg tests/regression/configs/topology_bhz_z2_topological.cfg`

Expected: Z2=0 for trivial (d=58Å), Z2=1 for topological (d=70Å)

**Step 2: Commit verification**

None needed — this is a read-only verification step.

---

## Task 2: Add Zeeman to QW Mode (Phase 2)

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:430` (end of ZB8bandQW)
- Test: `tests/regression/configs/landau_InAs.cfg` (already exists)

**Step 1: Add Zeeman block to ZB8bandQW**

Read `src/physics/hamiltonianConstructor.f90` around line 430 (end of ZB8bandQW, before optional `g` derivative argument).

Add after the Hamiltonian matrix elements are set:

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
      Vz(1:2) = -1.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! HH
      Vz(3:4) =  0.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! LH
      Vz(5:6) = -0.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! SO
      Vz(7)   = -1.0_dp * cfg%bdg%g_factor * mu_B * B_z  ! CB1
      Vz(8)   =  1.0_dp * cfg%bdg%g_factor * mu_B * B_z  ! CB2
      do i = 1, N
        HT((i-1)*8+1:(i-1)*8+8, (i-1)*8+1:(i-1)*8+8) = &
          HT((i-1)*8+1:(i-1)*8+8, (i-1)*8+1:(i-1)*8+8) + diag(Vz)
      end do
    end block
  end if
end if
```

**Step 2: Verify g_factor is accessible**

Confirm `cfg%bdg%g_factor` exists in `bdg_config` type in `defs.f90`. It should already exist from prior commit fbab6fe.

Run: `grep -n "g_factor" src/core/defs.f90 | grep bdg`

Expected: line with `g_factor` in `bdg_config` derived type

**Step 3: Build and test**

Run: `cmake --build build && ctest -R "test_qw|test_hamiltonian" --output-on-failure`

Expected: All QW/hamiltonian tests pass

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "feat: add Zeeman splitting to QW mode for B_perp"
```

---

## Task 3: Material-Generic Landau Level Verification (Phase 3)

**Files:**
- Modify: `scripts/verify_landau_levels.py`
- Config: `tests/regression/configs/landau_InAs.cfg` (already exists)

**Step 1: Update verify_landau_levels.py**

Read `scripts/verify_landau_levels.py` and update to:

1. Use `meff` from material database if available, else compute `m* = 1/gamma1`
2. Compute cyclotron frequency: `hbar_omega = hbar * e * B / (m* * m0)` in meV
3. Plot analytical Landau levels: E_n = CB_edge + (n + 0.5) * hbar_omega
4. Support InAs (meff=0.026), GaAs (meff=0.067), InP (meff=0.077)

```python
MATERIAL_DB = {
    'InAs': {'meff': 0.026, 'gamma1': 19.0},
    'GaAs': {'meff': 0.067, 'gamma1': 6.95},
    'InP':  {'meff': 0.077, 'gamma1': 5.08},
}

def compute_cyclotron_energy(m_eff, B):
    """Returns hbar*omega_c in meV."""
    hbar_J = 1.054571817e-34
    e_C = 1.602176634e-19
    m_kg = m_eff * 9.109e-31
    omega_c = e_C * B / m_kg
    hbar_omega_J = hbar_J * omega_c
    return hbar_omega_J * 6.241509e18  # convert to meV
```

**Step 2: Test with InAs at B=5T**

Run: `python scripts/verify_landau_levels.py`

Expected output:
- InAs at B=5T: hbar_omega_c ≈ 22.26 meV, E_0 ≈ 11.13 meV above CB
- GaAs at B=5T: hbar_omega_c ≈ 8.67 meV, E_0 ≈ 4.34 meV above CB

**Step 3: Commit**

```bash
git add scripts/verify_landau_levels.py
git commit -m "feat: make Landau level verification material-generic"
```

---

## Task 4: BdG Majorana Phase Diagram (Phase 4)

**Files:**
- Execute: `scripts/sweep_rashra_bdg.py`
- Output: `docs/lecture/figures/rashba_majorana_phase_diagram.png`

**Step 1: Quick stability test**

Create temp directory and run single BdG config:

```bash
cd /tmp && mkdir -p test_bdg && cd test_bdg && mkdir output
cat > input.cfg << 'EOF'
waveVector: kz
confinement: 2
wire_nx: 21; wire_ny: 21; wire_width: 50.0
material1: GaAs
numcb: 4; numvb: 4
topology: T
mode: bdg
compute_z2: F
compute_ldos: F
bdg T
mu: 0.0005
delta_0: 0.0003
g_factor: 2.0
EOF
/data/8bandkp-fdm/build/src/topologicalAnalysis
echo "Exit code: $?"
```

Expected: Exit code 0, `output/topology_result.dat` exists

**Step 2: Execute sweep script**

Run: `python scripts/sweep_rashra_bdg.py`

Expected: Sweeps B from 0 to 10T in 0.5T steps, generates `docs/lecture/figures/rashra_majorana_phase_diagram.png`

**Step 3: Verify figure**

Check that `docs/lecture/figures/rashra_majorana_phase_diagram.png` exists and is non-empty.

Expected: Plot shows gap closure near B_crit ≈ 5.0 T

**Step 4: Commit**

```bash
git add docs/lecture/figures/rashba_majorana_phase_diagram.png docs/lecture/figures/rashra_majorana_phase_diagram.txt
git commit -m "feat: add BdG Majorana phase diagram"
```

---

## Task 5: QW Z2 Fu-Kane — Deferred (Phase 5)

No implementation. Stub `compute_z2_fukane` in `topological_analysis.f90:169-176` remains as documentation.

Decision: QW mode is primarily for optical spectra, not topological analysis. Z2 via Fu-Kane requires eigenvector parity analysis at 4 TRIM points, significant implementation work. Deferred to future.

---

## Summary

| Task | Status | Notes |
|------|--------|-------|
| 1 | Verify Wire Peierls | Run existing test |
| 2 | QW Zeeman | ~15 lines Fortran |
| 3 | Bulk Landau markers | Python post-processing |
| 4 | BdG phase diagram | Execute sweep script |
| 5 | QW Z2 Fu-Kane | Deferred |

**Verification Commands:**
```bash
# Run BHZ Z2 test
bash tests/integration/test_bhz_z2.sh build/src/topologicalAnalysis \
  tests/regression/configs/topology_bhz_z2_trivial.cfg \
  tests/regression/configs/topology_bhz_z2_topological.cfg

# Build and test
cmake --build build && ctest --output-on-failure
```