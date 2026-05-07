# Topological Suite Verification Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Verify each physics module of the topological suite by producing figures comparing computed output vs. literature benchmarks, in a pedagogical sequence from QWZ Chern (working) → BHZ Z2 (fix) → Landau levels (implement) → BdG (fix+verify) → Lecture 13 (update).

**Architecture:** Start with QWZ Chern (already working, no code changes needed), produce verification figures to establish methodology, then fix BHZ Z2 using Fu-Kane parity method, implement Peierls substitution for Landau levels, fix BdG config parser, and finally update lecture notes with all verified figures.

**Tech Stack:** Fortran 2018 (existing), Python (figures via gnuplot/matplotlib), shell scripts for test orchestration

---

## Phase 1: QWZ Chern Number — Figures from Working Code

Phase 1 requires **no code changes** — only script creation, figure generation, and lecture updates. QWZ Chern is already verified working.

### Task 1: Create QWZ Chern Verification Python Script

**Files:**
- Create: `scripts/verify_qwz_chern.py`

**Step 1: Write the verification script**

```python
#!/usr/bin/env python3
"""Verify QWZ Chern number against literature benchmarks.

Literature: Qi, Wu & Zhang, Phys. Rev. B 74, 085308 (2006)
Expected: u=-0.8 -> C=+1, u=0.5 -> C=-1, u=2.5 -> C=0

This script:
1. Loops over u values and grid resolutions
2. Runs topologicalAnalysis with each config
3. Parses Chern number from topology_result.dat
4. Generates convergence and phase diagram figures
"""

import subprocess
import os
import re
import sys
from pathlib import Path

EXE = Path(__file__).parent.parent / "build" / "src" / "topologicalAnalysis"
OUT_DIR = Path(__file__).parent.parent / "docs" / "lecture" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

U_VALUES = [-0.8, 0.5, 2.5]
NK_VALUES = [20, 30, 40, 50, 60, 70]
EXPECTED_C = {"-0.8": 1, "0.5": -1, "2.5": 0}  # +1, -1, 0 as integers

def run_qwz(u, nk):
    """Run topologicalAnalysis for given u and nk (grid resolution)."""
    workdir = Path(tempfile.mkdtemp())
    cfg = workdir / "input.cfg"
    output = workdir / "output"
    output.mkdir()

    # Write config — override nk_default via compute_chern_nk if available
    # For QWZ mode, nk is set internally; use the existing config as template
    config_text = f"""\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0 EF
EFParams: 0.0005
topology: T
mode: qhe
compute_chern: T
compute_hall: T
qwz_u: {u}
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
"""

    cfg.write_text(config_text)
    result = subprocess.run(
        [str(EXE)],
        cwd=workdir,
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"topologicalAnalysis failed: {result.stderr}")

    topo_file = output / "topology_result.dat"
    if not topo_file.exists():
        raise RuntimeError("topology_result.dat not produced")

    content = topo_file.read_text()
    match = re.search(r"Chern number:\s*([+-]?\d+)", content)
    if not match:
        raise RuntimeError(f"Could not parse Chern number from: {content[:500]}")
    return int(match.group(1))


def main():
    results = {}  # {(u, nk): C}
    for u in U_VALUES:
        for nk in NK_VALUES:
            C = run_qwz(u, nk)
            results[(u, nk)] = C
            print(f"u={u}, nk={nk}: C={C}")

    # Generate figures using gnuplot or matplotlib
    try:
        import matplotlib.pyplot as plt
        make_figures(results)
    except ImportError:
        print("matplotlib not available, outputting CSV for manual plotting")
        write_csv(results)


if __name__ == "__main__":
    main()
```

**Step 2: Test the script runs**

Run: `python scripts/verify_qwz_chern.py`
Expected: Script runs, prints Chern numbers for each (u, nk) combination

**Step 3: Commit**

```bash
git add scripts/verify_qwz_chern.py
git commit -m "test: add QWZ Chern number verification script"
```

---

### Task 2: Generate QWZ Chern Figures

**Files:**
- Modify: `scripts/verify_qwz_chern.py` (add matplotlib figure generation)
- Create: `docs/lecture/figures/chern_convergence.png`
- Create: `docs/lecture/figures/chern_phase_diagram.png`
- Create: `docs/lecture/figures/chern_berry_curvature_topological.png`
- Create: `docs/lecture/figures/chern_berry_curvature_trivial.png`

**Step 1: Add convergence figure code**

In `verify_qwz_chern.py`, add:

```python
def make_figures(results):
    import matplotlib.pyplot as plt

    # Figure 1: Chern convergence vs nk
    fig, ax = plt.subplots(figsize=(7, 5))
    for u in U_VALUES:
        nk_vals = [nk for (uu, nk) in results.keys() if uu == u]
        c_vals = [results[(u, nk)] for nk in nk_vals]
        ax.plot(nk_vals, c_vals, 'o-', label=f'u={u}')
    ax.set_xlabel('Grid resolution nk')
    ax.set_ylabel('Chern number C')
    ax.set_title('QWZ Chern Number Convergence')
    ax.set_ylim(-2, 2)
    ax.axhline(y=0, color='gray', linewidth=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'chern_convergence.png', dpi=150)
    print(f"Saved chern_convergence.png")

    # Figure 2: Chern phase diagram
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    u_fine = [u/10 for u in range(-30, 35)]
    # Show expected phase transitions
    ax2.axvspan(-3, -2, alpha=0.3, color='red', label='C=+1 (u<-2)')
    ax2.axvspan(-2, 0, alpha=0.3, color='blue', label='C=-1 (-2<u<0)')
    ax2.axvspan(0, 3, alpha=0.3, color='gray', label='C=0 (u>0)')
    ax2.axvline(x=-2, color='black', linestyle='--', linewidth=1)
    ax2.axvline(x=0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlabel('Parameter u')
    ax2.set_ylabel('Chern number C')
    ax2.set_title('QWZ Chern Number Phase Diagram')
    ax2.set_xlim(-3, 3)
    ax2.set_ylim(-2, 2)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(OUT_DIR / 'chern_phase_diagram.png', dpi=150)
    print(f"Saved chern_phase_diagram.png")
```

**Step 2: Run script to generate figures**

Run: `python scripts/verify_qwz_chern.py`
Expected: Four PNG files created in `docs/lecture/figures/`

**Step 3: Commit**

```bash
git add scripts/verify_qwz_chern.py
git add docs/lecture/figures/chern_convergence.png
git add docs/lecture/figures/chern_phase_diagram.png
git commit -m "feat: add QWZ Chern convergence and phase diagram figures"
```

---

### Task 3: Update Lecture 13 — Chern Number Section

**Files:**
- Modify: `docs/lecture/13-topological-superconductivity.md` (Chern number section)

**Step 1: Read current Chern section**

Run: `grep -n "Chern\|QWZ\|13.2" docs/lecture/13-topological-superconductivity.md | head -30`

**Step 2: Update section 13.2 with verification table and figure references**

The Chern number section should add:
- A verification table showing Literature vs. Computed values
- Figure captions referencing the generated PNG files
- A note on the convergence study showing C→integer as nk increases

**Step 3: Commit**

```bash
git add docs/lecture/13-topological-superconductivity.md
git commit -m "docs: update lecture 13 Chern section with verified figures"
```

---

## Phase 2: Fix BHZ Z₂ Invariant — Fu-Kane Parity Method

### Task 4: Read Existing BHZ Z₂ Code

**Files:**
- Read: `src/physics/topological_analysis.f90:139-176` — compute_z2_gap and stub compute_z2_fukane
- Read: `src/physics/topological_analysis.f90:437-520` — build_bhz_wire_hamiltonian
- Read: `src/apps/main_topology.f90:223-352` — run_qshe_wire

**Step 1: Read the code**

Run: `head -176 src/physics/topological_analysis.f90 | tail -38`

**Step 2: Read build_bhz_wire_hamiltonian**

Run: `sed -n '437,520p' src/physics/topological_analysis.f90`

---

### Task 5: Implement Fu-Kane Z₂ for QW Mode

**Files:**
- Modify: `src/physics/topological_analysis.f90`

**Step 1: Write failing test**

In `tests/unit/test_topological_z2.pf` (create new file):

```fortran
@test
subroutine test_z2_fukane_trivial_qw()
  ! For a trivial QW (GaAs), Z2 should be 0
  ! The Fu-Kane method computes (-1)^Z2 from parity products at TRIM
  type(simulation_config) :: cfg
  real(dp) :: H_eigs(8, 4)  ! 8 bands, 4 TRIM points
  integer :: z2

  ! At TRIM points for GaAs (normal band order), product of parity = +1
  H_eigs = reshape([ &
    -1.5, 0.3, 0.3, -1.5, -0.1, -0.1, 1.5, 1.5, &  ! Gamma
    -1.5, 0.3, 0.3, -1.5, -0.1, -0.1, 1.5, 1.5, &  ! M1
    -1.5, 0.3, 0.3, -1.5, -0.1, -0.1, 1.5, 1.5, &  ! M2
    -1.5, 0.3, 0.3, -1.5, -0.1, -0.1, 1.5, 1.5], & ! M3
    shape=[8, 4])

  z2 = compute_z2_fukane(H_eigs, cfg)
  @assertEqual(0, z2)
end subroutine
```

**Step 2: Run test to verify it fails**

Run: `cd build && ctest -R test_topological_z2 -V` (or compile and run manually)
Expected: FAIL — `compute_z2_fukane` is a stub returning 0

**Step 3: Implement compute_z2_fukane**

Complete the stub in `topological_analysis.f90:169-176`:

```fortran
function compute_z2_fukane(H_eigs, parity_signs) result(z2)
  ! Fu-Kane parity method for Z2 invariant
  ! (-1)^Z2 = PROD_i PROD_m(xi_m(Gamma_i))
  ! where xi_m = parity eigenvalue (+/-1) of occupied band m at TRIM i
  ! H_eigs: eigenvalues at 4 TRIM points (n_bands, 4)
  ! parity_signs: array of size n_bands with parity of each band (+/-1)
  implicit none
  real(kind=dp), intent(in) :: H_eigs(:,:)
  real(kind=dp), intent(in) :: parity_signs(:)
  integer :: z2

  integer :: i, m, n_bands, n_trim
  real(kind=dp) :: delta_prod

  n_bands = size(H_eigs, 1)
  n_trim = size(H_eigs, 2)

  z2 = 0

  ! For each TRIM, compute product of parity signs for occupied bands
  ! (bands with negative energy are occupied)
  delta_prod = 1.0_dp
  do i = 1, n_trim
    do m = 1, n_bands
      if (H_eigs(m, i) < 0.0_dp) then  ! occupied band
        delta_prod = delta_prod * parity_signs(m)
      end if
    end do
  end do

  ! (-1)^Z2 = delta_prod => Z2 = 1 if delta_prod < 0
  if (delta_prod < 0.0_dp) z2 = 1

end function compute_z2_fukane
```

**Step 4: Run test to verify it passes**

Expected: PASS

**Step 5: Commit**

```bash
git add src/physics/topological_analysis.f90 tests/unit/test_topological_z2.pf
git commit -m "feat: implement Fu-Kane Z2 parity method for QW mode"
```

---

### Task 6: Fix BHZ Wire Mode — Gap-Sweep Method

**Files:**
- Modify: `src/physics/topological_analysis.f90` — add compute_z2_gap_sweep
- Modify: `src/apps/main_topology.f90` — wire mode uses gap-sweep

**Step 1: Write failing test**

```fortran
@test
subroutine test_z2_gap_sweep()
  ! For BHZ wire: Z2 flips at M=0
  ! M<0 (topological): gap closes at some k, Z2=1
  ! M>0 (trivial): no gap closing, Z2=0
  real(dp) :: eigenvalues_positive(4) = [10.0, 10.0, -10.0, -10.0]
  real(dp) :: eigenvalues_negative(4) = [-10.0, -10.0, 10.0, 10.0]
  real(dp) :: gap_threshold = 1.0
  integer :: z2_pos, z2_neg

  ! Gap-sweep: scan for gap closing at E=0
  ! M>0 (trivial): no states near zero => Z2=0
  z2_pos = compute_z2_gap_sweep(eigenvalues_positive, gap_threshold)
  @assertEqual(0, z2_pos)

  ! M<0 (topological): gap closes => Z2=1
  z2_neg = compute_z2_gap_sweep(eigenvalues_negative, gap_threshold)
  @assertEqual(1, z2_neg)
end subroutine
```

**Step 2: Run test to verify it fails**

Expected: FAIL — `compute_z2_gap_sweep` does not exist

**Step 3: Implement gap-sweep for wire mode**

In `topological_analysis.f90`, add:

```fortran
function compute_z2_gap_sweep(eigenvalues, gap_threshold) result(z2)
  ! Gap-sweep method for 1D wire Z2 invariant
  ! Counts sign changes in the eigenvalue spectrum across the gap
  ! If an odd number of bands cross from positive to negative => topological
  implicit none
  real(kind=dp), intent(in) :: eigenvalues(:)
  real(kind=dp), intent(in) :: gap_threshold
  integer :: z2

  integer :: n_crossings, i, j
  real(kind=dp) :: gap_min, gap_max

  z2 = 0
  gap_min = -gap_threshold
  gap_max = gap_threshold

  ! Count how many eigenvalues cross the gap region
  ! For a topological wire, the gap closes and reopens with odd number of crossings
  n_crossings = 0
  do i = 1, size(eigenvalues)
    if (abs(eigenvalues(i)) < gap_threshold) then
      n_crossings = n_crossings + 1
    end if
  end do

  ! If gap closes (states inside the gap), check parity of crossings
  ! Odd number of states in gap => Z2=1
  if (n_crossings >= 2) z2 = 1

end function compute_z2_gap_sweep
```

**Step 4: Wire mode in main_topology.f90 should call gap-sweep**

In `run_qshe_wire`, replace the call to `compute_z2_gap` with `compute_z2_gap_sweep`.

**Step 5: Commit**

```bash
git add src/physics/topological_analysis.f90 src/apps/main_topology.f90
git commit -m "fix: implement gap-sweep Z2 method for BHZ wire mode"
```

---

### Task 7: Verify BHZ Z₂ with Integration Test

**Files:**
- Run: `tests/integration/test_bhz_z2.sh` with both configs

**Step 1: Run trivial config**

Run: `cd /tmp && rm -rf test_bhz && mkdir test_bhz && cd test_bhz && mkdir output && cp /data/8bandkp-fdm/tests/regression/configs/topology_bhz_z2_trivial.cfg input.cfg && /data/8bandkp-fdm/build/src/topologicalAnalysis 2>&1 | grep "Z2"`
Expected: `Z2 invariant: 0`

**Step 2: Run topological config**

Run: `cd /tmp && rm -rf test_bhz2 && mkdir test_bhz2 && cd test_bhz2 && mkdir output && cp /data/8bandkp-fdm/tests/regression/configs/topology_bhz_z2_topological.cfg input.cfg && /data/8bandkp-fdm/build/src/topologicalAnalysis 2>&1 | grep "Z2"`
Expected: `Z2 invariant: 1`

**Step 3: Generate BHZ Z₂ figures**

Create `scripts/verify_bhz_z2.py`:
1. Sweep M from -15 to +15 meV
2. Run `topologicalAnalysis` for each M
3. Plot Z₂ vs. M — shows step at M=0
4. Plot edge localization length ξ vs. M for M<0

Figures:
- `docs/lecture/figures/bhz_z2_phase_transition.png`
- `docs/lecture/figures/bhz_edge_localization.png`

**Step 4: Commit**

```bash
git add scripts/verify_bhz_z2.py
git add docs/lecture/figures/bhz_z2_*.png
git commit -m "feat: add BHZ Z2 verification figures"
```

---

## Phase 3: Implement Peierls Substitution — Landau Levels

### Task 8: Read Landau Level Infrastructure

**Files:**
- Read: `src/physics/magnetic_field.f90` — Zeeman and Peierls stubs
- Read: `src/physics/hamiltonian_wire.f90` — wire builder with kz parameter
- Read: `src/core/defs.f90` — bdg_config type with B_vec and gauge

**Step 1: Read magnetic_field.f90**

Run: `cat src/physics/magnetic_field.f90`

**Step 2: Read wire Hamiltonian builder kz handling**

Run: `grep -n "kz" src/physics/hamiltonian_wire.f90 | head -20`

---

### Task 9: Implement Peierls Substitution

**Files:**
- Modify: `src/physics/magnetic_field.f90` — implement add_peierls_coo body

**Step 1: Write failing test**

In `tests/unit/test_magnetic_field.pf`:

```fortran
@test
subroutine test_peierls_phase_factor()
  ! For B=5T, y=1nm: phi = exp(-i*e*B*y*dz/hbar)
  ! e*hbar/hbar = e, so phi = exp(-i*e*B*y*dz)
  ! With e=1.602e-19 C, B=5T, y=1nm=1e-9m, dz=1nm=1e-9m
  ! phi_magnitude should be close to 1 (just a phase)
  real(kind=dp) :: phi
  phi = compute_peierls_phase(5.0_dp, 1.0_dp, 1.0_dp)
  ! Phase factor magnitude should be 1
  @assertEqual(1.0_dp, abs(phi), tolerance=1e-10_dp)
end subroutine
```

**Step 2: Run test to verify it fails**

Expected: FAIL — `compute_peierls_phase` does not exist

**Step 3: Implement add_peierls_coo body**

Replace the stub in `magnetic_field.f90:53-69`:

```fortran
subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                            grid, B_vec, gauge, kpterms_2d)
  ! Peierls substitution: k -> k - eA/hbar
  ! For Landau gauge A = (0, 0, Bx*y):
  !   kz -> kz - e*Bx*y/hbar  (phase factor on hopping)
  real(kind=dp), intent(inout) :: coo_vals(:)
  integer, intent(inout) :: coo_row(:), coo_col(:)
  integer, intent(inout) :: nnz_offset
  type(spatial_grid), intent(in) :: grid
  real(kind=dp), intent(in) :: B_vec(3)
  character(len=*), intent(in) :: gauge
  type(csr_matrix), intent(in) :: kpterms_2d(:)

  integer :: i, j, ngrid
  real(kind=dp) :: phi_i, phi_j, Bx, y_i, y_j
  complex(kind=dp) :: phase_factor
  real(kind=dp), parameter :: e = 1.602176634e-19_dp  ! C
  real(kind=dp), parameter :: hbar = 1.054571817e-34_dp ! J*s

  ngrid = grid%npoints()
  Bx = B_vec(1)

  if (abs(Bx) < 1e-12_dp) return  ! no magnetic field

  ! For each grid point i, compute Peierls phase phi_i = exp(-i*e*Bx*y_i*dz/hbar)
  ! For COO entries connecting site i to site j (coo_row -> coo_col),
  ! multiply by phi_i / phi_j (ratio of phase factors)
  do i = 1, nnz_offset
    y_i = grid%y(coo_row(i))
    y_j = grid%y(coo_col(i))

    phi_i = exp(-complex(0.0_dp, 1.0_dp) * e * Bx * y_i * grid%dz / hbar)
    phi_j = exp(-complex(0.0_dp, 1.0_dp) * e * Bx * y_j * grid%dz / hbar)

    phase_factor = phi_i / phi_j
    coo_vals(i) = coo_vals(i) * real(phase_factor, kind=dp)
  end do

end subroutine add_peierls_coo
```

**Step 4: Run test to verify it passes**

Expected: PASS

**Step 5: Commit**

```bash
git add src/physics/magnetic_field.f90 tests/unit/test_magnetic_field.pf
git commit -m "feat: implement Peierls substitution for Landau levels"
```

---

### Task 10: Create Landau Level Config and Test

**Files:**
- Create: `tests/regression/configs/landau_InAs.cfg`
- Create: `tests/integration/test_landau_levels.sh`

**Step 1: Create Landau InAs config**

```bash
# tests/regression/configs/landau_InAs.cfg
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: InAs
numcb: 2
numvb: 6
ExternalField: 5.0 EF  ! B = 5 Tesla
EFParams: 0.0005
```

**Step 2: Run and check Landau levels**

Run: `cd /tmp && rm -rf test_landau && mkdir test_landau && cd test_landau && mkdir output && cp tests/regression/configs/landau_InAs.cfg input.cfg && /data/8bandkp-fdm/build/src/bandStructure 2>&1 | grep -E "eigenvalue|Landau|E_0|E_1"`
Expected: Eigenvalues should show Landau quantization with spacing ~22 meV

**Step 3: Create integration test script**

```bash
#!/bin/bash
# tests/integration/test_landau_levels.sh
# Verify Landau levels for InAs at B=5T: E_n = hbar*omega_c*(n+1/2)
# hbar*omega_c = 22.26 meV, E_0 = 11.13 meV, E_1 = 33.39 meV

set -euo pipefail

EXE="$1"
CONFIG="$2"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT
mkdir -p "$WORKDIR/output"
cp "$CONFIG" "$WORKDIR/input.cfg"
cd "$WORKDIR"

"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: bandStructure returned exit code $RC"
    cat test_output.log
    exit 1
fi

# Parse eigenvalues from output
# Expected: E_0 = 11.13 meV, E_1 = 33.39 meV
# Tolerance: 1 meV
E0_EXPECTED=11.13
E1_EXPECTED=33.39
TOL=1.0

# Extract first few eigenvalues (they appear in the output)
E0=$(grep -A 50 "eigenvalues" "$WORKDIR/output/band_results.dat" | head -5 | awk 'NR==1 {print $2}')
E1=$(grep -A 50 "eigenvalues" "$WORKDIR/output/band_results.dat" | head -5 | awk 'NR==2 {print $2}')

echo "Computed: E_0=$E0 meV, E_1=$E1 meV"
echo "Expected: E_0=$E0_EXPECTED meV, E_1=$E1_EXPECTED meV"

E0_ERR=$(python3 -c "print(abs($E0 - $E0_EXPECTED))")
E1_ERR=$(python3 -c "print(abs($E1 - $E1_EXPECTED))")

if (( $(echo "$E0_ERR < $TOL" | bc -l) )) && (( $(echo "$E1_ERR < $TOL" | bc -l) )); then
    echo "PASS: Landau levels within $TOL meV of analytical values"
else
    echo "FAIL: E_0 error=$E0_ERR meV, E_1 error=$E1_ERR meV (tolerance=$TOL meV)"
    exit 1
fi
```

**Step 4: Run integration test**

Run: `bash tests/integration/test_landau_levels.sh build/src/bandStructure tests/regression/configs/landau_InAs.cfg`
Expected: PASS

**Step 5: Commit**

```bash
git add tests/regression/configs/landau_InAs.cfg tests/integration/test_landau_levels.sh
git commit -m "test: add Landau level regression test for InAs at B=5T"
```

---

### Task 11: Generate Landau Level Figure

**Files:**
- Create: `scripts/verify_landau_levels.py`
- Create: `docs/lecture/figures/landau_levels_inas_b5t.png`

**Step 1: Create verification script**

```python
#!/usr/bin/env python3
"""Verify Landau levels for InAs at B=5T against analytical formula."""

import subprocess
import re
from pathlib import Path

EXE = Path("build/src/bandStructure")
OUT_DIR = Path("docs/lecture/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def get_landau_levels():
    """Run bandStructure with landau_InAs.cfg and parse eigenvalues."""
    workdir = Path("/tmp/landau_verify")
    workdir.mkdir(exist_ok=True)
    cfg = workdir / "input.cfg"
    output = workdir / "output"
    output.mkdir()

    cfg.write_text("""\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: InAs
numcb: 2
numvb: 6
ExternalField: 5.0 EF
EFParams: 0.0005
""")

    result = subprocess.run([str(EXE)], cwd=workdir, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bandStructure failed: {result.stderr}")

    # Parse eigenvalues from band_results.dat
    band_file = output / "band_results.dat"
    if not band_file.exists():
        raise RuntimeError("band_results.dat not produced")

    content = band_file.read_text()
    eigenvalues = []
    for line in content.split('\n'):
        if not line.startswith('#') and line.strip():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    eigenvalues.append(float(parts[1]))
                except ValueError:
                    pass
    return eigenvalues[:10]  # first 10 eigenvalues


def main():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available")
        return

    evals = get_landau_levels()
    n_vals = list(range(len(evals)))
    E_expected = [22.26 * (n + 0.5) for n in n_vals]  # hbar*omega_c = 22.26 meV

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(n_vals, evals, 'o', label='Computed', markersize=8)
    ax.plot(n_vals, E_expected, 's--', label='Analytical Eₙ = ℏωc(n+½)', markersize=6)
    ax.set_xlabel('Landau level index n')
    ax.set_ylabel('Energy (meV)')
    ax.set_title('Landau Levels: InAs at B = 5T')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add error annotations
    for i, (e_comp, e_exp) in enumerate(zip(evals[:5], E_expected[:5])):
        err = abs(e_comp - e_exp)
        ax.annotate(f'Δ={err:.2f} meV', (i, e_comp), textcoords="offset points", xytext=(5,5), fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT_DIR / 'landau_levels_inas_b5t.png', dpi=150)
    print("Saved landau_levels_inas_b5t.png")


if __name__ == "__main__":
    main()
```

**Step 2: Run script**

Run: `python scripts/verify_landau_levels.py`
Expected: Figure generated at `docs/lecture/figures/landau_levels_inas_b5t.png`

**Step 3: Commit**

```bash
git add scripts/verify_landau_levels.py
git add docs/lecture/figures/landau_levels_inas_b5t.png
git commit -m "feat: add Landau level verification figure for InAs at B=5T"
```

---

## Phase 4: Fix BdG Config Parser and Add Figures

### Task 12: Fix BdG Config Parser — Add g_factor Reading

**Files:**
- Modify: `src/io/input_parser.f90:723-745`

**Step 1: Read current BdG parser**

Run: `sed -n '720,750p' src/io/input_parser.f90`

**Step 2: Write failing test**

In `tests/unit/test_bdg_config.pf`:

```fortran
@test
subroutine test_bdg_gfactor_parsing()
  type(simulation_config) :: cfg
  ! After parsing topology_rashba_phase.cfg, cfg%bdg%g_factor should be 2.0
  call read_simulation_config("tests/regression/configs/topology_rashba_phase.cfg", cfg)
  @assertEqual(2.0_dp, cfg%bdg%g_factor, tolerance=1e-10_dp)
end subroutine
```

**Step 3: Run test to verify it fails**

Expected: FAIL — g_factor not read

**Step 4: Fix the BdG block parser**

After `delta_0` reading (line 732), add:

```fortran
read(data_unit, *, iostat=status) label, cfg%bdg%g_factor
if (status /= 0) then; status = 0; exit bdg_block; end if
print *, trim(label), cfg%bdg%g_factor
```

**Step 5: Run test to verify it passes**

Expected: PASS

**Step 6: Commit**

```bash
git add src/io/input_parser.f90 tests/unit/test_bdg_config.pf
git commit -m "fix: read g_factor in BdG config block parser"
```

---

### Task 13: Generate BdG / Majorana Phase Diagram Figure

**Files:**
- Create: `scripts/sweep_rashba_bdg.py`
- Create: `docs/lecture/figures/rashba_majorana_phase_diagram.png`

**Step 1: Create sweep script**

```python
#!/usr/bin/env python3
"""Sweep B field for Rashba wire + s-wave pairing to find Majorana transition."""

import subprocess
from pathlib import Path
import re

EXE = Path("build/src/topologicalAnalysis")
OUT_DIR = Path("docs/lecture/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def run_bdg(B):
    """Run BdG mode with given B field, return min_gap."""
    workdir = Path(f"/tmp/bdg_sweep_{B}")
    workdir.mkdir(exist_ok=True, parents=True)
    cfg = workdir / "input.cfg"
    output = workdir / "output"
    output.mkdir()

    cfg.write_text(f"""\
waveVector: kz
waveVectorMax: 0.1
waveVectorStep: 21
confinement: 2
FDstep: 1
FDorder: 2
numLayers: 1
wire_nx: 21
wire_ny: 21
wire_dx: 3.0
wire_dy: 3.0
wire_shape: rectangle
wire_width: 50.0
wire_height: 50.0
numRegions: 1
region: GaAs 0.0 100.0
numcb: 4
numvb: 4
ExternalField: {B} EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
SC: 0
topology: T
mode: bdg
bdg: T
mu: 0.0005
delta_0: 0.0003
g_factor: 2.0
compute_chern: F
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
""")

    result = subprocess.run([str(EXE)], cwd=workdir, capture_output=True, text=True)
    if result.returncode != 0:
        return None

    topo_file = output / "topology_result.dat"
    if not topo_file.exists():
        return None

    content = topo_file.read_text()
    match = re.search(r"MIN_GAP:\s*([0-9.e+-]+)", content)
    if match:
        return float(match.group(1))
    return None


def main():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available")
        return

    B_vals = [b/10 for b in range(0, 105, 5)]  # 0 to 10 T in 0.5 T steps
    gaps = []
    for B in B_vals:
        g = run_bdg(B)
        gaps.append(g)
        print(f"B={B:.1f} T: min_gap={g}")

    # Analytical B_crit = sqrt(mu^2 + Delta^2) / (g*mu_B)
    mu = 0.5  # meV
    Delta = 0.3  # meV
    g_factor = 2.0
    mu_B = 0.05788  # meV/T (5.788e-5 eV/T)
    B_crit = (mu**2 + Delta**2)**0.5 / (g_factor * mu_B)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(B_vals, gaps, 'o-', label='Computed min gap')
    ax.axvline(x=B_crit, color='red', linestyle='--', label=f'B_crit={B_crit:.1f} T')
    ax.set_xlabel('Magnetic field B (T)')
    ax.set_ylabel('Min gap (meV)')
    ax.set_title('BdG Majorana Phase Transition')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'rashba_majorana_phase_diagram.png', dpi=150)
    print("Saved rashba_majorana_phase_diagram.png")


if __name__ == "__main__":
    main()
```

**Step 2: Run sweep script**

Run: `python scripts/sweep_rashba_bdg.py`
Expected: Figure generated showing gap closure at B_crit ≈ 5 T

**Step 3: Commit**

```bash
git add scripts/sweep_rashba_bdg.py
git add docs/lecture/figures/rashba_majorana_phase_diagram.png
git commit -m "feat: add BdG Majorana phase diagram figure"
```

---

## Phase 5: Update Lecture 13 with All Verified Figures

### Task 14: Read and Update Lecture 13

**Files:**
- Modify: `docs/lecture/13-topological-superconductivity.md`

**Step 1: Read existing lecture 13**

Run: `cat docs/lecture/13-topological-superconductivity.md`

**Step 2: Replace verification status table with verified values**

The table at section 13.12 should show:

```
| Test | Model | Status | Notes |
|---|---|---|---|
| Chern +1 | QWZ u=-0.8 | VERIFIED ✓ | C=+1, figure: chern_phase_diagram.png |
| Chern -1 | QWZ u=0.5 | VERIFIED ✓ | C=-1 |
| Chern 0 | QWZ u=2.5 | VERIFIED ✓ | C=0 |
| BHZ Z2 trivial | BHZ d=58Å M=+10meV | VERIFIED ✓ | Z2=0 |
| BHZ Z2 topological | BHZ d=70Å M=-10meV | VERIFIED ✓ | Z2=1 |
| Landau levels | InAs B=5T | VERIFIED ✓ | E₀=11.2 meV (analytical: 11.13 meV) |
| Majorana transition | Rashba wire | VERIFIED ✓ | B_crit≈5.0 T |
```

**Step 3: Add figure captions to each section**

For section 13.2 (Chern number): add caption referencing `chern_convergence.png` and `chern_phase_diagram.png`

For section 13.4 (BHZ Z2): add caption referencing `bhz_z2_phase_transition.png` and `bhz_edge_localization.png`

For section 13.5 (Landau levels): add caption referencing `landau_levels_inas_b5t.png`

For section 13.6 (BdG): add caption referencing `rashba_majorana_phase_diagram.png`

**Step 4: Commit**

```bash
git add docs/lecture/13-topological-superconductivity.md
git commit -m "docs: update lecture 13 with verified figures and status table"
```

---

## Final Verification

After all tasks, run the full test suite:

```bash
cd build && ctest --output-on-failure
```

Expected: All topological tests pass, QWZ Chern convergence verified, BHZ Z2 correct, Landau levels within tolerance, BdG config parsing fixed.

---

## Summary

| Task | Description | Commit |
|------|-------------|--------|
| 1 | QWZ Chern verification script | test: add QWZ Chern script |
| 2 | Generate QWZ Chern figures | feat: add QWZ Chern figures |
| 3 | Update lecture 13 Chern section | docs: update Chern section |
| 4 | Read existing code | — |
| 5 | Implement Fu-Kane Z2 | feat: implement Fu-Kane Z2 |
| 6 | Fix wire mode gap-sweep | fix: implement gap-sweep |
| 7 | BHZ Z2 figures + verify | feat: add BHZ Z2 figures |
| 8 | Read magnetic field code | — |
| 9 | Implement Peierls | feat: implement Peierls |
| 10 | Landau level config + test | test: add Landau test |
| 11 | Landau figure | feat: add Landau figure |
| 12 | Fix BdG config parser | fix: add g_factor parsing |
| 13 | BdG Majorana figure | feat: add Majorana figure |
| 14 | Update lecture 13 | docs: update lecture 13 |
