# U8 — BdG Wire Path: Corrected μ + Window Routing + Non-Zero-Minigap Regression

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Lock in the confirmed BdG all-zero fix on the wire path (μ at the conduction subband edge + transverse B → open→close→reopen) and remove the auto-window fallback that masked it — with a regression test and a validation guard.

**Architecture:** Three changes, each TDD-covered. (1) Remove the Gershgorin auto-window retry in `run_bdg_wire` and route its window through `apply_solver_window` (ADR 0005); on 0-found, warn + return the `min_gap=-1` sentinel instead of tail states. (2) Add a `validate_semantic` guard rejecting explicitly Gershgorin-scale BdG windows. (3) A regression config + verifier asserting the open→close→reopen curve and the absence of the auto-window retry. `eval_wire_bdg_gap_app` (the gap-sweep evaluator) is deferred to U10 — it has no retry bug and is untested here.

**Tech Stack:** Fortran 2018 (`-std=f2018`), CMake/Ninja, MKL FEAST, pFUnit; Python verifiers + bash wrappers under `tests/integration/`; ctest.

**Spec:** `docs/superpowers/specs/2026-06-21-u8-bdg-window-routing-design.md`

**Branch:** `feat/bdg-u8-window-routing` (already created; spec committed at `dcdea33`).

---

## File Structure

- **Create** `tests/regression/configs/wire_inas_gaas_bdg_topological.toml` — core/shell wire, μ=0.6601 eV (conduction edge), transverse B, g=15, ±5 meV window. The confirmed-working lineage.
- **Create** `tests/integration/verify_wire_bdg_topological.py` — runs `topologicalAnalysis` at Bx∈{0,2.8,5} T (asserts open→close→reopen) and at μ=0.0 in-gap (asserts NO auto-window retry). Covers T1 + T2.
- **Create** `tests/integration/test_wire_bdg_topological.sh` — bash wrapper (delegates to the verifier).
- **Create** `tests/integration/test_topology_validate_rejects.sh` — rejection test for the BdG window guard (T3).
- **Modify** `tests/CMakeLists.txt` — register the two new ctests.
- **Modify** `src/apps/main_topology.f90` — `run_bdg_wire`: import `apply_solver_window`, route window, remove fallback.
- **Modify** `src/core/defs.f90` — `validate_semantic`: BdG solver-window guard.
- **Modify** `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` — correct the U8 μ prescription.

---

## Task 1: Create the regression config

**Files:**
- Create: `tests/regression/configs/wire_inas_gaas_bdg_topological.toml`

- [ ] **Step 1: Write the config**

This mirrors the exact config that produced the confirmed open→close→reopen (B_crit≈2.8 T) on 2026-06-21. `B_vec` is `[0,0,0]` here; the verifier overrides it to `[Bx,0,0]` per run.

```toml
# U8 regression: topological open->close->reopen on the core/shell InAs/GaAs wire.
# mu=0.6601 eV sits just inside the conduction subband edge (edge = +0.6586 eV,
# located via a DENSE-FULL normal solve). transverse B_vec=[Bx,0,0] activates
# Peierls orbital coupling + Zeeman (B along the wire gives Zeeman only).
# Confirmed: minigap 2.9 meV (B=0) -> 0.019 meV (Bx=2.8T) -> 5.5 meV (Bx=6T).
confinement = "wire"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "kz"
max = 0.0
nsteps = 2

[bands]
num_cb = 4
num_vb = 8

[wire]
nx = 13
ny = 13
dx = 5.0
dy = 5.0

[wire.geometry]
shape = "rectangle"
width = 65.0
height = 65.0

[[region]]
material = "GaAs"
inner = 25.0
outer = 65.0

[[region]]
material = "InAs"
inner = 0.0
outer = 25.0

[strain]
reference = "GaAs"
solver = "pardiso"
piezoelectric = false

[bdg]
mu = 0.6601
delta_0 = 0.0002
kz = 0.0
B_vec = [0.0, 0.0, 0.0]
g_factor = 15.0

[topology]
mode = "bdg"
compute_chern = false
compute_hall = false
qwz_u = 0.0
compute_z2 = false
extract_edge_states = false
edge_E_window = 0.01
compute_ldos = false

[solver]
method = "FEAST"
mode = "ENERGY"
emin = -0.005
emax = 0.005
m0 = 128
```

- [ ] **Step 2: Commit**

```bash
git add tests/regression/configs/wire_inas_gaas_bdg_topological.toml
git commit -m "test(config): add U8 wire BdG topological regression config"
```

---

## Task 2: T1+T2 verifier + wrapper + ctest registration (RED)

**Files:**
- Create: `tests/integration/verify_wire_bdg_topological.py`
- Create: `tests/integration/test_wire_bdg_topological.sh`
- Modify: `tests/CMakeLists.txt` (register `regression_wire_bdg_topological`)

- [ ] **Step 1: Write the verifier**

`verify_wire_bdg_topological.py` — parses stdout for `Min gap:` (stable contract printed by `run_bdg_wire`); overrides `B_vec` and `mu` per run via regex on the config text.

```python
#!/usr/bin/env python3
# COVERAGE: observable=minigap geometry=wire material=InAs ref=topological_transition
"""U8 regression: wire BdG open->close->reopen + no auto-window fallback.

Part 1 (lock-in): at mu=0.6601 eV (conduction subband edge, core/shell wire),
transverse B_vec=[Bx,0,0] for Bx in {0, 2.8, 5} T must show the topological
open->close->reopen: minigap(2.8) < 0.5*minigap(0) and minigap(5) > minigap(2.8).

Part 2 (fallback removal, U8): a mu-in-gap run (mu=0.0) must NOT use the
Gershgorin auto-window retry (run_bdg_wire no longer falls back).

Args: <topologicalAnalysis_exe> <config_file>
"""
import os
import re
import sys
import subprocess
from pathlib import Path

EXE = str(Path(sys.argv[1]).resolve())
CONFIG = Path(sys.argv[2])
OMP = "4"
RE_GAP = re.compile(r"Min gap:\s+([-\d.eE+]+)")
RETRY = "Retrying with auto-computed energy window"
WARN_NONE = "found no eigenvalues in the search window"


def run(text, tag, timeout=900):
    d = Path(f"run_{tag}")
    d.mkdir(parents=True, exist_ok=True)
    (d / "input.toml").write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS=OMP)
    p = subprocess.run([EXE], cwd=d, env=env, capture_output=True,
                       text=True, timeout=timeout)
    (d / "run.log").write_text(p.stdout + "\n--STDERR--\n" + p.stderr)
    return p


def minigap(stdout):
    m = RE_GAP.search(stdout)
    return float(m.group(1)) if m else None


def cfg_with(bx, mu=None):
    t = CONFIG.read_text()
    t = re.sub(r"B_vec\s*=\s*\[[-\d.,\s]+\]",
               f"B_vec = [{bx}, 0.0, 0.0]", t, count=1)
    if mu is not None:
        t = re.sub(r"^mu\s*=\s*[-\d.eE+]+\s*$",
                   f"mu = {mu}", t, count=1, flags=re.MULTILINE)
    return t


def main():
    # --- Part 1: open->close->reopen at mu=0.6601 ---
    gaps = {}
    for bx in (0.0, 2.8, 5.0):
        p = run(cfg_with(bx), f"bx{bx}")
        if p.returncode != 0:
            print(f"FAIL: Bx={bx} exited {p.returncode}")
            return 1
        g = minigap(p.stdout)
        if g is None or g < 0:
            print(f"FAIL: no minigap at Bx={bx} (got {g})")
            return 1
        gaps[bx] = g
    print("minigap(meV): " + "  ".join(
        f"Bx={b}:{gaps[b]*1000:.3f}" for b in sorted(gaps)))
    g0, gc, g5 = gaps[0.0], gaps[2.8], gaps[5.0]
    if not (gc < 0.5 * g0 and g5 > gc):
        print(f"FAIL: no open->close->reopen "
              f"(g0={g0*1000:.3f} gc={gc*1000:.3f} g5={g5*1000:.3f} meV)")
        return 1
    print(f"PASS: open->close->reopen "
          f"(g0={g0*1000:.3f} >= gc={gc*1000:.3f} <= g5={g5*1000:.3f} meV)")

    # --- Part 2: mu-in-gap must NOT auto-window-retry ---
    p = run(cfg_with(0.0, mu=0.0), "mu_gap")
    if RETRY in p.stdout:
        print("FAIL: auto-window fallback still used for mu-in-gap "
              "(U8 fallback removal not applied)")
        return 1
    if WARN_NONE not in p.stdout:
        print("FAIL: expected 'found no eigenvalues' warning for mu-in-gap")
        return 1
    print("PASS: mu-in-gap warns without auto-window fallback")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Write the wrapper**

`test_wire_bdg_topological.sh`:

```bash
#!/bin/bash
# U8 regression: wire BdG open->close->reopen + no auto-window fallback.
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
HERE="$(cd "$(dirname "$0")" && pwd)"
python3 "$HERE/verify_wire_bdg_topological.py" "$EXE" "$CONFIG"
```

Make it executable: `chmod +x tests/integration/test_wire_bdg_topological.sh`.

- [ ] **Step 3: Register the ctest**

In `tests/CMakeLists.txt`, append after the existing
`regression_wire_bdg_strain_shift` block (around L797):

```cmake
# U8: wire BdG all-zero fix regression. mu at the conduction subband edge
# (0.6601 eV, core/shell wire) + transverse B_vec=[Bx,0,0] gives the
# topological open->close->reopen (B_crit~2.8T). Also asserts the auto-window
# fallback is NOT used on the BdG path (U8 removes it).
add_test(
    NAME regression_wire_bdg_topological
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_bdg_topological.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_bdg_topological.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    regression_wire_bdg_topological
    PROPERTIES LABELS "regression;strain-validation" TIMEOUT 900
)
```

- [ ] **Step 4: Build**

```bash
cmake --build build
```
Expected: builds cleanly (no new code yet; only test files added).

- [ ] **Step 5: Run the new test — verify RED**

```bash
OMP_NUM_THREADS=4 ctest --test-dir build -R regression_wire_bdg_topological --output-on-failure
```
Expected: **FAIL**. Part 1 (open→close→reopen) passes on current code; Part 2 fails because the current `run_bdg_wire` still prints `Retrying with auto-computed energy window` for the μ-in-gap case. (This RED run is slow — ~5 min — because the μ-in-gap retry diagonalizes the ±78 eV tail. That is exactly the behavior U8 removes.)

- [ ] **Step 6: Commit**

```bash
git add tests/integration/verify_wire_bdg_topological.py \
        tests/integration/test_wire_bdg_topological.sh tests/CMakeLists.txt
git commit -m "test(bdg): add U8 wire topological regression (RED — auto-window fallback)"
```

---

## Task 3: run_bdg_wire — route window + remove fallback (GREEN)

**Files:**
- Modify: `src/apps/main_topology.f90` (import L13-15; window L502-510; fallback L517-528)

- [ ] **Step 1: Import `apply_solver_window`**

In `src/apps/main_topology.f90` L13-15, change the `use eigensolver` statement:

```fortran
  use eigensolver, only: eigensolver_base, make_eigensolver, eigensolver_config, &
    & eigensolver_result, eigensolver_result_free, auto_compute_energy_window, &
    & apply_solver_window, EIGEN_MODE_ENERGY, EIGEN_MODE_FULL
```

(`auto_compute_energy_window` is kept — it remains used elsewhere; removing it is out of scope.)

- [ ] **Step 2: Route the window through `apply_solver_window`**

In `run_bdg_wire`, replace the window-set block (L502-510):

```fortran
    if (cfg%solver%emin /= 0.0_dp .or. cfg%solver%emax /= 0.0_dp) then
      emin_local = cfg%solver%emin
      emax_local = cfg%solver%emax
    else
      emin_local = -max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
      emax_local =  max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
    end if
    eigen_cfg_local%emin = emin_local
    eigen_cfg_local%emax = emax_local
```

with:

```fortran
    if (cfg%solver%emin /= 0.0_dp .or. cfg%solver%emax /= 0.0_dp) then
      emin_local = cfg%solver%emin
      emax_local = cfg%solver%emax
    else
      emin_local = -max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
      emax_local =  max(50.0_dp * cfg%bdg%delta_0, 0.01_dp)
    end if
    ! Route the physics-sized window through the window authority (ADR 0005 /
    ! KTD6). apply_solver_window honors a nonzero user override verbatim, so
    ! the physics window is preserved while window selection has one home.
    ! NEVER the Gershgorin auto-window for BdG (samples the FD-Nyquist tail;
    ! see docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md).
    call apply_solver_window(H_bdg_csr, emin_local, emax_local, &
                             eigen_cfg_local%emin, eigen_cfg_local%emax)
```

- [ ] **Step 3: Remove the auto-window fallback**

In `run_bdg_wire`, replace the two `nev_found == 0` blocks (L517-528):

```fortran
    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: ', eigen_solver_local%backend_name(), &
        ' found no eigenvalues in the search window'
      print *, '  Retrying with auto-computed energy window...'
      call auto_compute_energy_window(H_bdg_csr, eigen_cfg_local%emin, eigen_cfg_local%emax)
      print *, '  Auto window: [', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'
      call eigen_solver_local%solve_sparse(H_bdg_csr, eigen_cfg_local, eigen_res_local)
    end if

    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: auto-window retry also found no eigenvalues'
      result%min_gap = -1.0_dp
```

with a single no-retry block:

```fortran
    ! U8: no auto-window fallback on the BdG path. If the physics window found
    ! nothing, mu is in the band gap or the window is mis-sized -- fail loudly
    ! with the sentinel rather than silently returning FD-Nyquist tail states.
    if (eigen_res_local%nev_found == 0) then
      print *, 'Warning: ', eigen_solver_local%backend_name(), &
        ' found no eigenvalues in the search window'
      print *, '  mu=', cfg%bdg%mu, ' eV window=[', eigen_cfg_local%emin, ',', eigen_cfg_local%emax, ']'
      print *, '  (mu likely in the band gap, or the window is mis-sized; no auto-window retry on BdG)'
      result%min_gap = -1.0_dp
```

The existing `else` branch (allocate `eigvals_bdg`, `result%min_gap = 2*bdg_zero_energy_gap(...)`, write eigenvalues, Majorana check) is unchanged and now pairs with this single `if`.

- [ ] **Step 4: Build**

```bash
cmake --build build
```
Expected: clean build.

- [ ] **Step 5: Run the U8 test — verify GREEN**

```bash
OMP_NUM_THREADS=4 ctest --test-dir build -R regression_wire_bdg_topological --output-on-failure
```
Expected: **PASS**. Both parts now pass — the curve (unchanged physics) and the no-retry check (fallback removed). The μ-in-gap run is now fast.

- [ ] **Step 6: Confirm no regression in the sibling strain test**

```bash
OMP_NUM_THREADS=4 ctest --test-dir build -R regression_wire_bdg_strain_shift --output-on-failure
```
Expected: **PASS** (its ±0.5 eV window still works; routing is behavior-preserving for a user override).

- [ ] **Step 7: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "fix(bdg): route wire BdG window through apply_solver_window; remove auto-window fallback (U8)"
```

---

## Task 4: validate_semantic guard + T3 rejection test (RED → GREEN)

**Files:**
- Modify: `src/core/defs.f90` (`validate_semantic`, BdG block ~L944-949)
- Create: `tests/integration/test_topology_validate_rejects.sh`
- Modify: `tests/CMakeLists.txt` (register `validation_rejects_bad_topology`)

> **Approval note:** `defs.f90` is flagged approval-gated in CLAUDE.md (derived-type changes). This adds a check to an existing validator — no type changes. Surface it in the PR description.

- [ ] **Step 1: Write the rejection test (RED)**

`tests/integration/test_topology_validate_rejects.sh`:

```bash
#!/bin/bash
# U8: BdG solver-window validate guard. A BdG config with a Gershgorin-scale
# [solver] window must error stop at validate_semantic.
# Args: <topologicalAnalysis_exe>
set -uo pipefail
EXE="$1"
WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT
PASS=0
FAIL=0

run_test() {
  local name="$1"
  local pat="$2"
  "$EXE" > "$WORKDIR/test_output.log" 2>&1 && RC=0 || RC=$?
  if [ "$RC" -ne 0 ] && grep -qi "$pat" "$WORKDIR/test_output.log"; then
    PASS=$((PASS + 1))
  else
    FAIL=$((FAIL + 1))
    echo "FAIL: $name — exit=$RC, pattern='$pat'"
    grep -i 'error\|STOP' "$WORKDIR/test_output.log" | head -3 | sed 's/^/  /'
  fi
}

cd "$WORKDIR"

# T3: BdG mode + Gershgorin-scale solver window -> reject
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
fd_step = 1
[bands]
num_cb = 4
num_vb = 8
[wire]
nx = 9
ny = 9
dx = 5.0
dy = 5.0
[wire.geometry]
shape = "rectangle"
width = 45.0
height = 45.0
[[region]]
material = "InAs"
inner = 0.0
outer = 45.0
[bdg]
mu = 0.66
delta_0 = 0.0002
[topology]
mode = "bdg"
[solver]
method = "FEAST"
mode = "ENERGY"
emin = -70.0
emax = 70.0
EOF
run_test "T3_bdg_wide_window" "BdG solver window"

if [ "$FAIL" -ne 0 ]; then
  echo "FAIL: topology validate rejects ($FAIL failed)"
  exit 1
fi
echo "PASS: topology validate rejects ($PASS case(s))"
exit 0
```

`chmod +x tests/integration/test_topology_validate_rejects.sh`.

- [ ] **Step 2: Register the ctest**

In `tests/CMakeLists.txt`, append near the existing `validation_rejects_bad_configs` registration (around L1254):

```cmake
# U8: BdG solver-window validate guard rejection test.
add_test(
    NAME validation_rejects_bad_topology
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_validate_rejects.sh
        $<TARGET_FILE:topologicalAnalysis>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(validation_rejects_bad_topology PROPERTIES LABELS "validation" TIMEOUT 60)
```

- [ ] **Step 3: Build + run — verify RED**

```bash
cmake --build build
OMP_NUM_THREADS=4 ctest --test-dir build -R validation_rejects_bad_topology --output-on-failure
```
Expected: **FAIL** — the wide-window BdG config is accepted today (the guard does not exist yet), so the run proceeds past validation.

- [ ] **Step 4: Implement the guard**

In `src/core/defs.f90` `validate_semantic`, inside the BdG block, immediately after the S4 confinement check (its `end if` at ~L948, before the block's closing `end if` at ~L949), insert. (Indentation is 8 spaces — the level of S3/S4, which are inside the 6-space `if (trim(cfg%topo%mode) == 'bdg') then` block.)

```fortran
        ! U8: BdG solver window, if explicitly set, must be physics-sized.
        ! Rejects Gershgorin-scale windows that sample the FD-Nyquist tail
        ! (see docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md).
        ! Bound 1 eV is generous (physics BdG windows are <=~50 meV) but rejects
        ! the ±tens-of-eV auto/Gershgorin scale; does not reject the existing
        ! strain-BdG regression config (±0.5 eV).
        if (cfg%solver%emin /= 0.0_dp .or. cfg%solver%emax /= 0.0_dp) then
          if (max(abs(cfg%solver%emin), abs(cfg%solver%emax)) > 1.0_dp) then
            error stop 'validate_semantic: BdG solver window too wide ' // &
              '(max(|emin|,|emax|) > 1 eV); use a physics-sized window around E=0'
          end if
        end if
```

- [ ] **Step 5: Build + run — verify GREEN**

```bash
cmake --build build
OMP_NUM_THREADS=4 ctest --test-dir build -R validation_rejects_bad_topology --output-on-failure
```
Expected: **PASS** — the wide-window config now `error stop`s with the "BdG solver window" message.

- [ ] **Step 6: Confirm the U8 + strain regressions still pass (guard didn't over-reject)**

```bash
OMP_NUM_THREADS=4 ctest --test-dir build -R "regression_wire_bdg_topological|regression_wire_bdg_strain_shift" --output-on-failure
```
Expected: both **PASS** (T1 config ±5 meV and strain config ±0.5 eV are both under 1 eV).

- [ ] **Step 7: Commit**

```bash
git add src/core/defs.f90 tests/integration/test_topology_validate_rejects.sh tests/CMakeLists.txt
git commit -m "feat(validate): reject Gershgorin-scale BdG solver windows (U8 guard)"
```

---

## Task 5: Correct the U8 text in the validation plan

**Files:**
- Modify: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` (U8 section, ~L299-318)

- [ ] **Step 1: Correct the μ prescription**

In the U8 section's Approach paragraph (~L309), replace the sentence containing "μ ≈ EC=+0.719 eV or EV=−0.8 eV" with:

> The headline fix is μ + window: place μ at the conduction subband edge in the solver's EV=0/EC=Eg convention — located directly via a DENSE-FULL normal solve (`bandStructure`, read `output/eigenvalues.dat`; for the core/shell InAs/GaAs wire the edge is +0.659 eV) — and use transverse `B_vec=[Bx,0,0]` so Peierls orbital coupling activates (B along the wire gives diagonal Zeeman only). Confirmed 2026-06-21: μ=0.6601 eV + g=15 + transverse B gives open→close→reopen at B_crit≈2.8 T (gap closes to 0.019 meV ≈ 0.1·δ₀); μ-shift confirms B_crit tracks μ_eff.

Also update the U8 result/verification bullets to reference `regression_wire_bdg_topological` and the removed auto-window fallback.

- [ ] **Step 2: Commit**

```bash
git add docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md
git commit -m "docs(plan): correct U8 mu prescription (conduction edge + transverse B)"
```

---

## Task 6: Full ctest suite + final verification

- [ ] **Step 1: Run the full suite (thread-capped)**

```bash
OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 --output-on-failure
```
Expected: all tests PASS (was 124/124; now +2 new = 126), including the two new U8 tests. Watch for the ctest `-j` oversubscription gotcha (cap OMP to nproc/N).

- [ ] **Step 2: Sanity-check the branch diff**

```bash
git log --oneline main..HEAD
git diff --stat main..HEAD
```
Expected: the spec commit + 5 task commits; only U8 files (no sweep-in of the pre-existing skills churn).

- [ ] **Step 3: Report**

Summarize: T1 (open→close→reopen) green, T2 (no auto-window fallback) green, T3 (guard) green, full suite green, plan text corrected. Note `defs.f90` in the PR description (approval-gated file touched for a validator check only).

---

## Out of scope (deferred)

- `eval_wire_bdg_gap_app` window routing + error-stop→sentinel → **U10** (it has no auto-window retry; tested with the gap-sweep there).
- Third KTD6 bypass `compute_spectral_function_wire` → **U9**.
- Dead `[bdg] gauge` knob + principled `g_factor` → own small task.
- Pure per-point evaluator extraction → **U2**; (B, μ) phase diagram → **U10/U11**.
