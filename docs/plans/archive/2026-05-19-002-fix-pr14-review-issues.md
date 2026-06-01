---
title: Fix PR #14 Review Issues
type: fix
status: active
date: 2026-05-19
origin: PR review of feat/kdotpy-cross-validation
---

# Fix PR #14 Review Issues Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Address all review findings from the comprehensive PR review of the cross-code validation pipeline (PR #14), covering doc comment fixes, error handling improvements, comparison logic fixes, test coverage gaps, and physics parameter verification.

**Architecture:** Four batches of fixes, ordered by risk and dependency: (1) quick edits to comments and simple logic, (2) error handling refactor across all test files, (3) new test coverage for tensile strain, (4) physics parameter research and verification. Each batch commits independently.

**Tech Stack:** Fortran 2018, Python 3.14 (validation scripts), kdotpy library

---

## File Structure

### Files Modified

| File | Responsibility | Tasks |
|------|---------------|-------|
| `src/physics/strain_solver.f90` | Fix stale Q_eps doc comment, fix ZB8bandGeneralized reference | T1 |
| `src/core/parameters.f90` | Separate HgTe/CdTe comment, verify physics parameters | T1, T9 |
| `validation/shared/comparison.py` | Fix empty-array max_delta, update docstring | T2 |
| `validation/shared/fortran_runner.py` | Warn on skipped parser lines | T3 |
| `validation/shared/kdotpy_runner.py` | Wrap kdotpy calls with context | T4 |
| `validation/run_all.py` | Fix truncation, improve timeout handling, warn on all-skip | T5 |
| `validation/gfactor/test_gfactor_qw.py` | Fix bare pass, remove spurious kdotpy check | T6 |
| `validation/strain/test_strain_qw.py` | Remove dead functions, DRY up runners | T7 |
| `validation/strain/test_strain_bandedge.py` | DRY up runners, add tensile config | T7, T8 |
| `validation/landau/test_landau_bulk.py` | Remove spurious kdotpy check | T6 |
| `validation/wire/test_wire_subbands.py` | Fix return None pattern | T10 |
| `validation/selfconsistent/test_sc_qw.py` | Fix return None pattern | T10 |
| All 12 `validation/*/test_*.py` | Replace `return None` with proper exceptions (C4) | T10 |

---

## Task 1: Fix Stale Doc Comments (C1, S4, S5)

**Files:**
- Modify: `src/physics/strain_solver.f90:723`
- Modify: `src/physics/strain_solver.f90:715`
- Modify: `src/core/parameters.f90:734-737`

- [ ] **Step 1: Fix stale Q_eps formula in strain_solver.f90 doc comment**

In `src/physics/strain_solver.f90`, line 723, change the doc comment to include the minus sign:

```
! where P_eps = -av * Tr(eps), Q_eps = -(b/2) * (ezz - 0.5*(exx+eyy)).
```

This matches the implementation at line 841 (`Q_eps = -params%b_dp * 0.5_dp * ...`) and the CLAUDE.md convention.

- [ ] **Step 2: Fix ZB8bandGeneralized reference**

In `src/physics/strain_solver.f90`, line 715, change:

```
! Hamiltonian by ZB8bandQW or ZB8bandGeneralized.
```

to:

```
! Hamiltonian by ZB8bandQW or ZB8bandBulk (hamiltonianConstructor.f90).
```

There is no subroutine named `ZB8bandGeneralized` in the codebase.

- [ ] **Step 3: Separate inverted-gap comment to HgTe only**

In `src/core/parameters.f90`, move the "Inverted-gap semimetal" line from the shared block comment (line 736) into the HgTe case block. Change lines 734-737:

```fortran
          ! II-VI materials (Pfeuffer-Jeschke PhD thesis, U. Wuerzburg, 2000;
          ! Novik et al., PRB 72, 035321, 2005; Becker et al., PRB 62, 10353, 2000).
          case ("HgTe")
```

to:

```fortran
          ! II-VI materials (Pfeuffer-Jeschke PhD thesis, U. Wuerzburg, 2000;
          ! Novik et al., PRB 72, 035321, 2005; Becker et al., PRB 62, 10353, 2000).
          case ("HgTe")
            ! Inverted-gap semimetal: Gamma_6 below Gamma_8.
```

And remove line 736 (`! Inverted-gap semimetal: Gamma_6 below Gamma_8.`) from the shared block.

- [ ] **Step 4: Verify build still compiles**

Run: `cmake --build build`
Expected: clean build, no warnings

- [ ] **Step 5: Commit**

```bash
git add src/physics/strain_solver.f90 src/core/parameters.f90
git commit -m "fix(docs): stale Q_eps comment, ZB8bandGeneralized reference, HgTe/CdTe comment"
```

---

## Task 2: Fix Empty-Array Comparison Logic (C2)

**Files:**
- Modify: `validation/shared/comparison.py:50-56`

- [ ] **Step 1: Fix `compare_eigenvalues` for empty arrays**

In `validation/shared/comparison.py`, replace lines 50-56:

```python
    if len(ours_meV) == 0 and len(theirs_meV) == 0:
        return {
            "passed": False,
            "max_delta_meV": 0.0,
            "per_band": [],
            "error": "No eigenvalues to compare (both arrays empty)",
        }
```

with:

```python
    if len(ours_meV) == 0 or len(theirs_meV) == 0:
        return {
            "passed": False,
            "max_delta_meV": float("inf"),
            "per_band": [],
            "error": (f"Empty eigenvalue array: ours={len(ours_meV)}, "
                      f"kdotpy={len(theirs_meV)}"),
        }
```

This ensures `max_delta_meV` is `inf` (not 0.0) when no comparison actually happened, preventing false confidence in downstream consumers.

- [ ] **Step 2: Update the docstring**

Change line 33:

```python
    Both arrays must have the same length.
```

to:

```python
    Arrays of different lengths are allowed; the overlapping subset is compared
    for diagnostics, but `passed` will be False if lengths differ. Empty arrays
    produce `max_delta_meV=inf`.
```

- [ ] **Step 3: Run a quick smoke test**

Run: `cd /data/8bandkp-fdm && python3 -c "from validation.shared.comparison import compare_eigenvalues; r = compare_eigenvalues([], []); assert r['max_delta_meV'] == float('inf'); assert not r['passed']; print('PASS')"`

- [ ] **Step 4: Commit**

```bash
git add validation/shared/comparison.py
git commit -m "fix(validation): return inf max_delta for empty eigenvalue arrays"
```

---

## Task 3: Warn on Skipped Parser Lines (I5)

**Files:**
- Modify: `validation/shared/fortran_runner.py:272-292`

- [ ] **Step 1: Add skipped-line tracking to `_parse_eigenvalues`**

In `validation/shared/fortran_runner.py`, replace the `_parse_eigenvalues` function (lines 272-292):

```python
def _parse_eigenvalues(filepath):
    """Parse eigenvalues.dat, returning list of (|k|, [eigenvalues]).

    Each line: |k| eval_1 eval_2 ... eval_N
    Comment lines start with '#'.
    Eigenvalues are sorted ascending within each k-point.

    Returns:
        list of (float, list[float])
    """
    results = []
    skipped = 0
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                vals = [float(x) for x in line.split()]
            except ValueError:
                skipped += 1
                continue
            if len(vals) >= 2:
                k = vals[0]
                evals = sorted(vals[1:])
                results.append((k, evals))
            else:
                skipped += 1
    if skipped > 0:
        print(f"  Warning: skipped {skipped} malformed lines in {filepath}")
    return results
```

- [ ] **Step 2: Commit**

```bash
git add validation/shared/fortran_runner.py
git commit -m "fix(validation): warn on skipped malformed lines in eigenvalue parser"
```

---

## Task 4: Wrap kdotpy Calls with Context (I8)

**Files:**
- Modify: `validation/shared/kdotpy_runner.py:24-48, 64-98`

- [ ] **Step 1: Add error wrapping to `run_bulk`**

In `validation/shared/kdotpy_runner.py`, replace the loop in `run_bulk` (lines 43-48):

```python
    results = []
    for kx, ky, kz in k_points:
        k = Vector(kx, ky, kz)
        dd = hbulk(k, pp)
        results.append(sorted(dd.eival))
    return results
```

with:

```python
    results = []
    for kx, ky, kz in k_points:
        k = Vector(kx, ky, kz)
        try:
            dd = hbulk(k, pp)
        except Exception as e:
            raise RuntimeError(
                f"kdotpy hbulk failed for material, k=({kx},{ky},{kz}): {e}"
            ) from e
        results.append(sorted(dd.eival))
    return results
```

- [ ] **Step 2: Add error wrapping to `run_qw`**

Replace the loop in `run_qw` (lines 93-98):

```python
    results = []
    for kx, ky in k_points:
        k = Vector(kx, ky, 0.0)
        dd = hz(k, pp, neig=neig, energy=energy)
        results.append(sorted(dd.eival))
    return results
```

with:

```python
    results = []
    for kx, ky in k_points:
        k = Vector(kx, ky, 0.0)
        try:
            dd = hz(k, pp, neig=neig, energy=energy)
        except Exception as e:
            raise RuntimeError(
                f"kdotpy hz failed for QW, k=({kx},{ky}): {e}"
            ) from e
        results.append(sorted(dd.eival))
    return results
```

- [ ] **Step 3: Commit**

```bash
git add validation/shared/kdotpy_runner.py
git commit -m "fix(validation): wrap kdotpy calls with contextual error messages"
```

---

## Task 5: Fix run_all.py Truncation and Timeout Handling (I7)

**Files:**
- Modify: `validation/run_all.py:56-70, 111, 136`

- [ ] **Step 1: Preserve line boundaries in truncation**

In `validation/run_all.py`, replace lines 68-69:

```python
        "output": result.stdout[-2000:] if result.stdout else "",
        "error": result.stderr[-500:] if result.stderr else "",
```

with:

```python
        "output": '\n'.join(result.stdout.split('\n')[-50:]) if result.stdout else "",
        "error": '\n'.join(result.stderr.split('\n')[-20:]) if result.stderr else "",
```

This preserves line boundaries instead of cutting mid-line.

- [ ] **Step 2: Capture partial output on timeout**

Replace lines 56-58:

```python
    except subprocess.TimeoutExpired:
        return {"name": name, "script": script, "status": "TIMEOUT",
                "output": f"Exceeded {timeout}s timeout", "error": ""}
```

with:

```python
    except subprocess.TimeoutExpired as e:
        partial = ""
        if hasattr(e, 'stdout') and e.stdout:
            lines = e.stdout.decode('utf-8', errors='replace').split('\n')
            partial = '\n'.join(lines[-10:])
        return {"name": name, "script": script, "status": "TIMEOUT",
                "output": f"Exceeded {timeout}s timeout\n{partial}", "error": ""}
```

- [ ] **Step 3: Add all-skip warning**

After line 111 (`print(f"SUMMARY: {n_pass} passed, {n_fail} failed, {n_skip} skipped")`), add:

```python
    if n_pass == 0 and n_skip > 0:
        print("\nWARNING: No tests passed. All tests were skipped.")
        print("  Check: source validation/kdotpy_env/bin/activate")
```

- [ ] **Step 4: Commit**

```bash
git add validation/run_all.py
git commit -m "fix(validation): preserve line boundaries in output, capture timeout partials"
```

---

## Task 6: Fix parse_gfactor and Remove Spurious kdotpy Checks (I6, S3)

**Files:**
- Modify: `validation/gfactor/test_gfactor_qw.py:67-87, 97-103`
- Modify: `validation/landau/test_landau_bulk.py:96-102`

- [ ] **Step 1: Fix bare `except ValueError: pass` in parse_gfactor**

In `validation/gfactor/test_gfactor_qw.py`, replace lines 80-86:

```python
        elif current and line:
            try:
                vals = [float(x) for x in line.split()]
                if vals:
                    result[current] = vals
            except ValueError:
                pass
```

with:

```python
        elif current and line:
            try:
                vals = [float(x) for x in line.split()]
                if vals:
                    result[current] = vals
            except ValueError:
                unparseable.append(line)
```

Add `unparseable = []` before the loop (after `current = None`) and add after the loop, before the return:

```python
    if unparseable and not result:
        raise RuntimeError(
            f"Could not parse any g-factor values. "
            f"Unparseable lines: {unparseable[:5]}"
        )
```

- [ ] **Step 2: Remove spurious kdotpy availability check in test_gfactor_qw.py**

Replace lines 97-103:

```python
    # Check kdotpy availability
    try:
        from kdotpy.config import initialize_config
        initialize_config()
    except ImportError:
        print("SKIP: kdotpy not available (activate kdotpy_env first)")
        return False
```

with:

```python
    # This test validates against analytical references, not kdotpy.
    # No kdotpy availability check needed.
```

- [ ] **Step 3: Remove spurious kdotpy availability check in test_landau_bulk.py**

Replace lines 96-102:

```python
    # Check kdotpy availability (used for analytical comparison)
    try:
        from kdotpy.config import initialize_config
        initialize_config()
    except ImportError:
        print("SKIP: kdotpy not available (activate kdotpy_env first)")
        return False
```

with:

```python
    # This test validates against analytical formulas, not kdotpy.
    # No kdotpy availability check needed.
```

- [ ] **Step 4: Commit**

```bash
git add validation/gfactor/test_gfactor_qw.py validation/landau/test_landau_bulk.py
git commit -m "fix(validation): improve parse_gfactor error handling, remove spurious kdotpy checks"
```

---

## Task 7: Remove Dead Code and DRY Up Strain Runners (I10, S6)

**Files:**
- Modify: `validation/strain/test_strain_qw.py:180-214` (remove dead kdotpy functions)
- Modify: `validation/strain/test_strain_qw.py:41-165` (DRY up duplicated runners)
- Modify: `validation/strain/test_strain_bandedge.py:39-135` (DRY up duplicated runners)

- [ ] **Step 1: Remove dead kdotpy runner functions from test_strain_qw.py**

Delete lines 180-214 in `validation/strain/test_strain_qw.py` (the `run_kdotpy_qw_strained` and `run_kdotpy_qw_unstrained` functions). These are defined but never called.

- [ ] **Step 2: Update docstring in test_strain_qw.py**

Change lines 4-6:

```python
Also runs kdotpy as a secondary check, but notes that kdotpy's QW strain
model differs due to how substrate strain interacts with band offsets in
the 3-layer model.
```

to:

```python
kdotpy QW strain comparison is deferred — its substrate model interacts with
band offsets differently than our code.
```

- [ ] **Step 3: DRY up test_strain_qw.py — extract shared runner**

Replace the two duplicated functions (`run_fortran_qw_strained` and `run_fortran_qw_unstrained`, lines 41-165) with a single parameterized function:

```python
def _run_fortran_qw(barrier, well, l_well_ang, l_barrier_ang, total_ang,
                    fdstep=201, substrate=None, work_dir=None, timeout=120):
    """Run Fortran bandStructure for a QW, optionally with strain."""
    exe = os.path.join(BUILD_DIR, "src", "bandStructure")
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")

    cleanup = work_dir is None
    if work_dir is None:
        prefix = "strain_qw_" if substrate else "unstr_qw_"
        work_dir = tempfile.mkdtemp(prefix=prefix)

    try:
        half_total = total_ang / 2.0
        well_start = -l_well_ang / 2.0
        well_end = l_well_ang / 2.0

        lines = [
            "waveVector: k0",
            "waveVectorMax: 0",
            "waveVectorStep: 1",
            "confinement: 1",
            f"FDstep: {fdstep}",
            "FDorder: 2",
            "numLayers: 2",
            f"material1: {barrier} {-half_total:.1f} {half_total:.1f} 0",
            f"material2: {well} {well_start:.1f} {well_end:.1f} 0",
            "numcb: 6",
            "numvb: 12",
            "ExternalField: 0  EF",
            "EFParams: 0.0",
        ]
        if substrate:
            lines.append("strain: T")
            lines.append(f"strain_ref: {substrate}")

        with open(os.path.join(work_dir, "input.cfg"), 'w') as f:
            f.write('\n'.join(lines) + '\n')
        os.makedirs(os.path.join(work_dir, "output"), exist_ok=True)

        result = subprocess.run(
            [exe], cwd=work_dir, capture_output=True, text=True, timeout=timeout
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"bandStructure failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

        eig_path = os.path.join(work_dir, "output", "eigenvalues.dat")
        if not os.path.exists(eig_path):
            raise RuntimeError(f"No eigenvalues.dat produced in {work_dir}")

        rows = []
        with open(eig_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                rows.append([float(x) for x in line.split()])
        return rows[0][1:] if rows else None
    finally:
        if cleanup:
            shutil.rmtree(work_dir, ignore_errors=True)
```

Then update the call sites in `test_strain_qw()` (around lines 258-265) from:

```python
f_strained = run_fortran_qw_strained(...)
f_unstrained = run_fortran_qw_unstrained(...)
```

to:

```python
f_strained = _run_fortran_qw(..., substrate=substrate)
f_unstrained = _run_fortran_qw(..., substrate=None)
```

And update the None check (line 272) from `if f_strained is None or f_unstrained is None:` to handle `FileNotFoundError` and `RuntimeError` from the new function.

- [ ] **Step 4: DRY up test_strain_bandedge.py — extract shared runner**

Replace the two duplicated functions (`run_fortran_strained` and `run_fortran_unstrained`, lines 39-135) with:

```python
def _run_fortran_bulk(material, build_dir, substrate_a0=None):
    """Run Fortran bandStructure for bulk, optionally with strain substrate."""
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")

    workdir = tempfile.mkdtemp(prefix="strain_" if substrate_a0 else "unstr_")
    try:
        lines = [
            "waveVector: k0",
            "waveVectorMax: 0",
            "waveVectorStep: 1",
            "confinement:  0",
            "FDstep: 1",
            "FDorder: 2",
            "numLayers:  1",
            f"material1: {material}",
            "numcb: 4",
            "numvb: 4",
            "ExternalField: 0  EF",
            "EFParams: 0.0",
        ]
        if substrate_a0 is not None:
            lines.append(f"strainSubstrate: {substrate_a0}")

        with open(os.path.join(workdir, "input.cfg"), 'w') as f:
            f.write('\n'.join(lines) + '\n')
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)

        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"bandStructure failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

        eig_path = os.path.join(workdir, "output", "eigenvalues.dat")
        if not os.path.exists(eig_path):
            raise RuntimeError(f"No eigenvalues.dat produced in {workdir}")

        rows = []
        with open(eig_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                rows.append([float(x) for x in line.split()])
        return rows[0][1:] if rows else None
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
```

Update call sites similarly (lines 184-185):

```python
f_unstrained = _run_fortran_bulk(material, BUILD_DIR)
f_strained = _run_fortran_bulk(material, BUILD_DIR, substrate_a0=substrate_a0)
```

- [ ] **Step 5: Verify strain tests still work**

Run: `source validation/kdotpy_env/bin/activate && python3 validation/strain/test_strain_bandedge.py && python3 validation/strain/test_strain_qw.py`

- [ ] **Step 6: Commit**

```bash
git add validation/strain/test_strain_qw.py validation/strain/test_strain_bandedge.py
git commit -m "refactor(validation): DRY up duplicated strain runners, remove dead kdotpy functions"
```

---

## Task 8: Add Tensile Strain Test Coverage (C3)

**Files:**
- Modify: `validation/strain/test_strain_bandedge.py:29-36`

- [ ] **Step 1: Add tensile strain configuration to STRAIN_CONFIGS**

In `validation/strain/test_strain_bandedge.py`, add a second config to the `STRAIN_CONFIGS` list (after line 35):

```python
STRAIN_CONFIGS = [
    {
        "name": "InAs on GaAs substrate (compressive)",
        "material": "InAs",
        "substrate": "GaAs",
        "substrate_a0": 5.65325,
    },
    {
        "name": "GaAs on InP substrate (tensile)",
        "material": "GaAs",
        "substrate": "InP",
        "substrate_a0": 5.8687,  # InP lattice constant in Angstrom
    },
]
```

InP lattice constant 5.8687 A from Vurgaftman 2001. GaAs on InP gives tensile strain (a_sub > a_well), which reverses the HH/LH splitting direction compared to compressive strain.

Note: InP must already be in `FORTRAN_MATERIALS` in `param_mapper.py` for the kdotpy comparison. If not, add it there too.

- [ ] **Step 2: Verify InP is in param_mapper**

Check `validation/shared/param_mapper.py` for an InP entry. If missing, add it with parameters from `parameters.f90`.

- [ ] **Step 3: Run the strain test**

Run: `source validation/kdotpy_env/bin/activate && python3 validation/strain/test_strain_bandedge.py`

Expected: both configs pass (compressive and tensile).

- [ ] **Step 4: Commit**

```bash
git add validation/strain/test_strain_bandedge.py
git commit -m "test(validation): add tensile strain config (GaAs/InP) for Q_eps sign coverage"
```

---

## Task 9: Physics Parameter Verification (C5, I1, I2)

**Files:**
- Modify: `src/core/parameters.f90` (if verification confirms values need changing)
- Modify: `validation/shared/param_mapper.py` (if convention adjustments needed)

This task requires reference verification before making changes. Do NOT modify parameters without confirming against published references (CLAUDE.md boundary rule).

- [ ] **Step 1: Research HgTe deltaSO value**

Check kdotpy's material database for HgTe deltaSO:

Run: `python3 -c "source validation/kdotpy_env/bin/activate; from kdotpy.materials import allMaterials; m = allMaterials.get('HgTe'); print(m.param.get('delta_so') if m else 'NOT FOUND')"`

Compare against our value of 1.003 eV. The Novik 2005 table commonly uses 1.08 eV. If kdotpy uses 1.08 and our value differs, either:
(a) Update our value to match Novik (requires CLAUDE.md boundary approval), or
(b) Document the discrepancy and adjust the validation tolerance for HgTe

- [ ] **Step 2: Research CdTe gamma convention**

Check whether kdotpy uses bare or modified Luttinger parameters for CdTe:

Run: `python3 -c "source validation/kdotpy_env/bin/activate; from kdotpy.materials import allMaterials; m = allMaterials.get('CdTe'); print('gamma1:', m.param.get('gamma1')); print('gamma2:', m.param.get('gamma2')); print('gamma3:', m.param.get('gamma3'))"`

If kdotpy uses modified gammas (typically smaller values like gamma1~1.47 for CdTe), our bare gammas (gamma1=5.0) need conversion in `param_mapper.py` — NOT a change to `parameters.f90` (our code uses bare gammas correctly).

- [ ] **Step 3: Verify CdTe elastic constants**

Check literature values for CdTe: C11, C12, C44 from Landolt-Bornstein III/41B. Current code uses HgTe values (532, 368, 201) for both materials. CdTe literature values are typically C11~537, C12~374, C44~202 (close but not identical).

If the values are genuinely different, update CdTe in `parameters.f90` lines 773-775 with the correct reference.

- [ ] **Step 4: Apply fixes or document discrepancies**

Based on research results:
- If values need updating in `parameters.f90`: make the change with the reference citation
- If values are correct but param_mapper needs adjustment: update `param_mapper.py`
- If values are close enough (< 1% impact on validation): document the discrepancy in a comment

- [ ] **Step 5: Re-run validation if parameters changed**

Run: `source validation/kdotpy_env/bin/activate && python3 validation/bulk/test_bulk_k0.py`

- [ ] **Step 6: Commit**

```bash
git add src/core/parameters.f90 validation/shared/param_mapper.py
git commit -m "fix(physics): verify HgTe/CdTe parameters against Novik 2005 and kdotpy"
```

---

## Task 10: Replace `return None` with Proper Exceptions (C4)

This is the most widespread change — all test files that use `return None` for both "executable not found" and "no output produced" need to distinguish these failure modes.

**Files:**
- Modify: `validation/gfactor/test_gfactor_qw.py:47`
- Modify: `validation/wire/test_wire_subbands.py:40-41, 57-58`
- Modify: `validation/landau/test_landau_bulk.py:49-50, 66-67`
- Modify: `validation/selfconsistent/test_sc_qw.py:39-40, 57-58`
- Modify: `validation/gfactor/test_gfactor_qw.py:67-87` (run_gfactor already done in T6)
- Modify: Callers in all test files that check `if result is None:`

- [ ] **Step 1: Update test_gfactor_qw.py `run_gfactor`**

In `validation/gfactor/test_gfactor_qw.py`, replace lines 47-48:

```python
    if not os.path.isfile(exe) or not os.path.isfile(config_path):
        return None
```

with:

```python
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")
```

Update the caller (around line 124) to handle `FileNotFoundError`:

```python
        except (RuntimeError, FileNotFoundError) as e:
            print(f"  FAIL: {e}")
            all_pass = False
```

- [ ] **Step 2: Update test_wire_subbands.py `run_fortran_wire`**

In `validation/wire/test_wire_subbands.py`, replace lines 40-41:

```python
    if not os.path.isfile(exe) or not os.path.isfile(config_path):
        return None
```

with:

```python
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")
```

Replace lines 57-58:

```python
        if not os.path.exists(eig_path):
            return None
```

with:

```python
        if not os.path.exists(eig_path):
            raise RuntimeError(f"No eigenvalues.dat produced in {workdir}")
```

Update the caller (around line 98) to catch `FileNotFoundError`:

```python
        except (RuntimeError, FileNotFoundError) as e:
            print(f"  FAIL: {e}")
            all_pass = False
```

And remove the `if rows is None:` check (replaced by exceptions).

- [ ] **Step 3: Update test_landau_bulk.py `run_fortran_landau`**

Same pattern as Step 2. Replace lines 49-50 and 66-67 with `raise FileNotFoundError` and `raise RuntimeError`.

- [ ] **Step 4: Update test_sc_qw.py `run_fortran_sc`**

Same pattern. Replace lines 39-40 and 57-58 with `raise FileNotFoundError` and `raise RuntimeError`.

- [ ] **Step 5: Verify all tests still pass**

Run: `source validation/kdotpy_env/bin/activate && python3 validation/run_all.py`

Expected: all 12 tests pass (or skip due to kdotpy env).

- [ ] **Step 6: Commit**

```bash
git add validation/gfactor/test_gfactor_qw.py validation/wire/test_wire_subbands.py \
        validation/landau/test_landau_bulk.py validation/selfconsistent/test_sc_qw.py
git commit -m "fix(validation): distinguish exe-not-found from no-output with proper exceptions"
```

---

## Task 11: Tighten SC Tolerance and Add KeyError to Exception Handlers (I9, review finding)

**Files:**
- Modify: `validation/selfconsistent/test_sc_qw.py:23, 30`
- Modify: `validation/bulk/test_bulk_k0.py:95` (add KeyError)

- [ ] **Step 1: Tighten SC QW tolerance**

In `validation/selfconsistent/test_sc_qw.py`, change line 23:

```python
TOL_CB1_MEV = 10.0
```

to:

```python
TOL_CB1_MEV = 5.0
```

Also update the docstring (line 9) and print statement (line 80).

- [ ] **Step 2: Add KeyError to exception handler in test_bulk_k0.py**

In `validation/bulk/test_bulk_k0.py`, add `KeyError` to the except tuple at line 95:

```python
except (ImportError, FileNotFoundError, RuntimeError, ValueError,
        OSError, KeyError) as e:
```

This prevents a typo in a material name from crashing the entire test.

- [ ] **Step 3: Verify SC test still passes with tighter tolerance**

Run: `python3 validation/selfconsistent/test_sc_qw.py`

If it fails, check the actual CB1 value and adjust `expected_cb1_meV` if needed.

- [ ] **Step 4: Commit**

```bash
git add validation/selfconsistent/test_sc_qw.py validation/bulk/test_bulk_k0.py
git commit -m "fix(validation): tighten SC tolerance, add KeyError to bulk exception handler"
```

---

## Task 12: Update CLAUDE.md Cross-Validation Section

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update test count in CLAUDE.md**

In `CLAUDE.md`, the cross-validation section mentions "12 cross-code validation tests". Verify this is still accurate after all changes. The tests are still 12, but 3 of them (wire, gfactor, landau) are now documented as single-code. Update the description if needed.

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update cross-validation description after review fixes"
```

---

## Summary

| Task | Issue(s) | Risk | Est. Time |
|------|----------|------|-----------|
| T1 | C1, S4, S5 | Low (doc changes) | 5 min |
| T2 | C2 | Low (shared module) | 5 min |
| T3 | I5 | Low (shared module) | 5 min |
| T4 | I8 | Low (shared module) | 5 min |
| T5 | I7 | Low (pipeline runner) | 5 min |
| T6 | I6, S3 | Low (two test files) | 5 min |
| T7 | I10, S6 | Medium (refactor) | 15 min |
| T8 | C3 | Medium (new test) | 10 min |
| T9 | C5, I1, I2 | High (physics params) | 15 min |
| T10 | C4 | Medium (all test files) | 15 min |
| T11 | I9, review | Low | 5 min |
| T12 | docs | Low | 5 min |

**Total: ~95 min**

### Deferred to Follow-Up

- S1: Unit tests for `param_mapper.py` error paths
- S2: Smoke test for `regenerate_references.py`
- I3: Automated cross-check of `param_mapper.py` against `parameters.f90`
- I4: Relabel single-code tests (cosmetic, does not affect functionality)
