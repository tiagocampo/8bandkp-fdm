---
title: "Standard-star benchmark suite: three logic errors in aggregation, timeout handling, and layer config"
module: tests/integration
date: 2026-05-09
problem_type: logic_error
component: testing_framework
severity: medium
tags:
  - standard-star
  - benchmarks
  - markdown
  - subprocess
  - timeout
  - config
  - layer-overlap
  - confinement-init
symptoms:
  - "aggregate_star_benchmarks.py produces malformed markdown with colliding pipe delimiters when wrapping inner benchmark tables in a 2-column outer table"
  - "run_executable() in star_helpers.py crashes all 7 verify_star scripts with raw Python subprocess.TimeoutExpired traceback instead of clean FAIL message"
  - "qw_gaas_algaas.cfg produces flat AlGaAs profile instead of GaAs QW — layer 3 overwrites well due to last-layer-wins semantics with overlapping z-ranges"
root_cause: logic_error
resolution_type: code_fix
related_components:
  - confinement_init.f90
  - star_helpers.py
  - aggregate_star_benchmarks.py
---

# Standard-star benchmark suite: three logic errors in aggregation, timeout handling, and layer config

## Problem

Three independent logic errors in the standard-star benchmark test suite (S1-S7) produced malformed markdown output, crashed test scripts on executable hangs, and silently computed incorrect quantum-well band structures due to overlapping material layers. All found via code review (session history), all fixed in commit `efdb0c5`.

## Symptoms

- **Aggregation output illegible**: `aggregate_star_benchmarks.py` wrapped inner benchmark table rows (already pipe-delimited) inside a 2-column outer table, producing colliding `|` delimiters that no markdown renderer could parse.
- **Test scripts crash on timeout**: If a Fortran executable hung, `subprocess.TimeoutExpired` propagated uncaught from `run_executable()`, crashing the calling `verify_star_*.py` script with a raw Python traceback instead of a clean FAIL.
- **Silent physics error**: The GaAs quantum well in `qw_gaas_algaas.cfg` was never actually present in the simulation. Layer 3 (`Al30Ga70As -200 200 0`) covered the entire domain and overwrote the GaAs well set by layer 2 in `[-50, 50]`.

## What Didn't Work

All three bugs were directly identified during code review. The overlapping-layer bug was non-obvious — the 3-layer config *looks* correct at a glance (barrier, well, barrier), but the z-ranges reveal the overlap. Session history shows the same pattern was discovered 5 days earlier in `generate_all_figures.py` (fig_wavefunctions_qw) but not propagated to the config files.

## Solution

### Bug 1: Broken markdown table in aggregate_star_benchmarks.py

Before (wraps inner rows in outer table — pipe collision):

```python
print("| Star | Benchmark Table |")
print("|------|----------------|")
for star_id, rows in all_rows:
    for row in rows:
        print(f"| {star_id} | {row} |")
```

After (per-star sections, no nesting):

```python
for star_id, rows in all_rows:
    print(f"## {star_id}\n")
    for row in rows:
        print(row)
    print()
```

### Bug 2: Unhandled subprocess timeout in star_helpers.py

Before (exception propagates to caller):

```python
result = subprocess.run(
    [exe_path], cwd=work_dir, capture_output=True,
    text=True, timeout=timeout,
)
return result.returncode, output_dir
```

After (catch and return structured failure):

```python
try:
    result = subprocess.run(
        [exe_path], cwd=work_dir, capture_output=True,
        text=True, timeout=timeout,
    )
    return result.returncode, output_dir
except subprocess.TimeoutExpired:
    return -1, output_dir
```

Callers already treat nonzero returncode (including -1) as failure — no downstream changes needed.

### Bug 3: Overlapping material layers in qw_gaas_algaas.cfg

Before (3 layers — layer 3 overwrites well):

```
numLayers:  3
material1: Al30Ga70As -200 200 0
material2: GaAs -50 50 0
material3: Al30Ga70As -200 200 0
```

After (2 layers — barrier covers full domain, well overwrites center):

```
numLayers:  2
material1: Al30Ga70As -200 200 0
material2: GaAs -50 50 0
```

## Why This Works

**Bug 1**: Inner benchmark rows already contain pipe-delimited columns (`| Material | Observable | ... |`). Wrapping them as a cell value in an outer table creates colliding `|` delimiters. Per-star `##` headings with raw rows avoids nesting entirely.

**Bug 2**: `subprocess.run(..., timeout=300)` raises `TimeoutExpired` rather than returning a result with nonzero returncode. The `try/except` catches this and returns `(-1, output_dir)`, which callers handle through existing `if rc != 0` paths.

**Bug 3**: `confinement_init.f90` line 117 uses last-layer-wins assignment: "later layers overwrite earlier ones." With 3 layers, layer 3 (`Al30Ga70As` spanning `-200` to `200`) overwrites the GaAs well that layer 2 set in `[-50, 50]`. With 2 layers, the barrier covers the full domain first, then the well overwrites just the central region — matching the documented 2-layer pattern in the code comment.

## Prevention

1. **Never nest pre-formatted markdown tables inside other tables.** When aggregating tabular output, use section headings (`##`) and emit raw content. A CI check validating markdown output from test scripts would catch this.

2. **Every `subprocess.run` with `timeout` must handle `TimeoutExpired`.** The `run_executable` wrapper now catches it centrally. Audit any direct `subprocess.run` calls to ensure they go through this wrapper.

3. **Validate QW configs for overlapping layers.** A config with N layers where layer N spans the full domain after a narrower layer N-1 is almost certainly wrong. Four additional SC configs have this same pattern and need fixing:
   - `sc_qcse_gaas_algaas_ef.cfg`
   - `sc_mod_doped_gaas_algaas.cfg`
   - `sc_gaas_alas_qw_ef.cfg`
   - `sc_qcse_gaas_algaas.cfg`

   The last-layer-wins semantics are documented in `confinement_init.f90:117-119` and in `docs/lecture/07-self-consistent-sp.md:617`, but the pitfall keeps recurring (found twice in 5 days across two independent code paths).

## Related Issues

- Plan: `docs/plans/2026-05-09-003-fix-standard-star-review-findings-plan.md` — the fix plan that addressed 24 review findings
- Parent framework: `docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md` — verification ladder extended by standard-star suite
- Sibling: `docs/solutions/test-failures/csr-test-infra-stop1-oob-tautology-fixes-2026-05-09.md` — parallel test-infra fixes in same work stream
