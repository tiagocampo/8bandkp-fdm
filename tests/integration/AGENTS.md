# Integration Tests

Owns all full-executable integration tests: regression/smoke (shell), verification ladder (Python), standard-star benchmarks (S1–S7), convergence tests (U4–U8), and coverage tooling. Does NOT own: unit tests (`tests/unit/`), regression golden data (`tests/regression/data/`), cross-code validation (`validation/`).

## Entry Points

- **ctest labels**: `ctest --test-dir build -L regression`, `-L verification`, `-L standard-star`, `-L strain-validation`, `-L coverage`
- **Direct Python**: `python3 tests/integration/verify_8band_rung1_bulk_k0.py`
- **Lecture scripts** call back into `star_helpers.py` here: `sys.path.insert(0, str(REPO / "tests" / "integration"))`

## Test Categories

### A. Shell Regression Tests (`test_*.sh`) — 32 files
Run Fortran executables, compare output against golden reference data in `tests/regression/data/`. Pattern: `mktemp → copy config → run exe → grep/awk compare → exit 0/1`.

Subdivided by physics domain: bulk (2), QW (7), SC (5), wire (7), g-factor (2), topology (7), validation (1), scattering (1).

### B. 8-Band Verification Ladder (Rungs 1–8) — 8 files
`verify_8band_rung{1..8}_*.py` + `verify_strain_rung{5,6}_*.py`. Systematic structural validation from bulk k=0 through wire, strain, g-factor, and optics.

### C. Standard-Star Benchmarks (S1–S7) — 7 files
`verify_star_*.py`. Publication-ready benchmark comparisons against literature values (Vurgaftman, Kane, Roth, Winkler). Each produces a markdown benchmark table row.

### D. Convergence Tests (U4–U8) — 6 files
`test_*_convergence.py`. Multi-resolution runs with Richardson extrapolation, GCI uncertainty quantification, and convergence rate assertions.

### E. Inline Verification Scripts — 11 files
`verify_bulk_benchmarks.py`, `verify_qw_benchmarks.py`, `verify_sc_benchmark.py`, etc. Called by `test_*.sh` for quantitative checks.

### F. Coverage Tooling — 2 files
`coverage_matrix.py` (scans COVERAGE annotations → markdown report), `aggregate_star_benchmarks.py` (runs S1–S7, concatenates benchmark tables).

## Shared Infrastructure

### `star_helpers.py` (430 lines) — Core shared library
- **Executable runners**: `run_exe(build_dir, name, config, workdir)`, `run_executable(exe_path, ...)`
- **Parsers**: `parse_eigenvalues()`, `parse_gfactor()`, `parse_absorption()`, `parse_topology_result()`
- **Physics**: `roth_gfactor()`, `bir_pikus_biaxial_001()`, `extract_effective_mass()`
- **Validation**: `compare_value(actual, expected, tol, name, unit)` → (passed, delta, row_dict)
- **Constants**: `HBAR_EV_S`, `HBAR2_OVER_2M0 = 3.80998`, `E_CHARGE`, `M0_KG`
- **Tolerances**: `TOL_EXACT=1e-12`, `TOL_ANALYTICAL=0.01`, `TOL_NUMERICAL=0.05`

### `convergence_helpers.py` (416 lines) — Convergence analysis
- `richardson_extrapolate()`, `richardson_extrapolate_3grid()`, `compute_gci()`
- `extract_convergence_rates()`, `check_monotonic()`, `make_convergence_report()`

### `coverage_matrix.py` (266 lines) — Coverage tracking
- Reads `validation_universe.yml`, scans `verify_*.py` and `test_*.sh` for `# COVERAGE:` annotations
- Regex: `r"^#\s*COVERAGE:\s*(.+)$"`

## COVERAGE Annotation Format

```
# COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001
```

Fields: `observable` (required), `geometry` (required), `material`, `ref`, `tier`. Place annotations in `verify_*.py` files (not in `.sh` wrappers). New tests should include them so the coverage matrix tracks physics coverage.

## Tolerance Tiers

| Tier | Value | Use case |
|------|-------|----------|
| `TOL_EXACT` | 1e-12 | Machine-precision, FD-order independence |
| `TOL_ANALYTICAL` | 0.01 (1%) | Analytical formula comparison (Kane, Roth, Bir-Pikus) |
| `TOL_NUMERICAL` | 0.05 (5%) | FD-numerical results, strained systems |
| Golden regression | 1e-8 | Bit-identical output comparison |
| SC eigenvalues | 1e-5 | Iterative convergence paths |

## Patterns

### Shell test boilerplate
```bash
set -euo pipefail
EXE="$1"; CONFIG="$2"
WORKDIR=$(mktemp -d); trap "rm -rf $WORKDIR" EXIT
/bin/cp "$CONFIG" "$WORKDIR/input.toml"
mkdir -p "$WORKDIR/output"; cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
```

### Python test pattern
```python
rc, output_dir = run_exe(BUILD_DIR, "bandStructure", config_path, tmpdir, timeout=120)
eigs = parse_eigenvalues(output_dir / "eigenvalues.dat")
passed, delta, _ = compare_value(actual, expected, TOL_ANALYTICAL, "name", "eV")
```

### Lecture–test correspondence
Each `scripts/lecture_NN_*.py` has matching integration tests in this directory covering the same physics domain. Lectures are pedagogical; tests are automated regression.

## Anti-patterns

- Never create configs at runtime in shell tests — use `tests/regression/configs/`
- Never skip `# COVERAGE:` annotations in new `verify_*.py` files
- Never hardcode executable paths — use `$1` (shell) or `run_exe()` (Python)
- Never put COVERAGE annotations in `.sh` wrappers — put them in the verifier Python script

## Related Context

- Lecture scripts: `scripts/AGENTS.md`
- Regression configs: `tests/regression/configs/`
- Golden reference data: `tests/regression/data/`
- Cross-code validation: `validation/`
- Coverage universe: `tests/integration/validation_universe.yml`
