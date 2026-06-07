# Lecture Scripts

Owns 15 pedagogical lecture-companion scripts (L00–L14) that serve dual roles as integration tests and educational demonstrations. Also owns figure generation (`plotting/`) and the legacy config converter. Does NOT own: integration test shell scripts (`tests/integration/`), regression configs (`tests/regression/configs/`), unit tests (`tests/unit/`).

## Entry Points

- **Individual**: `python3 scripts/lecture_01_bulk.py`
- **All figures**: `python3 scripts/plotting/generate_all_figures.py` (90+ publication figures)
- **Config conversion**: `python3 scripts/convert_cfg_to_toml.py <file.cfg>`

## Lecture Scripts (L00–L14)

| Script | Executable | What it validates |
|--------|-----------|-------------------|
| `lecture_00_quickstart.py` | bandStructure | Demo only — k=0 spectrum, E(k) dispersion |
| `lecture_01_bulk.py` | bandStructure | GaAs eigenvalues, CB m* (Kane), InAs Eg |
| `lecture_02_qw.py` | bandStructure | Subband count, anticrossing splitting |
| `lecture_03_wavefunctions.py` | bandStructure | Band character >90%, normalization |
| `lecture_04_strain.py` | bandStructure | Bir-Pikus HH-LH splitting, CB hydrostatic shift |
| `lecture_05_gfactor.py` | gfactorCalculation | Roth g-factor, Landau levels, wire anisotropy |
| `lecture_06_optical.py` | opticalProperties | Absorption onset, TE/TM, ISBT, gain, spin-resolved |
| `lecture_07_scsp.py` | bandStructure | SC convergence, DIIS vs linear mixing |
| `lecture_08_wire.py` | bandStructure | CSR assembly, dense-sparse consistency |
| `lecture_09_numerical.py` | bandStructure | FD-order independence (2–8), Richardson extrapolation |
| `lecture_10_qcse.py` | bandStructure | QCSE Stark red-shift at E=-700 kV/cm |
| `lecture_11_convergence.py` | bandStructure | Grid convergence sweep, Richardson, log-log plot |
| `lecture_12_extending.py` | bandStructure | Architecture walkthrough — no assertions |
| `lecture_13_topological.py` | bandStructure, gfactorCalculation | Chern, Z2, BdG Majorana, Landau, spectral, Hall |
| `lecture_14_excitons_scattering.py` | bandStructure | Exciton binding energy, LO-phonon scattering rates |

## Shared Infrastructure

All lecture scripts import from `tests/integration/star_helpers.py`:
```python
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import run_exe, parse_eigenvalues, compare_value, ...
```

**Provided by star_helpers**: `run_exe()`, `parse_eigenvalues()`, `parse_gfactor()`, `parse_absorption()`, `compare_value()`, `roth_gfactor()`, `bir_pikus_biaxial_001()`, `extract_effective_mass()`, tolerance constants.

## Patterns

### Standard boilerplate (every lecture script)
```python
REPO = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO / "build"
CONFIGS_DIR = REPO / "tests" / "regression" / "configs"
FIGURES_DIR = REPO / "docs" / "lecture" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
```

### Tempdir isolation
Each script uses `tempfile.TemporaryDirectory()`. Config copied to `<tmpdir>/input.toml`, executable runs there, results read from `<tmpdir>/output/`.

### Validation pattern
```python
passed, delta, _ = compare_value(actual, expected, TOL_*, "name", "eV")
if not passed: failures.append(...)
sys.exit(1 if failures else 0)
```

### Figure output
- All figures: `matplotlib.use("Agg")` (headless)
- Saved to `docs/lecture/figures/lecture_XX_*.png` at 300 dpi
- `fig.tight_layout()` always called before save

### Config sourcing
Most scripts reference pre-existing TOML configs from `tests/regression/configs/`. Exceptions: L07 (SC) and L11 (convergence) generate configs inline with template strings.

## Plotting Subdirectory

| File | Lines | Purpose |
|------|------:|---------|
| `plotting/generate_all_figures.py` | ~6800 | Master figure generator: 90+ publication figures. Independent parser infrastructure. CLI: `--skip-build`, `--only FIG ...` |
| `plotting/generate_benchmark_figures.py` | ~180 | 4 benchmark figures from pre-existing data (no Fortran execution) |

## Other Scripts

| File | Lines | Purpose |
|------|------:|---------|
| `generate_all_figures.py` | ~140 | Legacy orchestrator: runs 4 topology verification scripts, generates 4 figures |
| `convert_cfg_to_toml.py` | ~800 | Converts legacy `.cfg` to TOML. Single file or `--batch` directory mode |

## Anti-patterns

- Never hardcode executable paths — use `BUILD_DIR / "src" / "name"` via `run_exe()`
- Never generate TOML configs inline except for convergence sweeps (L07, L11 pattern)
- Never skip `sys.exit(1)` on failures — CI depends on exit codes
- Never use `plt.show()` — always `matplotlib.use("Agg")` and save to file

## Related Context

- Shared infrastructure: `tests/integration/star_helpers.py` and `tests/integration/convergence_helpers.py`
- Integration tests: `tests/integration/AGENTS.md`
- Regression configs: `tests/regression/configs/`
- Output figures: `docs/lecture/figures/` and `docs/figures/`
