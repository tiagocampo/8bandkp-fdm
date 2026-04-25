# Docs Physics Revamp Design

## Goal

Revamp the lecture suite into `research-grade internal documentation` backed by
traceable simulation output, reproducible figure generation, explicit external
benchmarks, and stronger regression coverage for the underlying physics code.

This is not a prose-only rewrite. The docs become the visible surface of a
code-validation effort:

1. inspect the current solver and post-processing paths,
2. reproduce as many relevant nextnano and literature examples as practical,
3. fix code paths that fail those reproductions,
4. regenerate figures from validated configs and outputs,
5. rewrite the chapters so every strong claim is matched to evidence.

The revamp must be especially strict for the sparse `quantum wire` path, where
the current documentation is the least mature and the hidden numerical risk is
highest.

## Product Definition

The target artifact is the existing lecture tree in `docs/lecture/00-12`, kept
mostly intact but upgraded to a higher evidence standard. Moderate restructuring
is allowed where accuracy demands it, but the default is to preserve the
current chapter boundaries and improve them in place.

The governing documentation rule is:

- If a figure or claim cannot be validated cleanly, remove it or downgrade it.

This means:

- no quantitative caption without a traceable config and output path,
- no validation table that quotes values the current code does not reproduce,
- no schematic figure presented as if it were computed output,
- no strong wire physics claims without sparse-path validation.

## Scope

### In scope

- `docs/lecture/*.md`
- `docs/figures/*.png`
- `scripts/plotting/generate_all_figures.py`
- benchmark and regression configs under `tests/regression/configs/`
- integration and unit tests under `tests/integration/` and `tests/unit/`
- core physics and output paths in:
  - `src/physics/hamiltonianConstructor.f90`
  - `src/physics/strain_solver.f90`
  - `src/physics/optical_spectra.f90`
  - `src/physics/scattering.f90`
  - `src/physics/exciton.f90`
  - `src/physics/sc_loop.f90`
  - `src/physics/gfactor_functions.f90`
  - `src/io/outputFunctions.f90`
  - `src/math/geometry.f90`
  - `src/math/sparse_matrices.f90`
  - `src/math/eigensolver.f90`

### Out of scope

- public-facing packaging and styling work,
- speculative physics features not needed for benchmark reproduction,
- broad architecture rewrites not justified by a benchmark failure,
- preserving weak figures only for coverage.

## Quality Bar

A chapter is considered complete only when it passes all four gates:

### 1. Reproducibility

Every quantitative figure maps to:

- a known config in `tests/regression/configs/`,
- a known executable path,
- a known output file set,
- a deterministic generator or chapter-local parsing path.

### 2. Physics Defensibility

Each chapter needs at least one external anchor where practical:

- nextnano tutorial reproduction,
- analytic limit,
- or literature comparison.

Agreement does not need to be exact everywhere, but the core physical
signatures must match and any tolerance must be explained.

### 3. Solver Confidence

If a figure was wrong because of a code or parsing defect, the revamp must add a
test or regression asset that would fail if the defect returned.

### 4. Documentation Discipline

Each chapter must clearly distinguish:

- validated quantitative results,
- qualitative interpretation,
- provisional or future-work material.

## External Benchmark Strategy

The nextnano++ tutorial suite is the primary benchmark framework because it
already covers the same physics families as the lecture set:

- bulk k.p dispersion,
- quantum well dispersion,
- strain in QWs,
- optical absorption and polarization,
- excitons,
- scattering times,
- QCSE,
- quantum wires,
- magnetic effects,
- numerics.

The revamp should build a benchmark matrix that maps each lecture chapter to one
or more external anchors:

| Chapter | Internal target | External anchor |
|---|---|---|
| `01` | bulk GaAs, strained bulk GaAs | nextnano bulk GaAs tutorial |
| `02` | GaAs/AlGaAs and broken-gap QW dispersion | nextnano QW dispersion tutorials |
| `03` | wavefunction composition and state indexing | internal output consistency + QW tutorials |
| `04` | pseudomorphic strain and HH/LH splitting | nextnano strain/QW references |
| `05` | conduction-band g-factor and Zeeman splitting | internal perturbation theory + magnetic-effect references |
| `06` | TE/TM absorption, excitons, scattering | nextnano absorption/gain/exciton tutorials + Bastard/Ferreira |
| `07` | SP convergence and band bending | Snider/Tan + nextnano-style SP benchmark |
| `08` | wire confinement, strain, anisotropy | nextnano wire references + internal sparse-path checks |
| `09` | FD and solver scaling | internal convergence/timing only |
| `10` | QCSE field response | nextnano QCSE tutorials + internal Stark checks |
| `11` | convergence | internal manufactured/reference studies |
| `12` | extension guidance | codebase reality, not external physics |

## Workstreams

### Workstream A: Benchmark Matrix and Failure Ledger

Create a living control document that records:

- chapter-to-benchmark mapping,
- current figure defects,
- missing outputs,
- known solver suspicions,
- acceptance tolerances,
- completion state by chapter.

This is the source of truth for what remains untrusted.

### Workstream B: Figure Provenance and Post-Processing Hardening

Repair the generation layer so figures do not lie about:

- units,
- loaded state indices,
- polarization labels,
- field values,
- parsed columns,
- exciton/sharp-peak semantics.

This workstream fixes known falsehoods quickly, but it must not hide solver
defects behind prettier plotting.

### Workstream C: Solver Validation and Repair

For each benchmark miss, determine whether the problem lives in:

- Hamiltonian assembly,
- strain terms,
- optical matrix elements,
- scattering implementation,
- electric-field handling,
- sparse wire assembly,
- eigensolver stability,
- or output formatting/parsing.

Code changes must be benchmark-driven, not style-driven.

### Workstream D: Quantum Wire Hardening

Treat wire mode as a dedicated subproject. Before keeping strong wire claims,
validate:

- geometry/material masks,
- sparse assembly and indexing,
- 2D strain field mapping,
- eigensolver behavior,
- optical and g-factor post-processing on sparse outputs,
- consistency between generated figures and wire outputs.

### Workstream E: Chapter Rebuild

After a benchmark path is trusted, rewrite the corresponding chapter to include:

- the validated example and config,
- the relevant outputs,
- the regenerated figures,
- the benchmark summary,
- and explicit limitations where the code still falls short.

## Milestones

### Milestone 0: Audit Freeze

Capture the current bad state without trying to hide it. Record all known figure
and prose failures before touching the code.

### Milestone 1: Known-False Figure Cleanup

Fix obviously wrong unit conversions, state-index mistakes, missing figure
embeds, and broken parser assumptions.

### Milestone 2: Core Benchmark Reproduction

Reproduce the strongest external anchors for bulk, QW, strain, optics, SP, and
QCSE. Use misses to drive code repair.

### Milestone 3: Quantum Wire Hardening

Benchmark and repair the sparse wire path. The wire chapter is rebuilt only
after this milestone, not before.

### Milestone 4: Lecture Rebuild

Rewrite the chapters around validated outputs and remove weak material.

### Milestone 5: Regressionization

Encode the repaired behavior into unit, integration, and regression checks so
the docs do not rot again.

## Priority Order

The recommended priority sequence is:

1. Control documents and provenance inventory.
2. Fix known-false figures:
   - `gfactor_zeeman.png`
   - `qcse_stark_shift.png`
   - `qw_wavefunctions.png`
   - `perband_density.png`
   - Chapter 6 absorption/exciton/scattering figures
3. Reproduce bulk, QW, strain, SP, and QCSE external anchors.
4. Harden the wire code and wire chapter.
5. Rebuild all lecture chapters to the new evidence standard.

## Success Criteria

The revamp succeeds when:

- every surviving quantitative figure is traceable and regenerated,
- all strong claims in `docs/lecture/` are benchmark-backed or explicitly
  qualified,
- the code reproduces the maximum practical subset of relevant nextnano examples,
- wire-mode claims are limited to what the sparse path can defend,
- the main physics regressions are encoded in tests.

## Deliverables

- a benchmark matrix and failure ledger under `docs/plans/`,
- corrected configs under `tests/regression/configs/`,
- repaired tests under `tests/integration/` and `tests/unit/`,
- regenerated figures under `docs/figures/`,
- updated lecture chapters under `docs/lecture/`,
- a final validation summary documenting what was reproduced and what remains
  out of scope.
