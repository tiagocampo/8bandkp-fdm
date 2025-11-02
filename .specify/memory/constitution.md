<!--
Sync Impact Report:
Version change: 1.0.0 → 1.1.0 (MINOR: expanded scientific methodology guidance)
Modified principles:
  - Scientific Accuracy (expanded to include numerical perturbation methods)
  - Validation-First Development (enhanced with generic g-factor calculation approach)
Added sections:
  - Generic Numerical Perturbation Methods (new subsection under Scientific Standards)
Removed sections: N/A
Templates requiring updates:
  ✅ plan-template.md (Constitution Check section updated with numerical perturbation gates)
  ✅ spec-template.md (no changes needed)
  ✅ tasks-template.md (no changes needed)
Follow-up TODOs: None
-->

# 8bandkp-fdm Constitution

## Core Principles

### I. Scientific Accuracy (NON-NEGOTIABLE)
Every calculation MUST be traceable to peer-reviewed literature and validated against known results. Material parameters MUST be sourced from authoritative references (Vurgaftman et al., 2001). Numerical methods MUST follow established finite difference formulations with documented convergence properties. All physical constants and approximations MUST be explicitly documented with citations. Theoretical approaches MUST distinguish between analytical Löwdin partitioning approximations and exact numerical diagonalization methods, with clear documentation of method choice and limitations.

### II. Reproducible Results
All calculations MUST produce identical results across different platforms and compiler versions. Input files MUST be self-contained with complete parameter specifications. Input selection MUST be explicit and deterministic: executables MUST accept a configuration filename argument (with a documented default of `input.cfg` if omitted). Output artifacts MUST be written under a dedicated `outputs/` directory, isolated per run (e.g., `outputs/<timestamp-or-run-id>/`), and include sufficient metadata for verification. Build process MUST be deterministic and documented in Makefile with clear dependency management.

### III. Modular Architecture
Core physics modules (Hamiltonian construction, finite differences, g-factor calculations) MUST be independently testable and documented. Each module MUST have a single, well-defined responsibility. Interfaces between modules MUST be minimal and stable. New features MUST be added as separate modules without modifying existing core functionality. Generic numerical perturbation implementations MUST separate Hamiltonian construction from physics calculation to enable reuse across different band structures.

### IV. Validation-First Development

Every new feature MUST include validation against analytical solutions or published benchmarks before implementation. Test cases MUST cover both bulk and quantum well configurations. Performance benchmarks MUST be established and maintained. Regression testing MUST catch numerical accuracy degradation. G-factor calculations MUST validate against both analytical formulas (where available) and published numerical results, demonstrating that generic numerical methods capture band mixing, anisotropy, and confinement effects automatically.

### V. Documentation Standards

All scientific methods MUST be documented with mathematical formulations and references. Code MUST include inline comments explaining physical significance, not just implementation details. User documentation MUST include example input files and expected output ranges. API documentation MUST specify input/output units and physical meaning. Implementation choices between analytical and numerical approaches MUST be justified with clear trade-off analysis.

## Scientific Standards

### Generic Numerical Perturbation Methods

For g-factor calculations and similar perturbation theory applications, numerical methods MUST be preferred over analytical approximations for generality and accuracy. Generic implementations MUST compute physical properties as numerical derivatives of eigenvalues with respect to perturbation parameters (e.g., magnetic field). Hamiltonian construction MUST be separated from physics calculation to enable reuse across different band structures. Numerical perturbation MUST automatically capture band mixing, anisotropy, and confinement effects without requiring new analytical derivations. Central finite difference schemes MUST be used for derivative calculations with appropriate step size validation. Results MUST be validated against analytical formulas in limiting cases to ensure numerical accuracy.

### Material Parameter Management

Material parameters MUST be centralized in parameters.f90 with clear source citations. New materials MUST include complete parameter sets (band gaps, effective masses, Kane parameters, etc.). Parameter validation MUST check physical reasonableness (positive masses, appropriate band gaps). Deprecated parameters MUST be clearly marked with migration guidance.

### Numerical Method Requirements

Finite difference schemes MUST use consistent discretization across all spatial derivatives. Convergence studies MUST be performed for new geometries or parameter ranges. Matrix conditioning MUST be monitored and reported for ill-conditioned systems. Memory usage MUST be optimized for large quantum well calculations.

## Development Workflow

### Code Review Process

All Fortran code changes MUST be reviewed by domain expert (Tiago de Campos) for scientific accuracy. Mathematical formulations MUST be verified against literature references. Performance impact MUST be assessed for large-scale calculations. Documentation updates MUST accompany code changes.

### Testing Requirements

Unit tests MUST validate individual module functionality with known analytical solutions. Integration tests MUST verify end-to-end calculation accuracy. Regression tests MUST catch numerical precision changes. Build tests MUST verify compilation across different Fortran compilers and library versions.

### Release Standards

Each release MUST include complete input/output examples for bulk and quantum well cases. Performance benchmarks MUST be documented with timing and memory usage. Scientific accuracy MUST be validated against published results. Documentation MUST be updated to reflect any API changes.

## Governance

Constitution supersedes all other development practices. Amendments require scientific justification, impact assessment, and approval by project maintainer. All contributions must verify compliance with scientific accuracy and reproducibility standards. Use agents.md for runtime development guidance and collaboration protocols.

**Version**: 1.1.0 | **Ratified**: 2025-01-27 | **Last Amended**: 2025-11-02