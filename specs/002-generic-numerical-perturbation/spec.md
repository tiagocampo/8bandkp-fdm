# Feature Specification: Generic Numerical Perturbation Method for G-Factor Calculations

**Feature Branch**: `002-generic-numerical-perturbation`
**Created**: 2025-11-02
**Status**: Draft
**Input**: User description: "Implement generic numerical perturbation method for g-factor calculations in k·p solver using finite difference method, moving from analytical Löwdin partitioning to robust numerical differentiation approach that works for any band structure"

## Clarifications

### Session 2025-11-02

- Q: What are the specific numerical precision targets and validation tolerances for different material complexities? → A: High precision: 0.01% for simple bulk materials, 0.1% for complex quantum wells

## Requirements *(mandatory)*

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Generic G-Factor Calculation for Any Band (Priority: P1)

Research physicists need to calculate Landé g-factors for arbitrary electronic states in semiconductor heterostructures without being limited to specific bands or requiring new analytical derivations for each material system.

**Why this priority**: This is the core functionality that enables the solver to work with any band structure, making it truly generic and reusable across different research scenarios.

**Independent Test**: Can be fully tested by implementing the numerical perturbation method and calculating g-factors for multiple band types (conduction, heavy-hole, light-hole) in a single quantum well structure, demonstrating that the same code works for all bands without modification.

**Acceptance Scenarios**:

1. **Given** a quantum well structure with defined material parameters, **When** a researcher requests g-factor calculation for any electronic state, **Then** the system must compute the full g-tensor using numerical differentiation without requiring band-specific analytical formulas.

2. **Given** a target electronic state (any eigenstate index), **When** the numerical perturbation calculation runs, **Then** the system must automatically capture band mixing, anisotropy, and confinement effects present in the full Hamiltonian solution.

---

### User Story 2 - Band-Agnostic Physics Engine (Priority: P2)

Material scientists working with novel semiconductor alloys need to explore g-factor properties without being constrained by predefined band models or limited to well-studied material systems.

**Why this priority**: This extends the solver's utility to cutting-edge research where new materials and band structures are constantly being discovered.

**Independent Test**: Can be fully tested by creating a custom 6-band or 14-band k·p model and verifying that the same generic perturbation routine works correctly without code modifications.

**Acceptance Scenarios**:

1. **Given** any k·p model (4-band, 6-band, 8-band, 14-band), **When** the Hamiltonian is constructed, **Then** the generic perturbation method must calculate g-factors without requiring analytical formula derivations for each model.

2. **Given** a complex quantum well with multiple material layers, **When** g-factors are calculated, **Then** the results must automatically include effects of interface mixing, strain, and electric fields present in the structure.

---

### User Story 3 - Validation Against Known Results (Priority: P3)

Researchers need confidence that the numerical perturbation method produces accurate results that match established analytical solutions and published experimental data.

**Why this priority**: Essential for scientific credibility and adoption of the new method in the research community.

**Independent Test**: Can be fully tested by comparing numerical results with analytical formulas for simple cases (like bulk GaAs conduction band) and published experimental data for well-characterized quantum wells.

**Acceptance Scenarios**:

1. **Given** a bulk semiconductor with known analytical g-factor formula, **When** the numerical method is applied, **Then** the results must match the analytical solution within numerical precision tolerance.

2. **Given** published experimental g-factor data for standard quantum well systems, **When** calculations are performed, **Then** the numerical results must agree within experimental uncertainty margins.

---

### Edge Cases

- What happens when the target state is highly degenerate or nearly degenerate with other states?
- How does the system handle materials with extremely small band gaps where numerical instabilities may occur?
- What happens when applying very small magnetic field perturbations that approach numerical precision limits?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST implement numerical differentiation of eigenvalues with respect to magnetic field using central finite difference schemes
- **FR-002**: System MUST separate Hamiltonian construction from physics calculation to enable reuse across different band structures
- **FR-003**: System MUST automatically capture band mixing, anisotropy, and confinement effects without requiring new analytical derivations
- **FR-004**: System MUST calculate full g-tensor (gx, gy, gz components) for any selected electronic state
- **FR-005**: System MUST validate results against analytical formulas in limiting cases to ensure numerical accuracy
- **FR-006**: System MUST work with any k·p model band configuration (4-band, 6-band, 8-band, 14-band, etc.)
- **FR-007**: System MUST handle complex non-Hermitian Hamiltonians that arise from magnetic field perturbations
- **FR-008**: System MUST use appropriate step size validation for numerical derivatives to balance accuracy and stability

### Key Entities

- **Hamiltonian Builder**: Module that constructs k·p Hamiltonian matrices for any band structure and material configuration
- **Numerical Perturbation Engine**: Core algorithm that computes g-factors through numerical differentiation of eigenvalues
- **Validation Framework**: System for comparing numerical results with analytical solutions and published data
- **G-Factor Calculator**: High-level interface that accepts target state parameters and returns g-tensor results

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Researchers can calculate g-factors for any electronic state in less than 30 seconds for standard quantum well structures
- **SC-002**: Numerical g-factor results match analytical formulas within 0.01% for bulk semiconductor test cases and within 0.1% for complex quantum well structures
- **SC-003**: System supports k·p models with up to 14 bands without requiring code modifications
- **SC-004**: Generic method correctly captures band mixing effects that are missed by analytical approximations
- **SC-005**: Results agree with published experimental data within experimental uncertainty margins for standard test cases
- **SC-006**: System can handle quantum well structures with up to 10 material layers and 1000 grid points