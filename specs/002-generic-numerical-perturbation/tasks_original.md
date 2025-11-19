---

description: "Task list for Generic Numerical Perturbation Method implementation"
---

# Tasks: Generic Numerical Perturbation Method

**Input**: Design documents from `/specs/002-generic-numerical-perturbation/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/

**Tests**: Fortran unit tests with validation against analytical solutions and published experimental data

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Path Conventions

- **Fortran project**: `src/`, `tests/` at repository root
- Paths follow the modular Fortran architecture from plan.md

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and basic structure

- [X] T001 Update Makefile to include new numerical perturbation modules and ARPACK-ng dependency
- [X] T002 Create test directories structure per implementation plan
- [X] T003 [P] Create placeholder files for new modules (numerical_perturbation.f90, perturbation_solver.f90, validation_framework.f90)
- [X] T004 [P] Add .gitkeep files to test output directories

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core infrastructure that MUST be complete before ANY user story can be implemented

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [X] T005 Extend defs.f90 with new physical constants and precision parameters for numerical perturbation
- [X] T006 [P] Add ARPACK-ng interface wrapper in src/mkl_spblas.f90 for sparse eigenvalue problems
- [X] T007 [P] Create workspace management utilities in src/utils.f90 for memory optimization
- [X] T008 Implement error handling framework for numerical perturbation operations
- [X] T009 Setup base configuration types for perturbation parameters

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Generic G-Factor Calculation for Any Band (Priority: P1) 🎯 MVP

**Goal**: Implement core numerical perturbation engine that works for any band structure without analytical formulas

**Independent Test**: Calculate g-factors for multiple band types (conduction, heavy-hole, light-hole) in a single quantum well structure using the same code, demonstrating band-agnostic operation

### Tests for User Story 1

- [ ] T010 [P] [US1] Create unit test framework in tests/unit/test_numerical_perturbation.f90
- [ ] T011 [P] [US1] Create validation test cases in tests/unit/test_validation.f90 for bulk semiconductors
- [ ] T012 [P] [US1] Create integration tests in tests/unit/test_integration.f90 for quantum well structures

### Implementation for User Story 1

- [ ] T013 [US1] Implement PerturbationConfig derived type in src/numerical_perturbation.f90
- [ ] T014 [P] [US1] Implement GFactorTensor derived type in src/numerical_perturbation.f90 (depends on T013)
- [ ] T015 [P] [US1] Implement HamiltonianPerturbation derived type in src/numerical_perturbation.f90 (depends on T013)
- [ ] T016 [P] [US1] Implement ConvergenceResult derived type in src/numerical_perturbation.f90 (depends on T013)
- [ ] T017 [US1] Implement central finite difference algorithm in src/numerical_perturbation.f90
- [ ] T018 [US1] Implement adaptive step size validation function in src/numerical_perturbation.f90
- [ ] T019 [US1] Implement high-level calculate_gfactor_tensor interface in src/numerical_perturbation.f90
- [ ] T020 [US1] Integrate numerical perturbation with existing hamiltonianConstructor.f90
- [ ] T021 [US1] Update main_gfactor.f90 to use numerical perturbation method when specified
- [ ] T022 [US1] Add input parameter parsing for numerical perturbation configuration
- [ ] T023 [US1] Implement band mixing detection and quantification algorithms in src/numerical_perturbation.f90 to validate automatic capture requirement from FR-003
- [ ] T024 [US1] Implement degenerate state handling and block-diagonalization for energy differences < 1e-6 eV in src/numerical_perturbation.f90
- [ ] T025 [US1] Add numerical precision adjustment for small band gap materials (<0.01 eV) in src/numerical_perturbation.f90
- [ ] T026 [US1] Implement higher-order finite difference schemes for small magnetic field perturbations (δB < 1e-8 T) in src/numerical_perturbation.f90

**Checkpoint**: At this point, User Story 1 should be fully functional and testable independently

---

## Phase 4: User Story 2 - Band-Agnostic Physics Engine (Priority: P2)

**Goal**: Extend system to work with any k·p model (4-band, 6-band, 8-band, 14-band) without code modifications

**Independent Test**: Create custom 6-band or 14-band k·p model and verify the same generic perturbation routine works correctly without modifications

### Implementation for User Story 2

- [ ] T027 [P] [US2] Create generic band configuration interface in src/perturbation_solver.f90
- [ ] T028 [P] [US2] Implement adaptive solver selection in src/perturbation_solver.f90 (dense vs sparse)
- [ ] T029 [P] [US2] Implement non-Hermitian eigenvalue solver for magnetic field perturbations in src/perturbation_solver.f90
- [ ] T030 [P] [US2] Add automatic Hermitian property detection in src/perturbation_solver.f90
- [ ] T031 [P] [US2] Implement band-agnostic matrix operations in src/perturbation_solver.f90
- [ ] T032 [P] [US2] Extend material parameter system for arbitrary band configurations in src/parameters.f90
- [ ] T033 [US2] Add support for complex material interfaces and mixing effects in src/parameters.f90 with documented interface models and material parameter interpolation following Vurgaftman et al. (2001) constitutional requirements
- [ ] T034 [US2] Implement automatic dimension detection and memory management

**Checkpoint**: At this point, User Stories 1 AND 2 should both work independently

---

## Phase 5: User Story 3 - Validation Against Known Results (Priority: P3)

**Goal**: Build comprehensive validation framework to ensure numerical results match analytical solutions and published data

**Independent Test**: Compare numerical results with analytical formulas for bulk GaAs conduction band and published experimental data for standard quantum wells

### Implementation for User Story 3

- [ ] T035 [P] [US3] Create ValidationFramework module in src/validation_framework.f90
- [ ] T036 [P] [US3] Implement analytical solution database for bulk semiconductors in src/validation_framework.f90
- [ ] T037 [P] [US3] Implement experimental data comparison module in src/validation_framework.f90
- [ ] T038 [P] [US3] Create validation test suite runner in src/validation_framework.f90
- [ ] T039 [P] [US3] Implement precision testing and tolerance verification in src/validation_framework.f90
- [ ] T040 [P] [US3] Add validation report generation in src/validation_framework.f90
- [ ] T041 [P] [US3] Create benchmark test cases in tests/examples/validation_comparison.example
- [ ] T042 [P] [US3] Implement convergence analysis tools in src/validation_framework.f90

**Checkpoint**: All user stories should now be independently functional

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Improvements that affect multiple user stories

- [ ] T043 [P] Update documentation in docs/ with numerical perturbation method details
- [ ] T044 [P] Create example input files in examples/ for numerical perturbation usage
- [ ] T045 [P] Extend plotting scripts in scripts/ to handle g-factor tensor visualization
- [ ] T046 Performance optimization across all perturbation calculations
- [ ] T047 Add comprehensive error handling and user-friendly error messages
- [ ] T048 Update quickstart.md with numerical perturbation examples
- [ ] T049 Run validation test suite and ensure all tests pass
- [ ] T050 Performance benchmarking to meet <30 second calculation target

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3-5)**: All depend on Foundational phase completion
  - User stories can proceed in parallel (if staffed) or sequentially in priority order (P1 → P2 → P3)
- **Polish (Phase 6)**: Depends on all desired user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Can start after Foundational (Phase 2) - No dependencies on other stories
- **User Story 2 (P2)**: Can start after Foundational (Phase 2) - Builds on US1 but independently testable
- **User Story 3 (P3)**: Can start after Foundational (Phase 2) - Validates both US1 and US2 results

### Within Each User Story

- Unit tests MUST be written before implementation (TDD approach for scientific code)
- Configuration types before algorithms
- Core algorithms before high-level interfaces
- Integration with existing codebase after core implementation
- Story complete before moving to next priority

### Parallel Opportunities

- All Setup tasks marked [P] can run in parallel
- All Foundational tasks marked [P] can run in parallel (within Phase 2)
- Once Foundational phase completes, all user stories can start in parallel (if team capacity allows)
- All derived types within a story marked [P] can run in parallel
- Different user stories can be worked on in parallel by different team members

---

## Parallel Example: User Story 1

```bash
# Launch all unit tests for User Story 1 together:
Task: "Create unit test framework in tests/unit/test_numerical_perturbation.f90"
Task: "Create validation test cases in tests/unit/test_validation.f90"
Task: "Create integration tests in tests/unit/test_integration.f90"

# Launch all derived types for User Story 1 together:
Task: "Implement PerturbationConfig derived type in src/numerical_perturbation.f90"
Task: "Implement GFactorTensor derived type in src/numerical_perturbation.f90"
Task: "Implement HamiltonianPerturbation derived type in src/numerical_perturbation.f90"
Task: "Implement ConvergenceResult derived type in src/numerical_perturbation.f90"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup
2. Complete Phase 2: Foundational (CRITICAL - blocks all stories)
3. Complete Phase 3: User Story 1
4. **STOP and VALIDATE**: Test User Story 1 independently with bulk and quantum well examples
5. Verify numerical results match analytical solutions within 0.01% tolerance

### Incremental Delivery

1. Complete Setup + Foundational → Foundation ready
2. Add User Story 1 → Test independently → Validate against analytical solutions (MVP!)
3. Add User Story 2 → Test independently → Validate with custom band models
4. Add User Story 3 → Test independently → Validate against published experimental data
5. Each story adds scientific capability without breaking previous functionality

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together
2. Once Foundational is done:
   - Developer A: User Story 1 (core numerical perturbation)
   - Developer B: User Story 2 (band-agnostic extensions)
   - Developer C: User Story 3 (validation framework)
3. Stories complete and integrate independently
4. Final integration and performance optimization

---

## Success Criteria Validation

### User Story 1 Success Metrics
- Calculate g-factors for conduction, heavy-hole, light-hole states in same quantum well
- Results match analytical formulas within 0.01% for bulk test cases
- Calculation time <30 seconds for standard quantum well structures
- Memory usage <1GB for typical problems

### User Story 2 Success Metrics
- Support 4-band, 6-band, 8-band, 14-band models without code changes
- Automatic handling of non-Hermitian Hamiltonians
- Correct capture of band mixing and interface effects
- Scalable to 1000 grid points and 10 material layers

### User Story 3 Success Metrics
- Validation test suite with 100% pass rate
- Comparison with published experimental data within uncertainty margins
- Convergence testing with automatic step size optimization
- Performance benchmarks meeting all success criteria

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story should be independently completable and testable
- Scientific accuracy validation is critical for all numerical methods
- Performance targets (<30 seconds, <1GB memory) must be verified
- Follow Fortran best practices and existing code style conventions
- Maintain backward compatibility with existing input/output formats