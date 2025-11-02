# Tasks: Build Verification and Plotting

**Input**: Design documents from `/specs/001-build-verification-plotting/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), data-model.md, contracts/

**Tests**: Manual validation tasks included for result verification. No automated test framework required for this scientific computing project.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Path Conventions

- **Single project**: `src/`, `scripts/`, `docs/`, `examples/` at repository root
- Paths shown below assume single Fortran project structure

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and basic structure

- [x] T001 Create project structure per implementation plan
- [x] T002 [P] Create scripts directory for plotting scripts
- [x] T003 [P] Create docs directory for documentation
- [x] T004 [P] Create examples directory for input files
- [x] T005 [P] Verify existing Makefile configuration in Makefile
- [x] T006 [P] Verify existing source files in src/ directory

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core infrastructure that MUST be complete before ANY user story can be implemented

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [x] T007 Verify Fortran compiler availability and version
- [x] T008 Verify BLAS/LAPACK library availability and linking
- [x] T009 Verify FFTW3 library availability and linking
- [x] T010 Verify Make build system functionality
- [ ] T011 [P] Create comprehensive build documentation in docs/BUILD.md
- [ ] T012 [P] Document dependency requirements and installation instructions
- [ ] T013 [P] Create example input files with documented parameters in examples/

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Build Process Verification (Priority: P1) 🎯 MVP

**Goal**: Verify that the 8bandkp-fdm codebase builds correctly and produces reliable executables

**Independent Test**: Run `make all` and verify that both `bandStructure` and `gfactorCalculation` executables are created successfully, then run basic calculations with provided example input files

### Implementation for User Story 1

- [x] T014 [US1] Test build process with make all command
- [x] T015 [US1] Verify bandStructure executable creation and functionality
- [x] T016 [US1] Verify gfactorCalculation executable creation and functionality
- [x] T017 [US1] Test bandStructure with bulk.example input file
- [x] T018 [US1] Test bandStructure with quantumwell.example input file
- [x] T019 [US1] Test gfactorCalculation with gfactor.example input file
- [x] T020 [US1] Verify output files are generated in correct format
- [ ] T021 [US1] Document build process and troubleshooting in docs/BUILD.md
- [ ] T022 [US1] Add error handling and clear error messages for common build issues
- [ ] T023 [US1] Create build verification script for automated testing

**Checkpoint**: At this point, User Story 1 should be fully functional and testable independently

---

## Phase 4: User Story 2 - Result Verification Framework (Priority: P2)

**Goal**: Verify that bulk and quantum well calculations produce physically accurate results that match expected scientific behavior

**Independent Test**: Run calculations with known material parameters and compare results against analytical solutions or published benchmarks for simple cases

### Implementation for User Story 2

- [x] T024 [US2] Create validation test case for bulk InAs60Sb40 calculation
- [x] T025 [US2] Create validation test case for GaSb/InAs/AlSb quantum well calculation
- [x] T026 [US2] Run bulk calculation and verify band structure results match expected energy levels
- [x] T027 [US2] Run quantum well calculation and verify confinement effects produce quantized levels
- [x] T028 [US2] Run g-factor calculation and verify values fall within physically reasonable ranges
- [x] T029 [US2] Document expected results and validation criteria in docs/VALIDATION.md
- [x] T030 [US2] Create result validation script for automated verification
- [x] T031 [US2] Add numerical accuracy checks and convergence validation
- [x] T032 [US2] Document material parameter sources and references
- [x] T033 [US2] Create benchmark results database for future regression testing

**Checkpoint**: At this point, User Stories 1 AND 2 should both work independently

---

## Phase 5: User Story 3 - Standardized Plotting Scripts (Priority: P3)

**Goal**: Create reusable gnuplot scripts to visualize calculation results for verification and analysis purposes

**Independent Test**: Run calculations to generate data files, then use gnuplot scripts to create publication-quality plots that clearly show the expected physical behavior

### Implementation for User Story 3

- [x] T034 [P] [US3] Create band structure plotting script in scripts/plot_band_structure.gp
- [x] T035 [P] [US3] Create quantum well plotting script in scripts/plot_quantum_well.gp
- [x] T036 [P] [US3] Create g-factor plotting script in scripts/plot_gfactor.gp
- [x] T037 [US3] Test band structure plotting script with generated data
- [x] T038 [US3] Test quantum well plotting script with generated data
- [x] T039 [US3] Test g-factor plotting script with generated data
- [x] T040 [US3] Verify plots are publication-quality with proper labels and units
- [x] T041 [US3] Document plotting script usage and parameters in docs/PLOTTING.md
- [x] T042 [US3] Create plotting automation script for batch processing
- [x] T043 [US3] Add plot customization options and styling templates

**Checkpoint**: All user stories should now be independently functional

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Improvements that affect multiple user stories

- [ ] T044 [P] Update README.md with comprehensive usage instructions
- [ ] T045 [P] Create quickstart guide in docs/QUICKSTART.md
- [ ] T046 [P] Add performance benchmarking and timing measurements
- [ ] T047 [P] Create comprehensive error handling and troubleshooting guide
- [ ] T048 [P] Add cross-platform compatibility testing and documentation
- [ ] T049 [P] Create release notes and version documentation
- [ ] T050 [P] Add code comments explaining physical significance in source files
- [ ] T051 [P] Create user manual with complete examples and use cases
- [ ] T052 [P] Add automated build and validation pipeline
- [ ] T053 [P] Create contribution guidelines and development workflow documentation

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3+)**: All depend on Foundational phase completion
  - User stories can then proceed in parallel (if staffed)
  - Or sequentially in priority order (P1 → P2 → P3)
- **Polish (Final Phase)**: Depends on all desired user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Can start after Foundational (Phase 2) - No dependencies on other stories
- **User Story 2 (P2)**: Can start after Foundational (Phase 2) - Depends on US1 for working executables
- **User Story 3 (P3)**: Can start after Foundational (Phase 2) - Depends on US1 and US2 for data generation

### Within Each User Story

- Build verification before calculation testing
- Calculation testing before result validation
- Result validation before plotting script creation
- Core implementation before documentation
- Story complete before moving to next priority

### Parallel Opportunities

- All Setup tasks marked [P] can run in parallel
- All Foundational tasks marked [P] can run in parallel (within Phase 2)
- Once Foundational phase completes, User Stories 2 and 3 can start in parallel after US1
- All plotting script creation tasks marked [P] can run in parallel
- Documentation tasks marked [P] can run in parallel
- Different user stories can be worked on in parallel by different team members

---

## Parallel Example: User Story 3

```bash
# Launch all plotting script creation tasks together:
Task: "Create band structure plotting script in scripts/plot_band_structure.gp"
Task: "Create quantum well plotting script in scripts/plot_quantum_well.gp"
Task: "Create g-factor plotting script in scripts/plot_gfactor.gp"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup
2. Complete Phase 2: Foundational (CRITICAL - blocks all stories)
3. Complete Phase 3: User Story 1
4. **STOP and VALIDATE**: Test User Story 1 independently
5. Deploy/demo if ready

### Incremental Delivery

1. Complete Setup + Foundational → Foundation ready
2. Add User Story 1 → Test independently → Deploy/Demo (MVP!)
3. Add User Story 2 → Test independently → Deploy/Demo
4. Add User Story 3 → Test independently → Deploy/Demo
5. Each story adds value without breaking previous stories

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together
2. Once Foundational is done:
   - Developer A: User Story 1 (build verification)
   - Developer B: User Story 2 (result validation) - after US1 complete
   - Developer C: User Story 3 (plotting scripts) - after US1 complete
3. Stories complete and integrate independently

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story should be independently completable and testable
- Manual validation tasks included for scientific accuracy verification
- Commit after each task or logical group
- Stop at any checkpoint to validate story independently
- Avoid: vague tasks, same file conflicts, cross-story dependencies that break independence
