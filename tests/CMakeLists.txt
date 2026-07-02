# tests/CMakeLists.txt — pFUnit unit tests + regression tests

add_subdirectory(support)

# pFUnit installs as PFUNIT (uppercase) — use CMAKE_PREFIX_PATH or -DPFUNIT_DIR
find_package(PFUNIT QUIET)

if(PFUNIT_FOUND)
    # --- pFUnit unit tests ---
    add_pfunit_ctest(test_defs
        TEST_SOURCES unit/test_defs.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_finitedifferences
        TEST_SOURCES unit/test_finitedifferences.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_utils
        TEST_SOURCES unit/test_utils.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_parameters
        TEST_SOURCES unit/test_parameters.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_hamiltonian
        TEST_SOURCES unit/test_hamiltonian.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_poisson
        TEST_SOURCES unit/test_poisson.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_charge_density
        TEST_SOURCES unit/test_charge_density.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_sc_loop
        TEST_SOURCES unit/test_sc_loop.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_geometry
        TEST_SOURCES unit/test_geometry.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_hamiltonian_2d
        TEST_SOURCES unit/test_hamiltonian_2d.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_eigensolver
        TEST_SOURCES unit/test_eigensolver.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_strain_solver
        TEST_SOURCES unit/test_strain_solver.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_hamiltonian_blocks
        TEST_SOURCES unit/test_hamiltonian_blocks.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_confinement_assembly
        TEST_SOURCES unit/test_confinement_assembly.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_csr_spmv
        TEST_SOURCES unit/test_csr_spmv.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_csr_helpers
        TEST_SOURCES unit/test_csr_helpers.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_csr_structural
        TEST_SOURCES unit/test_csr_structural.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_krylov_helpers
        TEST_SOURCES unit/test_krylov_helpers.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_krylov_snapshots
        TEST_SOURCES unit/test_krylov_snapshots.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_optical
        TEST_SOURCES unit/test_optical.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_optical_qw
        TEST_SOURCES unit/test_optical_qw.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_topology_types
        TEST_SOURCES unit/test_topology_types.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_topology_parser
        TEST_SOURCES unit/test_topology_parser.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_magnetic_field
        TEST_SOURCES unit/test_magnetic_field.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_chern_number
        TEST_SOURCES unit/test_chern_number.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_bdg_hamiltonian
        TEST_SOURCES unit/test_bdg_hamiltonian.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Issue 02 (U5, Layer A of U4): PHS oracle pinning the BdG hole-block
    # canonical convention (-conjg(H0(-k)) per ADR 0007). Four-field matrix
    # (no_Z/Z x no_P/P) at generic k + cross-builder spectral identity.
    # RED on pre-Issue-03 code; GREEN post-Issue-03.
    add_pfunit_ctest(test_bdg_phs
        TEST_SOURCES unit/test_bdg_phs.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_bdg_config
        TEST_SOURCES unit/test_bdg_config.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_bdg_evaluator
        TEST_SOURCES unit/test_bdg_evaluator.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Task 1.5 (P1 BdG stabilization, ADR 0008 §2): cross-module pinning
    # of pairing_sign/pairing_partner as PUBLIC exports from bdg_hamiltonian.
    # RED until Task 1.6 adds them to the public list.
    add_pfunit_ctest(test_pairing_sign_cross_module
        TEST_SOURCES unit/test_pairing_sign_cross_module.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_phase_diagram
        TEST_SOURCES unit/test_phase_diagram.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_green_functions
        TEST_SOURCES unit/test_green_functions.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_edge_states
        TEST_SOURCES unit/test_edge_states.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_landau
        TEST_SOURCES unit/test_landau.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_z2_invariant
        TEST_SOURCES unit/test_z2_invariant.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_fd_convergence
        TEST_SOURCES unit/test_fd_convergence.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_simulation_setup
        TEST_SOURCES unit/test_simulation_setup.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Issue #02: drives the solver-config derivation entry point directly
    # across every (geometry x method x mode) combination without running a
    # full executable (the case untestable through simulation_setup_init).
    add_pfunit_ctest(test_solver_derivation
        TEST_SOURCES unit/test_solver_derivation.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_optical_spectra
        TEST_SOURCES unit/test_optical_spectra.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    add_pfunit_ctest(test_trk_sum_rule
        TEST_SOURCES unit/test_trk_sum_rule.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Issue 01: Pfaffian math (Parlett-Reid tridiagonalization)
    add_pfunit_ctest(test_pfaffian
        TEST_SOURCES unit/test_pfaffian.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Issue 01: Kitaev-chain Majorana number harness
    add_pfunit_ctest(test_kitaev_majorana
        TEST_SOURCES unit/test_kitaev_majorana.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Issue 04 (U6): Sticlet Majorana polarization as MZM discriminator.
    # Tests the pure function majorana_polarization + polarization_result_t
    # added to src/physics/topological_analysis.f90. Spin-sector signs s_σ
    # derived from KTD7 Nambu ordering per ADR 0007.
    add_pfunit_ctest(test_majorana_polarization
        TEST_SOURCES unit/test_majorana_polarization.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )

    # Issue 07 (U10): slim projected Pfaffian witness (S1 empirical + S2
    # analytical bands 7-8). Tests the pure wrapper wire_pfaffian_witness
    # added to src/physics/topological_analysis.f90.
    add_pfunit_ctest(test_wire_pfaffian_witness
        TEST_SOURCES unit/test_wire_pfaffian_witness.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )
else()
    message(STATUS "pFUnit not found — unit tests disabled. Install pFUnit and set -DPFUNIT_DIR to enable.")
endif()

# Issue 04 Fix Round 1 (Important 1): standalone Fortran driver that
# intentionally triggers `error stop` on size mismatch in
# majorana_polarization. pFUnit 4.x cannot catch Fortran `error stop`
# (per ADR 0002 + test_parameters.pf precedent), so the size-mismatch
# branch is verified by spawning this driver from a shell integration
# test (test_majorana_polarization_size_mismatch.sh) and checking the
# exit code + error message. Compiled unconditionally (no pFUnit needed)
# so the error-stop path is always covered.
add_executable(test_majorana_polarization_size_mismatch
    integration/test_majorana_polarization_size_mismatch.f90
)
target_link_libraries(test_majorana_polarization_size_mismatch PRIVATE 8bandkp_common)
add_test(
    NAME validation_majorana_polarization_size_mismatch
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_majorana_polarization_size_mismatch.sh
        $<TARGET_FILE:test_majorana_polarization_size_mismatch>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(validation_majorana_polarization_size_mismatch PROPERTIES LABELS "unit" TIMEOUT 30)

# --- Regression tests ---
# Use shell scripts from tests/integration/ to avoid CMake shell escaping issues

# Bulk GaAs band structure
add_test(
    NAME regression_bulk_gaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_bulk_bandstructure.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/bulk_gaas_kx.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/bulk_gaas_kx
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# QW AlSbW/GaSbW/InAsW band structure
add_test(
    NAME regression_qw_alsbw_gasbw_inasw
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_bandstructure.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_alsbw_gasbw_inasw.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/qw_alsb_gasb_inas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_qw_character_and_output_dir
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_character_and_output_dir.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_alsbw_gasbw_inasw.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_state_character.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_qw_absorption_polarization
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_absorption_polarization.sh
        $<TARGET_FILE:opticalProperties>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_gaas_algaas_absorption.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_absorption_polarization.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_qw_absorption_strained_polarization
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_absorption_polarization.sh
        $<TARGET_FILE:opticalProperties>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_ingaas_algaas_strained_absorption.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_absorption_polarization.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_scattering_transition_sanity
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_scattering_transition_sanity.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_gaas_algaas_qcse_scattering.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_scattering_transition_sanity.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor CB
add_test(
    NAME regression_gfactor_cb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_qw_cb.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_qw_cb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_gfactor_no_optics
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor_no_optics.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_qw_no_optics.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor Bulk GaAs VB
add_test(
    NAME regression_gfactor_bulk_gaas_vb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_bulk_gaas_vb.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_bulk_gaas_vb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor Bulk GaAsW CB
add_test(
    NAME regression_gfactor_bulk_gaasw_cb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_bulk_gaasw_cb.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_bulk_gaasw_cb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor QW VB (AlSbW/GaSbW/InAsW)
add_test(
    NAME regression_gfactor_qw_vb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_qw_vb.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_qw_vb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Bulk InAs band structure
add_test(
    NAME regression_bulk_inas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_bulk_bandstructure.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/bulk_inas_kx.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/bulk_inas_kx
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Landau levels in InAs bulk at B=5T (Zeeman + Peierls)
add_test(
    NAME regression_landau_inas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_landau_levels.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/landau_InAs.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/landau_InAs
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Landau mode (confinement=3) bulk regression tests
add_test(
    NAME regression_landau_bulk_gaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_landau_bulk.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/landau_bulk_GaAs.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/landau_bulk_GaAs
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_landau_bulk_inas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_landau_bulk.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/landau_bulk_InAs.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/landau_bulk_InAs
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_landau_bulk_inas_bsweep
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_landau_bulk.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/landau_bulk_InAs_Bsweep.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/landau_bulk_InAs_Bsweep
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Landau analytical validation (CB ground state vs E_n = E_C + hbar*omega_c/2)
add_test(
    NAME regression_landau_analytical
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_landau_analytical.py
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/landau_bulk_InAs.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_bulk_gaas regression_qw_alsbw_gasbw_inasw
    regression_qw_character_and_output_dir regression_qw_absorption_polarization
    regression_qw_absorption_strained_polarization
    regression_scattering_transition_sanity regression_gfactor_cb
    regression_gfactor_no_optics regression_gfactor_bulk_gaas_vb
    regression_gfactor_bulk_gaasw_cb regression_gfactor_qw_vb
    regression_bulk_inas
    regression_landau_inas regression_landau_bulk_gaas regression_landau_bulk_inas
    regression_landau_bulk_inas_bsweep regression_landau_analytical
    PROPERTIES LABELS "regression"
)

# SC GaAs/AlAs QW self-consistent calculation
add_test(
    NAME regression_sc_gaas_alas_qw
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_loop.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_gaas_alas_qw.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_gaas_alas_qw
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_sc_gaas_alas_qw
    PROPERTIES LABELS "regression"
)

# SC InAs/AlSb QW self-consistent (narrow-gap Type-II, Winkler parameters, 77K)
add_test(
    NAME regression_sc_qw_inas_alsb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_qw_inas_alsb.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qw_inas_alsb.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qw_inas_alsb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_sc_qw_inas_alsb
    PROPERTIES LABELS "regression"
)

add_test(
    NAME regression_qw_optical_matrix_selection
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_optical_matrix_selection.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_gaas_algaas_optics_selection.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_qw_optical_matrix_selection
    PROPERTIES LABELS "regression"
)

# QW optical properties via opticalProperties executable (commutator-based velocity)
add_test(
    NAME regression_qw_optics_commutator
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_optics_commutator.sh
        $<TARGET_FILE:opticalProperties>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_optics_commutator.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/qw_optics_commutator
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_qw_optics_commutator
    PROPERTIES LABELS "regression"
)

# QCSE GaAs/AlGaAs QW without electric field (Harrison benchmark)
add_test(
    NAME regression_qcse_gaas_algaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qcse_no_ef.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qcse_gaas_algaas.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# QCSE GaAs/AlGaAs QW with electric field -70 kV/cm (Stark shift validation)
add_test(
    NAME regression_qcse_gaas_algaas_ef
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qcse_with_ef.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qcse_gaas_algaas_ef.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas_ef
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_qcse_gaas_algaas regression_qcse_gaas_algaas_ef
    PROPERTIES LABELS "regression"
)

# Stark shift verification (compares reference eigenvalues, no simulation)
add_test(
    NAME regression_stark_shift
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_stark_shift.py
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qcse_gaas_algaas_ef.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas/eigenvalues.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas_ef/eigenvalues.dat
)

# Modulation-doped AlAs/GaAs/AlAs QW self-consistent (charge neutrality mode)
add_test(
    NAME regression_sc_mod_doped_gaas_algaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_mod_doped.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_mod_doped_gaas_algaas.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_mod_doped_gaas_algaas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_sc_mod_doped_gaas_algaas
    PROPERTIES LABELS "regression"
)

# Bulk doped GaAs self-consistent (uniform n-type, no barriers)
add_test(
    NAME regression_sc_bulk_gaas_doped
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_bulk_doped.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_bulk_gaas_doped.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_bulk_gaas_doped
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_sc_bulk_gaas_doped
    PROPERTIES LABELS "regression"
)

# Wire GaAs rectangular wire band structure
add_test(
    NAME regression_wire_gaas_rectangle
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_bandstructure.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_rectangle.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/wire_gaas_rectangle
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_gaas_rectangle
    PROPERTIES LABELS "regression"
)

# Wire InAs/GaAs core-shell with strain (Stier & Bimberg benchmark)
add_test(
    NAME regression_wire_inas_gaas_core_shell
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_core_shell.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_strain.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/wire_inas_gaas_core_shell
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_inas_gaas_core_shell
    PROPERTIES LABELS "regression"
)

# Wire InSb g-factor (extreme g-factor regime)
add_test(
    NAME regression_wire_insb_gfactor
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_insb_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_insb_gfactor.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/wire_insb_gfactor
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_insb_gfactor
    PROPERTIES LABELS "regression"
)

add_test(
    NAME regression_wire_optical_selection
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_optical_selection.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_optical_window.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_wire_optical_selection.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_optical_selection
    PROPERTIES LABELS "regression"
)

add_test(
    NAME regression_wire_parser_optional_sections
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_parser_optional_sections.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_rectangle.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_strain.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_parser_optional_sections
    PROPERTIES LABELS "regression"
)

add_test(
    NAME regression_wire_dense_sparse_consistency
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_dense_sparse_consistency.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_dense_sparse_dense.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_dense_sparse_sparse.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_dense_sparse_consistency
    PROPERTIES LABELS "regression"
)

# Wire hexagonal geometry
add_test(
    NAME regression_wire_gaas_hexagon
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_hexagon.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_hexagon.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/wire_gaas_hexagon
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_gaas_hexagon
    PROPERTIES LABELS "regression"
)

# SC wire convergence
add_test(
    NAME regression_sc_wire_gaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_wire.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_wire_gaas.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_wire_gaas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_sc_wire_gaas
    PROPERTIES LABELS "regression"
)

# Topological analysis QHE (Chern number via QWZ model)
add_test(
    NAME regression_topology_qhe
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_qhe.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_qhe_qwz.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_qhe
    PROPERTIES LABELS "regression"
)

# QWZ Chern number regression (three u values: -0.8, 0.5, 2.5)
add_test(
    NAME regression_topology_qwz_chern
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_chern_qwz.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_qwz_chern_u-0.8.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_qwz_chern_u0.5.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_qwz_chern_u2.5.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_qwz_chern
    PROPERTIES LABELS "regression"
)

# BHZ Z2 invariant regression (trivial d=58A and topological d=70A)
add_test(
    NAME regression_topology_bhz_z2
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_bhz_z2.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_bhz_z2_trivial.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_bhz_z2_topological.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_bhz_z2
    PROPERTIES LABELS "regression"
)

# QW Fu-Kane parity Z2 regression
add_test(
    NAME regression_topology_qw_fukane_z2
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_fukane_z2.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_qw_fukane.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_qw_fukane_z2
    PROPERTIES LABELS "regression"
)

# QWZ conductance regression
add_test(
    NAME regression_topology_conductance_qwz
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_conductance_qwz.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_conductance_qwz.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_conductance_qwz
    PROPERTIES LABELS "regression"
)

# Topological spectral mode regressions
add_test(
    NAME regression_topology_spectral_bulk
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_spectral.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_spectral_bulk.toml
        1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_topology_spectral_qw
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_spectral.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_spectral_qw.toml
        1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_topology_spectral_wire
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_spectral.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_spectral_wire.toml
        1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_spectral_bulk regression_topology_spectral_qw
    regression_topology_spectral_wire
    PROPERTIES LABELS "regression"
)

# Issue 06 / Unit U9: BdG LDOS + A(k,E) + Nambu-resolved LDOS on the
# topological wire. Verifies bdq_spectral mode produces all three
# output files and that the Majorana zero-energy LDOS peak exists.
add_test(
    NAME regression_bdq_spectral_wire
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_bdg_spectral.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_bdq_spectral_wire.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    regression_bdq_spectral_wire
    PROPERTIES LABELS "regression" TIMEOUT 600
)

# Wire BdG strain-omission regression (#04):
# Asserts the wire BdG spectrum SHIFTS when [strain] is applied. Before #04
# the wire BdG driver skipped compute_strain, so strained/unstrained runs
# were identical (the latent physics defect). Physically they must differ.
add_test(
    NAME regression_wire_bdg_strain_shift
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_bdg_strain_shift.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_strain_bdg.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    regression_wire_bdg_strain_shift
    PROPERTIES LABELS "regression;strain-validation" TIMEOUT 900
)

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

# Issue 07 (U10): 2D wire BdG minigap colormap regression. B swept,
# mu fixed at conduction-band edge; asserts the colormap is non-flat
# (gap varies across the grid) and matches the golden reference. The
# gap_map / z2_map columns are routed through wire_pfaffian_witness_sweep
# (Pfaffian, not heuristic) per the U10 acceptance criteria.
add_test(
    NAME regression_wire_bdg_topological_2d
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/test_wire_bdg_topological_2d.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_bdg_topological_2d.toml
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/wire_inas_gaas_bdg_topological_2d_phase.dat
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    regression_wire_bdg_topological_2d
    PROPERTIES LABELS "regression" TIMEOUT 900
)

# Issue 05 / Unit U7: dense-QW BdG rung end-to-end certification.
# Drives a B-sweep on the dense 16N x 16N BdG QW Hamiltonian and asserts
# the canonical class-D open->close->reopen pattern at B_crit (covers AE2
# + AE4). Uses the existing run_bdg_qw executable path + Issue 00's
# evaluator output (min_gap); the pure-function Pfaffian and
# polarization seams are consumed by Issue 01 / 04 pFUnit tests and by
# Issue 07's slim Pfaffian witness.
add_test(
    NAME regression_dense_qw_bdg_rung
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_dense_qw_bdg_rung.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_dense_qw_bdg.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    regression_dense_qw_bdg_rung
    PROPERTIES LABELS "regression" TIMEOUT 900
)

# Topological sweep mode regressions
add_test(
    NAME regression_topology_sweep_bhz
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_sweep.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_sweep_bhz.toml
        bhz_analytic
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Issue 08 / Unit U12: Lecture 13 acceptance gate. Runs the 5-section
# lecture script + asserts 4-witness B_crit agreement + asserts the
# false-PASS line is gone + asserts the Issue 07 regression is green.
# Regression lock against the all-zero B_crit false-PASS failure mode
# (per docs/solutions/patterns/2026-06-27-lecture-table-simulation-extracted.md).
add_test(
    NAME lecture_13_acceptance_gate
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_lecture_13_acceptance_gate.sh
        ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    lecture_13_acceptance_gate
    PROPERTIES LABELS "regression" TIMEOUT 1200
)

add_test(
    NAME regression_topology_sweep_qw
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_sweep.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_sweep_qw.toml
        qw_fukane
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_sweep_bhz regression_topology_sweep_qw
    PROPERTIES LABELS "regression"
)

# Rashba wire phase boundary regression (mu=0.5 meV, Delta=0.3 meV)
add_test(
    NAME regression_topology_rashba_phase
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_rashba_phase_boundary.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_rashba_phase.toml
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_topology_rashba_phase
    PROPERTIES LABELS "regression"
)

# --- 8-Band Verification Ladder (4 rungs) ---
# Rung 1: Bulk k=0 structural validation (R1-R5)
add_test(
    NAME verification_rung1_bulk_k0
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_8band_rung1_bulk_k0.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Rung 2: Bulk dispersion effective mass validation (R6-R9)
add_test(
    NAME verification_rung2_dispersion
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_8band_rung2_dispersion.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Rung 3: QW subband energy validation (R10, R12, R13)
add_test(
    NAME verification_rung3_qw
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_8band_rung3_qw.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Rung 4: Wire internal consistency validation (R14-R16)
add_test(
    NAME verification_rung4_wire
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_8band_rung4_wire.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Rung 7: g-factor validation vs Roth formula (R7.1, R7.2)
add_test(
    NAME verification_rung7_gfactor
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_8band_rung7_gfactor.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Rung 8: Optical observables validation (R8.1, R8.2, R8.3)
add_test(
    NAME verification_rung8_optical
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_8band_rung8_optical.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_rung1_bulk_k0 verification_rung2_dispersion
    verification_rung3_qw verification_rung4_wire
    verification_rung7_gfactor verification_rung8_optical
    PROPERTIES LABELS "verification"
)

# --- Strain Verification (Rungs 5-6) ---
add_test(
    NAME verification_rung5_strain_bulk
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_rung5_bulk.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME verification_rung6_strain_qw
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_rung6_qw.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_rung5_strain_bulk verification_rung6_strain_qw
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- Solver Smart Defaults Verification ---
add_test(
    NAME verification_solver_defaults
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_solver_defaults.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_solver_defaults
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- QW Dense vs Sparse (FEAST) Solver Agreement ---
add_test(
    NAME verification_qw_sparse_solver
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_sparse_solver.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_qw_sparse_solver
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- Dense vs Sparse eigenvalue agreement across all geometries ---
# Extends the dense-vs-sparse rung to bulk/QW/wire/Landau. QW and wire
# are the genuine dense-vs-FEAST cross-check; bulk and Landau verify
# mode-dispatch consistency (FEAST invalid for 8x8 bulk; Landau FEAST
# enablement is issue #10).
add_test(
    NAME verification_dense_sparse_all_geometries
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_dense_sparse_all_geometries.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_dense_sparse_all_geometries
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- QW dispersion-aware energy-window envelope (issue #03) ---
# A wide-k QW FEAST sweep with AUTO window (emin = emax = 0). The envelope
# authority (apply_solver_window) unions the Gershgorin bounds at the sweep
# endpoints (k = 0 and k = k_max) into ONE stable window; asserts FEAST
# delivers the requested eigenvalue count at every k-point, including the
# dispersion extremum. At HEAD (pre-#03) the envelope block was unreachable
# and the sweep fell through to the type-default [-1, +1] window, which
# under-delivered at k_max.
add_test(
    NAME verification_qw_window_envelope
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_window_envelope.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_qw_window_envelope
    PROPERTIES LABELS "verification" TIMEOUT 300
)

# --- QW g-factor FEAST-enable regression (issue #08) ---
# Runs the gfactorCalculation executable on the SAME GaSbW/InAsW/AlSbW QW
# twice - method=DENSE+FULL and method=FEAST+ENERGY (AUTO window) - and
# asserts (1) the g-factor tensor matches within 1e-6 relative and (2)
# BOTH backends return the FULL 8N spectrum (8*fd_step eigenfunctions).
# Before #08 the QW+FEAST combo was rejected by a stopgap guard; this test
# gates the removal of that guard. A truncated FEAST spectrum would move
# the Lowdin g-factor, so the agreement check catches silent truncation.
add_test(
    NAME verification_qw_gfactor_dense_feast
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_gfactor_dense_feast.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_qw_gfactor_dense_feast
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- QW optics FEAST-enable regression (issue #09) ---
# Runs the opticalProperties executable on the SAME Al30Ga70As/GaAs QW
# twice - method=DENSE (AUTO->DENSE+INDEX, the golden default) and
# method=FEAST+ENERGY (AUTO sweep-envelope window) - and asserts (1) the
# FEAST path is actually taken (dispatch gate: a FEAST/envelope
# diagnostic appears in stdout; at HEAD method=FEAST was silently
# overridden to DENSE) and (2) the output spectra (absorption_TE/TM,
# gain_TE/TM, spontaneous_TE/TM, ISBT - whichever the config enables)
# match within 1e-6 relative. Divergent optics would mean the FEAST
# window truncated the spectrum or the band indexing shifted between
# INDEX and ENERGY mode.
add_test(
    NAME verification_qw_optics_dense_feast
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_optics_dense_feast.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_qw_optics_dense_feast
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- QW band-structure FEAST-enable regression (review #4) ---
# Runs bandStructure on the SAME AlGaAs/GaAs QW twice - method=DENSE
# (AUTO->DENSE+INDEX, the golden default) and method=FEAST+ENERGY (AUTO
# sweep-envelope window) - and asserts the gap-straddling eigenvalue bands
# agree within tolerance. Pre-#4 the QW FEAST CSR path returned the deepest
# valence bands (lowest-M copy); this locks the corrected [il,iuu] slice.
add_test(
    NAME verification_qw_bandstructure_dense_feast
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_bandstructure_dense_feast.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(
    verification_qw_bandstructure_dense_feast
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- Landau FEAST-enable regression (issue #10) ---
# Runs the bandStructure executable on the SAME InAs Landau config twice -
# method=DENSE (AUTO->DENSE+INDEX, the golden default) and method=FEAST+ENERGY
# (AUTO sweep-envelope window) - and asserts (1) the FEAST path is actually
# taken (dispatch gate: a FEAST/envelope diagnostic appears in stdout; at HEAD
# method=FEAST hit a broken result-handling path that assumed INDEX sizing)
# and (2) the eigenvalues.dat output matches within 1e-6 relative. Divergent
# eigenvalues would mean the FEAST window truncated the spectrum, the
# result-handling extracted the wrong slice, or the band indexing shifted
# between INDEX and ENERGY mode. The config is the golden landau_InAs.toml
# shape (InAs, nx=100, Bz=5T, ky-sweep nsteps=10) - a genuine in-plane
# k-dispersion sweep at FIXED B (the B-sweep fan diagram is a separate
# application-layer block, untouched by issue #10 per ADR 0003).
add_test(
    NAME verification_landau_dense_feast
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_landau_dense_feast.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_landau_dense_feast
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- Wire Strain Validation ---
add_test(
    NAME strain_validation_wire
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_wire_profile.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    strain_validation_wire
    PROPERTIES LABELS "strain-validation" TIMEOUT 600
)

# Wire strain quantitative validation (wide-core InAs/GaAs vs Bir-Pikus analytical)
add_test(
    NAME strain_validation_wire_quantitative
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_wire_quantitative.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    strain_validation_wire_quantitative
    PROPERTIES LABELS "strain-validation;verification" TIMEOUT 600
)

# --- Standard-Star Benchmarks (S1-S7) ---
# Physics validation against published literature observables.
# Each test runs the Fortran executables and compares against analytical
# predictions or regression reference values.

add_test(
    NAME star_s1_gaas_bulk
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_gaas_bulk.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME star_s2_inas_bulk
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_inas_bulk.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME star_s3_insb_bulk
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_insb_bulk.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME star_s4_gaas_algaas_qw
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_gaas_algaas_qw.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME star_s5_inas_gasb_qw
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_inas_gasb_qw.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME star_s6_inas_gaas_strained_qw
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_inas_gaas_qw.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME star_s7_inas_wire
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_inas_wire.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    star_s1_gaas_bulk star_s2_inas_bulk star_s3_insb_bulk
    star_s4_gaas_algaas_qw star_s5_inas_gasb_qw
    star_s6_inas_gaas_strained_qw star_s7_inas_wire
    PROPERTIES LABELS "standard-star;verification"
              TIMEOUT 600
)

# --- Phase 21: Publishable Validation Benchmarks ---

# ISBT oscillator strength benchmark (2 materials × 3 widths vs infinite-well)
add_test(
    NAME verification_isbt_benchmark
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_isbt_benchmark.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# SC Schrödinger-Poisson benchmark (charge neutrality, Fermi level, subband shift, potential)
add_test(
    NAME verification_sc_benchmark
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_sc_benchmark.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_isbt_benchmark
    PROPERTIES LABELS "verification;standard-star"
              TIMEOUT 1200
)

set_tests_properties(
    verification_sc_benchmark
    PROPERTIES LABELS "verification;standard-star"
              TIMEOUT 900
)

# Validation coverage matrix — cross-checks universe file vs test annotations
add_test(
    NAME coverage_matrix
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/coverage_matrix.py
        --source-dir ${CMAKE_SOURCE_DIR}
        --output ${CMAKE_BINARY_DIR}/coverage_report.md
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coverage_matrix PROPERTIES LABELS "coverage" TIMEOUT 120)

# --- Convergence Tests (Richardson Extrapolation on Physics Observables) ---
# On-demand only: ctest -L convergence

# Exciton smoke test (R11a)
add_test(
    NAME convergence_exciton_smoke
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_exciton_smoke.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(convergence_exciton_smoke PROPERTIES LABELS "convergence" TIMEOUT 120)

# QW grid convergence (U4): S4 GaAs/AlGaAs, S5 InAs/GaSb, S6 InAs/GaAs strained
add_test(
    NAME convergence_qw_grid
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_grid_convergence.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(convergence_qw_grid PROPERTIES LABELS "convergence" TIMEOUT 1200)

# QW FD order convergence (U5): S4-S6 at orders {2,4,6,8,10}
add_test(
    NAME convergence_qw_order
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_order_convergence.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(convergence_qw_order PROPERTIES LABELS "convergence" TIMEOUT 1200)

# Wire grid convergence (U6): S7 InAs wire
add_test(
    NAME convergence_wire
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_wire_convergence.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(convergence_wire PROPERTIES LABELS "convergence" TIMEOUT 3600)

# SC convergence benchmark (U7): GaAs/AlAs QW SC loop
add_test(
    NAME convergence_sc
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_convergence.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(convergence_sc PROPERTIES LABELS "convergence" TIMEOUT 3600)

# Exciton convergence benchmark (U8): GaAs/AlGaAs QW exciton
add_test(
    NAME convergence_exciton
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_exciton_convergence.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(convergence_exciton PROPERTIES LABELS "convergence" TIMEOUT 1200)

# Validation rejection tests: verify error-stop branches reject bad configs
# Covers error-stop paths in validate()/validate_semantic() that pFUnit cannot test.
add_test(
    NAME validation_rejects_bad_configs
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_validate_rejects_bad_configs.sh
        $<TARGET_FILE:bandStructure>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(validation_rejects_bad_configs PROPERTIES LABELS "validation" TIMEOUT 60)

# U8: BdG solver-window validate guard rejection test.
add_test(
    NAME validation_rejects_bad_topology
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_validate_rejects.sh
        $<TARGET_FILE:topologicalAnalysis>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(validation_rejects_bad_topology PROPERTIES LABELS "validation" TIMEOUT 60)
