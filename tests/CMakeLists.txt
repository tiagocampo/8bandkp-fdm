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

    add_pfunit_ctest(test_bdg_config
        TEST_SOURCES unit/test_bdg_config.pf
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

    add_pfunit_ctest(test_optical_spectra
        TEST_SOURCES unit/test_optical_spectra.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )
else()
    message(STATUS "pFUnit not found — unit tests disabled. Install pFUnit and set -DPFUNIT_DIR to enable.")
endif()

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

# Topological sweep mode regressions
add_test(
    NAME regression_topology_sweep_bhz
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_topology_sweep.sh
        $<TARGET_FILE:topologicalAnalysis>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/topology_sweep_bhz.toml
        bhz_analytic
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
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

set_tests_properties(
    verification_rung1_bulk_k0 verification_rung2_dispersion
    verification_rung3_qw verification_rung4_wire
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
