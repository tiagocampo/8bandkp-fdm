# tests/CMakeLists.txt — pFUnit unit tests + regression tests

# pFUnit installs as PFUNIT (uppercase) — use CMAKE_PREFIX_PATH or -DPFUNIT_DIR
find_package(PFUNIT QUIET)

if(PFUNIT_FOUND)
    # --- pFUnit unit tests ---
    add_pfunit_ctest(test_defs
        TEST_SOURCES unit/test_defs.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_finitedifferences
        TEST_SOURCES unit/test_finitedifferences.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_utils
        TEST_SOURCES unit/test_utils.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_parameters
        TEST_SOURCES unit/test_parameters.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_hamiltonian
        TEST_SOURCES unit/test_hamiltonian.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_poisson
        TEST_SOURCES unit/test_poisson.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_charge_density
        TEST_SOURCES unit/test_charge_density.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_sc_loop
        TEST_SOURCES unit/test_sc_loop.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_geometry
        TEST_SOURCES unit/test_geometry.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_hamiltonian_2d
        TEST_SOURCES unit/test_hamiltonian_2d.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_eigensolver
        TEST_SOURCES unit/test_eigensolver.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_strain_solver
        TEST_SOURCES unit/test_strain_solver.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_csr_spmv
        TEST_SOURCES unit/test_csr_spmv.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_optical
        TEST_SOURCES unit/test_optical.pf
        LINK_LIBRARIES 8bandkp_common
        LABELS "unit"
    )

    add_pfunit_ctest(test_optical_qw
        TEST_SOURCES unit/test_optical_qw.pf
        LINK_LIBRARIES 8bandkp_common
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/bulk_gaas_kx.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/bulk_gaas_kx
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# QW AlSbW/GaSbW/InAsW band structure
add_test(
    NAME regression_qw_alsbw_gasbw_inasw
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_bandstructure.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_alsbw_gasbw_inasw.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/qw_alsb_gasb_inas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_qw_character_and_output_dir
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_character_and_output_dir.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_alsbw_gasbw_inasw.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_state_character.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_qw_absorption_polarization
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_absorption_polarization.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_gaas_algaas_absorption.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_qw_absorption_polarization.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_scattering_transition_sanity
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_scattering_transition_sanity.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/qw_gaas_algaas_qcse_scattering.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_scattering_transition_sanity.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# g-factor CB
add_test(
    NAME regression_gfactor_cb
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_qw_cb.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/gfactor_qw_cb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME regression_gfactor_no_optics
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_gfactor_no_optics.sh
        $<TARGET_FILE:gfactorCalculation>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/gfactor_qw_no_optics.cfg
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Bulk InAs band structure
add_test(
    NAME regression_bulk_inas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_bulk_bandstructure.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/bulk_inas_kx.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/bulk_inas_kx
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_bulk_gaas regression_qw_alsbw_gasbw_inasw
    regression_qw_character_and_output_dir regression_qw_absorption_polarization
    regression_scattering_transition_sanity regression_gfactor_cb
    regression_gfactor_no_optics regression_bulk_inas
    PROPERTIES LABELS "regression"
)

# SC GaAs/AlAs QW self-consistent calculation
add_test(
    NAME regression_sc_gaas_alas_qw
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_loop.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_gaas_alas_qw.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qw_inas_alsb.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qw_inas_alsb
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_sc_qw_inas_alsb
    PROPERTIES LABELS "regression"
)

# QCSE GaAs/AlGaAs QW without electric field (Harrison benchmark)
add_test(
    NAME regression_qcse_gaas_algaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qcse_no_ef.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qcse_gaas_algaas.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/compare_output.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# QCSE GaAs/AlGaAs QW with electric field -70 kV/cm (Stark shift validation)
add_test(
    NAME regression_qcse_gaas_algaas_ef
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qcse_with_ef.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qcse_gaas_algaas_ef.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_qcse_gaas_algaas_ef.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas/eigenvalues.dat
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/data/sc_qcse_gaas_algaas_ef/eigenvalues.dat
)

# Modulation-doped AlAs/GaAs/AlAs QW self-consistent (charge neutrality mode)
add_test(
    NAME regression_sc_mod_doped_gaas_algaas
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_sc_mod_doped.sh
        $<TARGET_FILE:bandStructure>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_mod_doped_gaas_algaas.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/sc_bulk_gaas_doped.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_rectangle.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_strain.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_insb_gfactor.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_optical_window.cfg
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
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_gaas_rectangle.cfg
        ${CMAKE_CURRENT_SOURCE_DIR}/regression/configs/wire_inas_gaas_strain.cfg
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    regression_wire_parser_optional_sections
    PROPERTIES LABELS "regression"
)
