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
    regression_gfactor_cb regression_bulk_inas
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
