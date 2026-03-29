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
