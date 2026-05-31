# tests/support/CMakeLists.txt — shared test-support library

add_library(8bandkp_test_support STATIC
    csr_test_helpers.f90
    krylov_helpers.f90
    krylov_reference_data.f90
    trk_helpers.f90
)
target_link_libraries(8bandkp_test_support PUBLIC 8bandkp_common)

# Reference data generator — links full physics stack
add_executable(generate_krylov_reference
    generate_krylov_reference.f90
)
target_link_libraries(generate_krylov_reference PRIVATE 8bandkp_common 8bandkp_test_support)

# Custom target to regenerate reference data in-place
# Usage: cmake --build build --target regenerate_krylov_references
add_custom_target(regenerate_krylov_references
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regenerate_references.sh
    COMMENT "Regenerating Krylov reference data"
)
