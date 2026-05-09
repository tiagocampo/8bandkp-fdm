# tests/support/CMakeLists.txt — shared test-support library

add_library(8bandkp_test_support STATIC
    csr_test_helpers.f90
    krylov_helpers.f90
    krylov_reference_data.f90
)
target_link_libraries(8bandkp_test_support PUBLIC 8bandkp_common)
