# Common source files for both executables
set(COMMON_SOURCES
    core/defs.f90
    core/parameters.f90
    core/utils.f90
    math/mkl_spblas.f90
    math/mkl_sparse_handle.f90
    math/finitedifferences.f90
    io/outputFunctions.f90
    io/input_parser.f90
    physics/hamiltonianConstructor.f90
)

# Create library of common code
add_library(8bandkp_common STATIC ${COMMON_SOURCES})

# Link MKL and FFTW3
target_link_libraries(8bandkp_common PUBLIC MKL::MKL FFTW3::FFTW3)

# Module include directory for consumers
target_include_directories(8bandkp_common PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Special compile flags for mkl_spblas.f90 (vendor MKL interface)
set_source_files_properties(math/mkl_spblas.f90
    PROPERTIES
    COMPILE_FLAGS "-w -fcray-pointer"
)

# Preprocessor for files with #ifdef/#include
set_source_files_properties(
    core/defs.f90
    core/parameters.f90
    core/utils.f90
    math/finitedifferences.f90
    io/outputFunctions.f90
    io/input_parser.f90
    physics/hamiltonianConstructor.f90
    physics/gfactor_functions.f90
    apps/main.f90
    apps/main_gfactor.f90
    PROPERTIES
    Fortran_PREPROCESS ON
)

# Main band structure executable
add_executable(bandStructure apps/main.f90)
target_link_libraries(bandStructure PRIVATE 8bandkp_common)

# G-factor calculation executable
add_executable(gfactorCalculation
    physics/gfactor_functions.f90
    apps/main_gfactor.f90
)
target_link_libraries(gfactorCalculation PRIVATE 8bandkp_common)

# Set module output directory for all targets
set_target_properties(8bandkp_common bandStructure gfactorCalculation
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
