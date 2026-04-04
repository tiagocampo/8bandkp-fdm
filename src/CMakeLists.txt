# Common source files for both executables
set(COMMON_SOURCES
    core/defs.f90
    core/parameters.f90
    core/utils.f90
    math/sparse_matrices.f90
    math/geometry.f90
    math/finitedifferences.f90
    math/eigensolver.f90
    math/linalg.f90
    io/outputFunctions.f90
    io/input_parser.f90
    physics/hamiltonianConstructor.f90
    physics/charge_density.f90
    physics/poisson.f90
    physics/sc_loop.f90
    physics/strain_solver.f90
    physics/gfactor_functions.f90
)

# Create library of common code
add_library(8bandkp_common STATIC ${COMMON_SOURCES})

# Link MKL and FFTW3
target_link_libraries(8bandkp_common PUBLIC MKL::MKL FFTW3::FFTW3)

# Module include directory for consumers
target_include_directories(8bandkp_common PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Enable MKL FEAST eigensolver (available when building with Intel MKL)
target_compile_definitions(8bandkp_common PUBLIC USE_MKL_FEAST)

# Preprocessor for files with #ifdef/#include
set_source_files_properties(
    core/defs.f90
    core/parameters.f90
    core/utils.f90
    math/sparse_matrices.f90
    math/geometry.f90
    math/finitedifferences.f90
    math/eigensolver.f90
    math/linalg.f90
    io/outputFunctions.f90
    io/input_parser.f90
    physics/hamiltonianConstructor.f90
    physics/charge_density.f90
    physics/poisson.f90
    physics/sc_loop.f90
    physics/strain_solver.f90
    physics/gfactor_functions.f90
    apps/main.f90
    apps/main_gfactor.f90
    PROPERTIES
    Fortran_PREPROCESS ON
)

# Main band structure executable
add_executable(bandStructure apps/main.f90)
target_link_libraries(bandStructure PRIVATE 8bandkp_common OpenMP::OpenMP_Fortran)

# G-factor calculation executable
add_executable(gfactorCalculation
    apps/main_gfactor.f90
)
target_link_libraries(gfactorCalculation PRIVATE 8bandkp_common OpenMP::OpenMP_Fortran)

# Set module output directory for all targets
set_target_properties(8bandkp_common bandStructure gfactorCalculation
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
