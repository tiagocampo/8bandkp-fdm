# Common source files for both executables
set(COMMON_SOURCES
    core/defs.f90
    core/parameters.f90
    core/utils.f90
    math/mkl_spblas.f90
    math/mkl_sparse_handle.f90
    math/finitedifferences.f90
    io/outputFunctions.f90
    physics/hamiltonianConstructor.f90
)

# Sources specific to gfactor calculation
set(GFACTOR_SOURCES
    physics/gfactor_functions.f90
    apps/main_gfactor.f90
)

# Create library of common code
add_library(8bandkp_common STATIC ${COMMON_SOURCES})
target_link_libraries(8bandkp_common PUBLIC MKL::MKL)

# Special compile flags for mkl_spblas.f90
set_source_files_properties(math/mkl_spblas.f90
    PROPERTIES 
    COMPILE_FLAGS "-w -fcray-pointer"
)

# Main band structure executable
add_executable(bandStructure apps/main.f90)
target_link_libraries(bandStructure PRIVATE 8bandkp_common)

# G-factor calculation executable
add_executable(gfactorCalculation ${GFACTOR_SOURCES})
target_link_libraries(gfactorCalculation PRIVATE 8bandkp_common)

# Set module output directory
set_target_properties(8bandkp_common bandStructure gfactorCalculation
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Dependencies between source files
set_source_files_properties(math/mkl_sparse_handle.f90 PROPERTIES DEPENDS math/mkl_spblas.f90)
set_source_files_properties(core/parameters.f90 PROPERTIES DEPENDS core/defs.f90)
set_source_files_properties(core/utils.f90 PROPERTIES DEPENDS "core/defs.f90;math/mkl_spblas.f90")
set_source_files_properties(math/finitedifferences.f90 PROPERTIES DEPENDS "core/defs.f90;core/utils.f90")
set_source_files_properties(io/outputFunctions.f90 PROPERTIES DEPENDS core/defs.f90)
set_source_files_properties(physics/hamiltonianConstructor.f90 PROPERTIES DEPENDS "core/defs.f90;math/finitedifferences.f90;math/mkl_spblas.f90;core/utils.f90")
set_source_files_properties(physics/gfactor_functions.f90 PROPERTIES DEPENDS "core/defs.f90;physics/hamiltonianConstructor.f90;math/mkl_spblas.f90")
set_source_files_properties(apps/main_gfactor.f90 PROPERTIES DEPENDS "physics/gfactor_functions.f90") 