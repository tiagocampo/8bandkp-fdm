# Common source files for both executables
set(COMMON_SOURCES
    defs.f90
    mkl_spblas.f90
    mkl_sparse_handle.f90
    parameters.f90
    utils.f90
    finitedifferences.f90
    outputFunctions.f90
    hamiltonianConstructor.f90
)

# Sources specific to gfactor calculation
set(GFACTOR_SOURCES
    gfactor_functions.f90
    main_gfactor.f90
)

# Create library of common code
add_library(8bandkp_common STATIC ${COMMON_SOURCES})
target_link_libraries(8bandkp_common PUBLIC MKL::MKL)

# Special compile flags for mkl_spblas.f90
set_source_files_properties(mkl_spblas.f90
    PROPERTIES 
    COMPILE_FLAGS "-w -fcray-pointer"
)

# Main band structure executable
add_executable(bandStructure main.f90)
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
set_source_files_properties(mkl_sparse_handle.f90 PROPERTIES DEPENDS mkl_spblas.f90)
set_source_files_properties(parameters.f90 PROPERTIES DEPENDS defs.f90)
set_source_files_properties(utils.f90 PROPERTIES DEPENDS "defs.f90;mkl_spblas.f90")
set_source_files_properties(finitedifferences.f90 PROPERTIES DEPENDS "defs.f90;utils.f90")
set_source_files_properties(outputFunctions.f90 PROPERTIES DEPENDS defs.f90)
set_source_files_properties(hamiltonianConstructor.f90 PROPERTIES DEPENDS "defs.f90;finitedifferences.f90;mkl_spblas.f90;utils.f90")
set_source_files_properties(gfactor_functions.f90 PROPERTIES DEPENDS "defs.f90;hamiltonianConstructor.f90;mkl_spblas.f90")
set_source_files_properties(main_gfactor.f90 PROPERTIES DEPENDS "gfactor_functions.f90") 