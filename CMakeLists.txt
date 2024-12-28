cmake_minimum_required(VERSION 3.12)
project(8bandkp-fdm Fortran)

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options
option(BUILD_TESTING "Build tests" OFF)
option(ENABLE_OPENMP "Enable OpenMP support" ON)
option(USE_MKL "Use Intel MKL" ON)

# Compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fimplicit-none -fno-second-underscore")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -fbounds-check -fcheck=all -Wall")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -mtune=native -march=native -flto -funroll-loops")
endif()

# OpenMP
if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# MKL
if(USE_MKL)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
endif()

# Create build directory for modules
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Add subdirectories
add_subdirectory(src)

# Install targets
install(TARGETS bandStructure gfactorCalculation
        RUNTIME DESTINATION bin)

# Install example files
install(FILES bulk.example gfactor.example quantumwell.example
        DESTINATION share/${PROJECT_NAME}/examples)

# Package configuration
include(CPack)
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "8-band k·p method implementation using finite differences")
set(CPACK_PACKAGE_VENDOR "Your Institution") 