cmake_minimum_required(VERSION 3.15)
project(8bandkp-fdm Fortran)

# --- Find dependencies ---
# Intel MKL (oneAPI provides MKLConfig.cmake)
# Defaults: LP64 interface + sequential threading (matching original Makefile)
set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface: lp64 or ilp64")
set(MKL_THREADING "sequential" CACHE STRING "MKL threading: sequential, intel_thread, etc.")
find_package(MKL CONFIG REQUIRED)

# FFTW3
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(FFTW3 REQUIRED)

# --- Fortran compiler flags (matching original Makefile) ---
set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -mtune=native -march=native -fno-second-underscore -ffree-form -fimplicit-none -m64 -O2 -flto -funroll-loops")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -fbounds-check -fcheck=all -fopenmp")

# --- Centralize .mod files ---
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# --- Build the project ---
add_subdirectory(src)

# --- Testing ---
option(BUILD_TESTING "Build tests (requires pFUnit for unit tests)" OFF)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
