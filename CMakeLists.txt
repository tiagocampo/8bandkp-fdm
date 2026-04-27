cmake_minimum_required(VERSION 3.15)
project(8bandkp-fdm Fortran)

# --- Find dependencies ---
# Intel MKL (oneAPI provides MKLConfig.cmake)
# Defaults: LP64 interface + intel_thread for internal MKL parallelism
# Note: CMPLR_ROOT must point to Intel compiler for libiomp5 (e.g. /opt/intel/oneapi/compiler/2025.0)
set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface: lp64 or ilp64")
set(MKL_THREADING "sequential" CACHE STRING "MKL threading: sequential, intel_thread")
find_package(MKL CONFIG REQUIRED)

# FFTW3
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(FFTW3 REQUIRED)

# OpenMP for parallel regions
find_package(OpenMP REQUIRED)

# ARPACK-NG (optional sparse eigensolver fallback)
find_package(ARPACK QUIET)

# BLAS/LAPACK backend selector (future portability)
set(LINALG_BACKEND "MKL" CACHE STRING "BLAS/LAPACK backend: MKL, OPENBLAS, AOCL")

# --- Fortran compiler flags (gfortran) ---
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -fno-second-underscore -ffree-form -fimplicit-none -m64 -fno-strict-aliasing")
    # NOTE: -ffast-math was replaced by its individual safe components.
    # Kept: -fno-signed-zeros (~2%), -freciprocal-math (~3%), -fno-trapping-math
    # Dropped: -ffinite-math-only (removes NaN/Inf guards), -funsafe-math-optimizations
    # Added: -fno-strict-aliasing (needed for MKL PARDISO which is polymorphic real/complex)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=haswell -mtune=haswell -ftree-vectorize -fstack-arrays -fno-signed-zeros -freciprocal-math -fno-trapping-math -fno-strict-aliasing -flto=auto -fuse-linker-plugin")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -fbounds-check -fcheck=all")
endif()

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

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
