# Mathematical Operations Module

This directory contains the mathematical operations and numerical methods used in the 8bandkp-fdm project.

## Contents

* `mkl_spblas.f90`: MKL sparse BLAS interface
  - Provides interface to Intel MKL sparse BLAS routines
  - Handles sparse matrix operations
  - Optimizes memory usage for large matrices

* `mkl_sparse_handle.f90`: MKL sparse matrix handling
  - Manages sparse matrix data structures
  - Provides utilities for sparse matrix manipulation
  - Implements sparse matrix conversions

* `finitedifferences.f90`: Finite difference method implementation
  - Implements finite difference discretization
  - Provides stencil operations for derivatives
  - Handles boundary conditions

## Dependencies

* Requires Intel MKL library for sparse matrix operations
* Uses core module for basic types and utilities
* Can be configured to use standard LAPACK/BLAS instead of MKL 