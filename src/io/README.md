# Input/Output Module

This directory contains the input/output operations for the 8bandkp-fdm project.

## Contents

* `outputFunctions.f90`: File I/O and data output
  - Handles configuration file reading
  - Manages output file creation and writing
  - Implements data formatting and organization
  - Provides utilities for file unit management

## Features

* Configuration file parsing
* Eigenvalue output
* Eigenfunction output
* Band structure data output
* g-factor results output
* Potential profile output
* Wavefunction visualization data

## Output Files

* `eigenvalues.dat`: Energy eigenvalues
* `eigenfunctions_k_*_ev_*.dat`: Wavefunctions
* `parts.dat`: Band character analysis
* `fort.101`: Potential profile
* Additional files for g-factor calculations

## Dependencies

* Uses core module for types and parameters
* Minimal external dependencies
* Standard Fortran I/O features 