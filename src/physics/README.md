# Physics Module

This directory contains the physics calculations and models for the 8bandkp-fdm project.

## Contents

* `hamiltonianConstructor.f90`: k·p Hamiltonian construction
  - Implements 8-band k·p Hamiltonian
  - Handles bulk and quantum well systems
  - Includes strain and external field effects
  - Supports multiple material layers

* `gfactor_functions.f90`: g-factor calculations
  - Implements g-factor calculations using Lowdin partitioning
  - Handles spin-orbit coupling effects
  - Calculates Zeeman splitting
  - Supports both bulk and quantum well systems

## Features

* Full 8-band k·p model implementation
* Support for heterostructures
* External field effects
* Strain effects
* Spin-orbit coupling
* g-factor calculations

## Dependencies

* Uses core module for types and parameters
* Requires math module for numerical operations
* Interfaces with MKL for matrix operations 