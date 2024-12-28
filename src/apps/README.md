# Applications Module

This directory contains the main applications of the 8bandkp-fdm project.

## Contents

* `main.f90`: Band structure calculation program
  - Implements the main band structure calculation workflow
  - Handles input parameter processing
  - Manages memory allocation and cleanup
  - Coordinates calculation steps
  - Outputs results to files

* `main_gfactor.f90`: g-factor calculation program
  - Implements the g-factor calculation workflow
  - Processes input parameters
  - Manages memory for g-factor calculations
  - Coordinates calculation steps
  - Outputs g-factor results

## Usage

### Band Structure Calculation
```bash
./bandStructure
# Reads input.cfg
# Outputs:
#   - eigenvalues.dat
#   - eigenfunctions_k_*_ev_*.dat
#   - parts.dat
#   - fort.101 (for quantum wells)
```

### G-factor Calculation
```bash
./gfactorCalculation
# Reads input.cfg
# Outputs g-factor specific results
```

## Dependencies

* Uses all other modules:
  - core: for types and parameters
  - math: for numerical operations
  - physics: for Hamiltonian and g-factor calculations
  - io: for file operations
* Requires input configuration file
* Generates output in the working directory 