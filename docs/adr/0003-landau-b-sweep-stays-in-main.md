# ADR 0003: Landau B-sweep stays in main.f90, not simulation_setup

## Status

Accepted

## Context

C1 (architecture cleanup) adds Landau as a fourth confinement mode to `simulation_setup`, absorbing the init, Hamiltonian build, and k-point solve into the existing `simulation_setup_init` / `setup_build_H` / `setup_solve_kpoint_serial` interface. This collapses ~188 lines of Landau-specific code in `main.f90` down to ~20 lines for the standard k-sweep path.

However, Landau mode also has a unique **B-sweep** (~81 lines) that loops over magnetic field values instead of k-points, rebuilds the Hamiltonian with updated B-field components, and writes a fan diagram to `output/landau_fan.dat`. This sweep does not fit the existing `setup_solve_kpoint_serial(wavevector)` interface.

Three options were considered:

- **(A)** Leave B-sweep in main.f90; simulation_setup handles init/build/k-solve only
- **(B)** Add `setup_solve_Bsweep` method to simulation_setup
- **(C)** Generalize to `setup_solve_sweep(setup, cfg, sweep_params)` polymorphic interface

## Decision

**Option A**: the B-sweep stays in `main.f90`. simulation_setup handles the generic pipeline (init, build, k-solve); `main.f90` handles the Landau-specific B-sweep and fan diagram output.

## Why

The B-sweep is 81 lines of specialized physics output (magnetic field parameterization, fan diagram formatting). Options B and C would move this complexity into `simulation_setup` without reducing it — they'd just relocate it. The method would be Landau-specific (violating the generic intent of simulation_setup) or over-generalized for a single caller. The real win of C1 is collapsing init + workspace query + k-sweep; the B-sweep is a separate concern that naturally belongs in the application layer.
