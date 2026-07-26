# Application Entry Points

## Purpose
Owns the entry point programs for the four executables: bandStructure, gfactorCalculation, opticalProperties, and topologicalAnalysis.

## Entry Points
- `main.f90` - Band structure
- `main_gfactor.f90` - G-factor
- `main_optics.f90` - Optical spectra
- `main_topology.f90` - Topological analysis

## Contracts & Invariants
- Reads `input.toml` and delegates to `simulation_setup.f90`.
- Outputs to `output/` directory.

## Patterns
- Build targets defined in main `CMakeLists.txt` and built into `build/src/`.
- Run using `./build/src/<exe>`.

## Related Context
- Configuration parsing: `src/core/AGENTS.md`
- Core orchestraion: `src/core/AGENTS.md`
