# Core Infrastructure

## Purpose
Owns fundamental constants, simulation orchestration, and configuration parsing. Defines the `simulation_config` derived type.

## Entry Points
- `defs.f90` - Kinds, global parameters, derived types
- `input_parser.f90` - TOML input parsing
- `simulation_setup.f90` - Global orchestration

## Contracts & Invariants
- No deps on other src/ subdirectories
- All config validation consolidated here
- Use `error stop` for all fatal config errors

## Patterns
- Always check `input.toml` against documented sections in `docs/reference/input-reference.md`.
- Validate config via `validate()` and `validate_semantic()` before any computational step.

## Related Context
- Primary config source: `./input.toml`
- Downstream dependencies: `src/math/` and `src/physics/`
