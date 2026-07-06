#!/usr/bin/env bash
# TDD-red: ctest wrapper for slim Pfaffian witness projection test.
# Pins the slim-Pfaffian row-index bug (C-1) at
# src/physics/topological_analysis.f90:1722, :1827. Fails until A.2 lands
# (multi-site scan band-major fix).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
python3 tests/integration/test_slim_pfaffian_witness_projection.py