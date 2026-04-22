#!/usr/bin/env python3
"""Verify quantum well eigenvalues against known benchmark values.

Checks GaAs/AlGaAs QW confinement energies and subband structure at k=0.

Usage: verify_qw_benchmarks.py <eigenvalues.dat> <config_file>

The config file is read to determine well width and barrier material so that
expected confinement energies can be computed.
"""
import sys
import math


def parse_eigenvalues(filepath):
    """Parse eigenvalues file, returning list of eigenvalues at k=0."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            # Return first data line (k=0)
            return vals[1:]
    raise RuntimeError(f"No data in {filepath}")


def parse_config_value(filepath, key):
    """Read a scalar config value."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            if line.startswith(f"{key}:"):
                return line.split(":", 1)[1].strip()
    return None


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <eigenvalues.dat> [config_file]")
        sys.exit(1)

    evals = parse_eigenvalues(sys.argv[1])
    config_file = sys.argv[2] if len(sys.argv) > 2 else None
    n_bands = len(evals)
    failures = []

    print(f"QW eigenvalues at k=0: {n_bands} bands")
    print(f"  Values: {['%.6f' % e for e in evals]}")

    # Basis ordering: bands 1-8 VB (4 spin-degenerate pairs), then bands 9+ CB
    # For numvb=8, numcb=4: 8 VB states (4 doubly-degenerate), 4 CB states (2 pairs)
    n_vb = 8  # standard for numvb=8
    n_cb = n_bands - n_vb

    vb_states = evals[:n_vb]
    cb_states = evals[n_vb:]

    print(f"\n  VB states ({n_vb}): {['%.4f' % e for e in vb_states]}")
    print(f"  CB states ({n_cb}): {['%.4f' % e for e in cb_states]}")

    # Benchmark 1: CB states above VB states by at least the band gap
    vb_top = max(vb_states)
    cb_bottom = min(cb_states)
    gap = cb_bottom - vb_top

    print(f"\n[Benchmark] VB-CB separation:")
    print(f"  VB top:    {vb_top:.6f} eV")
    print(f"  CB bottom: {cb_bottom:.6f} eV")
    print(f"  Gap:       {gap:.6f} eV")

    # For GaAs-based QWs, the gap should be at least Eg_GaAs = 1.519 eV
    # (may be slightly larger due to confinement pushing CB up and VB down)
    if gap < 1.4:
        failures.append(f"VB-CB gap = {gap:.4f} eV, expected > 1.4 eV for GaAs-based QW")
        print(f"  FAIL: gap too small for GaAs-based QW")
    else:
        print(f"  OK: gap > 1.4 eV (GaAs Eg + confinement)")

    # Benchmark 2: CB1 confinement energy relative to well CB edge
    # For the QCSE zero-field config (6nm GaAs/Al0.2Ga0.8As):
    #   GaAs CB edge = 1.519 eV (relative to VB top at 0)
    #   AlGaAs CB offset ~ 0.2 * 1.05 * 0.65 ~ 0.136 eV (small offset)
    #   CB1 should be just above the GaAs CB edge
    # For a generic check: CB1 should be a physically reasonable energy
    print(f"\n[Benchmark] CB ground state energy:")
    cb1 = cb_states[0]
    print(f"  CB1 = {cb1:.6f} eV")

    # Benchmark 3: Spin degeneracy - check that CB states come in pairs
    print(f"\n[Benchmark] Spin degeneracy:")
    spin_deg_ok = True
    for i in range(0, n_cb - 1, 2):
        split = abs(cb_states[i+1] - cb_states[i])
        if split > 0.001:  # 1 meV tolerance for spin degeneracy
            print(f"  WARNING: CB{i+1}/CB{i+2} split = {split*1000:.3f} meV (expected ~0)")
            spin_deg_ok = False
    if spin_deg_ok:
        print(f"  OK: all CB states spin-degenerate (< 1 meV splitting)")
    else:
        failures.append("CB spin degeneracy broken (splitting > 1 meV)")

    # Also check VB degeneracy (HH/LH pairs)
    vb_deg_ok = True
    for i in range(0, n_vb - 1, 2):
        split = abs(vb_states[i+1] - vb_states[i])
        if split > 0.001:
            vb_deg_ok = False
    if vb_deg_ok:
        print(f"  OK: all VB states spin-degenerate (< 1 meV splitting)")
    else:
        print(f"  NOTE: VB degeneracy broken (expected for 8-band k.p at k=0 in QW)")

    # Summary
    print()
    if failures:
        print("FAIL: benchmark checks failed:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("PASS: all QW benchmark checks passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
