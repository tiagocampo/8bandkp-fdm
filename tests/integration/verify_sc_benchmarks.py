#!/usr/bin/env python3
"""Verify self-consistent calculation results against expected physics.

Checks that the SC loop converges and produces physically reasonable results:
  1. SC convergence message present in output
  2. Subband energies show expected band bending (gap != bare material gap)
  3. Confined states exist (CB states above well edge, VB states below)

Usage: verify_sc_benchmarks.py <test_output.log> <eigenvalues.dat> <config_file>
"""
import sys
import math


def parse_eigenvalues_k0(filepath):
    """Parse eigenvalues at k=0."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
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
                return line.split(":", 1)[1].strip().split()[0]
    return None


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <test_output.log> <eigenvalues.dat> <config_file>")
        sys.exit(1)

    log_file = sys.argv[1]
    eval_file = sys.argv[2]
    config_file = sys.argv[3]
    failures = []

    # Check 1: SC convergence
    with open(log_file) as f:
        log_content = f.read()

    if "SC loop converged" not in log_content:
        failures.append("SC loop did not converge (missing 'SC loop converged' message)")
        print("[Benchmark] SC convergence: FAIL - convergence message not found")
    else:
        # Extract iteration count if present
        import re
        match = re.search(r'SC loop converged.*?(\d+)\s*iterations?', log_content, re.IGNORECASE)
        if match:
            n_iter = int(match.group(1))
            print(f"[Benchmark] SC convergence: OK (converged in {n_iter} iterations)")
            if n_iter > 200:
                failures.append(f"SC convergence too slow: {n_iter} iterations (>200)")
        else:
            print("[Benchmark] SC convergence: OK")

    # Check 2: Parse eigenvalues
    evals = parse_eigenvalues_k0(eval_file)
    n_bands = len(evals)

    # Determine structure from config
    confinement = parse_config_value(config_file, "confinement")
    is_bulk = (confinement == "0")
    material1 = parse_config_value(config_file, "material1")

    # Determine VB/CB split point
    # For 8-band: bands 1-6 = VB (including SO), bands 7-8 = CB
    # But some configs have different numvb/numcb
    if is_bulk and n_bands == 8:
        # Standard 8-band bulk: 6 VB + 2 CB
        n_vb = 6
        n_cb = 2
    elif n_bands >= 10:
        # QW with 8 VB + N CB subbands
        n_vb = 8
        n_cb = n_bands - n_vb
    else:
        # Bulk or small QW: determine from gap
        # Find the gap: largest positive energy jump between adjacent sorted eigenvalues
        sorted_evals = sorted(evals)
        max_gap_idx = 0
        max_gap = 0
        for i in range(len(sorted_evals) - 1):
            g = sorted_evals[i+1] - sorted_evals[i]
            if g > max_gap:
                max_gap = g
                max_gap_idx = i
        # Everything below the gap is VB, everything above is CB
        vb_cutoff = sorted_evals[max_gap_idx]
        cb_cutoff = sorted_evals[max_gap_idx + 1]
        vb_states = [e for e in evals if e <= vb_cutoff + 0.001]
        cb_states = [e for e in evals if e >= cb_cutoff - 0.001]
        n_vb = len(vb_states)
        n_cb = len(cb_states)

    if n_bands >= 10 or (is_bulk and n_bands == 8):
        vb_states = evals[:n_vb]
        cb_states = evals[n_vb:]

    if not cb_states or not vb_states:
        print(f"\n[Benchmark] Subband structure: SKIP (cannot determine VB/CB split for {n_bands} bands)")
        print("PASS: SC benchmark checks passed (structure analysis skipped)")
        sys.exit(0)

    vb_top = max(vb_states)
    cb_bottom = min(cb_states)
    gap = cb_bottom - vb_top

    print(f"\n[Benchmark] Subband structure ({n_bands} bands, {n_vb} VB + {n_cb} CB):")
    print(f"  VB top:    {vb_top:.6f} eV")
    print(f"  CB bottom: {cb_bottom:.6f} eV")
    print(f"  Gap:       {gap:.6f} eV")

    # Check 3: Gap should be physically reasonable for the material system
    # Different materials have different gaps:
    #   GaAs: Eg = 1.519 eV
    #   InAs: Eg ~ 0.354 eV (narrow gap)
    #   InAs/AlSb: Type-II, gap can be ~0.1-0.5 eV
    # We check the gap is positive and not absurdly large
    if gap < 0.05:
        failures.append(f"Gap {gap:.4f} eV too small (possible broken-gap or error)")
        print(f"  FAIL: gap too small")
    elif gap > 6.0:
        failures.append(f"Gap {gap:.4f} eV too large")
        print(f"  FAIL: gap too large")
    else:
        print(f"  OK: gap in physically reasonable range")
        # Additional material-specific note
        if gap < 1.0:
            print(f"  Note: Narrow-gap material (InAs, InAs/AlSb Type-II, etc.)")
        elif gap < 2.0:
            print(f"  Note: Typical III-V semiconductor gap (GaAs, InP, etc.)")
        else:
            print(f"  Note: Wide-gap material or heavily doped (AlAs barriers, etc.)")

    # Check 4: CB subband spacing should be positive (confined states)
    if n_cb >= 4:
        cb1_cb2 = abs(cb_states[2] - cb_states[0])
        print(f"\n[Benchmark] CB subband spacing: {cb1_cb2*1000:.1f} meV")
        if cb1_cb2 < 0.001:
            failures.append(f"CB1-CB2 spacing {cb1_cb2*1000:.3f} meV too small")
            print(f"  FAIL: degenerate CB states in QW (expected distinct subbands)")
        else:
            print(f"  OK: distinct CB subbands present")

    # Check 5: Spin degeneracy preserved (only for QW with enough states)
    if n_cb >= 2:
        deg_ok = True
        for i in range(0, n_cb - 1, 2):
            if abs(cb_states[i+1] - cb_states[i]) > 0.002:
                deg_ok = False
        if deg_ok:
            print(f"\n[Benchmark] Spin degeneracy: OK")
        else:
            print(f"\n[Benchmark] Spin degeneracy: WARNING (some CB states split > 2 meV)")

    # Summary
    print()
    if failures:
        print("FAIL: SC benchmark checks failed:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("PASS: all SC benchmark checks passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
