#!/usr/bin/env python3
"""Verify Stark shift magnitude from QCSE reference eigenvalue files.

Usage: verify_stark_shift.py <no_field_eigenvalues> <with_field_eigenvalues>

Checks that CB1 shifts in the expected direction and magnitude under electric field.
"""
import sys


def parse_eigenvalues(filepath):
    """Parse first data line, return list of floats."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            return [float(t) for t in line.split()]
    raise RuntimeError(f"No data in {filepath}")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <no_field.dat> <with_field.dat>")
        sys.exit(1)

    no_field = parse_eigenvalues(sys.argv[1])
    with_field = parse_eigenvalues(sys.argv[2])

    # CB1 is the 10th column (index 9): 8 VB eigenvalues, then CB1
    cb1_no = no_field[9]
    cb1_ef = with_field[9]

    shift = cb1_ef - cb1_no

    print(f"CB1 no-field:   {cb1_no:.6f} eV")
    print(f"CB1 with-field: {cb1_ef:.6f} eV")
    print(f"CB1 shift:      {shift:.6f} eV (includes linear potential)")

    # CB1 should increase with negative field (linear potential ramp dominates)
    if shift <= 0.0:
        print(f"FAIL: CB1 shift is non-positive ({shift:.6f} eV)")
        sys.exit(1)

    # Shift should be in reasonable range (~1 eV from linear ramp for -70 kV/cm)
    if shift < 0.5 or shift > 2.0:
        print(f"FAIL: CB1 shift {shift:.4f} eV outside range [0.5, 2.0]")
        sys.exit(1)

    print("PASS: Stark shift direction and magnitude verified")


if __name__ == "__main__":
    main()
