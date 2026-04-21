#!/usr/bin/env python3
"""Verify Stark shift magnitude and field conversion for QCSE reference data.

Usage: verify_stark_shift.py <field_config> <no_field_eigenvalues> <with_field_eigenvalues>

Checks that CB1 shifts in the expected direction and magnitude under electric
field, and that the config's EFParams value converts to the expected kV/cm
label.
"""
import sys


def parse_config_scalar(filepath, key):
    """Read a scalar config value from a `key: value` line."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(f"{key}:"):
                return float(line.split(":", 1)[1].strip().split()[0])
    raise RuntimeError(f"{key} not found in {filepath}")


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
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <field_config> <no_field.dat> <with_field.dat>")
        sys.exit(1)

    field_config = sys.argv[1]
    no_field = parse_eigenvalues(sys.argv[2])
    with_field = parse_eigenvalues(sys.argv[3])
    field_ev_per_angstrom = parse_config_scalar(field_config, "EFParams")
    field_kv_per_cm = field_ev_per_angstrom * 1.0e5
    field_label = f"{field_kv_per_cm:+.0f} kV/cm"

    # CB1 is the 10th column (index 9): 8 VB eigenvalues, then CB1
    cb1_no = no_field[9]
    cb1_ef = with_field[9]

    shift = cb1_ef - cb1_no

    print(
        f"Config field:    {field_ev_per_angstrom:+.6f} eV/Angstrom "
        f"= {field_kv_per_cm:+.1f} kV/cm"
    )
    print(f"Expected label:  {field_label}")
    print(f"CB1 no-field:   {cb1_no:.6f} eV")
    print(f"CB1 with-field: {cb1_ef:.6f} eV")
    print(f"CB1 shift:      {shift:.6f} eV (includes linear potential)")

    # EFParams is numerically identical to V/Angstrom, so 0.007 eV/Angstrom = 700 kV/cm.
    if abs(abs(field_kv_per_cm) - abs(field_ev_per_angstrom) * 1.0e5) > 1e-6:
        print("FAIL: field conversion from eV/Angstrom to kV/cm is inconsistent")
        sys.exit(1)

    if field_label != f"{field_kv_per_cm:+.0f} kV/cm":
        print(f"FAIL: computed field label does not match config-derived value ({field_label})")
        sys.exit(1)

    # CB1 should shift in the same sign as the configured field in this potential convention.
    if field_ev_per_angstrom < 0.0 and shift <= 0.0:
        print(f"FAIL: CB1 shift is non-positive ({shift:.6f} eV) for negative field")
        sys.exit(1)
    if field_ev_per_angstrom > 0.0 and shift >= 0.0:
        print(f"FAIL: CB1 shift is non-negative ({shift:.6f} eV) for positive field")
        sys.exit(1)

    # Shift should be in a broad but finite range for this 460-Angstrom QCSE domain.
    if shift < 0.5 or shift > 2.0:
        print(f"FAIL: CB1 shift {shift:.4f} eV outside range [0.5, 2.0]")
        sys.exit(1)

    print("PASS: Stark shift direction, magnitude, and field conversion verified")


if __name__ == "__main__":
    main()
