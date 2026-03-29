#!/usr/bin/env python3
"""Compare numerical output files against reference data for regression testing.

Usage: compare_output.py <reference_file> <test_file> [--tolerance <tol>]

Compares corresponding numerical columns. Lines starting with '#' or non-numeric
content are skipped. Exit code 0 = pass, 1 = fail.
"""

import sys
import argparse
import math


def parse_numeric_lines(filepath):
    """Parse file, returning list of lists of floats (one per numeric line)."""
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            tokens = line.split()
            try:
                vals = [float(t) for t in tokens]
                rows.append(vals)
            except ValueError:
                continue
    return rows


def compare(ref_rows, test_rows, tolerance):
    """Compare two parsed datasets. Returns (pass, max_rel_error, details)."""
    if len(ref_rows) != len(test_rows):
        return False, 0.0, f"Row count mismatch: ref={len(ref_rows)} test={len(test_rows)}"

    max_rel_err = 0.0
    details = []
    for i, (ref_row, test_row) in enumerate(zip(ref_rows, test_rows)):
        if len(ref_row) != len(test_row):
            details.append(f"Line {i+1}: column count mismatch {len(ref_row)} vs {len(test_row)}")
            return False, max_rel_err, "\n".join(details)
        for j, (rv, tv) in enumerate(zip(ref_row, test_row)):
            if rv == 0.0 and tv == 0.0:
                continue
            if rv == 0.0:
                rel_err = abs(tv)
            else:
                rel_err = abs((tv - rv) / rv)
            if rel_err > tolerance:
                details.append(f"Line {i+1} col {j+1}: ref={rv:.10e} test={tv:.10e} rel_err={rel_err:.2e}")
            max_rel_err = max(max_rel_err, rel_err)

    if details:
        return False, max_rel_err, "\n".join(details)
    return True, max_rel_err, "All values match within tolerance"


def main():
    parser = argparse.ArgumentParser(description="Compare numerical output files")
    parser.add_argument("reference", help="Reference (golden) data file")
    parser.add_argument("test", help="Test output file to compare")
    parser.add_argument("--tolerance", type=float, default=1e-8, help="Relative tolerance (default: 1e-8)")
    args = parser.parse_args()

    ref_rows = parse_numeric_lines(args.reference)
    test_rows = parse_numeric_lines(args.test)

    passed, max_err, details = compare(ref_rows, test_rows, args.tolerance)

    if passed:
        print(f"PASS: max relative error = {max_err:.2e} (tolerance = {args.tolerance:.2e})")
        sys.exit(0)
    else:
        print(f"FAIL: {details}")
        sys.exit(1)


if __name__ == "__main__":
    main()
