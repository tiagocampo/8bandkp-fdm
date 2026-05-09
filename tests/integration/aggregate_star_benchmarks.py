#!/usr/bin/env python3
"""Aggregate standard-star benchmark results into a single markdown table.

Runs all 7 standard-star benchmark scripts (S1-S7) via subprocess and parses
their stdout benchmark table output. Concatenates parsed rows into a unified
markdown table.

Usage: aggregate_star_benchmarks.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/ executables)
  source_dir -- path to repo root

Output: markdown table to stdout, summary to stderr.
"""

import os
import subprocess
import sys


SCRIPTS = [
    ("S1", "verify_star_gaas_bulk.py"),
    ("S2", "verify_star_inas_bulk.py"),
    ("S3", "verify_star_insb_bulk.py"),
    ("S4", "verify_star_gaas_algaas_qw.py"),
    ("S5", "verify_star_inas_gasb_qw.py"),
    ("S6", "verify_star_inas_gaas_qw.py"),
    ("S7", "verify_star_inas_wire.py"),
]


def run_benchmark(script_path, build_dir, source_dir):
    """Run a benchmark script and return (star_id, stdout, returncode)."""
    result = subprocess.run(
        [sys.executable, script_path, build_dir, source_dir],
        capture_output=True, text=True, timeout=600,
    )
    return result.stdout, result.returncode


def parse_benchmark_rows(output):
    """Extract benchmark table rows from script stdout.

    Looks for lines matching the benchmark table format.
    Returns list of raw row strings.
    """
    rows = []
    in_table = False
    for line in output.splitlines():
        stripped = line.strip()
        if stripped.startswith("Material") and "Observable" in stripped:
            in_table = True
            rows.append(stripped)
            continue
        if in_table:
            if stripped.startswith("---") or not stripped:
                if rows:
                    rows.append(stripped)
                continue
            if any(s in stripped for s in ("PASS", "FAIL", "SKIP")):
                rows.append(stripped)
            elif stripped and not stripped.startswith("="):
                break
    return rows


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>", file=sys.stderr)
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])
    integration_dir = os.path.join(source_dir, "tests", "integration")

    all_rows = []
    n_pass = 0
    n_fail = 0
    n_skip = 0

    for star_id, script_name in SCRIPTS:
        script_path = os.path.join(integration_dir, script_name)
        if not os.path.isfile(script_path):
            print(f"SKIP: {script_name} not found", file=sys.stderr)
            n_skip += 1
            continue

        print(f"Running {star_id} ({script_name})...", file=sys.stderr)
        try:
            stdout, rc = run_benchmark(script_path, build_dir, source_dir)
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT after 600s", file=sys.stderr)
            n_fail += 1
            continue

        rows = parse_benchmark_rows(stdout)
        if rows:
            all_rows.append((star_id, rows))
        else:
            all_rows.append((star_id, [f"No benchmark table produced (rc={rc})"]))

        if rc == 0:
            n_pass += 1
        else:
            n_fail += 1

    print("# Standard-Star Benchmark Results\n")
    print("| Star | Benchmark Table |")
    print("|------|----------------|")
    for star_id, rows in all_rows:
        for row in rows:
            print(f"| {star_id} | {row} |")

    print(f"\nSummary: {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP",
          file=sys.stderr)
    sys.exit(1 if n_fail > 0 else 0)


if __name__ == "__main__":
    main()
