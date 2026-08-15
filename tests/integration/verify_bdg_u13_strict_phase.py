#!/usr/bin/env python3
# COVERAGE: observable=majorana_number geometry=wire material=InAs tier=verification
"""Verify the U13 strict wire-BdG phase-diagram contract."""

from pathlib import Path
import re
import sys


def _read_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 4:
            raise AssertionError(f"{path}: expected 4 columns, got {len(parts)}: {line}")
        rows.append((float(parts[0]), float(parts[1]), int(round(float(parts[2]))), float(parts[3])))
    return rows


def main():
    if len(sys.argv) != 4:
        print("usage: verify_bdg_u13_strict_phase.py PHASE CONFIG WITNESS")
        return 2

    phase_path, config_path, witness_path = map(Path, sys.argv[1:])
    config = config_path.read_text()
    match = re.search(r"^\s*nk_par\s*=\s*(\d+)\s*$", config, flags=re.MULTILINE)
    if match is None or int(match.group(1)) <= 1:
        print("FAIL: U13 phase fixture must configure nk_par > 1")
        return 1

    try:
        rows = _read_rows(phase_path)
    except (OSError, ValueError, AssertionError) as exc:
        print(f"FAIL: {exc}")
        return 1
    if not rows:
        print(f"FAIL: {phase_path} has no data rows")
        return 1

    z2_values = {row[2] for row in rows}
    if not z2_values <= {-1, 0, 1}:
        print(f"FAIL: strict z2 column has non-native values {sorted(z2_values)}")
        return 1
    if -1 not in z2_values or not ({0, 1} & z2_values):
        print(f"FAIL: strict z2 column is flat or lacks a closure/trivial region: {sorted(z2_values)}")
        return 1

    pattern = re.compile(
        r"^B=([\d.eE+-]+)\s+\|Pf\|=([\d.eE+-]+)\s+reason=(\d+)\s*$",
        flags=re.MULTILINE,
    )
    try:
        witness_text = witness_path.read_text()
    except OSError as exc:
        print(f"FAIL: cannot read {witness_path}: {exc}")
        return 1
    witness_rows = pattern.findall(witness_text)
    if not witness_rows:
        print(f"FAIL: {witness_path} has no strict witness rows with reason=...")
        return 1
    reasons = {int(reason) for _, _, reason in witness_rows}
    if not reasons <= {0, 1, 2, 3}:
        print(f"FAIL: disagreement_reason contains unknown values {sorted(reasons)}")
        return 1
    print(
        f"PASS: strict U13 phase rows={len(rows)} z2={sorted(z2_values)} "
        f"witness_rows={len(witness_rows)} reasons={sorted(reasons)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
