#!/usr/bin/env python3
# COVERAGE: observable=bdg_ldos geometry=wire material=InAsW
# COVERAGE: observable=bdg_spectral_function geometry=wire material=InAsW
"""Issue 06 / Unit U9: BdG LDOS + A(k,E) + Nambu-resolved LDOS verifier.

Runs `topologicalAnalysis` with `mode = bdq_spectral` and asserts that:
  1. All three output files exist (bdg_ldos.dat, bdg_ldos_nambu.dat,
     bdg_spectral.dat) and have non-empty content.
  2. bdg_ldos.dat has a peak at E=0 (the Majorana zero-energy mode)
     in the topological phase.
  3. bdg_spectral.dat shows non-negative A(k,E) (the spectral
     function is non-negative by construction).

Args: <topologicalAnalysis_exe> <config_file>
"""
import os
import re
import sys
import subprocess
from pathlib import Path

EXE = str(Path(sys.argv[1]).resolve())
CONFIG = Path(sys.argv[2])
OMP = "4"


def run(text, tag, timeout=600):
    d = Path(f"run_{tag}")
    d.mkdir(parents=True, exist_ok=True)
    (d / "input.toml").write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS=OMP)
    p = subprocess.run([EXE], cwd=d, env=env, capture_output=True,
                       text=True, timeout=timeout)
    (d / "run.log").write_text(p.stdout + "\n--STDERR--\n" + p.stderr)
    if p.returncode != 0:
        print(f"FAIL: topologicalAnalysis exited {p.returncode} for {tag}")
        print(p.stdout)
        print(p.stderr)
        sys.exit(1)
    return d


def read_2col(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return rows


def read_3col(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                pass
    return rows


def read_5col(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                rows.append((float(parts[0]), float(parts[1]),
                             float(parts[2]), float(parts[3]),
                             float(parts[4])))
            except ValueError:
                pass
    return rows


def main():
    config_text = CONFIG.read_text()
    workdir = run(config_text, "bdq_spectral")

    # 1. All three output files exist and have valid content
    ldos_path = workdir / "output" / "bdg_ldos.dat"
    nambu_path = workdir / "output" / "bdg_ldos_nambu.dat"
    spectral_path = workdir / "output" / "bdg_spectral.dat"

    for p in (ldos_path, nambu_path, spectral_path):
        if not p.exists() or p.stat().st_size == 0:
            print(f"FAIL: {p.name} missing or empty")
            sys.exit(1)

    ldos_rows = read_2col(ldos_path)
    if len(ldos_rows) < 3:
        print(f"FAIL: bdg_ldos.dat has too few rows ({len(ldos_rows)})")
        sys.exit(1)

    nambu_rows = read_3col(nambu_path)
    if len(nambu_rows) < 3:
        print(f"FAIL: bdg_ldos_nambu.dat has too few rows ({len(nambu_rows)})")
        sys.exit(1)

    spectral_rows = read_5col(spectral_path)
    if len(spectral_rows) < 3:
        print(f"FAIL: bdg_spectral.dat has too few rows ({len(spectral_rows)})")
        sys.exit(1)

    # 2. bdg_ldos.dat has a peak at E=0 (the Majorana mode is the largest
    #    LDOS feature for a BdG spectrum in the topological phase).
    peak_zero = 0.0
    peak_far = 0.0
    ldos_at_zero = [v for (_, v) in ldos_rows if abs(_) < 0.01]
    ldos_far = [v for (_, v) in ldos_rows if abs(_) > 0.05]
    if ldos_at_zero:
        peak_zero = max(ldos_at_zero)
    if ldos_far:
        peak_far = max(ldos_far)
    if ldos_at_zero and ldos_far:
        if peak_zero <= peak_far:
            print(f"FAIL: E=0 LDOS peak ({peak_zero}) not > off-peak LDOS ({peak_far})")
            sys.exit(1)

    # 3. bdg_spectral.dat: A(k,E) must be non-negative everywhere.
    for r in spectral_rows:
        if r[4] < 0.0:
            print(f"FAIL: A(k,E) < 0 at row {r}")
            sys.exit(1)

    print(f"PASS: bdq_spectral (Issue 06 / U9)")
    print(f"  bdg_ldos.dat:        {len(ldos_rows)} rows, peak_at_zero={peak_zero:.4e}")
    print(f"  bdg_ldos_nambu.dat:  {len(nambu_rows)} rows")
    print(f"  bdg_spectral.dat:    {len(spectral_rows)} rows, max A={max(r[4] for r in spectral_rows):.4e}")


if __name__ == "__main__":
    main()