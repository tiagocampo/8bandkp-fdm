#!/usr/bin/env python3
# COVERAGE: observable=minigap geometry=wire material=InAs ref=topological_transition
"""U8 regression: wire BdG open->close->reopen + no auto-window fallback.

Part 1 (lock-in): at mu=0.6601 eV (conduction subband edge, core/shell wire),
transverse B_vec=[Bx,0,0] for Bx in {0, 2.8, 5} T must show the topological
open->close->reopen: minigap(2.8) < 0.5*minigap(0) and minigap(5) > minigap(2.8).

Part 2 (fallback removal, U8): a mu-in-gap run (mu=0.0) must NOT use the
Gershgorin auto-window retry (run_bdg_wire no longer falls back).

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
RE_GAP = re.compile(r"Min gap:\s+([-\d.eE+]+)")
RETRY = "Retrying with auto-computed energy window"
WARN_NONE = "found no eigenvalues in the search window"


def run(text, tag, timeout=900):
    d = Path(f"run_{tag}")
    d.mkdir(parents=True, exist_ok=True)
    (d / "input.toml").write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS=OMP)
    p = subprocess.run([EXE], cwd=d, env=env, capture_output=True,
                       text=True, timeout=timeout)
    (d / "run.log").write_text(p.stdout + "\n--STDERR--\n" + p.stderr)
    return p


def minigap(stdout):
    m = RE_GAP.search(stdout)
    return float(m.group(1)) if m else None


def cfg_with(bx, mu=None):
    t = CONFIG.read_text()
    t = re.sub(r"B_vec\s*=\s*\[[-\d.,\s]+\]",
               f"B_vec = [{bx}, 0.0, 0.0]", t, count=1)
    if mu is not None:
        t = re.sub(r"^mu\s*=\s*[-\d.eE+]+\s*$",
                   f"mu = {mu}", t, count=1, flags=re.MULTILINE)
    return t


def main():
    # --- Part 1: open->close->reopen at mu=0.6601 ---
    gaps = {}
    for bx in (0.0, 2.8, 5.0):
        p = run(cfg_with(bx), f"bx{bx}")
        if p.returncode != 0:
            print(f"FAIL: Bx={bx} exited {p.returncode}")
            return 1
        g = minigap(p.stdout)
        if g is None or g < 0:
            print(f"FAIL: no minigap at Bx={bx} (got {g})")
            return 1
        gaps[bx] = g
    print("minigap(meV): " + "  ".join(
        f"Bx={b}:{gaps[b]*1000:.3f}" for b in sorted(gaps)))
    g0, gc, g5 = gaps[0.0], gaps[2.8], gaps[5.0]
    if not (gc < 0.5 * g0 and g5 > gc):
        print(f"FAIL: no open->close->reopen "
              f"(g0={g0*1000:.3f} gc={gc*1000:.3f} g5={g5*1000:.3f} meV)")
        return 1
    print(f"PASS: open->close->reopen "
          f"(g0={g0*1000:.3f} >= gc={gc*1000:.3f} <= g5={g5*1000:.3f} meV)")

    # --- Part 2: mu-in-gap must NOT auto-window-retry ---
    p = run(cfg_with(0.0, mu=0.0), "mu_gap")
    if RETRY in p.stdout:
        print("FAIL: auto-window fallback still used for mu-in-gap "
              "(U8 fallback removal not applied)")
        return 1
    if WARN_NONE not in p.stdout:
        print("FAIL: expected 'found no eigenvalues' warning for mu-in-gap")
        return 1
    print("PASS: mu-in-gap warns without auto-window fallback")
    return 0


if __name__ == "__main__":
    sys.exit(main())
