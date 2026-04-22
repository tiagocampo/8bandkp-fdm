#!/usr/bin/env python3
"""Verify wire optical transitions are selected from near-gap states.

Usage:
    verify_wire_optical_selection.py <stdout_log> <optical_transitions.dat>
"""

import re
import sys

import numpy as np


def parse_selected_indices(log_path):
    pattern = re.compile(
        r"Selecting VB indices\s+(\d+)\s*:\s*(\d+)\s+and CB indices\s+(\d+)\s*:\s*(\d+)"
    )
    with open(log_path) as f:
        text = f.read()
    match = pattern.search(text)
    if not match:
        raise RuntimeError("selection indices not found in stdout log")
    vb_start, vb_end, cb_start, cb_end = [int(g) for g in match.groups()]
    return vb_start, vb_end, cb_start, cb_end


def parse_transition_energies(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 2]


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <stdout_log> <optical_transitions.dat>")
        sys.exit(1)

    vb_start, vb_end, cb_start, cb_end = parse_selected_indices(sys.argv[1])
    energies = parse_transition_energies(sys.argv[2])

    print(
        f"Selected VB indices {vb_start}:{vb_end}, CB indices {cb_start}:{cb_end}, "
        f"min transition energy {energies.min():.6f} eV"
    )

    if vb_start <= 1:
        print("FAIL: wire selection still starts from the deepest valence state")
        sys.exit(1)

    if energies.min() > 1.0:
        print("FAIL: minimum wire optical transition energy is too large for a near-gap selection")
        sys.exit(1)

    print("PASS: wire optical transitions are selected from the near-gap manifold")


if __name__ == "__main__":
    main()
