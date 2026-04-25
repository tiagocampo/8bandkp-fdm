#!/usr/bin/env python3
"""Verify broken-gap QW state ordering near the valence/conduction boundary.

Usage:
    verify_qw_state_character.py <config> <eigenvalues.dat> <parts.dat>
"""

import sys

import numpy as np


def parse_config_int(filepath, key):
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(f"{key}:"):
                return int(float(line.split(":", 1)[1].strip().split()[0]))
    raise RuntimeError(f"{key} not found in {filepath}")


def load_first_k_eigenvalues(filepath):
    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[0, 1:]


def grouped_parts(row):
    row = row / row.sum()
    return {
        "HH": row[0] + row[3],
        "LH": row[1] + row[2],
        "SO": row[4] + row[5],
        "CB": row[6] + row[7],
    }


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <config> <eigenvalues.dat> <parts.dat>")
        sys.exit(1)

    config_path, eigenvalues_path, parts_path = sys.argv[1:]
    numvb = parse_config_int(config_path, "numvb")
    eig = load_first_k_eigenvalues(eigenvalues_path)
    parts = np.loadtxt(parts_path)

    if eig.shape[0] < numvb + 2 or parts.shape[0] < numvb + 2:
        print("FAIL: not enough eigenstates to validate the valence/conduction boundary")
        sys.exit(1)

    state1 = grouped_parts(parts[0])
    cb1 = grouped_parts(parts[numvb])
    cb1_partner = grouped_parts(parts[numvb + 1])

    print(f"numvb = {numvb}")
    print(f"state 1 energy      = {eig[0]:+.6f} eV, grouped parts = {state1}")
    print(f"state {numvb + 1} energy = {eig[numvb]:+.6f} eV, grouped parts = {cb1}")
    print(f"state {numvb + 2} energy = {eig[numvb + 1]:+.6f} eV, grouped parts = {cb1_partner}")

    if state1["HH"] < 0.70:
        print("FAIL: state 1 is not predominantly HH-like")
        sys.exit(1)

    if cb1["CB"] < 0.75:
        print(f"FAIL: state {numvb + 1} is not predominantly conduction-like")
        sys.exit(1)

    if abs(eig[numvb] - eig[numvb + 1]) > 1e-8:
        print("FAIL: first conduction-like Kramers pair is not degenerate at k=0")
        sys.exit(1)

    print("PASS: QW state ordering and conduction-state offset verified")


if __name__ == "__main__":
    main()
